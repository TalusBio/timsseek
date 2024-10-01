use log::{info, warn};
use std::collections::{BTreeMap, HashMap};
use std::ops::Deref;
use std::time::Instant;

use timsquery::models::aggregators::raw_peak_agg::multi_chromatogram_agg::{
    FinalizedMultiCMGStatsArrays, NaturalFinalizedMultiCMGStatsArrays,
};
use timsquery::models::aggregators::{self, MultiCMGStatsFactory};
use timsquery::models::indices::transposed_quad_index::QuadSplittedTransposedIndex;
use timsquery::queriable_tims_data::queriable_tims_data::query_multi_group;
use timsquery::traits::tolerance::{
    DefaultTolerance, MobilityTolerance, MzToleramce, QuadTolerance, RtTolerance,
};
use timsquery::ElutionGroup;
use timsseek::fragment_mass::elution_group_converter::SequenceToElutionGroupConverter;

use timsseek::digest::digestion::as_decoy_string;
use timsseek::digest::digestion::{Digest, DigestionEnd, DigestionParameters, DigestionPattern};
use timsseek::fragment_mass::fragment_mass_builder::FragmentMassBuilder;
use timsseek::fragment_mass::fragment_mass_builder::SafePosition;
use timsseek::protein::fasta::ProteinSequenceCollection;

use rayon::prelude::*;
use serde::ser::SerializeMap;
use serde::Serialize;

use std::collections::HashSet;

// timseek --fasta asdasdad --config asdad.json --out_dir outputs # should generate a
// 'results.sqlite' with 1 score per unique peptide+charge in the fasta file.
// "Score" is used loosely here, its the combination of score + RT
//
//
//

use csv::Writer;
use std::path::Path;

fn write_results_to_csv<P: AsRef<Path>>(
    results: &[IonSearchResults],
    out_path: P,
) -> Result<(), Box<dyn std::error::Error>> {
    let start = Instant::now();
    let mut writer = Writer::from_path(out_path.as_ref())?;

    // Write the headers
    writer.write_record(&[
        "sequence",
        "precursor_charge",
        "lazy_hyperscore",
        "lazy_hyperscore_vs_baseline",
        "apex_npeaks",
        "apex_intensity",
        "apex_rt_milliseconds",
        "avg_mobility_at_apex",
    ])?;

    for result in results {
        writer.serialize((
            result.sequence,
            result.elution_group.precursor_charge,
            result.score_data.lazy_hyperscore,
            result.score_data.lazy_hyperscore_vs_baseline,
            result.score_data.apex_npeaks,
            result.score_data.apex_intensity,
            result.score_data.apex_rt_miliseconds,
            result.score_data.avg_mobility_at_apex,
        ))?;
    }
    writer.flush()?;
    log::info!(
        "Writing took {:?} -> {:?}",
        start.elapsed(),
        out_path.as_ref()
    );
    Ok(())
}

#[derive(Debug, Serialize)]
struct ScoreData {
    lazy_hyperscore: f64,
    lazy_hyperscore_vs_baseline: f64,
    apex_npeaks: usize,
    apex_intensity: u64,
    apex_rt_miliseconds: u32,
    avg_mobility_at_apex: f64,

    #[serde(skip_serializing)]
    mzs_at_apex: BTreeMap<SafePosition, f64>,
}

impl From<NaturalFinalizedMultiCMGStatsArrays<SafePosition>> for ScoreData {
    fn from(x: NaturalFinalizedMultiCMGStatsArrays<SafePosition>) -> Self {
        let apex_index = x.apex_hyperscore_index;
        let lazy_hyperscore = x.lazy_hyperscore[apex_index];
        let lazy_hyperscore_vs_baseline = x.lazy_hyperscore_vs_baseline[apex_index];
        let apex_npeaks = x.npeaks[apex_index];
        let apex_intensity = x.summed_intensity[apex_index];
        let apex_rt_miliseconds = x.retention_time_miliseconds[apex_index];
        let avg_mobility_at_apex = x.average_mobility[apex_index];
        let mzs_at_apex = x
            .transition_mzs
            .iter()
            .map(|(lab, mzs_vec)| (lab.clone(), mzs_vec[apex_index]))
            .collect();

        ScoreData {
            lazy_hyperscore,
            lazy_hyperscore_vs_baseline,
            apex_npeaks,
            apex_intensity,
            apex_rt_miliseconds,
            avg_mobility_at_apex,
            mzs_at_apex,
        }
    }
}

#[derive(Debug, Serialize)]
struct IonSearchResults<'a> {
    sequence: &'a str,
    score_data: ScoreData,
    #[serde(skip_serializing)]
    elution_group: &'a ElutionGroup<SafePosition>,
}

impl<'a> IonSearchResults<'a> {
    fn new(
        sequence: &'a str,
        elution_group: &'a ElutionGroup<SafePosition>,
        finalized_scores: NaturalFinalizedMultiCMGStatsArrays<SafePosition>,
    ) -> Self {
        let score_data = ScoreData::from(finalized_scores);
        Self {
            sequence,
            elution_group,
            score_data,
        }
    }
}

fn process_chunk<'a>(
    eg_chunk: &'a [ElutionGroup<SafePosition>],
    seq_chunk: &'a [&str],
    index: &'a QuadSplittedTransposedIndex,
    factory: &'a MultiCMGStatsFactory<SafePosition>,
    tolerance: &'a DefaultTolerance,
) -> Vec<IonSearchResults<'a>> {
    let start = Instant::now();
    let res = query_multi_group(index, index, tolerance, &eg_chunk, &|x| factory.build(x));
    let elap_time = start.elapsed();
    info!("Querying + Aggregation took {:?}", elap_time);

    let start = Instant::now();

    let (out, hyperscores): (Vec<IonSearchResults>, Vec<f64>) = res
        .into_par_iter()
        .zip(seq_chunk.par_iter())
        .zip(eg_chunk.par_iter())
        .map(|((res_elem, seq_elem), eg_elem)| {
            let res = IonSearchResults::new(seq_elem, eg_elem, res_elem);
            let hyperscore = res.score_data.lazy_hyperscore;
            (res, hyperscore)
        })
        .unzip();

    // TODO: make this a compiler flag
    // let (out, hyperscores): (Vec<IonSearchResults>, Vec<f64>) = res
    //     .into_iter()
    //     .zip(seq_chunk.iter())
    //     .zip(eg_chunk.iter())
    //     .map(|((res_elem, seq_elem), eg_elem)| {
    //         let res = IonSearchResults::new(seq_elem, eg_elem, res_elem);
    //         let hyperscore = res.score_data.lazy_hyperscore;
    //         (res, hyperscore)
    //     })
    //     .unzip();

    let avg_hyperscore_vs_baseline = hyperscores.iter().sum::<f64>() / hyperscores.len() as f64;
    let elapsed = start.elapsed();
    log::info!(
        "Bundling took {:?} for {} elution_groups",
        elapsed,
        eg_chunk.len()
    );
    println!(
        "Avg hyperscore vs baseline: {:?}",
        avg_hyperscore_vs_baseline
    );
    out
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    // Initialize logging
    env_logger::init();

    // let fasta_location = "/Users/sebastianpaez/Downloads/UP000000625_83333.fasta";

    // Covid 19 fasta
    // let fasta_location = "/Users/sebastianpaez/Downloads/UP000464024_2697049.fasta";

    // Human
    let fasta_location = "/Users/sebastianpaez/git/ionmesh/benchmark/UP000005640_9606.fasta";

    let digestion_params = DigestionParameters {
        min_length: 6,
        max_length: 20,
        pattern: DigestionPattern::trypsin(),
        digestion_end: DigestionEnd::CTerm,
        max_missed_cleavages: 0,
    };
    let fasta_proteins = ProteinSequenceCollection::from_fasta_file(&fasta_location)?;
    let sequences: Vec<&str> = fasta_proteins
        .sequences
        .iter()
        .map(|x| x.sequence.as_str())
        .collect();

    let start = Instant::now();
    let digest_sequences: HashSet<&str> = digestion_params
        .digest_multiple(&sequences)
        .iter()
        .map(|x| x.sequence)
        .collect();

    let digest_sequences: Vec<&str> = digest_sequences.into_iter().collect();

    let elap_time = start.elapsed();
    let time_per_digest = elap_time / digest_sequences.len() as u32;
    info!("Digestion took {:?}", elap_time);
    info!("Time per digestion: {:?}", time_per_digest);
    info!(
        "Digests per second: {:?}",
        digest_sequences.len() as f32 / time_per_digest.as_secs_f32()
    );

    let def_converter = SequenceToElutionGroupConverter::default();

    // ====================
    let dotd_file_location =
        "/Users/sebastianpaez/git/ionmesh/benchmark/240402_PRTC_01_S1-A1_1_11342.d";

    let tolerance = DefaultTolerance {
        ms: MzToleramce::Ppm((50.0, 50.0)),
        rt: RtTolerance::None,
        mobility: MobilityTolerance::Pct((15.0, 15.0)),
        quad: QuadTolerance::Absolute((0.1, 0.1, 1)),
    };
    let index = QuadSplittedTransposedIndex::from_path(dotd_file_location)?;
    let factory = MultiCMGStatsFactory {
        converters: (index.mz_converter, index.im_converter),
        _phantom: std::marker::PhantomData::<SafePosition>,
    };

    let start = Instant::now();
    const CHUNK_SIZE: usize = 20000;
    let tot_elems = digest_sequences.len();
    let tot_chunks = tot_elems / CHUNK_SIZE;
    let mut chunk_num = 0;

    let target_dir = std::path::Path::new("./results/");
    if !target_dir.exists() {
        std::fs::create_dir(target_dir)?;
    }

    digest_sequences.chunks(CHUNK_SIZE).for_each(|seq_chunk| {
        log::info!("Chunk - Targets {}/{}", chunk_num, tot_chunks);
        let eg_chunk = def_converter.convert_sequences(seq_chunk);
        let out = process_chunk(&eg_chunk, seq_chunk, &index, &factory, &tolerance);
        println!("{:?}", out[0]);
        println!("Chunk -Targets {}/{}", chunk_num, tot_chunks);

        let start = Instant::now();
        let out_path = target_dir.join(format!("targets_chunk_{}.json", chunk_num));
        // TODO replace this unwrap with some nice error handling ...
        let mut out_file = std::fs::File::create(out_path.clone()).unwrap();
        serde_json::to_writer_pretty(&mut out_file, &out).unwrap();
        log::info!("Writing took {:?} -> {:?}", start.elapsed(), out_path);

        let out_path = target_dir.join(format!("targets_chunk_{}.csv", chunk_num));
        write_results_to_csv(&out, out_path).unwrap();
        log::info!("Chunk - Decoys {}/{}", chunk_num, tot_chunks);

        let decoys: Vec<String> = seq_chunk.iter().map(|x| as_decoy_string(x)).collect();
        let decoy_str_slc = decoys.iter().map(|x| x.as_str()).collect::<Vec<&str>>();
        let decoy_eg_chunk = def_converter.convert_sequences(&decoy_str_slc);
        let out = process_chunk(
            &decoy_eg_chunk,
            &decoy_str_slc,
            &index,
            &factory,
            &tolerance,
        );
        println!("{:?}", out[0]);
        println!("Chunk - Decoys {}/{}", chunk_num, tot_chunks);

        let start = Instant::now();
        let out_path = target_dir.join(format!("decoys_chunk_{}.json", chunk_num));
        // TODO replace this unwrap with some nice error handling ...
        let mut out_file = std::fs::File::create(out_path.clone()).unwrap();
        serde_json::to_writer_pretty(&mut out_file, &out).unwrap();
        log::info!("Writing took {:?} -> {:?}", start.elapsed(), out_path);

        let out_path = target_dir.join(format!("decoys_chunk_{}.csv", chunk_num));
        write_results_to_csv(&out, out_path).unwrap();

        chunk_num += 1;
    });
    let elap_time = start.elapsed();
    info!("Querying took {:?}", elap_time);

    // println!("{:?}", out);
    Ok(())
}
