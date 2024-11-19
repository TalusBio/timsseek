use log::info;
use rayon::prelude::*;
use serde::Serialize;
use std::collections::HashSet;
use std::collections::BTreeMap;
use std::time::Instant;
use timsquery::models::aggregators::raw_peak_agg::multi_chromatogram_agg::multi_chromatogram_agg::{NaturalFinalizedMultiCMGStatsArrays, ApexScores};
use timsquery::models::aggregators::MultiCMGStatsFactory;
use timsquery::models::indices::transposed_quad_index::QuadSplittedTransposedIndex;
use timsquery::queriable_tims_data::queriable_tims_data::query_multi_group;
use timsquery::traits::tolerance::{
    DefaultTolerance, MobilityTolerance, MzToleramce, QuadTolerance, RtTolerance,
};
use timsquery::ElutionGroup;
use timsseek::digest::digestion::as_decoy_string;
use timsseek::digest::digestion::{DigestionEnd, DigestionParameters, DigestionPattern};
use timsseek::errors::TimsSeekError;
use timsseek::fragment_mass::elution_group_converter::SequenceToElutionGroupConverter;
use timsseek::fragment_mass::fragment_mass_builder::SafePosition;
use timsseek::protein::fasta::ProteinSequenceCollection;

use csv::Writer;
use std::path::Path;

fn write_results_to_csv<P: AsRef<Path>>(
    results: &[IonSearchResults],
    out_path: P,
) -> std::result::Result<(), Box<dyn std::error::Error>> {
    let start = Instant::now();
    let mut writer = Writer::from_path(out_path.as_ref())?;

    // Write the headers
    writer.write_record(IonSearchResults::get_csv_labels())?;

    for result in results {
        let record = result.as_csv_record();
        writer.write_record(&record)?;
    }
    writer.flush()?;
    log::info!(
        "Writing took {:?} -> {:?}",
        start.elapsed(),
        out_path.as_ref()
    );
    Ok(())
}

type MzDeltaPair = (f64, f64);

#[derive(Debug, Serialize, Clone)]
struct ScoreData {
    lazy_hyperscore: f64,
    lazy_hyperscore_vs_baseline: f64,
    norm_hyperscore_vs_baseline: f64,
    lazier_score: f64,
    lazier_score_vs_baseline: f64,
    norm_lazyerscore_vs_baseline: f64,
    apex_npeaks: usize,
    apex_intensity: u64,
    apex_rt_miliseconds: u32,
    avg_mobility_at_apex: f64,

    #[serde(skip_serializing)]
    mzs_at_apex: BTreeMap<SafePosition, MzDeltaPair>,
}

// impl ScoreData {
//     fn new(
//         x: NaturalFinalizedMultiCMGStatsArrays<SafePosition>,
//         elution_group: &ElutionGroup<SafePosition>,
//     ) -> Self {
//         let apex_index = x.ms2_stats.apex_primary_score_index;
//         let lazy_hyperscore = x.ms2_stats.scores.lazy_hyperscore[apex_index];
//         let lazy_hyperscore_vs_baseline =
//             x.ms2_stats.scores.lazy_hyperscore_vs_baseline[apex_index];
//         let norm_hyperscore_vs_baseline =
//             x.ms2_stats.scores.norm_hyperscore_vs_baseline[apex_index];
//
//         let lazier_score = x.ms2_stats.scores.lazyerscore[apex_index];
//         let lazier_score_vs_baseline = x.ms2_stats.scores.lazyerscore_vs_baseline[apex_index];
//         let norm_lazyerscore_vs_baseline =
//             x.ms2_stats.scores.norm_lazyerscore_vs_baseline[apex_index];
//
//         let apex_npeaks = x.ms2_stats.scores.npeaks[apex_index];
//         let apex_intensity = x.ms2_stats.summed_intensity[apex_index];
//         let apex_rt_miliseconds = x.ms2_stats.retention_time_miliseconds[apex_index];
//         let avg_mobility_at_apex = x.ms2_stats.average_mobility[apex_index];
//
//         // let mut mzs_at_apex = BTreeMap::new();
//         // for (lab, theo_mz) in elution_group.fragment_mzs.iter() {
//         //     let obs_mz = match x.ms2_stats.transition_mzs.get(lab) {
//         //         Some(x) => x[apex_index],
//         //         None => f64::NAN,
//         //     }; //[lab][apex_index];
//         //     mzs_at_apex.insert(lab.clone(), (*theo_mz, *theo_mz - obs_mz));
//         // }
//
//         ScoreData {
//             lazy_hyperscore,
//             lazy_hyperscore_vs_baseline,
//             norm_hyperscore_vs_baseline,
//             lazier_score,
//             lazier_score_vs_baseline,
//             norm_lazyerscore_vs_baseline,
//             apex_npeaks,
//             apex_intensity,
//             apex_rt_miliseconds,
//             avg_mobility_at_apex,
//             mzs_at_apex,
//         }
//     }
// }

#[derive(Debug, Serialize, Clone)]
struct PrecursorData {
    pub charge: u8,
    pub mz: f64,
    pub mobility: f32,
    pub rt: f32,
}

#[derive(Debug, Serialize, Clone)]
struct IonSearchResults<'a> {
    sequence: &'a str,
    // score_data: ScoreData,
    score_data: ApexScores,
    precursor_data: PrecursorData,
}

impl<'a> IonSearchResults<'a> {
    fn new(
        sequence: &'a str,
        charge: u8,
        elution_group: &'a ElutionGroup<SafePosition>,
        finalized_scores: NaturalFinalizedMultiCMGStatsArrays<SafePosition>,
    ) -> Self {
        // let score_data = ScoreData::new(finalized_scores, elution_group);
        let score_data = finalized_scores.finalized_score().unwrap();
        let precursor_data = PrecursorData {
            charge,
            mz: elution_group.precursor_mzs[0],
            mobility: elution_group.mobility,
            rt: elution_group.rt_seconds,
        };

        Self {
            sequence,
            score_data,
            precursor_data,
        }
    }

    fn get_info_labels() -> [&'static str; 5] {
        [
            "sequence",
            "precursor_mz",
            "precursor_charge",
            "precursor_mobility_query",
            "precursor_rt_query",
        ]
    }

    fn get_ms2_scoring_labels() -> [&'static str; 11] {
        [
            // Combined
            "lazyerscore",
            "lazyerscore_vs_baseline",
            "norm_lazyerscore_vs_baseline",
            "cosine_similarity",
            "npeaks",
            "summed_transition_intensity",
            "rt_ms",
            // MS2 - Split
            "ms2_mz_errors",
            "ms2_mobility_errors",
            "ms2_intensity",
            "main_score",
        ]
    }

    fn get_ms1_scoring_labels() -> [&'static str; 5] {
        [
            "ms1_cosine_similarity",
            "ms1_summed_precursor_intensity",
            "ms1_mz_errors",
            "ms1_mobility_errors",
            "ms1_intensity",
        ]
    }

    fn get_scoring_labels() -> [&'static str; 16] {
        let mut out: [&'static str; 16] = [""; 16];
        let (id_sec, score_sec) = out.split_at_mut(5);
        id_sec.copy_from_slice(&Self::get_ms1_scoring_labels());
        score_sec.copy_from_slice(&Self::get_ms2_scoring_labels());
        out
    }

    fn get_csv_labels() -> [&'static str; 21] {
        let out = {
            let mut whole: [&'static str; 21] = [""; 21];
            let (id_sec, score_sec) = whole.split_at_mut(5);
            id_sec.copy_from_slice(&Self::get_info_labels());
            score_sec.copy_from_slice(&Self::get_scoring_labels());
            whole
        };
        out
    }

    fn get_csv_record_lab_sec(&self) -> [String; 5] {
        [
            self.sequence.to_string(),
            self.precursor_data.mz.to_string(),
            self.precursor_data.charge.to_string(),
            self.precursor_data.mobility.to_string(),
            self.precursor_data.rt.to_string(),
        ]
    }

    fn get_csv_record_ms1_score_sec(&self) -> [String; 5] {
        let fmt_mz_errors = format!("{:?}", self.score_data.ms1_scores.mz_errors.clone());
        let fmt_mobility_errors =
            format!("{:?}", self.score_data.ms1_scores.mobility_errors.clone());
        let fmt_intensity = format!("{:?}", self.score_data.ms1_scores.transition_intensities);

        [
            self.score_data.ms1_scores.cosine_similarity.to_string(),
            self.score_data.ms1_scores.summed_intensity.to_string(),
            fmt_mz_errors,
            fmt_mobility_errors,
            fmt_intensity,
        ]
    }

    fn get_csv_record_ms2_score_sec(&self) -> [String; 11] {
        let fmt_mz_errors = format!("{:?}", self.score_data.ms2_scores.mz_errors.clone());
        let fmt_mobility_errors =
            format!("{:?}", self.score_data.ms2_scores.mobility_errors.clone());
        let fmt_intensity = format!("{:?}", self.score_data.ms2_scores.transition_intensities);

        [
            self.score_data.ms2_scores.lazyerscore.to_string(),
            self.score_data
                .ms2_scores
                .lazyerscore_vs_baseline
                .to_string(),
            self.score_data
                .ms2_scores
                .norm_lazyerscore_vs_baseline
                .to_string(),
            self.score_data.ms2_scores.cosine_similarity.to_string(),
            self.score_data.ms2_scores.npeaks.to_string(),
            self.score_data.ms2_scores.summed_intensity.to_string(),
            self.score_data
                .ms2_scores
                .retention_time_miliseconds
                .to_string(),
            fmt_mz_errors,
            fmt_mobility_errors,
            fmt_intensity,
            self.score_data.main_score.to_string(),
        ]
    }

    fn as_csv_record(&self) -> [String; 21] {
        let mut out: [String; 21] = core::array::from_fn(|_| "".to_string());
        let lab_sec = self.get_csv_record_lab_sec();
        let mut offset = 0;
        for (i, x) in lab_sec.into_iter().enumerate() {
            out[i + offset] = x;
        }
        offset += 5;

        let ms1_sec = self.get_csv_record_ms1_score_sec();
        for (i, x) in ms1_sec.into_iter().enumerate() {
            out[i + offset] = x;
        }
        offset += 5;

        let ms2_sec = self.get_csv_record_ms2_score_sec();
        for (i, x) in ms2_sec.into_iter().enumerate() {
            out[i + offset] = x;
        }
        offset += 11;

        assert!(offset == 21);
        out
    }
}

fn process_chunk<'a>(
    eg_chunk: &'a [ElutionGroup<SafePosition>],
    crg_chunk: &'a [u8],
    seq_chunk: &'a [&str],
    index: &'a QuadSplittedTransposedIndex,
    factory: &'a MultiCMGStatsFactory<SafePosition>,
    tolerance: &'a DefaultTolerance,
) -> Vec<IonSearchResults<'a>> {
    let start = Instant::now();
    let res = query_multi_group(index, tolerance, eg_chunk, &|x| {
        factory.build_with_elution_group(x)
    });
    let elap_time = start.elapsed();
    info!("Querying + Aggregation took {:?}", elap_time);

    let start = Instant::now();

    let (out, main_scores): (Vec<IonSearchResults>, Vec<f64>) = res
        .into_par_iter()
        .zip(seq_chunk.par_iter().zip(crg_chunk.par_iter()))
        .zip(eg_chunk.par_iter())
        .map(|((res_elem, (seq_elem, crg_elem)), eg_elem)| {
            let res = IonSearchResults::new(seq_elem, *crg_elem, eg_elem, res_elem);
            let main_score = res.score_data.main_score;
            (res, main_score)
        })
        .unzip();

    let avg_main_scores = main_scores.iter().sum::<f64>() / main_scores.len() as f64;

    assert!(!avg_main_scores.is_nan());
    let elapsed = start.elapsed();
    log::info!(
        "Bundling took {:?} for {} elution_groups",
        elapsed,
        eg_chunk.len()
    );
    log::info!("Avg main score: {:?}", avg_main_scores);

    out
}

fn main() -> std::result::Result<(), TimsSeekError> {
    // Initialize logging
    env_logger::init();

    // let fasta_location = "/Users/sebastianpaez/Downloads/UP000000625_83333.fasta";

    // Covid 19 fasta
    // let fasta_location = "/Users/sebastianpaez/Downloads/UP000464024_2697049.fasta";

    // Human
    // let fasta_location = "/Users/sebastianpaez/git/ionmesh/benchmark/UP000005640_9606.fasta";

    // Only HeLa proteins fasta
    let fasta_location = "/Users/sebastianpaez/git/timsseek/data/HeLa_cannonical_proteins.fasta";

    let digestion_params = DigestionParameters {
        min_length: 6,
        max_length: 20,
        pattern: DigestionPattern::trypsin(),
        digestion_end: DigestionEnd::CTerm,
        max_missed_cleavages: 0,
    };
    info!(
        "Digesting {} with parameters: \n {:?}",
        fasta_location, digestion_params
    );
    let fasta_proteins = ProteinSequenceCollection::from_fasta_file(fasta_location)?;
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

    let digest_sequences: Vec<&str> = digest_sequences.iter().cloned().collect();

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
        ms: MzToleramce::Ppm((15.0, 15.0)),
        rt: RtTolerance::None,
        mobility: MobilityTolerance::Pct((10.0, 10.0)),
        quad: QuadTolerance::Absolute((0.1, 0.1)),
    };
    // let index = QuadSplittedTransposedIndex::from_path(dotd_file_location)?;
    let index = QuadSplittedTransposedIndex::from_path_centroided(dotd_file_location)?;

    let factory = MultiCMGStatsFactory {
        converters: (index.mz_converter, index.im_converter),
        _phantom: std::marker::PhantomData::<SafePosition>,
    };

    let start = Instant::now();
    const CHUNK_SIZE: usize = 5000;
    let tot_elems = digest_sequences.len();
    let tot_chunks = tot_elems / CHUNK_SIZE;
    let mut chunk_num = 0;

    let target_dir = std::path::Path::new("./results/");
    if !target_dir.exists() {
        std::fs::create_dir(target_dir)?;
    }

    digest_sequences.chunks(CHUNK_SIZE).for_each(|seq_chunk| {
        log::info!("Chunk - Targets {}/{}", chunk_num, tot_chunks);
        let (eg_chunk, crg_chunk) = def_converter.convert_sequences(seq_chunk).unwrap();
        let out = process_chunk(
            &eg_chunk, &crg_chunk, seq_chunk, &index, &factory, &tolerance,
        );
        let first_target = out[0].clone();
        println!("Chunk -Targets {}/{}", chunk_num, tot_chunks);

        let out_path = target_dir.join(format!("targets_chunk_{}.csv", chunk_num));
        write_results_to_csv(&out, out_path).unwrap();
        log::info!("Chunk - Decoys {}/{}", chunk_num, tot_chunks);

        let decoys: Vec<String> = seq_chunk.iter().map(|x| as_decoy_string(x)).collect();
        let enum_decoy_str_slc = decoys
            .iter()
            .enumerate()
            .map(|(i, x)| (i, x.as_str()))
            .filter(|(_i, x)| !digest_sequences.contains(x))
            .collect::<Vec<(usize, &str)>>();

        log::info!("Number of Decoys to process: {}", enum_decoy_str_slc.len());

        let decoy_str_slc = enum_decoy_str_slc
            .iter()
            .map(|x| x.1)
            .collect::<Vec<&str>>();
        let (decoy_eg_chunk, decoy_crg_chunk) = def_converter
            .convert_enumerated_sequences(&enum_decoy_str_slc)
            .unwrap();
        let out = process_chunk(
            &decoy_eg_chunk,
            &decoy_crg_chunk,
            &decoy_str_slc,
            &index,
            &factory,
            &tolerance,
        );
        let first_decoy = out[0].clone();
        println!("Chunk - Decoys {}/{}", chunk_num, tot_chunks);

        let out_path = target_dir.join(format!("decoys_chunk_{}.csv", chunk_num));
        write_results_to_csv(&out, out_path).unwrap();

        println!("First target in chunk: {:#?}", first_target);
        println!("First decoy in chunk: {:#?}", first_decoy);

        chunk_num += 1;
    });
    let elap_time = start.elapsed();
    info!("Querying took {:?}", elap_time);

    // println!("{:?}", out);
    Ok(())
}
