use log::info;
use rayon::prelude::*;
use std::collections::HashSet;
use std::path::Path;
use std::time::Instant;
use timsquery::models::aggregators::raw_peak_agg::multi_chromatogram_agg::multi_chromatogram_agg::{NaturalFinalizedMultiCMGStatsArrays, ApexScores};
use timsquery::models::aggregators::MultiCMGStatsFactory;
use timsquery::models::indices::transposed_quad_index::QuadSplittedTransposedIndex;
use timsquery::queriable_tims_data::queriable_tims_data::query_multi_group;
use timsquery::traits::tolerance::{
    DefaultTolerance, MobilityTolerance, MzToleramce, QuadTolerance, RtTolerance,
};
use timsquery::ElutionGroup;
use timsseek::digest::digestion::{as_decoy_string, DigestSlice};
use timsseek::digest::digestion::{DigestionEnd, DigestionParameters, DigestionPattern};
use timsseek::errors::TimsSeekError;
use timsseek::fragment_mass::elution_group_converter::SequenceToElutionGroupConverter;
use timsseek::fragment_mass::fragment_mass_builder::SafePosition;
use timsseek::protein::fasta::ProteinSequenceCollection;
use timsseek::scoring::search_results::{IonSearchResults, write_results_to_csv};
use timsseek::models::DecoyMarking;

fn process_chunk<'a>(
    eg_chunk: &'a [ElutionGroup<SafePosition>],
    crg_chunk: &'a [u8],
    seq_chunk: &'a [&str],
    index: &'a QuadSplittedTransposedIndex,
    factory: &'a MultiCMGStatsFactory<SafePosition>,
    tolerance: &'a DefaultTolerance,
    decoy: DecoyMarking,
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
            let res = IonSearchResults::new(seq_elem, *crg_elem, eg_elem, res_elem, decoy);
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

#[derive(Debug, Clone)]
struct NamedQueryChunk {
    names: Vec<DigestSlice>,
    decoy: DecoyMarking,
    charges: Vec<u8>,
    queries: Vec<ElutionGroup<SafePosition>>,
}

impl NamedQueryChunk {
    fn new(
        names: Vec<DigestSlice>,
        decoy: DecoyMarking,
        charges: Vec<u8>,
        queries: Vec<ElutionGroup<SafePosition>>,
    ) -> Self {
        Self {
            names,
            decoy,
            charges,
            queries,
        }
    }
}

struct DigestedSequenceIterator {
    digest_sequences: Vec<DigestSlice>,
    chunk_size: usize,
    max_iterations: usize,
    iteration_index: usize,
    converter: SequenceToElutionGroupConverter,
}

impl DigestedSequenceIterator {
    fn new(
        digest_sequences: Vec<DigestSlice>,
        chunk_size: usize,
        converter: SequenceToElutionGroupConverter,
    ) -> Self {
        let max_iterations = digest_sequences.len() / chunk_size;
        Self {
            digest_sequences,
            chunk_size,
            max_iterations,
            converter,
            iteration_index: 0,
        }
    }

    fn get_chunk_str(&self, chunk_index: usize) -> &[DigestSlice] {
        let start = chunk_index * self.chunk_size;
        let end = start + self.chunk_size;
        &self.digest_sequences[start..end]
    }

    fn target_chunk(&self, chunk_index: usize) -> NamedQueryChunk {
        let decoy = DecoyMarking::Target;
        let seqs = self.get_chunk_str(chunk_index);
        let (eg_seq, eg_chunk, charge_chunk) = self.converter.convert_sequences(seqs).unwrap();
        NamedQueryChunk::new(eg_seq, decoy, charge_chunk, eg_chunk)
    }

    fn decoy_chunk(&self, chunk_index: usize) -> NamedQueryChunk {
        let decoy = DecoyMarking::Decoy;
        let seqs = self.get_chunk_str(chunk_index);
        let decoys = seqs
            .iter()
            .map(|x| as_decoy_string(x))
            .enumerate()
            .filter(|(_i, x)| !self.digest_sequences.contains(&x.as_str()))
            .collect::<Vec<(usize, String)>>();
        let (eg_seq, eg_chunk, charge_chunk) = self
            .converter
            .convert_enumerated_sequences(&decoys)
            .unwrap();
        NamedQueryChunk::new(eg_seq, decoy, charge_chunk, eg_chunk)
    }
}

impl<'a> Iterator for DigestedSequenceIterator {
    type Item = NamedQueryChunk;

    fn next(&mut self) -> Option<Self::Item> {
        // If its an even iteration, we return the targets.
        // And if its an odd iteration, we return the decoys.
    }
}

fn main_loop<'a>(
    chunked_query_iterator: impl ExactSizeIterator<Item = NamedQueryChunk>,
    // def_converter: &SequenceToElutionGroupConverter,
    index: &'a QuadSplittedTransposedIndex,
    factory: &'a MultiCMGStatsFactory<SafePosition>,
    tolerance: &'a DefaultTolerance,
    out_path: &Path,
) -> std::result::Result<(), TimsSeekError> {
    let start = Instant::now();
    let tot_chunks = chunked_query_iterator.len();
    let mut chunk_num = 0;

    chunked_query_iterator.for_each(|chunk| {
        log::info!("Chunk - Targets {}/{}", chunk_num, tot_chunks);
        let out = process_chunk(
            &chunk.queries,
            &chunk.charges,
            &chunk.names,
            &index,
            &factory,
            &tolerance,
            chunk.decoy,
        );
        let first_target = out[0].clone();
        println!("Chunk -Targets {}/{}", chunk_num, tot_chunks);

        let out_path = out_path.join(format!("targets_chunk_{}.csv", chunk_num));
        write_results_to_csv(&out, out_path).unwrap();
        // log::info!("Chunk - Decoys {}/{}", chunk_num, tot_chunks);

        // let decoys: Vec<String> = seq_chunk.iter().map(|x| as_decoy_string(x)).collect();
        // let enum_decoy_str_slc = decoys
        //     .iter()
        //     .enumerate()
        //     .map(|(i, x)| (i, x.as_str()))
        //     .filter(|(_i, x)| !digest_sequences.contains(x))
        //     .collect::<Vec<(usize, &str)>>();

        // log::info!("Number of Decoys to process: {}", enum_decoy_str_slc.len());

        // let decoy_str_slc = enum_decoy_str_slc
        //     .iter()
        //     .map(|x| x.1)
        //     .collect::<Vec<&str>>();
        // let (decoy_eg_chunk, decoy_crg_chunk) = def_converter
        //     .convert_enumerated_sequences(&enum_decoy_str_slc)
        //     .unwrap();
        // let out = process_chunk(
        //     &decoy_eg_chunk,
        //     &decoy_crg_chunk,
        //     &decoy_str_slc,
        //     &index,
        //     &factory,
        //     &tolerance,
        // );
        // let first_decoy = out[0].clone();
        // println!("Chunk - Decoys {}/{}", chunk_num, tot_chunks);

        // let out_path = out_path.join(format!("decoys_chunk_{}.csv", chunk_num));
        // write_results_to_csv(&out, out_path).unwrap();

        println!("First target in chunk: {:#?}", first_target);
        // println!("First decoy in chunk: {:#?}", first_decoy);

        chunk_num += 1;
    });
    Ok(())
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
        .map(|x| Into::<String>::into(x))
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

    let target_dir = std::path::Path::new("./results/");
    if !target_dir.exists() {
        std::fs::create_dir(target_dir)?;
    }

    main_loop(
        &chunked_query_iterator,
        &index,
        &factory,
        &tolerance,
        &target_dir,
    )?;
    let elap_time = start.elapsed();
    info!("Querying took {:?}", elap_time);

    // println!("{:?}", out);
    Ok(())
}
