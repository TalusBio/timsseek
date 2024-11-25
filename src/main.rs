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
use timsseek::digest::digestion::{DigestionEnd, DigestionParameters, DigestionPattern};
use timsseek::errors::TimsSeekError;
use timsseek::fragment_mass::elution_group_converter::SequenceToElutionGroupConverter;
use timsseek::fragment_mass::fragment_mass_builder::SafePosition;
use timsseek::protein::fasta::ProteinSequenceCollection;
use timsseek::scoring::search_results::{IonSearchResults, write_results_to_csv};
use timsseek::models::{DigestSlice, deduplicate_digests, NamedQueryChunk};
use core::marker::Send;
use std::sync::Arc;
use rayon::prelude::*;
use timsseek::data_sources::speclib::Speclib;

fn process_chunk<'a>(
    queries: NamedQueryChunk,
    index: &'a QuadSplittedTransposedIndex,
    factory: &'a MultiCMGStatsFactory<SafePosition>,
    tolerance: &'a DefaultTolerance,
) -> Vec<IonSearchResults> {
    let start = Instant::now();
    let num_queries = queries.len();
    let res = query_multi_group(index, tolerance, &queries.queries, &|x| {
        factory.build_with_elution_group(x)
    });
    let elap_time = start.elapsed();
    info!("Querying + Aggregation took {:?}", elap_time);

    let start = Instant::now();

    let tmp: Vec<(IonSearchResults, f64)> = res
        .into_par_iter()
        .zip(queries.into_zip_par_iter())
        .map(|(res_elem, (eg_elem, (digest, charge_elem)))| {
            let decoy = digest.decoy;
            let res = IonSearchResults::new(digest.clone(), charge_elem, &eg_elem, res_elem, decoy);
            if res.is_err() {
                log::error!(
                    "Error creating Digest: {:#?} \nElutionGroup: {:#?}\n Error: {:?}",
                    digest,
                    eg_elem,
                    res,
                );
                return None;
            }
            let res = res.unwrap();
            let main_score = res.score_data.main_score;
            Some((res, main_score))
        })
        .flatten()
        .collect();

    if tmp.is_empty() {
        // TODO: Remove this and check the error elsewhere.
        panic!("No results found");
    }

    let (out, main_scores): (Vec<IonSearchResults>, Vec<f64>) = tmp.into_iter().unzip();

    let avg_main_scores = main_scores.iter().sum::<f64>() / main_scores.len() as f64;

    assert!(!avg_main_scores.is_nan());
    let elapsed = start.elapsed();
    log::info!(
        "Bundling took {:?} for {} elution_groups",
        elapsed,
        num_queries,
    );
    log::info!("Avg main score: {:?}", avg_main_scores);

    out
}

struct DigestedSequenceIterator {
    digest_sequences: Vec<DigestSlice>,
    chunk_size: usize,
    max_iterations: usize,
    iteration_index: usize,
    converter: SequenceToElutionGroupConverter,
    build_decoys: bool,
}

impl DigestedSequenceIterator {
    fn new(
        digest_sequences: Vec<DigestSlice>,
        chunk_size: usize,
        converter: SequenceToElutionGroupConverter,
        build_decoys: bool,
    ) -> Self {
        let max_iterations = digest_sequences.len() / chunk_size;
        Self {
            digest_sequences,
            chunk_size,
            max_iterations,
            converter,
            iteration_index: 0,
            build_decoys,
        }
    }

    fn get_chunk_digests(&self, chunk_index: usize) -> &[DigestSlice] {
        let start = chunk_index * self.chunk_size;
        let end = start + self.chunk_size;
        let end = if end > self.digest_sequences.len() {
            self.digest_sequences.len()
        } else {
            end
        };
        &self.digest_sequences[start..end]
    }

    fn get_chunk(&self, chunk_index: usize) -> NamedQueryChunk {
        let seqs = self.get_chunk_digests(chunk_index);
        let (eg_seq, eg_chunk, charge_chunk) = self.converter.convert_sequences(seqs).unwrap();
        let eg_seq = eg_seq.into_iter().cloned().collect();
        NamedQueryChunk::new(eg_seq, charge_chunk, eg_chunk)
    }

    fn get_decoy_chunk(&self, chunk_index: usize) -> NamedQueryChunk {
        let seqs = self.get_chunk_digests(chunk_index);
        let decoys = seqs
            .iter()
            .map(|x| x.as_decoy())
            .enumerate()
            .collect::<Vec<(usize, DigestSlice)>>();
        // NOTE: RN I am not checking if the decoy is also a target ... bc its hard ...
        // .filter(|(_i, x)| !self.digest_sequences.contains(&x.as_str()))

        let (eg_seq, eg_chunk, charge_chunk) = self
            .converter
            .convert_enumerated_sequences(&decoys)
            .unwrap();
        let eg_seq = eg_seq.into_iter().cloned().collect();
        NamedQueryChunk::new(eg_seq, charge_chunk, eg_chunk)
    }
}

impl Iterator for DigestedSequenceIterator {
    type Item = NamedQueryChunk;

    fn next(&mut self) -> Option<Self::Item> {
        // If its an even iteration, we return the targets.
        // And if its an odd iteration, we return the decoys.
        // IF the struct is requested to build decoys.
        let mut decoy_batch = false;
        let index_use;
        if self.build_decoys {
            index_use = self.iteration_index / 2;
            let decoy_index = self.iteration_index % 2;
            if decoy_index == 1 {
                decoy_batch = true;
            }
            self.iteration_index += 1;
        } else {
            index_use = self.iteration_index;
            self.iteration_index += 1;
        }

        let out = if decoy_batch {
            self.get_decoy_chunk(index_use)
        } else {
            self.get_chunk(index_use)
        };

        if out.is_empty() {
            None
        } else {
            Some(out)
        }
    }
}

impl ExactSizeIterator for DigestedSequenceIterator {
    fn len(&self) -> usize {
        let num_chunks = self.digest_sequences.len() / self.chunk_size;
        if self.build_decoys {
            num_chunks * 2
        } else {
            num_chunks
        }
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
    let tot_chunks = chunked_query_iterator.len();
    let mut chunk_num = 0;

    chunked_query_iterator.for_each(|chunk| {
        log::info!("Chunk {}/{}", chunk_num, tot_chunks);
        let out = process_chunk(chunk, &index, &factory, &tolerance);
        println!("Chunk -Targets {}/{}", chunk_num, tot_chunks);

        let out_path = out_path.join(format!("chunk_{}.csv", chunk_num));
        write_results_to_csv(&out, out_path).unwrap();

        let first_target = out[0].clone();
        println!("First query in chunk: {:#?}", first_target);
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
    let fasta_location = "/Users/sebastianpaez/git/ionmesh/benchmark/UP000005640_9606.fasta";
    let speclib_loc = "/Users/sebastianpaez/git/timsseek/speclib_build/FUUUUU.ndjson";

    // Only HeLa proteins fasta
    // let fasta_location = "/Users/sebastianpaez/git/timsseek/data/HeLa_cannonical_proteins.fasta";
    // let speclib_loc = "/Users/sebastianpaez/git/timsseek/speclib_build/FUUUUU_small.ndjson";

    let speclib = Speclib::from_ndjson_file(&Path::new(speclib_loc))?;

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
    let sequences: Vec<Arc<str>> = fasta_proteins
        .sequences
        .iter()
        .map(|x| x.sequence.clone())
        .collect();

    let start = Instant::now();
    let digest_sequences: Vec<DigestSlice> =
        deduplicate_digests(digestion_params.digest_multiple(&sequences));

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
    let speclib_target_dir = std::path::Path::new("./speclib_results/");
    if !speclib_target_dir.exists() {
        std::fs::create_dir(speclib_target_dir)?;
    }

    const CHUNK_SIZE: usize = 1000;
    const BUILD_DECOYS: bool = true;
    let chunked_query_iterator =
        DigestedSequenceIterator::new(digest_sequences, CHUNK_SIZE, def_converter, BUILD_DECOYS);

    let speclib_iter = speclib.as_iterator(CHUNK_SIZE);
    main_loop(
        speclib_iter,
        &index,
        &factory,
        &tolerance,
        &speclib_target_dir,
    )?;

    todo!();
    main_loop(
        chunked_query_iterator,
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
