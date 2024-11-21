use serde::Serialize;
use crate::digest::digestion::DigestSlice;
use crate::fragment_mass::fragment_mass_builder::SafePosition;
use timsquery::models::aggregators::raw_peak_agg::multi_chromatogram_agg::multi_chromatogram_agg::{NaturalFinalizedMultiCMGStatsArrays, ApexScores};
use timsquery::ElutionGroup;
use std::path::Path;
use csv::Writer;
use std::time::Instant;
use crate::models::DecoyMarking;

#[derive(Debug, Serialize, Clone)]
pub struct PrecursorData {
    pub charge: u8,
    pub mz: f64,
    pub mobility: f32,
    pub rt: f32,
}

#[derive(Debug, Serialize, Clone)]
pub struct IonSearchResults {
    pub sequence: DigestSlice,
    pub score_data: ApexScores,
    pub precursor_data: PrecursorData,
    pub decoy: DecoyMarking,
}

impl IonSearchResults {
    pub fn new(
        digest_sequence: DigestSlice,
        charge: u8,
        elution_group: &ElutionGroup<SafePosition>,
        finalized_scores: NaturalFinalizedMultiCMGStatsArrays<SafePosition>,
        decoy: DecoyMarking,
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
            sequence: digest_sequence,
            score_data,
            precursor_data,
            decoy,
        }
    }

    pub fn get_csv_labels() -> [&'static str; 22] {
        let out = {
            let mut whole: [&'static str; 22] = [""; 22];
            let (id_sec, score_sec) = whole.split_at_mut(6);
            id_sec.copy_from_slice(&Self::get_info_labels());
            score_sec.copy_from_slice(&Self::get_scoring_labels());
            whole
        };
        out
    }

    pub fn as_csv_record(&self) -> [String; 22] {
        let mut out: [String; 22] = core::array::from_fn(|_| "".to_string());
        let lab_sec = self.get_csv_record_lab_sec();
        let mut offset = 0;
        for x in lab_sec.into_iter() {
            out[offset] = x;
            offset += 1;
        }

        let ms1_sec = self.get_csv_record_ms1_score_sec();
        for x in ms1_sec.into_iter() {
            out[offset] = x;
            offset += 1;
        }

        let ms2_sec = self.get_csv_record_ms2_score_sec();
        for x in ms2_sec.into_iter() {
            out[offset] = x;
            offset += 1;
        }

        assert!(offset == 22);
        out
    }

    fn get_info_labels() -> [&'static str; 6] {
        [
            "sequence",
            "precursor_mz",
            "precursor_charge",
            "precursor_mobility_query",
            "precursor_rt_query",
            "decoy",
        ]
    }

    fn get_csv_record_lab_sec(&self) -> [String; 6] {
        [
            self.sequence.clone().into(),
            self.precursor_data.mz.to_string(),
            self.precursor_data.charge.to_string(),
            self.precursor_data.mobility.to_string(),
            self.precursor_data.rt.to_string(),
            self.decoy.as_str().to_string(),
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
}

pub fn write_results_to_csv<P: AsRef<Path>>(
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
