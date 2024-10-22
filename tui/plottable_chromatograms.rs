use ratatui::layout::{Constraint, Layout};
use ratatui::symbols;
use ratatui::{
    buffer::Buffer,
    layout::{Alignment, Rect},
    style::{Color, Style, Stylize},
    widgets::{block::Title, Axis, Block, Chart, Dataset, GraphType, LegendPosition, Widget},
};
use serde::Serialize;
use std::collections::BTreeMap;
use timsquery::models::aggregators::raw_peak_agg::multi_chromatogram_agg::NaturalFinalizedMultiCMGStatsArrays;
use timsseek::fragment_mass::fragment_mass_builder::SafePosition;

#[derive(Debug, Clone, Serialize)]
pub struct PlottableChromatograms {
    pub arrays: NaturalFinalizedMultiCMGStatsArrays<SafePosition>,
    pub sequence: String,
    pub plottable_state: PlottableChromatogramsState,
}

impl PlottableChromatograms {
    pub fn new(
        arrays: NaturalFinalizedMultiCMGStatsArrays<SafePosition>,
        sequence: String,
    ) -> Self {
        let pstate = get_plottable_state(&arrays, &sequence);
        Self {
            arrays,
            sequence,
            plottable_state: pstate,
        }
    }
}

#[derive(Debug, Clone, Serialize)]
pub struct PlottableChromatogramsState {
    pub inten_data_tuples: BTreeMap<String, Vec<(f64, f64)>>,
    pub hyperscore_data_tuples: BTreeMap<String, Vec<(f64, f64)>>,
    pub lazyscore_data_tuples: BTreeMap<String, Vec<(f64, f64)>>,
    pub normscores_data_tuples: BTreeMap<String, Vec<(f64, f64)>>,
    pub min_max_inten: MinMax,
    pub min_max_hyperscore: MinMax,
    pub min_max_lazyscore: MinMax,
    pub min_max_norm_scores: MinMax,
    pub x_labels: [String; 2],
    pub min_rt: f64,
    pub max_rt: f64,
    pub title: String,
}

const COLOR_CYCLE: [Color; 4] = [Color::Cyan, Color::Magenta, Color::LightCyan, Color::Gray];

#[derive(Debug, Clone, Serialize, Copy)]
pub struct MinMax(f64, f64);

impl MinMax {
    fn as_scinot_labs(&self) -> [String; 2] {
        [format!("{:.2}", self.0), format!("{:.2e}", self.1)]
    }

    fn as_bounds(&self) -> [f64; 2] {
        [self.0, self.1]
    }
}

fn shared_axis_datatuples<'a>(
    xs: &'a [f64],
    ys: &'a BTreeMap<String, &'a [f64]>,
    xmin: f64,
    xmax: f64,
) -> (BTreeMap<String, Vec<(f64, f64)>>, MinMax) {
    let mut data_tuples = BTreeMap::new();
    let mut max_val = f64::NEG_INFINITY;
    let mut min_val = f64::INFINITY;
    for (k, v) in ys.iter() {
        let data_tuple: Vec<(f64, f64)> = xs
            .iter()
            .zip(v.iter())
            .filter_map(|(x, y)| {
                if x < &xmin || x > &xmax {
                    return None;
                }

                Some((*x, *y))
            })
            .collect();

        let local_max = data_tuple
            .iter()
            .map(|(_x, y)| *y)
            .reduce(f64::max)
            .unwrap_or(0.0);
        if local_max > max_val {
            max_val = local_max;
        }

        let local_min = data_tuple
            .iter()
            .map(|(_x, y)| *y)
            .reduce(f64::min)
            .unwrap_or(f64::INFINITY);
        if local_min < min_val {
            min_val = local_min;
        }

        data_tuples.insert(k.clone(), data_tuple);
    }

    (data_tuples, MinMax(min_val, max_val))
}

fn datasets_from_datatuples<'a>(
    data_tuples: &'a BTreeMap<String, Vec<(f64, f64)>>,
    color_cycle: &[Color],
    graph_type: GraphType,
    marker_type: symbols::Marker,
) -> Vec<Dataset<'a>> {
    data_tuples
        .iter()
        .enumerate()
        .map(|(i, (k, v))| {
            let color = color_cycle[i % color_cycle.len()];
            Dataset::default()
                .name(k.clone())
                .data(v.as_slice())
                .marker(marker_type)
                .style(Style::default().fg(color))
                .graph_type(graph_type)
        })
        .collect()
}

const RT_WIDTH_SECONDS: f64 = 8.0 * 20.0;

fn get_plottable_state(
    arrays: &NaturalFinalizedMultiCMGStatsArrays<SafePosition>,
    sequence: &str,
) -> PlottableChromatogramsState {
    let apex = arrays.apex_primary_score_index;
    let hyperscore = arrays.lazy_hyperscore[apex];
    let hyperscore_vs_baseline = arrays.lazy_hyperscore_vs_baseline[apex];
    let norm_hyperscore_vs_baseline = arrays.norm_hyperscore_vs_baseline[apex];
    let lazyscore = arrays.lazyerscore[apex];
    let lazyscore_vs_baseline = arrays.lazyerscore_vs_baseline[apex];
    let norm_lazyscore_vs_baseline = arrays.norm_lazyerscore_vs_baseline[apex];
    let npeaks = arrays.npeaks[apex];
    let intensity = arrays.summed_intensity[apex];

    let title = format!(
        "{}, hyperscore={}/{}/{} lazyscore={}/{}/{} npeaks={} intensity={}",
        sequence,
        hyperscore,
        hyperscore_vs_baseline,
        norm_hyperscore_vs_baseline,
        lazyscore,
        lazyscore_vs_baseline,
        norm_lazyscore_vs_baseline,
        npeaks,
        intensity
    );

    let rts: Vec<f64> = arrays
        .retention_time_miliseconds
        .iter()
        .map(|x| *x as f64 / 1000.0)
        .collect();

    // Pretty sure there is a simpler way to do this but I am tired rn ...
    let apex_rt = rts[apex];
    let mut min_rt = apex_rt - (RT_WIDTH_SECONDS / 2.0);
    min_rt = min_rt.max(rts[0]);
    let mut max_rt = min_rt + RT_WIDTH_SECONDS;
    max_rt = max_rt.min(rts[rts.len() - 1]);
    min_rt = min_rt.min(max_rt - RT_WIDTH_SECONDS);

    // TODO do here a binary search to find the min and max rt
    // indices, then we can use that to convert A LOT less data.

    // OR ... make the conversion and store the converted data.
    // This would be worth it if I implement the zooming in-out-panning functionality.

    let x_labels: [String; 2] = [format!("{:.2}", min_rt), format!("{:.2}", max_rt)];
    let f64_intensities: BTreeMap<String, Vec<f64>> = arrays
        .transition_intensities
        .iter()
        .map(|(k, v)| (format!("{}", k), v.into_iter().map(|x| *x as f64).collect()))
        .collect();

    let mut f64_inten_slices: BTreeMap<String, &[f64]> = f64_intensities
        .iter()
        .map(|(k, v)| (k.clone(), v.as_slice()))
        .collect();

    let slice_intensity_sums: BTreeMap<String, f64> = f64_intensities
        .iter()
        .map(|(k, v)| (k.clone(), v.iter().sum::<f64>()))
        .collect();
    let max_intensity = slice_intensity_sums
        .values()
        .max_by(|x, y| x.partial_cmp(y).unwrap())
        .unwrap();

    // Drop transitions with intensity under 0.1% of the max intensity.
    f64_inten_slices.retain(|_k, v| v.iter().sum::<f64>() > 0.001 * max_intensity);

    let (inten_data_tuples, min_max_inten) =
        shared_axis_datatuples(&rts, &f64_inten_slices, min_rt, max_rt);

    let mut hyperscore_section_data = BTreeMap::new();
    hyperscore_section_data.insert("Hyperscore".to_string(), arrays.lazy_hyperscore.as_slice());
    hyperscore_section_data.insert(
        "Hyperscore vs Baseline".to_string(),
        arrays.lazy_hyperscore_vs_baseline.as_slice(),
    );

    let mut lazyscore_section_data = BTreeMap::new();
    lazyscore_section_data.insert("Lazyscore".to_string(), arrays.lazyerscore.as_slice());
    lazyscore_section_data.insert(
        "Lazyscore vs Baseline".to_string(),
        arrays.lazyerscore_vs_baseline.as_slice(),
    );

    let mut norm_section_data = BTreeMap::new();
    norm_section_data.insert(
        "Norm Lazyscore vs Baseline".to_string(),
        arrays.norm_lazyerscore_vs_baseline.as_slice(),
    );
    norm_section_data.insert(
        "Norm Hyperscore vs Baseline".to_string(),
        arrays.norm_hyperscore_vs_baseline.as_slice(),
    );

    let (lazyscore_data_tuples, min_max_lazyscore) =
        shared_axis_datatuples(&rts, &lazyscore_section_data, min_rt, max_rt);

    let (hyperscore_data_tuples, min_max_hyperscore) =
        shared_axis_datatuples(&rts, &hyperscore_section_data, min_rt, max_rt);

    let (norm_lazyscore_data_tuples, min_max_norm_lazyscore) =
        shared_axis_datatuples(&rts, &norm_section_data, min_rt, max_rt);

    PlottableChromatogramsState {
        inten_data_tuples,
        hyperscore_data_tuples,
        lazyscore_data_tuples,
        normscores_data_tuples: norm_lazyscore_data_tuples,
        min_max_inten,
        min_max_hyperscore,
        min_max_lazyscore,
        min_max_norm_scores: min_max_norm_lazyscore,
        x_labels,
        min_rt,
        max_rt,
        title,
    }
}

impl Widget for PlottableChromatogramsState {
    fn render(self, area: Rect, buf: &mut Buffer) {
        let [top1, top2, bottom1, bottom2] = Layout::vertical([
            Constraint::Percentage(25),
            Constraint::Percentage(25),
            Constraint::Percentage(25),
            Constraint::Percentage(25),
        ])
        .areas(area);

        let datasets_inten = datasets_from_datatuples(
            &self.inten_data_tuples,
            &COLOR_CYCLE,
            GraphType::Line,
            symbols::Marker::Braille,
        );
        let datasets_hyperscore = datasets_from_datatuples(
            &self.hyperscore_data_tuples,
            &COLOR_CYCLE,
            GraphType::Scatter,
            symbols::Marker::Dot,
        );
        let datasets_lazyscore = datasets_from_datatuples(
            &self.lazyscore_data_tuples,
            &COLOR_CYCLE,
            GraphType::Scatter,
            symbols::Marker::Dot,
        );
        let datasets_normscore = datasets_from_datatuples(
            &self.normscores_data_tuples,
            &COLOR_CYCLE,
            GraphType::Scatter,
            symbols::Marker::Dot,
        );

        let ylabels_end = format!("{:.2e}", self.min_max_inten.1);
        let top_chart = Chart::new(datasets_inten)
            .block(
                Block::bordered().title(
                    Title::default()
                        .content(self.title.bold())
                        .alignment(Alignment::Center),
                ),
            )
            .x_axis(
                Axis::default()
                    .bounds([self.min_rt, self.max_rt])
                    .labels(self.x_labels.clone()),
            )
            .y_axis(
                Axis::default()
                    .style(Style::default().gray())
                    .bounds([0.0, self.min_max_inten.1])
                    .labels(["0.0", ylabels_end.as_str()]),
            )
            .legend_position(Some(LegendPosition::TopRight));

        let bottom_chart = Chart::new(datasets_hyperscore)
            .block(
                Block::bordered().title(
                    Title::default()
                        .content("Hyperscore".bold())
                        .alignment(Alignment::Center),
                ),
            )
            .x_axis(
                Axis::default()
                    .bounds([self.min_rt, self.max_rt])
                    .labels(self.x_labels.clone()),
            )
            .y_axis(
                Axis::default()
                    .style(Style::default().gray())
                    .bounds(self.min_max_hyperscore.as_bounds())
                    .labels(self.min_max_hyperscore.as_scinot_labs()),
            );
        let bottom2_chart = Chart::new(datasets_normscore)
            .block(
                Block::bordered().title(
                    Title::default()
                        .content("Hyperscore".bold())
                        .alignment(Alignment::Center),
                ),
            )
            .x_axis(
                Axis::default()
                    .bounds([self.min_rt, self.max_rt])
                    .labels(self.x_labels.clone()),
            )
            .y_axis(
                Axis::default()
                    .style(Style::default().gray())
                    .bounds(self.min_max_norm_scores.as_bounds())
                    .labels(self.min_max_norm_scores.as_scinot_labs()),
            );

        let mid_chart = Chart::new(datasets_lazyscore)
            .block(
                Block::bordered().title(
                    Title::default()
                        .content("Lazyscore".bold())
                        .alignment(Alignment::Center),
                ),
            )
            .x_axis(
                Axis::default()
                    .bounds([self.min_rt, self.max_rt])
                    .labels(self.x_labels.clone()),
            )
            .y_axis(
                Axis::default()
                    .style(Style::default().gray())
                    .bounds(self.min_max_lazyscore.as_bounds())
                    .labels(self.min_max_lazyscore.as_scinot_labs()),
            );

        top_chart.render(top1, buf);
        mid_chart.render(top2, buf);
        bottom_chart.render(bottom1, buf);
        bottom2_chart.render(bottom2, buf);
    }
}
