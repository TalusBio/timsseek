use core::panic;
use crossterm::event::{self, poll, Event, KeyCode, KeyEvent, KeyEventKind};
use ratatui::layout::{Constraint, Layout};
use ratatui::{
    buffer::Buffer,
    layout::{Alignment, Rect},
    style::Stylize,
    text::{Line, Text},
    widgets::{block::Title, Block, Paragraph, Widget},
    DefaultTerminal, Frame,
};
use rustyms::error::CustomError;
use std::io;
use std::time::Duration;
use timsquery::models::aggregators::MultiCMGStatsFactory;
use timsquery::models::indices::transposed_quad_index::QuadSplittedTransposedIndex;
use timsquery::queriable_tims_data::queriable_tims_data::query_indexed;
use timsquery::traits::tolerance::{
    DefaultTolerance, MobilityTolerance, MzToleramce, QuadTolerance, RtTolerance,
};
use timsseek::fragment_mass::elution_group_converter::SequenceToElutionGroupConverter;
use timsseek::fragment_mass::fragment_mass_builder::SafePosition;

mod datatable;
mod plottable_chromatograms;

use datatable::{Data, TableInfo};
use timsseek;

#[derive(Debug, Default, PartialEq, Clone)]
pub enum AppState {
    #[default]
    Startup,
    LoadingIndex,
    LoadingSampleData,
    Ready,
    Rendered,
    Exiting,
}

#[derive(Debug, PartialEq, Clone)]
pub enum AppMessages {
    LoadIndex(String),

    // Q: Why doesnt this raise a compiler error if there is no `use` for Data?
    LoadData(Option<Data>),
    Warn(String),
    MoveUp,
    MoveDown,
    Quit,
}

// TODO: Change this redraw by doing 2 things ...
// 1. Make this a message
// 2. change the way I handle messages from being a single message that gets
//    dispatched to a queue that gets processed.
type ShouldRedraw = bool;

/// The main application state.
///
/// This struct is meant to hall all mutable state in the application and
/// its implementations are responsible for piping the events to the correct
/// way to manage it.
///
/// I am attempting to maintain the ELM architecture for the purpose of
/// simpliciry in the development. I am aware this might have some drawbacks
/// and if there is a compeling reason to do so I will gladly consider changing
/// it. (Having said so I would appreciate a discussion on this topic)
///
#[derive(Debug)]
pub struct App {
    index: Option<QuadSplittedTransposedIndex>,
    table_info: Option<TableInfo>,
    sample_data: Option<plottable_chromatograms::PlottableChromatograms>,
    state: AppState,
}

impl Default for App {
    fn default() -> Self {
        Self {
            index: None,
            sample_data: None,
            table_info: Some(TableInfo::default()),
            state: AppState::Startup,
        }
    }
}

struct LoadingBanner {
    pub data_path: String,
}

impl Widget for LoadingBanner {
    fn render(self, area: Rect, buf: &mut Buffer) {
        let title = Title::from(" Loading ".bold());
        let instructions = Text::from(Line::from(vec![
            " Loading index ".into(),
            self.data_path.bold(),
        ]));
        let block = Block::default().title(title.alignment(Alignment::Center));

        Paragraph::new(instructions)
            .centered()
            .block(block)
            .render(area, buf);
    }
}

impl App {
    fn load_index(&mut self) {
        let index = QuadSplittedTransposedIndex::from_path(
            "/Users/sebastianpaez/git/ionmesh/benchmark/240402_PRTC_01_S1-A1_1_11342.d",
        )
        .unwrap();
        self.index = Some(index);
    }

    fn load_data(&mut self, data: Option<Data>) -> Result<(), ()> {
        let index = match &self.index {
            Some(x) => x,
            None => return Err(()),
        };

        let sample_data = get_sample_data(&index, data).unwrap();

        let out_path = std::path::Path::new("./last_chromatograms.json");
        let mut out_file = std::fs::File::create(out_path).unwrap();
        serde_json::to_writer_pretty(&mut out_file, &sample_data).unwrap();

        self.sample_data = Some(sample_data);
        Ok(())
    }
}

impl App {
    /// runs the application's main loop until the user quits
    pub fn run(&mut self, terminal: &mut DefaultTerminal) -> io::Result<()> {
        while self.state != AppState::Exiting {
            terminal.draw(|frame| self.draw(frame))?;

            let (mut should_redraw, mut curr_message) = self.handle_events();

            loop {
                if should_redraw {
                    terminal.draw(|frame| self.draw(frame))?;
                }
                match curr_message {
                    Some(msg_contents) => {
                        (should_redraw, curr_message) = self.handle_state(msg_contents);
                    }
                    None => {
                        break;
                    }
                }
            }
        }
        Ok(())
    }

    pub fn handle_state(&mut self, msg: AppMessages) -> (ShouldRedraw, Option<AppMessages>) {
        match (&self.state, msg) {
            (_, AppMessages::Quit) => {
                self.state = AppState::Exiting;
                (true, None)
            }
            (AppState::Startup, _msg) => {
                self.state = AppState::LoadingIndex;
                (
                    true,
                    Some(AppMessages::LoadIndex(
                        "/Users/sebastianpaez/git/ionmesh/benchmark/240402_PRTC_01_S1-A1_1_11342.d"
                            .to_string(),
                    )),
                )
            }
            (AppState::LoadingIndex, _) => {
                self.load_index();
                self.state = AppState::LoadingSampleData;
                (true, Some(AppMessages::LoadData(None)))
            }
            (AppState::LoadingSampleData, AppMessages::LoadData(data)) => {
                match self.load_data(data) {
                    Err(_) => self.state = AppState::LoadingIndex,
                    Ok(_) => self.state = AppState::Ready,
                }

                (true, None)
            }
            (AppState::LoadingSampleData, _) => {
                panic!("Shouldnt be able to get here ... state: {:?}", self);
            }
            (AppState::Ready, msg) => match (&mut self.table_info, msg) {
                (Some(ref mut tab), AppMessages::MoveUp) => {
                    let out_data = tab.previous();
                    self.state = AppState::LoadingSampleData;
                    (true, Some(AppMessages::LoadData(Some(out_data))))
                }
                (Some(ref mut tab), AppMessages::MoveDown) => {
                    let out_data = tab.next();
                    self.state = AppState::LoadingSampleData;
                    (true, Some(AppMessages::LoadData(Some(out_data))))
                }
                (None, _) => (true, None),
                (Some(_), msg) => {
                    panic!(
                        "Shouldnt be able to get here ... msg: {:?} state: {:?}",
                        msg, self
                    );
                }
            },
            (AppState::Rendered, _) => {
                // pass
                (true, None)
            }

            (AppState::Exiting, _) => (true, None),
        }
    }

    fn draw(&mut self, frame: &mut Frame) {
        match &self.state {
            // AppState::Startup => self.draw_startup(frame),
            AppState::LoadingIndex => self.draw_loading_index(frame),
            // AppState::LoadingSampleData => self.draw_loading_sample_data(frame),
            AppState::Ready => self.draw_ready(frame),
            _ => {}
        }
    }

    fn draw_loading_index(&self, frame: &mut Frame) {
        frame.render_widget(
            LoadingBanner {
                data_path:
                    "/Users/sebastianpaez/git/ionmesh/benchmark/240402_PRTC_01_S1-A1_1_11342.d"
                        .to_string(),
            },
            frame.area(),
        );
    }

    fn draw_ready(&mut self, frame: &mut Frame) {
        let [left, right] =
            Layout::horizontal([Constraint::Percentage(30), Constraint::Percentage(70)])
                .areas(frame.area());

        if self.sample_data.is_some() {
            let plottable = self.sample_data.as_mut().unwrap();
            frame.render_widget(plottable.plottable_state.clone(), right);
        }
        if self.table_info.is_none() {
            frame.render_widget(Block::bordered().title("Left Block"), left);
        } else {
            let table = self.table_info.as_mut().unwrap();
            table.draw(frame, left);
        }
    }

    // TECHNICALLY this does not handle them ... more accurately it dispatches them.
    fn handle_events(&self) -> (ShouldRedraw, Option<AppMessages>) {
        match &self.state {
            AppState::Startup => self.handle_events_startup(),
            AppState::LoadingIndex => self.handle_events_loading_index(),
            AppState::LoadingSampleData => self.handle_events_loading_sample_data(),
            AppState::Ready => self.handle_events_ready(),
            _ => (false, None),
        }
    }

    // on startup we dont pay attenton to any events, we just want to
    // get to plotting our data asap.
    fn handle_events_startup(&self) -> (ShouldRedraw, Option<AppMessages>) {
        (true, Some(AppMessages::LoadIndex("Stuff".to_string())))
    }

    // on loading index we just want to wait for the index to load.
    // and move to data loading once that is done.
    fn handle_events_loading_index(&self) -> (ShouldRedraw, Option<AppMessages>) {
        if self.index.is_some() {
            (true, Some(AppMessages::LoadData(None)))
        } else {
            panic!(
                "Requesting data load when index is not loaded state: {:?}",
                self
            );
        }
    }

    fn handle_events_loading_sample_data(&self) -> (ShouldRedraw, Option<AppMessages>) {
        if self.sample_data.is_some() {
            (true, None)
        } else {
            panic!(
                "Requesting render when sample data is not loaded {:?}",
                self
            );
        }
    }

    fn handle_events_ready(&self) -> (ShouldRedraw, Option<AppMessages>) {
        let polled = poll(Duration::from_millis(100));

        if polled.is_err() {
            return (true, None);
        }

        match event::read() {
            Ok(Event::Key(key_event)) => {
                // it's important to check that the event is a key press event as
                // crossterm also emits key release and repeat events on Windows.
                if key_event.kind == KeyEventKind::Press {
                    self.handle_key_event(key_event)
                } else {
                    (true, None)
                }
            }
            Ok(_) => (true, None),
            Err(e) => {
                log::error!("Error reading event: {:?}", e);
                (
                    true,
                    Some(AppMessages::Warn("Error reading event".to_string())),
                )
            }
        }
    }

    fn handle_key_event(&self, key_event: KeyEvent) -> (ShouldRedraw, Option<AppMessages>) {
        match key_event.code {
            KeyCode::Char('q') | KeyCode::Esc => (false, Some(AppMessages::Quit)),
            KeyCode::Char('j') | KeyCode::Down => (true, Some(AppMessages::MoveDown)),
            KeyCode::Char('k') | KeyCode::Up => (true, Some(AppMessages::MoveUp)),
            _ => (false, None),
        }
    }
}

fn main() -> io::Result<()> {
    let mut terminal = ratatui::init();
    let app_result = App::default().run(&mut terminal);
    ratatui::restore();
    app_result
}

fn get_sample_data(
    index: &QuadSplittedTransposedIndex,
    data: Option<Data>,
) -> Result<plottable_chromatograms::PlottableChromatograms, CustomError> {
    // Really good peptide

    let (sample_peptide, sample_charge) = match data {
        Some(x) => (x.peptide.clone(), x.charge.parse::<u8>().unwrap()),
        None => ("VTIAQGGVLPNIQAVLLPK".to_string(), 2),
    };
    // False positive peptide
    // let sample_peptide = "SYFNANTNVHMFK";
    let tolerance = DefaultTolerance {
        ms: MzToleramce::Ppm((50.0, 50.0)),
        rt: RtTolerance::None,
        mobility: MobilityTolerance::Pct((20.0, 20.0)),
        quad: QuadTolerance::Absolute((0.1, 0.1, 1)),
    };
    let def_converter = SequenceToElutionGroupConverter::default();
    let factory = MultiCMGStatsFactory {
        converters: (index.mz_converter, index.im_converter),
        _phantom: std::marker::PhantomData::<SafePosition>,
    };
    let eg = def_converter.convert_sequence(&sample_peptide, 0)?;

    let out = query_indexed(index, &|x| factory.build(x), index, &tolerance, &eg[0]);
    let out = plottable_chromatograms::PlottableChromatograms::new(out, sample_peptide.to_string());

    Ok(out)
}
