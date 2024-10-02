use ratatui::layout::{Constraint, Layout, Margin};
use ratatui::{
    layout::Rect,
    style::{Color, Modifier, Style, Stylize},
    text::{Line, Text},
    widgets::{
        Block, BorderType, Cell, HighlightSpacing, Paragraph, Row, Scrollbar, ScrollbarOrientation,
        ScrollbarState, Table, TableState,
    },
    Frame,
};
use serde::{Deserialize, Serialize};

#[derive(Debug, Default, PartialEq, Clone, Serialize, Deserialize)]
pub struct Data {
    pub peptide: String,
    pub charge: String,
    pub note: String,
}

impl Data {
    const fn ref_array(&self) -> [&String; 3] {
        [&self.peptide, &self.charge, &self.note]
    }
}

#[derive(Debug, PartialEq, Clone)]
pub struct TableColors {
    buffer_bg: Color,
    header_bg: Color,
    header_fg: Color,
    row_fg: Color,
    selected_style_fg: Color,
    normal_row_color: Color,
    alt_row_color: Color,
    footer_border_color: Color,
}

impl Default for TableColors {
    fn default() -> Self {
        Self {
            buffer_bg: Color::Black,
            header_bg: Color::Black,
            header_fg: Color::White,
            row_fg: Color::White,
            selected_style_fg: Color::Gray,
            normal_row_color: Color::DarkGray,
            alt_row_color: Color::DarkGray,
            footer_border_color: Color::Cyan,
        }
    }
}

#[derive(Debug, PartialEq, Clone)]
pub struct TableInfo {
    pub items: Vec<Data>,
    pub selected: Option<usize>,
    pub scroll_state: ScrollbarState,
    pub state: TableState,
    pub colors: TableColors,
}

impl Default for TableInfo {
    fn default() -> Self {
        let items = vec![
            Data {
                peptide: "VTIAQGGVLPNIQAVLLPK".to_string(),
                charge: "2".to_string(),
                note: "Note".to_string(),
            },
            Data {
                peptide: "SYFNANTNVHMFK".to_string(),
                charge: "2".to_string(),
                note: "Note".to_string(),
            },
            Data {
                peptide: "MYPEPTIDEK".to_string(),
                charge: "2".to_string(),
                note: "Note".to_string(),
            },
        ];
        Self {
            items,
            selected: None,
            scroll_state: ScrollbarState::default(),
            state: TableState::default(),
            colors: TableColors::default(),
        }
    }
}

impl TableInfo {
    pub fn read_csv(path: &str) -> Result<Self, Box<dyn std::error::Error>> {
        let mut rdr = csv::Reader::from_path(path)?;
        let mut items = Vec::new();
        for result in rdr.deserialize() {
            let record: Data = result?;
            items.push(record);
        }
        Ok(Self {
            items,
            selected: None,
            scroll_state: ScrollbarState::default(),
            state: TableState::default(),
            colors: TableColors::default(),
        })
    }
}

const ITEM_HEIGHT: usize = 2;
const INFO_TEXT: &str = "(Esc) quit | (↑) move up | (↓) move down ";

impl TableInfo {
    pub fn next(&mut self) -> Data {
        let i = match self.state.selected() {
            Some(i) => {
                if i >= self.items.len() - 1 {
                    0
                } else {
                    i + 1
                }
            }
            None => 0,
        };
        self.state.select(Some(i));
        self.scroll_state = self.scroll_state.position(i * ITEM_HEIGHT);
        self.items[i].clone()
    }

    pub fn previous(&mut self) -> Data {
        let i = match self.state.selected() {
            Some(i) => {
                if i == 0 {
                    self.items.len() - 1
                } else {
                    i - 1
                }
            }
            None => 0,
        };
        self.state.select(Some(i));
        self.scroll_state = self.scroll_state.position(i * ITEM_HEIGHT);
        self.items[i].clone()
    }

    pub fn draw(&mut self, frame: &mut Frame, area: Rect) {
        let vert = &Layout::vertical([Constraint::Min(5), Constraint::Length(3)]);
        let rects = vert.split(area);

        self.render_table(frame, rects[0]);
        self.render_scrollbar(frame, rects[0]);
        self.render_footer(frame, rects[1]);
    }

    fn render_table(&mut self, frame: &mut Frame, area: Rect) {
        let header_style = Style::default()
            .fg(self.colors.header_fg)
            .bg(self.colors.header_bg);
        let selected_style = Style::default()
            .add_modifier(Modifier::REVERSED)
            .fg(self.colors.selected_style_fg);

        let header = ["Peptide", "Charge", "Note"]
            .into_iter()
            .map(Cell::from)
            .collect::<Row>()
            .style(header_style)
            .height(1);
        let rows = self.items.iter().enumerate().map(|(i, data)| {
            let color = match i % 2 {
                0 => self.colors.normal_row_color,
                _ => self.colors.alt_row_color,
            };
            let item = data.ref_array();
            item.into_iter()
                .map(|content| Cell::from(Text::from(format!("\n{content}\n"))))
                .collect::<Row>()
                .style(Style::new().fg(self.colors.row_fg).bg(color))
                .height(4)
        });
        let bar = " █ ";
        let t = Table::new(
            rows,
            [
                // + 1 is for padding.
                // TODO: Set this constraints at build time
                Constraint::Percentage(70),
                Constraint::Percentage(10),
                Constraint::Percentage(20),
            ],
        )
        .header(header)
        .highlight_style(selected_style)
        .highlight_symbol(Text::from(vec![
            "".into(),
            bar.into(),
            bar.into(),
            "".into(),
        ]))
        .bg(self.colors.buffer_bg)
        .highlight_spacing(HighlightSpacing::Always);
        frame.render_stateful_widget(t, area, &mut self.state);
    }

    fn render_scrollbar(&mut self, frame: &mut Frame, area: Rect) {
        frame.render_stateful_widget(
            Scrollbar::default()
                .orientation(ScrollbarOrientation::VerticalRight)
                .begin_symbol(None)
                .end_symbol(None),
            area.inner(Margin {
                vertical: 1,
                horizontal: 1,
            }),
            &mut self.scroll_state,
        );
    }

    fn render_footer(&self, frame: &mut Frame, area: Rect) {
        let info_footer = Paragraph::new(Line::from(INFO_TEXT))
            .style(
                Style::new()
                    .fg(self.colors.row_fg)
                    .bg(self.colors.buffer_bg),
            )
            .centered()
            .block(
                Block::bordered()
                    .border_type(BorderType::Double)
                    .border_style(Style::new().fg(self.colors.footer_border_color)),
            );
        frame.render_widget(info_footer, area);
    }
}
