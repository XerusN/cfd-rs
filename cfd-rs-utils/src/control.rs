use std::{
    fs::{self, DirEntry, FileType},
    io::Result,
    path::Path,
};

use log::warn;

#[derive(Clone, Debug, Default, PartialEq)]
pub enum Frequency {
    TimeStep(f64),
    Iteration(usize),
    #[default]
    None,
}

#[derive(Clone, Debug, PartialEq)]
pub struct OutputControl {
    initial_output: bool,
    final_output: bool,
    frequency: Frequency,
}

impl Default for OutputControl {
    fn default() -> Self {
        OutputControl {
            frequency: Frequency::default(),
            initial_output: true,
            final_output: false,
        }
    }
}

impl OutputControl {
    pub fn new(initial_output: bool, final_output: bool, frequency: Frequency) -> Self {
        OutputControl {
            frequency,
            initial_output,
            final_output,
        }
    }

    pub fn initial_output(&self) -> bool {
        self.initial_output
    }

    pub fn final_output(&self) -> bool {
        self.final_output
    }

    pub fn frequency(&self) -> &Frequency {
        &self.frequency
    }
}

pub fn clear_dump_folder(path: &Path) -> Result<()> {
    for entry in fs::read_dir(path)? {
        let entry = entry?;
        let file_type = entry.file_type()?;
        if file_type.is_symlink() {
            ();
        } else if file_type.is_dir() {
            clear_dump_folder(&entry.path())?;
        } else if file_type.is_file() {
            if entry.file_name() == ".gitignore" {
                ();
            } else {
                fs::remove_file(entry.path())?;
            }
        }
    }
    Ok(())
}
