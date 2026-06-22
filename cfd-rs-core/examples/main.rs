use std::path::Path;

// #[path = "poisson/test_1"]
// mod test
mod poisson;

fn examples_default_logger(log_path: &Path, level: log::LevelFilter) -> Result<(), fern::InitError> {
    fern::Dispatch::new()
        .format(|out, message, record| {
            out.finish(format_args!(
                "[{} {}] {}",
                record.level(),
                record.target(),
                message
            ))
        })
        .level(level)
        // .chain(std::io::stdout())
        .chain(std::fs::OpenOptions::new()
            .write(true)
            .create(true)
            .truncate(true)
            .open(log_path)?)
        .apply()?;
    Ok(())
}

fn main() {

}

