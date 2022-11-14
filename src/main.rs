use anyhow::{Context, Result};
use clap::{Parser, Subcommand};
use std::path::PathBuf;

pub mod atac_shift_bam;

#[derive(Parser)]
#[command(author, version, about, long_about = None)]
struct Cli {
    #[command(subcommand)]
    command: Option<Commands>,
}

#[derive(Subcommand)]
enum Commands {
    Shift {
        /// Bam file for processing
        #[arg(short, long)]
        bam: Option<PathBuf>,

        /// Output file name
        #[arg(short, long)]
        output: Option<PathBuf>,
    },
}

fn main() -> Result<()> {
    let cli = Cli::parse();

    match &cli.command {
        Some(Commands::Shift { bam, output }) => match (bam, output) {
            (Some(bam_file), Some(output_file)) => {
                atac_shift_bam::atac_shift_bam(bam_file, output_file).with_context(|| {
                    format!(
                        "Shifting reads failed for file `{}`",
                        bam_file.to_string_lossy()
                    )
                })?;
            }
            (Some(bam_file), None) => {
                let output_file = &PathBuf::from("shifted.bam");
                atac_shift_bam::atac_shift_bam(bam_file, output_file).with_context(|| {
                    format!(
                        "Shifting reads failed for file `{}`",
                        bam_file.to_string_lossy()
                    )
                })?
            }

            _ => {
                println!("Options not provided, will not run")
            }
        },
        _ => {
            println!("Subcommand not provided, will not run")
        }
    }
    Ok(())
}
