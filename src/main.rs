use anyhow::{Context, Result};
use clap::{Parser, Subcommand};
use std::path::{PathBuf};

pub mod atac_shift_bam;
pub mod subtract_regions;

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

    Subtract {
        /// Bed file for processing
        #[arg(short='r', long="regions")]
        regions: Option<PathBuf>,

        /// Bam file for processing
        #[arg(short='b', long="bam")]
        bam: Option<PathBuf>,

        /// Output file name
        #[arg(short, long)]
        output: Option<PathBuf>,

        /// Number of threads to use
        #[arg(short, long)]
        threads: Option<usize>,
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

        Some(Commands::Subtract {
            regions: bed,
            bam,
            output,
            threads,
        }) => {
            println!("Running subtract subcommand. Will subtract regions from BAM file.");

            match (bed, bam) {
                (Some(bed_file), Some(bam_file)) => {
                    let output = match output {
                        Some(output) => output.to_owned(),
                        None => PathBuf::from("subtracted.bam"),
                    };
                    let threads = match threads {
                        Some(threads) => *threads,
                        None => 1,
                    };

                    println!("BED file: {}", bed_file.to_string_lossy());
                    println!("BAM file: {}", bam_file.to_string_lossy());
                    println!("Output file: {}", output.to_string_lossy());
                    println!("Threads: {}", threads);
                    
                    subtract_regions::remove_regions_from_bam(
                        bed_file.to_path_buf(),
                        bam_file.to_path_buf(),
                        output,
                        threads,
                    )?;
                }
                _ => {
                    println!("Options not provided, will not run");
                    return Ok(());
                }
            }
        }

        _ => {
            println!("Subcommand not provided, will not run")
        }
    }
    Ok(())
}
