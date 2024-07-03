use ahash::HashMap;
use anyhow::{Context, Result};
use bstr::ByteSlice;
use noodles::bam::io::Writer;
use noodles::bed::record;
use noodles::{bam, bgzf, sam};
use std::fmt::format;
use std::num::NonZeroUsize;
use std::path::PathBuf;
use std::prelude::v1::*;
use serde::{Serialize, Deserialize};
use indicatif::{ProgressBar, ProgressIterator};
use sam::header::record::value::{map::ReferenceSequence, Map};


#[derive(Debug, Serialize, Deserialize)]
pub struct SplitStats {
    filename: String,
    n_unmapped_reads: u64,
    n_qcfail_reads: u64,
    n_duplicate_reads: u64,
    n_secondary_reads: u64,
    n_low_maq: u64,
    n_both_genomes: u64,
    n_exogenous: u64,
    n_endogenous: u64,
}

impl SplitStats{
    fn new(filename: String) -> Self {
        Self {
            filename,
            n_unmapped_reads: 0,
            n_qcfail_reads: 0,
            n_duplicate_reads: 0,
            n_secondary_reads: 0,
            n_low_maq: 0,
            n_both_genomes: 0,
            n_exogenous: 0,
            n_endogenous: 0,
        }
    }

    fn add_unmapped(&mut self) {
        self.n_unmapped_reads += 1;
    }

    fn add_qcfail(&mut self) {
        self.n_qcfail_reads += 1;
    }

    fn add_duplicate(&mut self) {
        self.n_duplicate_reads += 1;
    }

    fn add_secondary(&mut self) {
        self.n_secondary_reads += 1;
    }

    fn add_low_maq(&mut self) {
        self.n_low_maq += 1;
    }

    fn add_both_genomes(&mut self) {
        self.n_both_genomes += 1;
    }

    fn add_exogenous(&mut self) {
        self.n_exogenous += 1;
    }

    fn add_endogenous(&mut self) {
        self.n_endogenous += 1;
    }

    pub fn print(&self) {
        println!("Filename: {}", self.filename);
        println!("Unmapped reads: {}", self.n_unmapped_reads);
        println!("QC fail reads: {}", self.n_qcfail_reads);
        println!("Duplicate reads: {}", self.n_duplicate_reads);
        println!("Secondary reads: {}", self.n_secondary_reads);
        println!("Low mapping quality reads: {}", self.n_low_maq);
        println!("Both genomes reads: {}", self.n_both_genomes);
        println!("Exogenous reads: {}", self.n_exogenous);
        println!("Endogenous reads: {}", self.n_endogenous);
    }

}


pub struct SplitBam {
    bam_input: bam::io::Reader<noodles::bgzf::Reader<std::fs::File>>,
    bam_endogenous: bam::io::Writer<noodles::bgzf::Writer<std::fs::File>>,
    bam_exogenous: bam::io::Writer<noodles::bgzf::Writer<std::fs::File>>,
    bam_both_genomes: bam::io::Writer<noodles::bgzf::Writer<std::fs::File>>,
    bam_unmapped: bam::io::Writer<noodles::bgzf::Writer<std::fs::File>>,
}

struct BamHeaders {
    header_input: sam::Header,
    header_endogenous: sam::Header,
    header_exogenous: sam::Header,
    header_both_genomes: sam::Header,
    header_unmapped: sam::Header,
}

impl SplitBam {
    pub fn new(bam_input: PathBuf, output_prefix: PathBuf) -> Result<Self> {
        let bam_input = bam::io::reader::Builder::default().build_from_path(bam_input)?;
        let bam_endogenous = bam::io::writer::Builder::default()
            .build_from_path(output_prefix.with_extension("endogenous.bam"))?;
        let bam_exogenous = bam::io::writer::Builder::default()
            .build_from_path(output_prefix.with_extension("exogenous.bam"))?;
        let bam_both_genomes = bam::io::writer::Builder::default()
            .build_from_path(output_prefix.with_extension("both_genomes.bam"))?;
        let bam_unmapped = bam::io::writer::Builder::default()
            .build_from_path(output_prefix.with_extension("unmapped.bam"))?;

        Ok(Self {
            bam_input,
            bam_endogenous,
            bam_exogenous,
            bam_both_genomes,
            bam_unmapped,
        })
    }

    fn make_headers(&mut self, exogenous_prefix: &[u8]) -> Result<BamHeaders> {
        let header_input = self.bam_input.read_header()?;

        let reference_seqs = header_input.reference_sequences().clone();

        // Split reference sequences into endogenous and exogenous based on prefixes present.
        // Endogenous sequences have no prefix, exogenous sequences have a prefix.
        let mut reference_seqs_endogenous = sam::header::ReferenceSequences::new();
        let mut reference_seqs_exogenous = sam::header::ReferenceSequences::new();

        for (name, len) in reference_seqs.iter() {
            if name.starts_with(&exogenous_prefix) {
                reference_seqs_exogenous.insert(name.clone(), len.clone());
            } else {
                reference_seqs_endogenous.insert(name.clone(), len.clone());
            }
        }

        let header_endogenous = sam::Header::builder()
            .set_header(header_input.header().expect("No header present").clone())
            .set_reference_sequences(reference_seqs_endogenous)
            .build();

        let header_exogenous = sam::Header::builder()
            .set_header(header_input.header().expect("No header present").clone())
            .set_reference_sequences(reference_seqs_exogenous)
            .build();

        let header_both_genomes = sam::Header::builder()
            .set_header(header_input.header().expect("No header present").clone())
            .set_reference_sequences(reference_seqs.clone())
            .build();

        // let header_unmapped = sam::Header::builder()
        //     .set_header(header_input.header().expect("No header present").clone())
        //     .add_reference_sequence("unmapped",  Map::<ReferenceSequence>::new(NonZeroUsize::try_from(1e6 as usize)?)) // Provide a dummy reference sequence argument
        //     .build();
        
        let header_unmapped = sam::Header::builder()
            .set_header(header_input.header().expect("No header present").clone())
            .set_reference_sequences(reference_seqs)
            .build();

        Ok(BamHeaders {
            header_input,
            header_endogenous,
            header_exogenous,
            header_both_genomes,
            header_unmapped,
        })
    }

    fn write_headers(&mut self, headers: &BamHeaders) -> Result<()> {
        self.bam_endogenous.write_header(&headers.header_endogenous)?;
        self.bam_exogenous.write_header(&headers.header_exogenous)?;
        self.bam_both_genomes.write_header(&headers.header_both_genomes)?;
        self.bam_unmapped.write_header(&headers.header_unmapped)?;
        Ok(())
    }

    pub fn split(&mut self, exogenous_prefix: &[u8]) -> Result<SplitStats> {
        
        let headers = self.make_headers(exogenous_prefix)?;
        self.write_headers(&headers)?;
        let mut stats = SplitStats::new("SplitBam".to_string());


        for (ii, record) in self.bam_input.records().enumerate() {
            let record = record.expect(format!("Error reading record {}", ii).as_str());
            if ii % 1_000_000 == 0 {
                println!("Processed {} reads", ii);
            }
    
            if record.flags().is_unmapped() {
                self.bam_unmapped
                    .write_record(&headers.header_unmapped, &record)
                    .expect("Error writing record");
                stats.add_unmapped();
                continue;
            } else if record.flags().is_qc_fail() {
                self.bam_unmapped
                    .write_record(&headers.header_unmapped, &record)
                    .expect("Error writing record");
                stats.add_qcfail();
                continue;
            } else if record.flags().is_duplicate() {
                self.bam_unmapped
                    .write_record(&headers.header_unmapped, &record)
                    .expect("Error writing record");
                stats.add_duplicate();
                continue;
            } else if record.flags().is_secondary() {
                self.bam_unmapped
                    .write_record(&headers.header_unmapped, &record)
                    .expect("Error writing record");
                stats.add_secondary();
                continue;
            } else if record.mapping_quality().expect("No mapping quality").get() < 30 {
                self.bam_unmapped
                    .write_record(&headers.header_unmapped, &record)
                    .expect("Error writing record");
                stats.add_low_maq();
                continue;
            } else if !record.flags().is_mate_unmapped() {
                let r1_seq_id = record
                    .reference_sequence_id()
                    .expect("No reference sequence ID")
                    .expect("Failed to get reference sequence ID");
                let r1_seq_name = headers
                    .header_input
                    .reference_sequences()
                    .get_index(r1_seq_id)
                    .expect("Failed to get reference sequence name")
                    .0;
                let r2_seq_id = record
                    .mate_reference_sequence_id()
                    .expect("No mate reference sequence ID")
                    .expect("Failed to get mate reference sequence ID");
                let r2_seq_name = headers
                    .header_input
                    .reference_sequences()
                    .get_index(r2_seq_id)
                    .expect("Failed to get mate reference sequence name")
                    .0;

                if r1_seq_name.starts_with(exogenous_prefix)
                    && r2_seq_name.starts_with(exogenous_prefix)
                {
                    let res = self.bam_exogenous
                        .write_record(&headers.header_exogenous, &record);


                    match res {
                        Ok(_) => {},
                        Err(e) => {
                            println!("Error writing record: {:?}", e);
                        }
                    }

                    stats.add_exogenous();
                    continue;
                } else if r1_seq_name.starts_with(exogenous_prefix)
                    || r2_seq_name.starts_with(exogenous_prefix)
                {
                    self.bam_both_genomes
                        .write_record(&headers.header_both_genomes, &record)
                        .expect("Error writing record");

                    stats.add_both_genomes();
                    continue;
                } else {
                    self.bam_endogenous
                        .write_record(&headers.header_endogenous, &record)
                        .expect("Error writing record");
                    stats.add_endogenous();
                    continue;
                };
            } else if record.flags().is_mate_unmapped() {
                let r1_seq_id = record
                    .reference_sequence_id()
                    .expect("No reference sequence ID")
                    .expect("Failed to get reference sequence ID");
                let r1_seq_name = headers
                    .header_input
                    .reference_sequences()
                    .get_index(r1_seq_id)
                    .expect("Failed to get reference sequence name")
                    .0;

                if r1_seq_name.starts_with(exogenous_prefix) {
                    self.bam_exogenous
                        .write_record(&headers.header_exogenous, &record)
                        .expect("Error writing record");
                    stats.add_exogenous();
                    continue;
                } else {
                    self.bam_endogenous
                        .write_record(&headers.header_endogenous, &record)
                        .expect("Error writing record");
                    stats.add_endogenous();
                    continue;
                }
            }
        }
        Ok(stats)
    }

}
