#[macro_use]
extern crate clap;
extern crate flate2;
extern crate debruijn;
extern crate dna_io;


use clap::{App};

use debruijn::*;
use debruijn::dna_string::*;
use debruijn::kmer::*;

use std::collections::HashMap;


use std::io::BufReader;
use std::io::BufRead;
use std::fs::File;

use std::cmp::min;

static mut KMER_SIZE: usize = 0; // meant to change. will fail if it doesn't

fn main() {
    let params = load_params();
    let kmers: HashMap<u64, usize> = load_kmers(&params.kmer_files);
    bin_long_reads(kmers, params);
}

fn bin_long_reads(kmers: HashMap<u64, usize>,  params: Params) {
    let mut writers : Vec<dna_io::DnaWriter> = Vec::new();
    let reader = dna_io::DnaReader::from_path(&params.input_files[0]);
    for writer_index in 0..params.bins {
        let mut filename = params.output_prefix.clone();
        filename.push_str(&writer_index.to_string());
        filename.push_str(&reader.extension());
        writers.push(dna_io::DnaWriter::from_reader(&filename, &reader));
    }
    let mut filename = params.output_prefix.clone();
    filename.push_str(&"_unplaced_");
    filename.push_str(&reader.extension());
    writers.push(dna_io::DnaWriter::from_reader(&filename, &reader));
    let unplaced_index = writers.len()-1;

    
    for input in &params.input_files {
        let mut bin_hits: Vec<i32> = Vec::new();
        for _bin in 0..params.bins { bin_hits.push(0); }
        let mut reader = dna_io::DnaReader::from_path(input);
        
        for record in reader {
            for k in KmerX::kmers_from_ascii(&record.seq.as_bytes()) {
                let to_hash = get_rc_invariant_kmer(k);
                if to_hash % params.modimizer != params.mod_index { continue }
                match kmers.get(&get_rc_invariant_kmer(k)) {
                    Some(x) => {
                        bin_hits[*x] += 1;   
                    },
                    None => (),
                };
            }
            let mut best_bin: usize = 0;
            let mut best_bin_count: i32 = -5000;
            let mut second_best_bin_count: i32 = -5000;
            for bin in 0..params.bins {
                if bin_hits[bin] > best_bin_count {
                    second_best_bin_count = best_bin_count;
                    best_bin = bin;
                    best_bin_count = bin_hits[bin];
                } else if bin_hits[bin] > second_best_bin_count {
                    second_best_bin_count = bin_hits[bin];
                }
            }
            if best_bin_count - second_best_bin_count > params.difference_threshold {
                match writers[best_bin].write(&record) {
                    Ok(_) => (),
                    Err(error) => panic!(error),
                }
            } else {
                match writers[unplaced_index].write(&record) {
                    Ok(_) => (),
                    Err(error) => panic!(error),
                }
            }
            for bin in 0..params.bins { bin_hits[bin] = 0; }

        }
        
    }
}

fn load_kmers(kmers: &Vec<String>) -> HashMap<u64, usize> {
    let mut to_ret: HashMap<u64, usize> = HashMap::new();
    for (index, kmer_file) in kmers.iter().enumerate() {
        let f = match File::open(kmer_file) {
            Ok(f) => f,
            Err(err) => panic!("There was a problem opening the file: {:?}",err),
        };
        let mut file = BufReader::new(&f);
        for line in file.lines() {
            eprintln!("read line");
            let tmp = line.unwrap();
            let tokens: Vec<&str> = tmp.split_whitespace().collect();
            
            let dna = DnaString::from_dna_string(&tokens[0]);
            unsafe {
                if KMER_SIZE == 0 && dna.len() != 0 {
                    KMER_SIZE = dna.len();                
                } else if dna.len() != KMER_SIZE {
                    panic!("kmer sizes in files are not consistent, previous kmers of length {}, vs {}",KMER_SIZE, dna.to_string());
                }
            }
            let to_hash = get_rc_invariant(&dna, 0);
            if to_ret.contains_key(&to_hash) {
                panic!("kmer {} seen before, kmer sets must be unique (and reverse compliment unique) and mutually exclusive", dna.to_string());
            }
            to_ret.insert(to_hash, index);
        }           
    }
    to_ret
}

fn get_rc_invariant(dna: &DnaString, pos: usize) -> u64 {
    let kmer: KmerX = dna.get_kmer(pos);
    min(kmer.to_u64(), kmer.rc().to_u64())
}

fn get_rc_invariant_kmer(kmer: KmerX) -> u64 {
    min(kmer.to_u64(), kmer.rc().to_u64())
}

type KmerX = VarIntKmer<u64, KX>;
#[derive(Debug, Hash, Copy, Clone, Ord, PartialOrd, Eq, PartialEq)]
pub struct KX;

impl KmerSize for KX {
    fn K() -> usize {
        unsafe {
            KMER_SIZE
        }
    }
}

pub fn logaddexp(a: f64, b: f64, base: f64) -> f64 {
    let c: f64;
    let d: f64;

    // powf(a-b) doesn't do the right thing if a > b.
    if a > b {
        c = b;
        d = a;
    } else {
        c = a;
        d = b;
    }

    (base.powf(c - d) + 1.0).log(base) + d
}


struct Params {
    input_files: Vec<String>,
    kmer_files: Vec<String>,
    output_prefix: String,
    modimizer: u64,
    mod_index: u64,
    difference_threshold: i32,
    bins: usize,
}

fn load_params() -> Params {
    let yaml_params = load_yaml!("params.yml");
    let params = App::from_yaml(yaml_params).get_matches();
    let mut inputs: Vec<String> = Vec::new();
    for input in params.values_of("input").unwrap() {
        inputs.push(input.to_string());
    }
    assert!(inputs.len() > 0, "Input file list empty");
    let mut kmer_files: Vec<String> = Vec::new();
    for kmer_file in params.values_of("kmers").unwrap() {
        kmer_files.push(kmer_file.to_string());
    }
    let output_prefix = params.value_of("output_prefix").unwrap_or("longreads_bin_").to_string();
    let modimizer = params.value_of("modimizer").unwrap_or("1");
    let modimizer: u64 = modimizer.to_string().parse::<u64>().unwrap();
    let mod_index = params.value_of("mod_index").unwrap_or("0");
    let mod_index: u64 = mod_index.to_string().parse::<u64>().unwrap();
    let difference_threshold = params.value_of("p_threshold").unwrap_or("1");
    let difference_threshold: i32 = difference_threshold.to_string().parse::<i32>().unwrap();
    let bins: usize = kmer_files.len();
    Params {input_files: inputs, 
            kmer_files: kmer_files, 
            bins: bins,
            output_prefix: output_prefix, 
            modimizer: modimizer, 
            mod_index: mod_index, 
            difference_threshold: difference_threshold
        }
}
