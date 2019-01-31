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
use std::io::BufWriter;
use std::fs::File;
use flate2::Compression;
use flate2::write::ZlibEncoder;

use std::cmp::min;

static mut KMER_SIZE: usize = 0; // meant to change. will fail if it doesn't

fn main() {
    let (fastqs, kmer_files, output_prefix) = load_params();
    let bins: usize = kmer_files.len();
    let kmers: HashMap<u64, usize> = load_kmers(kmer_files);
    bin_long_reads(kmers, bins, fastqs, output_prefix.to_string());
}

fn bin_long_reads(kmers: HashMap<u64, usize>, bins: usize, fastqs: Vec<String>, output_prefix: String) {
    let mut writers : Vec<ZlibEncoder<BufWriter<File>>> = Vec::new();
    for writer_index in 0..bins {
        let mut filename = output_prefix.clone();
        filename.push_str(&writer_index.to_string());
        filename.push_str(".fastq.gz");
        let mut file = File::create(filename).expect("Unable to create file");
        let mut buf = BufWriter::new(file);
        let mut writer = ZlibEncoder::new(buf, Compression::default());
        writers.push(writer);
    }
    for fastq in &fastqs {
        let mut bin_hits: Vec<u32> = Vec::new();
        for _bin in 0..bins { bin_hits.push(0); }
        let mut reader = dna_io::DnaReader::from_path(fastq);
        'lineloop: loop {
            let record = match reader.next() {
                Some(x) => x,
                None => break 'lineloop,
            };
            for k in KmerX::kmers_from_ascii(&record.seq.as_bytes()) {
                match kmers.get(&get_rc_invariant_kmer(k)) {
                    Some(x) => {
                        bin_hits[*x] += 1;                
                    },
                    None => (),
                }
            }
            println!("{:?}", bin_hits);
            for bin in 0..bins { bin_hits[bin] = 0; }
        }
    }
}

fn load_kmers(kmers: Vec<String>) -> HashMap<u64, usize> {
    let mut to_ret: HashMap<u64, usize> = HashMap::new();
    for (index, kmer_file) in kmers.iter().enumerate() {
        println!("kmer file {}",kmer_file);
        let f = match File::open(kmer_file) {
            Ok(f) => f,
            Err(err) => panic!("There was a problem opening the file: {:?}",err),
        };
        let mut file = BufReader::new(&f);
        for line in file.lines() {
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

fn load_params() -> (Vec<String>, Vec<String>, String) {
    let yaml_params = load_yaml!("params.yml");
    let params = App::from_yaml(yaml_params).get_matches();
    let mut fastqs: Vec<String> = Vec::new();
    for fastq in params.values_of("input").unwrap() {
        fastqs.push(fastq.to_string());
    }
    let mut kmer_files: Vec<String> = Vec::new();
    for kmer_file in params.values_of("kmers").unwrap() {
        kmer_files.push(kmer_file.to_string());
    }
    let output_prefix = params.value_of("output_prefix").unwrap_or("longreads_bin_").to_string();
    (fastqs, kmer_files, output_prefix)
}
