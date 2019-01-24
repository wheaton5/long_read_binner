extern crate clap;
extern crate flate2;
extern crate debruijn;

use clap::{Arg, App};

use debruijn::*;
use debruijn::dna_string::*;
use debruijn::kmer::*;

//use std::cmp::min;

use std::collections::HashMap;

mod extra_kmers;
use extra_kmers::*;

use std::io::BufReader;
use std::io::BufRead;
use std::io::BufWriter;
use std::fs::File;
use flate2::Compression;
use flate2::write::ZlibEncoder;
use flate2::read::GzDecoder;

//use std::path::Path;


fn main() {
    let matches = App::new("long_read_binner")
        .author("Haynes Heaton <whheaton@gmail.com>")
        .about("Uses distinguishing kmers to bin long reads")
        .arg(Arg::with_name("kmers")
            .long("kmers")
            .short("k")
            .takes_value(true)
            .required(true)
            .multiple(true)
            .help("set of kmer files with 1 kmer per line. Will output a long read file for each kmer file plus an unplaced long read file."))
        .arg(Arg::with_name("output_prefix")
            .long("output_prefix")
            .short("o")
            .takes_value(true)
            .required(true)
            .help("prefix for output long read files"))
        .arg(Arg::with_name("fastqs")
            .long("fastqs")
            .short("f")
            .required(true)
            .takes_value(true)
            .multiple(true)
            .help("long read fastq.gz files"))
        .get_matches();
    let fastqs: Vec<_> = matches.values_of("fastqs").unwrap().collect();
    let kmer_files: Vec<_> = matches.values_of("kmers").unwrap().collect();
    let output_prefix = matches.value_of("output_prefix").unwrap_or("longreads_bin_");
    let bins: u32 = kmer_files.len() as u32;
    let (kmers, k_size): (HashMap<u64,u32>, usize) = load_kmers(kmer_files);
    bin_long_reads(kmers, bins, fastqs, output_prefix.to_string(), k_size);
    println!("Hello, world!");
}

fn bin_long_reads(kmers: HashMap<u64, u32>, bins: u32, fastqs: Vec<&str>, output_prefix: String, k_size: usize) {
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
        let file = match File::open(fastq) {
            Ok(file) => file,
            Err(error) => panic!("There was a problem opening the file: {:?}", error),
        };
        let gz = GzDecoder::new(file);
        for (line_number, line) in BufReader::new(gz).lines().enumerate() {
            if line_number % 4 == 1 {
                let dna = DnaString::from_dna_string(&line.unwrap());
                for kmer_start in 0..(dna.len() - k_size + 1) {
                    println!("do stuff there");
                }
            }
        }
    }
}

fn load_kmers(kmers: Vec<&str>) -> (HashMap<u64,u32>, usize) {
    let mut to_ret: HashMap<u64, u32> = HashMap::new();
    let mut kmer_size: usize = 0;
    for (index, kmer_file) in kmers.iter().enumerate() {
        let f = match File::open(kmer_file) {
            Ok(f) => f,
            Err(err) => panic!("There was a problem opening the file: {:?}",err),
        };
        let mut file = BufReader::new(&f);
        for line in file.lines() {
            let tmp = line.unwrap();
            let tokens: Vec<&str> = tmp.split_whitespace().collect();
            
            let dna = DnaString::from_dna_string(&tokens[0]);
            if kmer_size == 0 {
                kmer_size = dna.len();
            } else if dna.len() != kmer_size {
                panic!("kmer sizes in files are not consistent, previous kmers of length {}, vs {}",kmer_size,dna.to_string());
            }
            let to_hash = get_rc_invariant(&dna, kmer_size);
            if to_ret.contains_key(&to_hash) {
                panic!("kmer {} seen before, kmer sets must be unique (and reverse compliment unique) and mutually exclusive", dna.to_string());
            }
            to_ret.insert(to_hash, index as u32);
        }           
    }
    (to_ret, kmer_size)
}

