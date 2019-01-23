extern crate clap;
extern crate flate2;
extern crate debruijn;

use clap::{Arg, App};

use debruijn::*;
use debruijn::dna_string::*;
use debruijn::kmer::*;

//use std::cmp::min;

use std::collections::HashMap;

use std::io::BufReader;
use std::io::BufRead;
use std::fs::File;
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
    let kmers: Vec<_> = matches.values_of("kmers").unwrap().collect();
    let output_prefix = matches.values_of("output_prefix").unwrap();
    let kmers: HashMap<u64,u32> = load_kmers(kmers);
    println!("Hello, world!");
}

fn load_kmers(kmers: Vec<&str>) -> HashMap<u64,u32> {
    let mut to_ret: HashMap<u64, u32> = HashMap::new();
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
            let kmer: Kmer21 = dna.get_kmer(0);
            println!("{}",kmer.to_string());
            to_ret.insert(kmer.to_u64(), index as u32);
        }           
    }
    to_ret
}


type Kmer21 = VarIntKmer<u64, K21>;
#[derive(Debug, Hash, Copy, Clone, Ord, PartialOrd, Eq, PartialEq)]
pub struct K21;

impl KmerSize for K21 {
    fn K() -> usize {
        21
    }
}
