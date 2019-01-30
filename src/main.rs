extern crate clap;
extern crate flate2;
extern crate debruijn;
extern crate dna_io;
//extern crate num;

use clap::{Arg, App};

use debruijn::*;
use debruijn::dna_string::*;
use debruijn::kmer::*;

use std::collections::HashMap;

//mod extra_kmers;
//use extra_kmers::*;

use std::io::BufReader;
use std::io::BufRead;
use std::io::BufWriter;
use std::fs::File;
use flate2::Compression;
use flate2::write::ZlibEncoder;
use flate2::read::GzDecoder;

use std::cmp::min;

static mut KMER_SIZE: usize = 0; // meant to change. will fail if it doesn't

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
            .help("set of kmer files with 1 kmer per line. Will output a long read file for each kmer file plus an unplaced long read file. Kmers from each file must be the same length and mutually exclusive. Kmers of up to size 32 supported."))
        .arg(Arg::with_name("output_prefix")
            .long("output_prefix")
            .short("o")
            .takes_value(true)
            .required(true)
            .help("prefix for output long read files"))
        .arg(Arg::with_name("input")
            .long("input")
            .short("f")
            .required(true)
            .takes_value(true)
            .multiple(true)
            .help("long read files (fastq/fastq.gz, fasta/fasta.gz, sam, bam supported)"))
        .get_matches();
    let fastqs: Vec<_> = matches.values_of("fastqs").unwrap().collect();
    let kmer_files: Vec<_> = matches.values_of("kmers").unwrap().collect();
    let output_prefix = matches.value_of("output_prefix").unwrap_or("longreads_bin_");
    let bins: usize = kmer_files.len();
    let kmers: HashMap<u64, usize> = load_kmers(kmer_files);
    bin_long_reads(kmers, bins, fastqs, output_prefix.to_string());
}

fn bin_long_reads(kmers: HashMap<u64, usize>, bins: usize, fastqs: Vec<&str>, output_prefix: String) {
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
        let filetype = fastq.split(".").collect::<Vec<&str>>();
        let filetype = filetype[filetype.len()-1];
        println!("{}",filetype);
        let gz = GzDecoder::new(file);
        let mut read_buf: String = String::new();
        let mut bin_hits: Vec<u32> = Vec::new();
        for bin in 0..bins { bin_hits.push(0); }
        //let mut next_line_seq = false;
        let reader = dna_io::DnaReader::from_path(fastq);
        //for (line_number, line) in BufReader::new(gz).lines().enumerate() {
        //    match line_number % 4 {
        //        0 => {
        //            read_buf.clear();
        //            read_buf.push_str(&line.unwrap());
        //        },
        //        1 => {
        //            let dna = DnaString::from_dna_string(&line.unwrap());
        //            for kmer_start in 0..(dna.len() - k_size + 1) {
        //                let to_hash = get_rc_invariant(&dna, kmer_start, k_size);
        //                match kmers.get(&to_hash) {
        //                    Some(bin) => bin_hits[*bin] += 1,
        //                    None => (),
        //                }     
        //            }
        //        },
        //        2 => { read_buf.push_str(&line.unwrap()); },
        //        3 => { 
        //            read_buf.push_str(&line.unwrap());
                    // decide where to bin read

        //            for bin in 0..bins { bin_hits[bin] = 0; }
        //        },
        //        _ => (),
        //    }
        //}
    }
}

fn load_kmers(kmers: Vec<&str>) -> HashMap<u64, usize> {
    let mut to_ret: HashMap<u64, usize> = HashMap::new();
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
