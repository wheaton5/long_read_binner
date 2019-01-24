extern crate debruijn;
use debruijn::*;
use debruijn::dna_string::*;
use debruijn::kmer::*;
use std::cmp::min;

pub fn get_rc_invariant(dna: &DnaString, kmer_size: usize) -> u64 {
    match kmer_size {
        9 => {
            let kmer: Kmer9 = dna.get_kmer(0);
            min(kmer.to_u64(), kmer.rc().to_u64())
        }
        10 => {
            let kmer: Kmer10 = dna.get_kmer(0);
            min(kmer.to_u64(), kmer.rc().to_u64())
        }
        11 => {
            let kmer: Kmer11 = dna.get_kmer(0);
            min(kmer.to_u64(), kmer.rc().to_u64())
        }
        12 => {
            let kmer: Kmer12 = dna.get_kmer(0);
            min(kmer.to_u64(), kmer.rc().to_u64())
        }
        13 => {
            let kmer: Kmer13 = dna.get_kmer(0);
            min(kmer.to_u64(), kmer.rc().to_u64())
        }
        14 => {
            let kmer: Kmer14 = dna.get_kmer(0);
            min(kmer.to_u64(), kmer.rc().to_u64())
        }
        15 => {
            let kmer: Kmer15 = dna.get_kmer(0);
            min(kmer.to_u64(), kmer.rc().to_u64())
        }
        16 => {
            let kmer: Kmer16 = dna.get_kmer(0);
            min(kmer.to_u64(), kmer.rc().to_u64())
        }
        17 => {
            let kmer: Kmer17 = dna.get_kmer(0);
            min(kmer.to_u64(), kmer.rc().to_u64())
        }
        18 => {
            let kmer: Kmer18 = dna.get_kmer(0);
            min(kmer.to_u64(), kmer.rc().to_u64())
        }
        19 => {
            let kmer: Kmer19 = dna.get_kmer(0);
            min(kmer.to_u64(), kmer.rc().to_u64())
        }
        20 => {
            let kmer: Kmer20 = dna.get_kmer(0);
            min(kmer.to_u64(), kmer.rc().to_u64())
        }
        21 => {
            let kmer: Kmer21 = dna.get_kmer(0);
            min(kmer.to_u64(), kmer.rc().to_u64())
        }
        22 => {
            let kmer: Kmer22 = dna.get_kmer(0);
            min(kmer.to_u64(), kmer.rc().to_u64())
        }
        23 => {
            let kmer: Kmer23 = dna.get_kmer(0);
            min(kmer.to_u64(), kmer.rc().to_u64())
        }
        24 => {
            let kmer: Kmer24 = dna.get_kmer(0);
            min(kmer.to_u64(), kmer.rc().to_u64())
        }
        25 => {
            let kmer: Kmer25 = dna.get_kmer(0);
            min(kmer.to_u64(), kmer.rc().to_u64())
        }
        26 => {
            let kmer: Kmer26 = dna.get_kmer(0);
            min(kmer.to_u64(), kmer.rc().to_u64())
        }
        27 => {
            let kmer: Kmer27 = dna.get_kmer(0);
            min(kmer.to_u64(), kmer.rc().to_u64())
        }
        28 => {
            let kmer: Kmer28 = dna.get_kmer(0);
            min(kmer.to_u64(), kmer.rc().to_u64())
        }
        29 => {
            let kmer: Kmer29 = dna.get_kmer(0);
            min(kmer.to_u64(), kmer.rc().to_u64())
        }
        30 => {
            let kmer: Kmer30 = dna.get_kmer(0);
            min(kmer.to_u64(), kmer.rc().to_u64())
        }
        31 => {
            let kmer: Kmer31 = dna.get_kmer(0);
            min(kmer.to_u64(), kmer.rc().to_u64())
        }
        32 => {
            let kmer: Kmer32 = dna.get_kmer(0);
            min(kmer.to_u64(), kmer.rc().to_u64())
        }
        _ => panic!("only supporting kmer sizes 9 to 32")
    }
}

type Kmer9 = VarIntKmer<u64, K9>;
#[derive(Debug, Hash, Copy, Clone, Ord, PartialOrd, Eq, PartialEq)]
pub struct K9;
impl KmerSize for K9 {
    fn K() -> usize {
        9
    }
}
type Kmer10 = VarIntKmer<u64, K10>;
#[derive(Debug, Hash, Copy, Clone, Ord, PartialOrd, Eq, PartialEq)]
pub struct K10;
impl KmerSize for K10 {
    fn K() -> usize {
        10
    }
}
type Kmer11 = VarIntKmer<u64, K11>;
#[derive(Debug, Hash, Copy, Clone, Ord, PartialOrd, Eq, PartialEq)]
pub struct K11;
impl KmerSize for K11 {
    fn K() -> usize {
        11
    }
}
type Kmer12 = VarIntKmer<u64, K12>;
#[derive(Debug, Hash, Copy, Clone, Ord, PartialOrd, Eq, PartialEq)]
pub struct K12;
impl KmerSize for K12 {
    fn K() -> usize {
        12
    }
}
type Kmer13 = VarIntKmer<u64, K13>;
#[derive(Debug, Hash, Copy, Clone, Ord, PartialOrd, Eq, PartialEq)]
pub struct K13;
impl KmerSize for K13 {
    fn K() -> usize {
        13
    }
}
type Kmer14 = VarIntKmer<u64, K14>;
#[derive(Debug, Hash, Copy, Clone, Ord, PartialOrd, Eq, PartialEq)]
pub struct K14;
impl KmerSize for K14 {
    fn K() -> usize {
        14
    }
}
type Kmer15 = VarIntKmer<u64, K15>;
#[derive(Debug, Hash, Copy, Clone, Ord, PartialOrd, Eq, PartialEq)]
pub struct K15;
impl KmerSize for K15 {
    fn K() -> usize {
        15
    }
}
type Kmer16 = VarIntKmer<u64, K16>;
#[derive(Debug, Hash, Copy, Clone, Ord, PartialOrd, Eq, PartialEq)]
pub struct K16;
impl KmerSize for K16 {
    fn K() -> usize {
        16
    }
}
type Kmer17 = VarIntKmer<u64, K17>;
#[derive(Debug, Hash, Copy, Clone, Ord, PartialOrd, Eq, PartialEq)]
pub struct K17;
impl KmerSize for K17 {
    fn K() -> usize {
        17
    }
}
type Kmer18 = VarIntKmer<u64, K18>;
#[derive(Debug, Hash, Copy, Clone, Ord, PartialOrd, Eq, PartialEq)]
pub struct K18;
impl KmerSize for K18 {
    fn K() -> usize {
        18
    }
}
type Kmer19 = VarIntKmer<u64, K19>;
#[derive(Debug, Hash, Copy, Clone, Ord, PartialOrd, Eq, PartialEq)]
pub struct K19;
impl KmerSize for K19 {
    fn K() -> usize {
        19
    }
}
type Kmer20 = VarIntKmer<u64, K20>;
#[derive(Debug, Hash, Copy, Clone, Ord, PartialOrd, Eq, PartialEq)]
pub struct K20;
impl KmerSize for K20 {
    fn K() -> usize {
        20
    }
}
type Kmer21 = VarIntKmer<u64, K21>;
#[derive(Debug, Hash, Copy, Clone, Ord, PartialOrd, Eq, PartialEq)]
pub struct K21;
impl KmerSize for K21 {
    fn K() -> usize {
        21
    }
}
type Kmer22 = VarIntKmer<u64, K22>;
#[derive(Debug, Hash, Copy, Clone, Ord, PartialOrd, Eq, PartialEq)]
pub struct K22;
impl KmerSize for K22 {
    fn K() -> usize {
        22
    }
}
type Kmer23 = VarIntKmer<u64, K23>;
#[derive(Debug, Hash, Copy, Clone, Ord, PartialOrd, Eq, PartialEq)]
pub struct K23;
impl KmerSize for K23 {
    fn K() -> usize {
        23
    }
}
type Kmer24 = VarIntKmer<u64, K24>;
#[derive(Debug, Hash, Copy, Clone, Ord, PartialOrd, Eq, PartialEq)]
pub struct K24;
impl KmerSize for K24 {
    fn K() -> usize {
        24
    }
}
type Kmer25 = VarIntKmer<u64, K25>;
#[derive(Debug, Hash, Copy, Clone, Ord, PartialOrd, Eq, PartialEq)]
pub struct K25;
impl KmerSize for K25 {
    fn K() -> usize {
        25
    }
}
type Kmer26 = VarIntKmer<u64, K26>;
#[derive(Debug, Hash, Copy, Clone, Ord, PartialOrd, Eq, PartialEq)]
pub struct K26;
impl KmerSize for K26 {
    fn K() -> usize {
        26
    }
}
type Kmer27 = VarIntKmer<u64, K27>;
#[derive(Debug, Hash, Copy, Clone, Ord, PartialOrd, Eq, PartialEq)]
pub struct K27;
impl KmerSize for K27 {
    fn K() -> usize {
        27
    }
}
type Kmer28 = VarIntKmer<u64, K28>;
#[derive(Debug, Hash, Copy, Clone, Ord, PartialOrd, Eq, PartialEq)]
pub struct K28;
impl KmerSize for K28 {
    fn K() -> usize {
        28
    }
}
type Kmer29 = VarIntKmer<u64, K29>;
#[derive(Debug, Hash, Copy, Clone, Ord, PartialOrd, Eq, PartialEq)]
pub struct K29;
impl KmerSize for K29 {
    fn K() -> usize {
        29
    }
}
type Kmer30 = VarIntKmer<u64, K30>;
#[derive(Debug, Hash, Copy, Clone, Ord, PartialOrd, Eq, PartialEq)]
pub struct K30;
impl KmerSize for K30 {
    fn K() -> usize {
        30
    }
}
type Kmer31 = VarIntKmer<u64, K31>;
#[derive(Debug, Hash, Copy, Clone, Ord, PartialOrd, Eq, PartialEq)]
pub struct K31;
impl KmerSize for K31 {
    fn K() -> usize {
        31
    }
}
