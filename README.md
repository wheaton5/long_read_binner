# long_read_binner

Install requirements: rust ver 1.3 or later, clang
```
curl https://sh.rustup.rs -sSf | sh
echo 'export PATH=~/.cargo/bin:$PATH' >> ~/.bashrc
source ~/.bashrc
which cargo
```
If the build fails on the htslib dependency you might need xz. 
```
export CFLAGS='-I/path/to/xz/<version>/include'
or add that to your .bashrc and source it
```
Then you should be able to clone and install the project.
```
git clone https://github.com/wheaton5/long_read_binner
cd long_read_binner
cargo build --release
```
And the usage is
```
./target/release/long_read_binner -h
long_read_binner 1.0
Haynes Heaton <whheaton@gmail.com>
Uses distinguishing kmers to bin long reads

USAGE:
    long_read_binner [OPTIONS] --input <input>... --kmers <kmers>... --output_prefix <output_prefix>

FLAGS:
    -h, --help       Prints help information
    -V, --version    Prints version information

OPTIONS:
        --difference_threshold <difference_threshold>
            kmer hit difference threshold for assigning a read to a kmer bin

    -i, --input <input>...
            long read files (fastq/fastq.gz, fasta/fasta.gz, sam, bam supported)

    -k, --kmers <kmers>...
            set of kmer files with 1 kmer per line. Will output a long read file for each kmer file plus an unplaced
            long read file. Kmers from each file must be the same length and mutually exclusive. Kmers of up to size 32
            supported.
        --mod_index <mod_index>                          see modimizer help
    -m, --modimizer <modimizer>
            only look at kmers where kmer % modimizer == mod_index. Use same values that you did for
            distinguishing_kmers
    -o, --output_prefix <output_prefix>                  prefix for output long read files
```
