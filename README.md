# long_read_binner


```
./target/release/long_read_binner -h
long_read_binner
Haynes Heaton <whheaton@gmail.com>
Uses distinguishing kmers to bin long reads

USAGE:
    long_read_binner --fastqs <fastqs>... --kmers <kmers>... --output_prefix <output_prefix>

FLAGS:
    -h, --help       Prints help information
    -V, --version    Prints version information

OPTIONS:
    -f, --fastqs <fastqs>...               long read fastq.gz files
    -k, --kmers <kmers>...                 set of kmer files with 1 kmer per line. Will output a long read file for each
                                           kmer file plus an unplaced long read file.
    -o, --output_prefix <output_prefix>    prefix for output long read files
```
