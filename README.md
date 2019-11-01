# fasta-downsampler
[![Travis Build Status](https://img.shields.io/travis/oskarvid/fasta-downsampler.svg?logo=travis)](https://travis-ci.org/oskarvid/fasta-downsampler)  
Make smaller fasta files and generate indexes, a bed file and a dict file

N.B: The default is that it only works for chromosomes named "chr1, chr2, chr3 ..."  
N.B2: If one of your contigs is shorter than your desired length you will end up with duplicated contigs.

## Features
* Takes your large fasta file and makes a smaller one
* Indexes your new small fasta file with samtools
* Indexes your new small fasta file with bwa
* Makes a dictionary file with gatk 

### Instructions
It runs in a docker container, install docker using the appropriate instructions here: https://docs.docker.com/install/  

Get yourself a fasta file to downsample and run this:
```bash
./downsampler.sh -d 1000 -f input.fasta
```
You must have the fasta file placed in this directory and you have to give the arguments in the order above.
