# fasta-downsampler
Make smaller fasta files and generate indexes, a bed file and a dict file

N.B: The default is that it only works for chromosomes named "chr1, chr2, chr3 ..."

## Features
* Takes your large fasta file and makes a smaller one
* Indexes your new small fasta file with samtools
* Indexes your new small fasta file with bwa
* Makes a dictionary file with gatk 

### Instructions
It runs in a docker container, install docker using the appropriate instructions here: https://docs.docker.com/install/  

Get yourself a fasta file to downsample and run this:
```bash
./downsampler.sh 1000 input.fasta
```
You have the fasta file placed in this directory and you have to give the arguments in the order above.
