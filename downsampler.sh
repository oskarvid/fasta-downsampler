#!/bin/bash

# Uncomment for debugging
#set -eoxv xtrace

# Check that the desired size is given as first input
if [[ $1 != [0-9]*  ]]; then
	echo "You forgot to supply a desired new fasta file length"
	exit 0
fi
SIZE=$1

# Check that the fasta file is uncompressed
if [[ $2 != *.fasta ]]; then
	echo "You must use an uncompressed fasta file as input"
	exit 0
fi

# Remove previous output files if you want to make them smaller or bigger
echo "Removing previous output files"
rm -f outputs/* 2>/dev/null

# I tried skipping the temp.fasta file by piping these grep and sed commands but it was way slower than splitting it like this
# Piping took 1 min 12s while splitting took 22s
echo "Running grep"
grep -A $SIZE ">chr" $2 >> temp.fasta
echo "Running sed"
sed 's/LN:[0-9]*/LN:'"$SIZE"'/' temp.fasta > outputs/tiny.fasta

# There's no use for this file, remove it!
echo "Removing temporary file"
rm temp.fasta

# Docker to the rescue once again!
echo "Making fasta index file"
docker run --rm -ti -u $UID:1000 -v $(pwd):/data -w /data oskarv/snakemake-germline-tools:4.1.2.0 samtools faidx outputs/tiny.fasta

# Fancy oneliner for bed making from a fai file
echo "It's time to make a bed file with awk"
awk 'BEGIN {FS="\t"}; {print $1 FS "0" FS $2}' outputs/tiny.fasta.fai > outputs/tiny.bed

# Docker to the rescue once again once again!
echo "Executing docker to make a dictionary"
docker run --rm -ti -u $UID:1000 -v $(pwd):/data -w /data oskarv/snakemake-germline-tools:4.1.2.0 gatk CreateSequenceDictionary -R outputs/tiny.fasta

# Run bwa index
echo "Running bwa index"
docker run --rm -ti -u $UID:1000 -v $(pwd):/data -w /data oskarv/snakemake-germline-tools:4.1.2.0 bwa index -a bwtsw outputs/tiny.fasta

# Use pigz because the name is funny and also because it's parallelized! Foar moar coares!
if [[ $(which pigz) ]]; then
	echo "Running pigz, single threaded compression is for losers!"
	pigz -c outputs/tiny.fasta > outputs/tiny.fasta.gz
else
	echo "I haven't installed pigz so I will only use one thread for compression. This is why we can't have nice things."
	gzip -c outputs/tiny.fasta > outputs/tiny.fasta.gz
fi

echo "All done!"
