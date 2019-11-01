#!/bin/bash

# Uncomment for debugging
set -Eeu
#set -oxv xtrace

LANG=C

trap die SIGINT SIGKILL

source functions.sh

if [[ $# -ne 4 ]]; then
	usage
fi

while getopts 'd:f:h' flag; do
	case "${flag}" in
	d)
		d=${OPTARG}
		if [[ $d != [0-9]*  ]]; then
			err "You forgot to supply a desired new contig length"
			usage
		elif [[ $d == [0-9]* ]]; then
			SIZE=$d
		fi
		;;
	f)
		f=${OPTARG}
		if [[ -f $f && $f == *.fasta || $f == *.fna || $f == *.fa ]]; then
			FASTA=$f
		elif [[ ! -f $f || $f != *.fasta ]]; then
			err "You must use an uncompressed fasta file as input"
			usage
		fi
		;;
	h)
		usage
		;;
	?)
		usage
		;;
	esac
done

mkdir -p outputs
# Remove previous output files if you want to make them smaller or bigger
inf "Removing previous output files"
rm -f outputs/* 2>/dev/null

# I tried skipping the temp.fasta file by piping these grep and sed commands but it was way slower than splitting it like this
# Piping took 1 min 12s while splitting took 22s
inf "Running grep"
grep -A $SIZE -i ">chr" $FASTA >> temp.fasta
inf "Running sed"
sed 's/LN:[0-9]*/LN:'"$SIZE"'/' temp.fasta > outputs/tiny.fasta

# There's no use for this file, remove it!
inf "Removing temporary file"
rm temp.fasta

# Docker to the rescue once again!
inf "Making fasta index file"
docker run --rm -ti -u $UID:1000 -v $(pwd):/data -w /data oskarv/snakemake-germline-tools:4.1.2.0 samtools faidx outputs/tiny.fasta

# Fancy oneliner for bed making from a fai file
inf "It's time to make a bed file with awk"
awk 'BEGIN {FS="\t"}; {print $1 FS "0" FS $2}' outputs/tiny.fasta.fai > outputs/tiny.bed

# Docker to the rescue once again once again!
#inf "Executing docker to make a dictionary"
#docker run --rm -ti -u $UID:1000 -v $(pwd):/data -w /data oskarv/snakemake-germline-tools:4.1.2.0 gatk CreateSequenceDictionary -R outputs/tiny.fasta

# Run bwa index
#inf "Running bwa index"
#docker run --rm -ti -u $UID:1000 -v $(pwd):/data -w /data oskarv/snakemake-germline-tools:4.1.2.0 bwa index -a bwtsw outputs/tiny.fasta

# Use pigz because the name is funny and also because it's parallelized! Foar moar coares!
if [[ $(which pigz) ]]; then
	inf "Running pigz, single threaded compression is for losers!"
	pigz -c outputs/tiny.fasta > outputs/tiny.fasta.gz
else
	inf "I haven't installed pigz so I will only use one thread for compression. This is why we can't have nice things."
	gzip -c outputs/tiny.fasta > outputs/tiny.fasta.gz
fi

succ "All done!"
