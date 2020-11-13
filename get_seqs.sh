#!/bin/bash

echo -n "Compiling required software..."
pushd $(pwd)/software > /dev/null
make &> /dev/null
popd > /dev/null
pushd $(pwd)/software/snp-sites-dir > /dev/null
autoreconf -i -f &> /dev/null
./configure --prefix=$(pwd)/../ &> /dev/null
make &> /dev/null
make install &> /dev/null
popd > /dev/null
echo "done."

export PATH=$(pwd)/software/bin:$(pwd)/software:$PATH

# DATE_FMT is used as part of the directory name
# LAST_TIME is used in info.json and status.json
# They have two different format so we need to get them twice
DATE_FMT=$(date +%Y%m%d-%H%M%S)
LAST_TIME=$(python3 -c 'import datetime; print(str(datetime.datetime.today().replace(microsecond=0)))')

THREADS=8
ACC_FILENAME="accession.id.list"
REF_FILENAME="reference.fasta"
SEQ_FILENAME="seqs.fasta"
GFF_FILENAME="GCF_009858895.2_ASM985889v3_genomic.gff"
REF_ID="NC_045512"

# Check if we already have a list of genomes we downloaded
# previously.  If so, move it.
if [ -f $ACC_FILENAME ]; then
    mv $ACC_FILENAME $ACC_FILENAME.bkp
fi

# Get IDS of the complete genomes for SARS-CoV-2
# We filter out the reference genome since we need it in
# a different file
echo -n "Querying entrez..."
esearch -db nucleotide -query "SARS-CoV-2 [ORGN] AND complete genome [TITLE]" |
    efetch -format acc |
    grep -v $REF_ID > $ACC_FILENAME
echo "done."

# Check if we have some new nucleotide in the acc list
if [ -f $ACC_FILENAME.bkp ]; then
    diff $ACC_FILENAME $ACC_FILENAME.bkp > /dev/null

    if [ $? -eq 0 ]; then
        echo "No new complete genomes found in entrez."
        mv $ACC_FILENAME.bkp $ACC_FILENAME
        exit 1
    fi
fi

echo "Downloading complete genomes sequences..."
# Format the sequences
tr '\n' ',' < $ACC_FILENAME > $ACC_FILENAME.temp
sed -i 's/,$//' $ACC_FILENAME.temp
echo "There are" $(wc -l $ACC_FILENAME) "sequences to download"

# Download the sequences
rm -f $SEQ_FILENAME.temp
split -t "," -l 1000 -a 3 $ACC_FILENAME.temp $ACC_FILENAME.temp.chunk_
for f in $ACC_FILENAME.temp.chunk_*
do
	efetch -format fasta -db nucleotide -id $(cat $f) >> $SEQ_FILENAME.temp
	echo "Downloaded" $(grep -c "^>" $SEQ_FILENAME.temp) "sequences"
done
rm -f $ACC_FILENAME.temp.chunk_*

# Clean the sequences
clean_fasta.py clean $SEQ_FILENAME.temp > $SEQ_FILENAME
rm $ACC_FILENAME.temp $SEQ_FILENAME.temp
echo "done."

# Create the directory
echo -n "Creating the db..."
IDD=NCBI-SARS-CoV-2-$DATE_FMT
UNIQ_ID=vcf/$IDD
mkdir -p $UNIQ_ID

# Copy the info
cp refs/NC_045512.2/$REF_FILENAME $UNIQ_ID
cp $SEQ_FILENAME $UNIQ_ID
cp refs/NC_045512.2/$GFF_FILENAME $UNIQ_ID

# Compute the VCF
REF=$UNIQ_ID/$REF_FILENAME
SEQ=$UNIQ_ID/$SEQ_FILENAME
GFF=$UNIQ_ID/$GFF_FILENAME
VCF=$UNIQ_ID/run.cleaned.vcf
INFO=$UNIQ_ID/info.json
STAT=$UNIQ_ID/status.json
TMPDIR=/tmp/$UNIQ_ID/
MSA=$TMPDIR/output.msa

mkdir -p $TMPDIR

mafft --thread $THREADS --auto --reorder --mapout --keeplength --addfragments $SEQ $REF > $MSA.unfilled 2> /dev/null
fill_msa $MSA.unfilled > $MSA
snp-sites -rmcv -o $TMPDIR/run $MSA
format_vcf.py clean $TMPDIR/run.vcf | cut -f 1-9,11- > $TMPDIR/run.1.vcf
samtools faidx $REF
format_vcf.py freq $TMPDIR/run.1.vcf $REF.fai 2 > $VCF
gzip -9v $VCF

echo "done."

rm $UNIQ_ID/*.map
rm $UNIQ_ID/*.fai

echo -n "Creating info.json and status.json..."

cat <<EOF > $STAT
{
	"status": "Precomputed",
	"last_time": "$LAST_TIME",
	"output": {
		  "vcf": "/jobs/$VCF.gz",
		  "reference": "/jobs/$REF"
	}
}
EOF

cat <<EOF > $INFO
{
	"filename": "/jobs/$VCF.gz",
	"id": "$IDD",
	"description": "SARS-CoV-2 precomputed VCF and ref downloaded from NCBI on $(date +%Y-%m-%d).\n List of genomes included: $(cat $ACC_FILENAME | tr '\n' ',' | sed 's/,$//' | sed 's/,/, /g')",
	"submission_time": $(date +%s),
	"alias": "NCBI-SARS-CoV-2-$DATE_FMT",
	"reference": "/jobs/$REF",
	"gtf": "/jobs/$GFF",
	"internal_ref": {
		"snpEff": {
			"id": "NC_045512.2"
		},
		"reference": {
			"file": "NC_045512.2/reference.fasta"
		},
		"id": "NC_045512.2",
		"annotation": {
			"file": "NC_045512.2/GCF_009858895.2_ASM985889v3_genomic.gff"
		},
		"alias": "SARS-CoV-2, Severe acute respiratory syndrome coronavirus 2 isolate Wuhan-Hu-1"
	}
}
EOF

echo "done"

echo -n "Creating git commit..."
NO_OF_GENOMES=$(grep -c "^>" $SEQ_FILENAME)
git add $UNIQ_ID > /dev/null
cat <<EOF > commit-msg
[AUTO] SARS-CoV-2 NCBI - $(date +%Y-%m-%d) automatic update

This release comprises $NO_OF_GENOMES SARS-CoV-2 genomes downloaded from NCBI.

List of accession keys of each genome in this release:
EOF
grep "^>" $SEQ_FILENAME | cut -d' ' -f1 | tr -d '>' | tr '\n' ', ' | sed 's/,$//' | sed 's/,/, /g' | fold -s -w 120 >> commit-msg
git commit -F commit-msg

# git push
echo "done."
