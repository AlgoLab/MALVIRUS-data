#!/bin/sh

# DATE_FMT is used as part of the directory name
# LAST_TIME is used in info.json and status.json
# They have two different format so we need to get them twice
DATE_FMT=$(date +%Y%m%d-%H%M%S)
LAST_TIME=$(python3 -c 'import datetime; print(str(datetime.datetime.today().replace(microsecond=0)))')

THREADS=3
ACC_FILENAME="accession.id.list"
REF_FILENAME="reference.fasta"
SEQ_FILENAME="seqs.fasta"
GFF_FILENAME="sars-cov-2.gff"
GFF_URL="https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/009/858/895/GCF_009858895.2_ASM985889v3/GCF_009858895.2_ASM985889v3_genomic.gff.gz"
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

# Check if the reference exists already
if [ ! -f $REF_FILENAME ]; then
    echo -n "reference.fasta not found, downloading..."
    efetch -format fasta -db nucleotide -id $REF_ID > $REF_FILENAME
    echo "done."
fi

# Check if we have some new nucleotide in the acc list
if [ -f $ACC_FILENAME.bkp ]; then
    diff $ACC_FILENAME $ACC_FILENAME.bkp > /dev/null

    if [ $? -eq 0 ]; then
	echo "No new complete genomes found in entrez."
	mv $ACC_FILENAME.bkp $ACC_FILENAME
	exit 1
    fi
fi

echo -n "Downloading complete genomes sequences..."
# Format the sequences
tr '\n' ',' < $ACC_FILENAME > $ACC_FILENAME.temp
sed -i 's/,$//' $ACC_FILENAME.temp

# Download the sequences
efetch -format fasta -db nucleotide -id $(cat $ACC_FILENAME.temp) > $SEQ_FILENAME.temp

# Clean the sequences
clean_fasta.py clean $SEQ_FILENAME.temp > $SEQ_FILENAME
rm $ACC_FILENAME.temp $SEQ_FILENAME.temp
echo "done."

# Download the annotation
echo -n "Downloading the annotation..."
wget -q -O - $GFF_URL | gunzip -c > $GFF_FILENAME
echo "done."

# Create the directory
echo -n "Creating the db..."
IDD=NCBI-SARS-CoV-2-$DATE_FMT
UNIQ_ID=vcf/$IDD
mkdir -p $UNIQ_ID

# Copy the info
cp $REF_FILENAME $UNIQ_ID
cp $SEQ_FILENAME $UNIQ_ID
cp $GFF_FILENAME $UNIQ_ID

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
format_vcf.py freq $TMPDIR/run.1.vcf $REF.fai > $VCF

echo "done."

rm $UNIQ_ID/*.map
rm $UNIQ_ID/*.fai

echo -n "Creating info.json and status.json..."

cat <<EOF > $STAT
{
	"status": "Precomputed",
	"last_time": "$LAST_TIME",
	"output": {
		  "vcf": "/jobs/$VCF",
		  "reference": "/jobs/$REF"
	}
}
EOF

cat <<EOF > $INFO
{
	"filename": "/jobs//$VCF",
	"id": "$IDD",
	"description": "SARS-CoV-2 precomputed VCF and ref downloaded from NCBI on $(date +%Y-%m-%d).\n List of genomes included: $(cat $ACC_FILENAME | tr '\n' ',' | sed 's/,$//' | sed 's/,/, /g')",
	"submission_time": $(date +%s),
     	"alias": "NCBI-SARS-CoV-2-$DATE_FMT",
	"reference": "/jobs/$REF",
	"gtf": "/jobs/$GFF"
}
EOF

echo "done"

echo -n "Creating git commit..."
NO_OF_GENOMES=$(grep -c "^>" $SEQ_FILENAME)
git add $UNIQ_ID > /dev/null
git commit -m """[AUTO] SARS-CoV-2 NCBI - $(date +%Y-%m-%d) automatic update

This release comprises $NO_OF_GENOMES SARS-CoV-2 genomes downloaded from NCBI.

List of accession keys of each genome in this release:
$(grep "^>" $SEQ_FILENAME | cut -d' ' -f1 | tr -d '>' | tr '\n' ', ' | sed 's/,$//' | sed 's/,/, /g')
"""
# git push
echo "done."
