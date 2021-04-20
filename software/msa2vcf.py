#!/usr/bin/env python3

import sys
import collections
import re
import os
import os.path
from multiprocessing import Pool

from Bio import SeqIO

def log(*args):
    print("[MSA2VCF]", *args, file=sys.stderr)

def process_sequence(record):
    mylog = []
    def log(*args):
        mylog.append("[MSA2VCF] " + " ".join((str(el) for el in args)))
    unk = set(['n', 'r', 'y', 'k', 'm', 's', 'w', 'b', 'd', 'h', 'v'])
    good_nts = set(['a', 'c', 'g', 't']).union(unk)

    rep = re.compile("[/|]")
    def normalize(id):
        return rep.sub("_", id[:60])

    with open(os.path.join(out_dir, f"{normalize(record.id)}.vcf"), "w") as vcf:
        print('##fileformat=VCFv4.1', file=vcf)
        print(f'##contig=<ID={reference.id},length={len(reference.seq)}>', file=vcf)
        print(f'##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">', file=vcf)
        print(f'#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{record.id}', file=vcf)

        log(f"Processing sequence {record.id} length={len(record.seq)}nt")
        pos = 0
        pos_ref = 0
        # Skip initial blanks
        while pos < len(record.seq) and (aln_ref.seq[pos] == '-' or record.seq[pos] == '-'):
            pos_ref += 1 if aln_ref.seq[pos] != '-' else 0
            pos += 1
        last_pos = len(record.seq)
        while last_pos > pos and (aln_ref[last_pos-1] == '-' or record.seq[last_pos-1] == '-'):
            last_pos -= 1
        log("Aligned region:", pos, "-", last_pos)
        #log(aln_ref.seq[:pos+30])
        #log(record.seq[:pos+30])
        while pos < last_pos:
            while pos < last_pos and (record.seq[pos] in unk or aln_ref.seq[pos] in unk or record.seq[pos] == aln_ref.seq[pos]):
                #log("while", record.seq[pos], aln_ref.seq[pos], reference.seq[pos_ref])
                pos_ref += 1 if pos < last_pos and aln_ref.seq[pos] != '-' else 0
                pos += 1
            if pos >= last_pos:
                break
            #log("end_w", record.seq[pos], aln_ref.seq[pos], reference.seq[pos_ref])
            start_pos = pos
            start_pos_ref = pos_ref
            while pos < last_pos and ((record.seq[pos] in good_nts and aln_ref.seq[pos] in good_nts) or (record.seq[pos] == '-' or aln_ref.seq[pos] == '-')) and ((record.seq[pos] != aln_ref.seq[pos] and record.seq[pos] not in unk and aln_ref.seq[pos] not in unk) or (record.seq[pos] == '-' and aln_ref.seq[pos] == '-')):
                pos_ref += 1 if pos < last_pos and aln_ref.seq[pos] != '-' else 0
                pos += 1
            ref_seq = str(aln_ref.seq[start_pos:pos])
            p_ref_seq = ref_seq.replace("-", "")
            var_seq = str(record.seq[start_pos:pos])
            p_var_seq = var_seq.replace("-", "")
            if len(p_var_seq) == 0:
                start_pos_ref -= 1
                log("del", start_pos, pos, f'"{ref_seq}"->"{var_seq}"')
            elif len(p_ref_seq) == 0:
                start_pos_ref -= 1
                log("ins", start_pos, pos, f'"{ref_seq}"->"{var_seq}"')
            else:
                log("snp", start_pos, pos, f'"{ref_seq}"->"{var_seq}"')
            while len(p_ref_seq) == 0 or len(p_var_seq) == 0:
                start_pos -= 1
                ref_seq = str(aln_ref.seq[start_pos:pos])
                p_ref_seq = ref_seq.replace("-", "")
                var_seq = str(record.seq[start_pos:pos])
                p_var_seq = var_seq.replace("-", "")
            #log("   ", start_pos, pos, f'"{ref_seq}"->"{var_seq}"', reference.seq[start_pos_ref-5:start_pos_ref], reference.seq[start_pos_ref], reference.seq[start_pos_ref:start_pos_ref+5])
            print('\t'.join((reference.id, str(start_pos_ref+1), '.', p_ref_seq.upper(), p_var_seq.upper(), '.', '.', '.', 'GT', '1')), file=vcf)
            #sys.exit(1)
        log("Stopped at position", pos)
    return "\n".join(mylog)

def initializer(ireference, ialn_ref, iout_dir):
    global reference
    reference = ireference
    global aln_ref
    aln_ref = ialn_ref
    global out_dir
    out_dir = iout_dir

if __name__ == '__main__':
    log("Reading reference from", sys.argv[1])
    global reference
    reference = None
    with open(sys.argv[1], "r") as handle:
        reference = SeqIO.read(handle, "fasta")
        reference.seq = reference.seq.lower()
    log(f"Read sequence {reference.id} length={len(reference.seq)}nt")

    global records
    log("Reading MSA from", sys.argv[2])
    with open(sys.argv[2], "r") as handle:
        records = list(SeqIO.parse(handle, "fasta"))
    global aln_ref
    aln_ref = records[0]
    if reference.id != aln_ref.id:
        log("The first sequence of the MSA must be the reference")
        sys.exit(1)
    os.makedirs(sys.argv[3], exist_ok=True)
    with Pool(int(sys.argv[4]) if len(sys.argv) > 4 else 1, initializer=initializer, initargs=(reference, aln_ref, sys.argv[3])) as p:
        for l in p.map(process_sequence, records[1:]):
            print(l, file=sys.stderr)
