#!/usr/bin/env python3

import sys
from collections import defaultdict
from pysam import VariantFile

def add_freqs():
    vcf_path = sys.argv[2]
    fai_path = sys.argv[3]
    min_samples = int(sys.argv[4]) if len(sys.argv) == 5 else 0

    vcf = VariantFile(vcf_path, 'r', drop_samples = False)

    ref_name, ref_len = open(fai_path).readlines()[0].strip('\n').split('\t')[0:2]
    new_contig = f"##contig=<ID={ref_name},length={ref_len}>"

    vcf.header.add_line(new_contig)
    vcf.header.add_line("##INFO=<ID=AF,Number=A,Type=Float,Description=\"Estimated allele frequency in the range (0,1)\">")
    print('\n'.join(str(vcf.header).split('\n')[:-1]))

    removed_variants = 0
    removed_some_alt = 0
    tot_removed_alleles = 0
    tot_removed_genotypes = 0
    next_milestone = 1000
    for record in vcf:
        if record.pos <= 50 or record.pos >= int(ref_len)-50:
            continue
        if record.pos > next_milestone:
            print("Reached position", record.pos, file=sys.stderr)
            next_milestone += 1000

        n_gts = defaultdict(int)
        low_samples = defaultdict(list)
        for sample in record.samples.itervalues():
            curr_gt = sample.allele_indices[0]
            n_gts[curr_gt] += 1
            if curr_gt != 0 and n_gts[curr_gt] < min_samples:
                low_samples[curr_gt].append(sample)

        to_delete = set([gt for gt in n_gts if n_gts[gt] < min_samples])
        if len(to_delete) > 0:
            if len(n_gts) - len(to_delete) == 1:
                removed_variants += 1
                continue
            tot_removed_alleles += len(to_delete)
            removed_some_alt += 1

            for gt in to_delete:
                for sample in low_samples[gt]:
                    sample.allele_indices = (0,)
                n_gts[0] += n_gts[gt]
                tot_removed_genotypes += n_gts[gt]
                n_gts[gt] = 0

        tot_alleles = len(n_gts)
        tot_samples = sum(n_gts.values())
        for gt in n_gts:
            n_gts[gt] /= tot_samples
            n_gts[gt] = round(n_gts[gt], 6)
        alt_freqs = [n_gts[i] for i in n_gts if i > 0]
        record.chrom = ref_name
        record.info.__setitem__("AF", alt_freqs)
        print(record, end='', flush=False)
    print("removed_variants=", removed_variants, file=sys.stderr)
    print("removed_some_alt=", removed_some_alt, file=sys.stderr)
    print("tot_removed_alleles=", tot_removed_alleles, file=sys.stderr)
    print("tot_removed_genotypes=", tot_removed_genotypes, file=sys.stderr)

def clean_header():
    vcf_path = sys.argv[2]
    for line in open(vcf_path, "r"):
        if line.startswith("#C"):
            fields = line.strip("\n").split("\t")
            repeats = {}
            new_fields = []
            for field in fields:
                if field in repeats:
                    repeats[field] = repeats[field]+1
                    field = field + "_" + str(repeats[field])
                else:
                    repeats[field] = 1
                new_fields.append(field)
            line = "\t".join(new_fields) + "\n"
        print(line, end="")

if __name__ == "__main__":
    if sys.argv[1] == "clean":
        clean_header()
    elif sys.argv[1] == "freq":
        add_freqs()
