import sys
import os
import csv

def parse_eff(info_field):
    if "EFF=" not in info_field:
        return "OTHER"

    try:
        eff_section = info_field.split("EFF=")[1].split(";")[0]
        effects = [eff.split('(')[0] for eff in eff_section.split(',')]
        effect_set = set(effects)

        if effect_set == {"NON_SYNONYMOUS_CODING"}:
            return "NONSYN"
        elif effect_set == {"SYNONYMOUS_CODING"}:
            return "SYN"
        else:
            return "OTHER"
    except Exception:
        return "OTHER"

def filter_vcf(input_vcf, output_csv):
    output_fields = ["CHROM", "POS", "REF", "ALT", "FREQ", "CLASSIFICATION"]

    with open(input_vcf, 'r') as fin, open(output_csv, 'w', newline='') as fout:
        writer = csv.writer(fout, delimiter='\t')
        writer.writerow(output_fields)

        for line in fin:
            if line.startswith("##") or line.startswith("#CHROM"):
                continue

            fields = line.strip().split('\t')
            if len(fields) < 8:
                continue

            chrom, pos, _, ref, alt, _, _, info = fields[:8]

            if ref == alt:
                continue

            classification = parse_eff(info)

            writer.writerow([chrom, pos, ref, alt, "", classification])

    print(f"File created: {output_csv}")

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("python sel_mutations.py <input.vcf> <output.csv>")
        sys.exit(1)

    input_vcf, output_csv = sys.argv[1], sys.argv[2]

    if not os.path.isfile(input_vcf):
        sys.exit(1)

    filter_vcf(input_vcf, output_csv)
