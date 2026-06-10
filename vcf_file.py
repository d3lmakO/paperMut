import pandas as pd
import sys
import os

def tsv_to_vcf(input_tsv, output_vcf):

    col_names = ["sampleID", "chr", "position", "reference", "mutation"]
    df = pd.read_csv(input_tsv, sep='\t', header=0, names=col_names, dtype=str)

    df['position'] = pd.to_numeric(df['position'], errors='coerce')

    sample_name = df['sampleID'].iloc[0] if not df['sampleID'].empty else "sample"

    vcf_df = pd.DataFrame({
        '#CHROM': df['chr'],
        'POS': df['position'],
        'ID': '.',
        'REF': df['reference'],
        'ALT': df['mutation'],
        'QUAL': '.',
        'FILTER': '.',
        'INFO': '.',
        'FORMAT': 'GT',
        sample_name: '0/1'
    })

    with open(output_vcf, 'w') as f:
        f.write("##fileformat=VCFv4.2\n")
        f.write("##source=tsv_to_vcf_converter\n")
        f.write("##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n")
        f.write("#" + "\t".join(vcf_df.columns) + "\n")
        vcf_df.to_csv(f, sep='\t', index=False, header=False)

    print(f"File VCF: {output_vcf}")

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("python vcf_file.py <input.tsv> <output.vcf>")
        sys.exit(1)

    input_tsv, output_vcf = sys.argv[1], sys.argv[2]

    if not os.path.isfile(input_tsv):
        print(f"{input_tsv} doesn't exist!")
        sys.exit(1)

    tsv_to_vcf(input_tsv, output_vcf)
