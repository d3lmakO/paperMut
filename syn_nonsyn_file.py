import sys
import os
import pandas as pd

def merge_freq_info(main_csv, file1, file2, output_csv):
    df_main = pd.read_csv(main_csv, sep='\t', dtype=str)
    df_main['POS'] = pd.to_numeric(df_main['POS'], errors='coerce')

    for col in ['FREQ', 'PLOIDY']:
        if col in df_main.columns:
            df_main = df_main.drop(columns=[col])

    colnames = ["SAMPLEID", "CHROM", "POS", "REF", "ALT", "FREQ", "PLOIDY"]
    df1 = pd.read_csv(file1, sep='\t', header=None, names=colnames, dtype=str)
    df2 = pd.read_csv(file2, sep='\t', header=None, names=colnames, dtype=str)

    df_freq = pd.concat([df1, df2], ignore_index=True)
    df_freq['POS'] = pd.to_numeric(df_freq['POS'], errors='coerce')


    df_merged = df_main.merge(df_freq[['CHROM', 'POS', 'FREQ', 'PLOIDY']],
                              on=['CHROM', 'POS'], how='left')

    final_columns = ["CHROM", "POS", "REF", "ALT", "FREQ", "CLASSIFICATION", "PLOIDY"]
    df_merged = df_merged[final_columns]

    df_merged.to_csv(output_csv, sep='\t', index=False)
    print(f"Saved {output_csv}")

if __name__ == "__main__":
    if len(sys.argv) != 5:
        print("python syn_nonsyn_file.py <main.csv> <file1.tsv> <file2.tsv> <output.csv>")
        sys.exit(1)

    main_csv, file1, file2, output_csv = sys.argv[1:5]

    for path in [main_csv, file1, file2]:
        if not os.path.isfile(path):
            sys.exit(1)

    merge_freq_info(main_csv, file1, file2, output_csv)
