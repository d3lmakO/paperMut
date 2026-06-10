import pandas as pd
import sys
import os

def load_and_merge(file1, file2, output_path):
    col_names = ["sampleID", "chr", "position", "reference", "mutation"]

    df1 = pd.read_csv(file1, sep='\t', header=None, names=col_names, dtype=str)
    df2 = pd.read_csv(file2, sep='\t', header=None, names=col_names, dtype=str)
    df = pd.concat([df1, df2], ignore_index=True)
    df['chr'] = pd.to_numeric(df['chr'], errors='coerce')
    df['position'] = pd.to_numeric(df['position'], errors='coerce')
    df.sort_values(by=["chr", "position"], inplace=True)
    #df['chr'] = df['chr'].astype(str).apply(lambda x: x if x.startswith('chr') else f"chr{x}")
    df.to_csv(output_path, sep='\t', index=False)

    print("Merging complete.")

if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("python merge_dnds_files.py <file1.tsv> <file2.tsv> <output_file.csv>")
        sys.exit(1)

    file1 = sys.argv[1]
    file2 = sys.argv[2]
    output_path = sys.argv[3]

    if not os.path.isfile(file1) or not os.path.isfile(file2):
        print("One of the two files doesn't exist!")
        sys.exit(1)

    load_and_merge(file1, file2, output_path)
