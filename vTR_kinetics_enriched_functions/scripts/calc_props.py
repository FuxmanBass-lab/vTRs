import pandas as pd

def normalize_table(input_tsv: str, output_tsv: str):
    """
    Reads a TSV file with counts, normalizes each column (except the first 'Protein' column)
    to proportions by dividing by the column sum, and saves the result as a TSV.
    """
    # Load table
    df = pd.read_csv(input_tsv, sep="\t")

    # Keep protein column separately
    protein_col = df.iloc[:, 0]  # first column
    counts = df.iloc[:, 1:]      # all numeric columns

    # Compute column sums
    col_sums = counts.sum(axis=0)

    # Normalize to proportions
    proportions = counts.div(col_sums, axis=1)

    # Reattach protein column
    result = pd.concat([protein_col, proportions], axis=1)

    # Save as TSV
    result.to_csv(output_tsv, sep="\t", index=False)

if __name__ == "__main__":
    import sys
    if len(sys.argv) != 3:
        print("Usage: python calc_props.py <input.tsv> <output.tsv>")
        sys.exit(1)
    
    input_file = sys.argv[1]
    output_file = sys.argv[2]
    normalize_table(input_file, output_file)