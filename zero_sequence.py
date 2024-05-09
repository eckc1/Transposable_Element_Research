import pandas as pd
from Bio import SeqIO
import sys

def extract_sequences(fasta_files, csv_file, output_file):
    """This script extracts sequences from a fasta file given a CSV file with Queries as input."""
    # Load the CSV file with the correct delimiter
    #df = pd.read_csv(csv_file, delimiter='\t')
    df = pd.read_csv(csv_file, delimiter=',')
    print("Columns found in CSV:", df.columns)  # Debug print to check column names

    # Check if 'Query' column exists
    if 'Query' not in df.columns:
        print("Error: No 'Query' column found in the CSV file.")
        return

    # Extract queries
    queries = df['Query'].unique()
    print(f"Unique queries found: {len(queries)}")

    # Dictionary to store query sequences
    query_sequences = {query: None for query in queries}

    # Search for queries in each FASTA file
    found_queries = set()
    for fasta_file in fasta_files:
        print(f"Processing {fasta_file}...")
        for record in SeqIO.parse(fasta_file, 'fasta'):
            if record.id in query_sequences and query_sequences[record.id] is None:
                query_sequences[record.id] = str(record.seq)
                found_queries.add(record.id)

    print(f"Found sequences for {len(found_queries)} queries.")

    # Write output to a new FASTA file
    with open(output_file, 'w') as f:
        for query in queries:
            if query_sequences[query]:
                f.write(f">{query}\n{query_sequences[query]}\n")
            else:
                print(f"Warning: Sequence not found for query {query}")

if __name__ == "__main__":
    if len(sys.argv) < 7:
        print(
            "Usage: python script.py <output.fasta> <input1.fasta> <input2.fasta> <input3.fasta> <input4.fasta> <input.csv>")
        sys.exit(1)

    output_fasta = sys.argv[1]
    fasta_files = sys.argv[2:6]
    csv_file = sys.argv[6]
    #fasta_files = sys.argv[2]
    #csv_file = sys.argv[3]

    extract_sequences(fasta_files, csv_file, output_fasta)
