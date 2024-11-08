#!/usr/bin/python
import sqlite3
import os
import re
import pandas as pd


def create_resFinder_table(conn):
    c = conn.cursor()
    c.execute(
        """
        CREATE TABLE IF NOT EXISTS resFinder_results (
            sample_id TEXT,
            resistance_gene TEXT,
            identity TEXT,
            alignment_length_gene_length TEXT,
            coverage TEXT,
            position_in_reference TEXT,
            contig TEXT,
            position_in_contig TEXT,
            phenotype TEXT,
            accession_no TEXT,
            UNIQUE (sample_id, resistance_gene)
        )
        """
    )
    conn.commit()

def parse_resFinder_file(filepath):

    filename = os.path.basename(filepath)
    sample_id = re.sub("_ResFinder_results_tab.txt", "", filename)
    resFinder_data = []
    
    try:
        # Check if the file is empty
        if os.stat(filepath).st_size == 0:
            print(f"File {filepath} is empty. Adding empty sample record.")            
            return [(sample_id, "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA")]

        # Ensure the file is read as tab-separated
        df = pd.read_csv(filepath, sep="\t", comment="#", engine='python')
        print(f"Parsed DataFrame for file {filepath}:")
        print(df.head())  # Print the first few rows of the DataFrame
        
        # Check if the column names are as expected
        print("Columns in file:", df.columns)

    except pd.errors.EmptyDataError:
        print(f"No columns to parse from file {filepath}. Adding empty sample record.")
        return [(sample_id, "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA")]
  
    print(f"Processing data from file: {filename}")
    for index, row in df.iterrows():
        try:
            resFinder_record = (
                sample_id,
                row.get("Resistance gene", "NA").strip(),
                row.get("Identity", "NA"),
                row.get("Alignment Length/Gene Length", "NA"),
                row.get("Coverage", "NA"),
                row.get("Position in reference", "NA"),
                row.get("Contig", "NA"),
                row.get("Position in contig", "NA"),
                row.get("Phenotype", "NA"),
                row.get("Accession no.", "NA"),
            )
            print(f"Adding record: {resFinder_record}")  # Debugging: Show the record being added
            resFinder_data.append(resFinder_record)
        except ValueError as e:
            print(f"Error parsing line {index}: {e}")

    if not resFinder_data:  # If no records were added, add an empty sample record
        resFinder_data.append(
            (sample_id, "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA")
        )

    return resFinder_data


def insert_resFinder_data(db_name, resFinder_data):
    conn = sqlite3.connect(db_name)
    c = conn.cursor()

    sql_insert = """
        INSERT OR IGNORE INTO resFinder_results (sample_id, resistance_gene, identity, alignment_length_gene_length, coverage, position_in_reference, contig, position_in_contig, phenotype, accession_no) 
        VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
    """

    c.executemany(sql_insert, resFinder_data)
    conn.commit()
    conn.close()


# Create the database and table
conn = sqlite3.connect("Ecoli.db")
create_resFinder_table(conn)
conn.close()

# Directory containing result files
directory = "/Users/juanjovel/OneDrive/jj/UofC/data_analysis/sylviaCheckley/alyssaButters/eColi_genomics/SQLite_database/hybrid_data/resfinder"

# Parse and insert data for all resFinder result files
for filename in os.listdir(directory):
    if filename.endswith("_ResFinder_results_tab.txt"):
        filepath = os.path.join(directory, filename)
        print(f"Processing file: {filepath}")
        resFinder_data = parse_resFinder_file(filepath)
        if resFinder_data:
            insert_resFinder_data("Ecoli.db", resFinder_data)
            print(f"Inserted data into Ecoli.db from file: {filename}")
        else:
            print(f"No data to insert for file: {filename}")