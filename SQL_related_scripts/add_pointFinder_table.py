#!/usr/bin/python
import sqlite3
import os
import re
import pandas as pd

def create_pointFinder_table(conn):
    c = conn.cursor()
    c.execute(
        """
        CREATE TABLE IF NOT EXISTS pointFinder_results (
            sample_id TEXT,
            mutation TEXT,
            nucleotide_change TEXT,
            amino_acid_change TEXT,
            resistance TEXT,
            pmid TEXT,
            UNIQUE (sample_id, mutation)
        )
        """
    )
    conn.commit()

def parse_pointFinder_file(filepath):
    filename = os.path.basename(filepath)
    sample_id = re.sub("_pointFinder.txt", "", filename)
    pointFinder_data = []

    try:
        # Check if the file is empty
        if os.stat(filepath).st_size == 0:
            print(f"File {filepath} is empty. Adding empty sample record.")
            return [(sample_id, "NA", "NA", "NA", "NA", "NA")]

        df = pd.read_csv(filepath, sep="\t", comment="#")
        print(f"Parsed DataFrame for file {filepath}:")
        print(df.head())  # Print the first few rows of the DataFrame

    except pd.errors.EmptyDataError:
        print(f"No columns to parse from file {filepath}. Adding empty sample record.")
        return [(sample_id, "NA", "NA", "NA", "NA", "NA")]

    if not df.empty and "# No point mutation found" not in df.values:
        print(f"Processing data from file: {filename}")
        for index, row in df.iterrows():
            try:
                pointFinder_record = (
                    sample_id,
                    row.get("Mutation", "NA"),
                    row.get("Nucleotide change", "NA"),
                    row.get("Amino acid change", "NA"),
                    row.get("Resistance", "NA"),
                    row.get("PMID", "NA")
                )
                pointFinder_data.append(pointFinder_record)
                print(f"Parsed record: {pointFinder_record}")  # Print each parsed record
            except ValueError as e:
                print(f"Error parsing line: {index} - {e}")

    if not pointFinder_data:  # If no records were added, add an empty sample record
        pointFinder_data.append((sample_id, "NA", "NA", "NA", "NA", "NA"))

    return pointFinder_data

def insert_pointFinder_data(db_name, pointFinder_data):
    conn = sqlite3.connect(db_name)
    c = conn.cursor()

    sql_insert = """
        INSERT OR IGNORE INTO pointFinder_results (sample_id, mutation, nucleotide_change, amino_acid_change, resistance, pmid) 
        VALUES (?, ?, ?, ?, ?, ?)
    """

    for record in pointFinder_data:
        print(f"Inserting record: {record}")
        c.execute(sql_insert, record)
    
    conn.commit()
    conn.close()

# Create the database and table
conn = sqlite3.connect("Ecoli.db")
create_pointFinder_table(conn)
conn.close()

# Directory containing result files
directory = "/Users/juanjovel/OneDrive/jj/UofC/data_analysis/sylviaCheckley/alyssaButters/eColi_genomics/SQLite_database/hybrid_data/resfinder"

# Parse and insert data for all pointFinder result files
for filename in os.listdir(directory):
    if filename.endswith("_PointFinder.txt"):
        filepath = os.path.join(directory, filename)
        print(f"Processing file: {filepath}")
        pointFinder_data = parse_pointFinder_file(filepath)
        insert_pointFinder_data("Ecoli.db", pointFinder_data)
        print(f"Inserted data into Ecoli.db from file: {filename}")