#!/usr/bin/python
import sqlite3
import os
import re
import pandas as pd


def create_rgi_table(conn):
    c = conn.cursor()
    c.execute(
        """
        CREATE TABLE IF NOT EXISTS rgi_results (
            my_id INTEGER PRIMARY KEY AUTOINCREMENT,
            sample_id TEXT,
            ORF_id TEXT,
            contig TEXT,
            start INTEGER,
            end INTEGER,
            orientation TEXT,
            cut_off TEXT,
            pass_bitscore TEXT,
            best_hit_bitscore TEXT,
            best_hit_ARO TEXT,
            best_identities TEXT,
            ARO TEXT,
            model_type TEXT,
            SNPs_in_best_hit_ARO TEXT,
            other_SNPs TEXT,
            drug_class TEXT,
            resistance_mechanism TEXT,
            AMR_gene_family TEXT,
            predicted_DNA TEXT,
            predicted_protein TEXT,
            CARD_protein_sequence TEXT,
            percentage_length_of_reference_sequence TEXT,
            ID TEXT,
            model_id TEXT,
            nudged TEXT,
            note TEXT,
            UNIQUE (sample_id, start, end)
        )
        """
    )
    conn.commit()


def parse_rgi_file(filepath):
    try:
        # Check if the file is empty
        if os.stat(filepath).st_size == 0:
            print(f"File {filepath} is empty. Skipping.")
            return []

        df = pd.read_csv(filepath, sep="\t")
    except pd.errors.EmptyDataError:
        print(f"No columns to parse from file {filepath}. Skipping.")
        return []

    filename = os.path.basename(filepath)
    sample_id = re.sub(r"_RGI.tsv", "", filename)

    # Check if the second line matches "# No Integron found"
    rgi_data = []
    if len(df) > 1 and "# No RGI record found" not in df.iloc[1].values:
        print(f"Processing data from file: {filename}")
        for index, row in df.iterrows():
            try:
                rgi_record = (
                    sample_id,
                    row["ORF_ID"],
                    row["Contig"],
                    row["Start"],
                    row["Stop"],
                    row["Orientation"],
                    row["Cut_Off"],
                    row["Pass_Bitscore"],
                    row["Best_Hit_Bitscore"],
                    row["Best_Hit_ARO"],
                    row["Best_Identities"],
                    row["ARO"],
                    row["Model_type"],
                    row["SNPs_in_Best_Hit_ARO"],
                    row["Other_SNPs"],
                    row["Drug Class"],
                    row["Resistance Mechanism"],
                    row["AMR Gene Family"],
                    row["Predicted_DNA"],
                    row["Predicted_Protein"],
                    row["CARD_Protein_Sequence"],
                    row["Percentage Length of Reference Sequence"],
                    row["ID"],
                    row["Model_ID"],
                    row["Nudged"],
                    row["Note"],
                )
                rgi_data.append(rgi_record)
            except ValueError as e:
                print(f"Error parsing line: {index} - {e}")

    else:
        print(f"No rgi data found in file: {filename}")
    return rgi_data


def insert_rgi_data(db_name, rgi_data):
    conn = sqlite3.connect(db_name)
    c = conn.cursor()

    sql_insert = """
        INSERT OR IGNORE INTO rgi_results (sample_id, ORF_ID, contig, start, end, orientation, cut_off, pass_bitscore, 
        best_hit_bitscore, best_hit_ARO, best_identities, ARO, model_type, SNPs_in_best_hit_ARO, other_SNPs, drug_class, 
        resistance_mechanism, AMR_gene_family, predicted_DNA, predicted_protein, CARD_protein_sequence, 
        percentage_length_of_reference_sequence, ID, model_id, nudged, note) 
        VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)
    """

    c.executemany(sql_insert, rgi_data)
    conn.commit()
    conn.close()


# Create the database and table
conn = sqlite3.connect("Ecoli.db")
create_rgi_table(conn)
conn.close()

# Directory containing CSV files
directory = "/Users/juanjovel/OneDrive/jj/UofC/data_analysis/sylviaCheckley/alyssaButters/eColi_genomics/SQLite_database/rgi"

# Parse and insert data for all integronFinder result files
for filename in os.listdir(directory):
    if filename.endswith("_RGI.tsv"):
        filepath = os.path.join(directory, filename)
        print(f"Processing file: {filepath}")
        rgi_data = parse_rgi_file(filepath)
        if rgi_data:
            insert_rgi_data("Ecoli.db", rgi_data)
            print(f"Inserted data into Ecoli.db from file: {filename}")
        else:
            print(f"No data to insert for file: {filename}")