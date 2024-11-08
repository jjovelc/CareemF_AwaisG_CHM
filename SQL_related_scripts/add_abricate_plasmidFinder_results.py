#!/usr/bin/python
import sqlite3
import os
import sys
import re
import pandas as pd


def create_plasmFinder_table(conn):
    c = conn.cursor()
    c.execute(
        """
        CREATE TABLE IF NOT EXISTS abricate_plasmidfinder_results (
            id INTEGER PRIMARY KEY AUTOINCREMENT,
            sample_id TEXT,
            contig TEXT,
            start INTEGER,
            end INTEGER,
            strand TEXT,
            gene TEXT,
            coverage TEXT,
            coverage_map TEXT,
            gaps TEXT,
            percent_coverage REAL,
            percent_identity REAL,
            database TEXT,
            accession TEXT,
            product TEXT,
            resistance TEXT,

            UNIQUE (sample_id, contig, start, end, gene)
        )
    """
    )
    conn.commit()


def parse_plasmFinder_file(filepath):
    df = pd.read_csv(filepath, sep="\t")
    filename = os.path.basename(filepath)
    sample_id = re.sub("_plasmidfinder.tab", "", filename)

    plasmFinder_data = []
    for index, row in df.iterrows():
        try:
            plasmFinder_record = (
                sample_id,
                row["SEQUENCE"],
                int(row["START"]),
                int(row["END"]),
                row["STRAND"],
                row["GENE"],
                row["COVERAGE"],
                row["COVERAGE_MAP"],
                row["GAPS"],
                float(row["%COVERAGE"]),
                float(row["%IDENTITY"]),
                row["DATABASE"],
                row["ACCESSION"],
                row["PRODUCT"],
                row["RESISTANCE"],
            )
            plasmFinder_data.append(plasmFinder_record)
        except ValueError as e:
            print(f"Error parsing line: {index} - {e}")
    return plasmFinder_data


def insert_plasmFinder_data(db_name, plasmFinder_data):
    conn = sqlite3.connect(db_name)
    c = conn.cursor()

    sql_insert = """
        INSERT OR IGNORE INTO abricate_plasmidfinder_results (sample_id, contig, start, end, strand, gene, coverage, coverage_map, gaps, percent_coverage, percent_identity, database, accession, product, resistance) 
        VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
    """

    c.executemany(sql_insert, plasmFinder_data)
    conn.commit()
    conn.close()


# Create the database and table
conn = sqlite3.connect("Ecoli.db")
create_plasmFinder_table(conn)
conn.close()

# Directory containing CSV files
directory = "/Users/juanjovel/OneDrive/jj/UofC/data_analysis/sylviaCheckley/alyssaButters/eColi_genomics/SQLite_database/hybrid_data/abricate_plasmidfinder"

# Parse and insert data for all plasmidFinder result files
for filename in os.listdir(directory):
    if filename.endswith(".tab"):
        filepath = os.path.join(directory, filename)
        print(f"Processing file: {filepath}")
        plasmFinder_data = parse_plasmFinder_file(filepath)
        if plasmFinder_data:
            insert_plasmFinder_data("Ecoli.db", plasmFinder_data)
            print(f"Inserted data into Ecoli.db from file: {filename}")
        else:
            print(f"No data to insert for file: {filename}")