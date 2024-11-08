#!/usr/bin/python
import sqlite3
import os
import re
import pandas as pd


def create_intFinder_table(conn):
    c = conn.cursor()
    c.execute(
        """
        CREATE TABLE IF NOT EXISTS intFinder_results (
            id INTEGER PRIMARY KEY AUTOINCREMENT,
            sample_id TEXT,
            int_id TEXT,
            contig TEXT,
            start INTEGER,
            end INTEGER,
            strand TEXT,
            evalue TEXT,
            type_elt TEXT,
            annotation TEXT,
            model TEXT,
            type TEXT,
            default_value TEXT,
            distance_2attC TEXT,
            considered_topology TEXT,
            UNIQUE (sample_id, start, end)
        )
        """
    )
    conn.commit()


def parse_intFinder_file(filepath):
    try:
        # Check if the file is empty
        if os.stat(filepath).st_size == 0:
            print(f"File {filepath} is empty. Skipping.")
            return []

        df = pd.read_csv(filepath, sep="\t", comment="#")
    except pd.errors.EmptyDataError:
        print(f"No columns to parse from file {filepath}. Skipping.")
        return []

    filename = os.path.basename(filepath)
    sample_id = re.sub("_final.integrons", "", filename)
    

    # Check if the second line matches "# No Integron found"
    intFinder_data = []
    if len(df) > 1 and "# No Integron found" not in df.iloc[1].values:
        print(f"Processing data from file: {filename}")
        for index, row in df.iterrows():
            try:
                intFinder_record = (
                    sample_id,
                    row.get("ID_integron", "NA"),
                    row.get("ID_replicon", "NA"),
                    row.get("pos_beg", "NA"),
                    row.get("pos_end", "NA"),
                    row.get("strand", "NA"),
                    row.get("evalue", "NA"),
                    row.get("type_elt", "NA"),
                    row.get("annotation", "NA"),
                    row.get("model", "NA"),
                    row.get("type", "NA"),
                    row.get("default", "NA"),
                    row.get("distance_2attC", "NA"),
                    row.get("considered_topology", "NA"),
                )
                intFinder_data.append(intFinder_record)
            except ValueError as e:
                print(f"Error parsing line: {index} - {e}")
        
    else:
        print(f"No integron data found in file: {filename}")
    return intFinder_data


def insert_intFinder_data(db_name, intFinder_data):
    conn = sqlite3.connect(db_name)
    c = conn.cursor()

    sql_insert = """
        INSERT OR IGNORE INTO intFinder_results (sample_id, int_id, contig, start, end, strand, evalue, type_elt, annotation, model, type, default_value, distance_2attC, considered_topology) 
        VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
    """

    c.executemany(sql_insert, intFinder_data)
    conn.commit()
    conn.close()


# Create the database and table
conn = sqlite3.connect("Ecoli.db")
create_intFinder_table(conn)
conn.close()

# Directory containing CSV files
directory = "/Users/juanjovel/OneDrive/jj/UofC/data_analysis/sylviaCheckley/alyssaButters/eColi_genomics/SQLite_database/hybrid_data/integron_finder"

# Parse and insert data for all integronFinder result files
for filename in os.listdir(directory):
    if filename.endswith(".integrons"):
        filepath = os.path.join(directory, filename)
        print(f"Processing file: {filepath}")
        intFinder_data = parse_intFinder_file(filepath)
        if intFinder_data:
            insert_intFinder_data("Ecoli.db", intFinder_data)
            print(f"Inserted data into Ecoli.db from file: {filename}")
        else:
            print(f"No data to insert for file: {filename}")
