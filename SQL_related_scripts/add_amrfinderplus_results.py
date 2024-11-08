#!/usr/bin/python
import sqlite3
import os
import re
import pandas as pd

def create_amrp_table(conn):
    c = conn.cursor()
    c.execute('''
        CREATE TABLE IF NOT EXISTS amrp_results (
            id INTEGER PRIMARY KEY AUTOINCREMENT,
            sample_id TEXT,
            contig_id TEXT,
            start INTEGER,
            stop INTEGER,
            strand TEXT,
            gene_symbol TEXT,
            sequence_name TEXT,
            scope TEXT,
            element_type TEXT,
            element_subtype TEXT,
            class TEXT,
            subclass TEXT,
            method TEXT,
            target_length INTEGER,
            reference_sequence_length INTEGER,
            coverage REAL,
            identity REAL,
            alignment_length INTEGER,
            closest_accession TEXT,
            closest_name TEXT,
            hmm_id TEXT,
            hmm_description TEXT,
            UNIQUE (sample_id, contig_id, start, stop)
        )
    ''')
    conn.commit()

def parse_amrp_file(filepath):
    df = pd.read_csv(filepath, sep='\t')
    filename = os.path.basename(filepath)
    sample_id = filename
    
    amr_data = []
    for index, row in df.iterrows():
        try:
            amr_record = (
                sample_id,
                row['Contig id'],
                int(row['Start']),
                int(row['Stop']),
                row['Strand'],
                row['Gene symbol'],
                row['Sequence name'],
                row['Scope'],
                row['Element type'],
                row['Element subtype'],
                row['Class'],
                row['Subclass'],
                row['Method'],
                int(row['Target length']),
                int(row['Reference sequence length']),
                float(row['% Coverage of reference sequence']),
                float(row['% Identity to reference sequence']),
                int(row['Alignment length']),
                row['Accession of closest sequence'],
                row['Name of closest sequence'],
                row['HMM id'],
                row['HMM description']
            )
            amr_data.append(amr_record)
        except ValueError as e:
            print(f"Error parsing line: {index} - {e}")
    return amr_data

def insert_amrp_data(db_name, amr_data):
    conn = sqlite3.connect(db_name)
    c = conn.cursor()
    
    sql_insert = '''
        INSERT OR IGNORE INTO amrp_results (sample_id, contig_id, start, stop, strand, gene_symbol, sequence_name, scope, element_type, element_subtype, class, subclass, method, target_length, reference_sequence_length, coverage, identity, alignment_length, closest_accession, closest_name, hmm_id, hmm_description) 
        VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
    '''
    
    c.executemany(sql_insert, amr_data)
    conn.commit()
    conn.close()

# Create the database and table
conn = sqlite3.connect('Ecoli.db')
create_amrp_table(conn)
conn.close()

# Directory containing TSV files
directory = '/Users/juanjovel/OneDrive/jj/UofC/data_analysis/sylviaCheckley/alyssaButters/eColi_genomics/SQLite_database/hybrid_data/amrfinderplus'

# Parse and insert data for all amrp result files
for filename in os.listdir(directory):
    if filename.startswith('EC-'):
        filepath = os.path.join(directory, filename)
        print(f"Processing file: {filepath}")
        amr_data = parse_amrp_file(filepath)
        if amr_data:
            insert_amrp_data('Ecoli.db', amr_data)
            print(f"Inserted data into Ecoli.db from file: {filename}")
        else:
            print(f"No data to insert for file: {filename}")
