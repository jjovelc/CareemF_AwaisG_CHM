#!/usr/bin/python
import sqlite3
import os
import csv
import re

def create_mge_table(conn):
    c = conn.cursor()
    c.execute('''
        CREATE TABLE IF NOT EXISTS mge_results (
            id INTEGER PRIMARY KEY AUTOINCREMENT,
            sample_id TEXT,
            mge_no INTEGER,
            name TEXT,
            synonyms TEXT,
            prediction TEXT,
            type TEXT,
            allele_len INTEGER,
            depth REAL,
            e_value REAL,
            identity REAL,
            coverage REAL,
            gaps INTEGER,
            substitution INTEGER,
            contig TEXT,
            start INTEGER,
            end INTEGER,
            cigar TEXT,
            UNIQUE (sample_id, mge_no, name, start, end)
        )
    ''')
    conn.commit()

def parse_mge_file(file_path):
    filename = os.path.basename(file_path)
    sample_id = re.sub(r'_mge.csv', '', filename)

    with open(file_path, 'r') as file:
        lines = file.readlines()
        # Extract mge data
        mge_data = []
        for line in lines:
            if not line.startswith('#'):
                parts = line.strip().split(',')
                if len(parts) == 16 and parts[0] != "mge_no":  # Skip header line and ensure proper format
                    try:
                        # Handle missing values by providing default values if the fields are empty
                        mge_record = (
                            sample_id,
                            int(parts[0]),  # mge_no
                            parts[1],       # name
                            parts[2],       # synonyms
                            parts[3],       # prediction
                            parts[4],       # type
                            int(parts[5]) if parts[5] else None,  # allele_len, default to None if empty
                            float(parts[6]) if parts[6] else None, # depth, default to None if empty
                            float(parts[7]) if parts[7] else None, # e_value, default to None if empty
                            float(parts[8]),  # identity
                            float(parts[9]),  # coverage
                            int(parts[10]) if parts[10] else 0,   # gaps, default to 0 if empty
                            int(parts[11]) if parts[11] else 0,   # substitution, default to 0 if empty
                            parts[12],      # contig
                            int(parts[13]), # start
                            int(parts[14]), # end
                            parts[15]       # cigar
                        )
                        mge_data.append(mge_record)
                    except ValueError as e:
                        print(f"Error parsing line: {line} - {e}")
        return mge_data
    
def insert_mge_data(db_name, mge_data):
    conn = sqlite3.connect(db_name)
    c = conn.cursor()
    
    sql_insert = '''
        INSERT OR IGNORE INTO mge_results (sample_id, mge_no, name, synonyms, prediction, type, allele_len, depth, e_value, identity, coverage, gaps, substitution, contig, start, end, cigar) 
        VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
    '''
    
    c.executemany(sql_insert, mge_data)
    conn.commit()
    conn.close()

# Create the database and table
conn = sqlite3.connect('Ecoli.db')
create_mge_table(conn)
conn.close()

# Directory containing CSV files
directory = '/Users/juanjovel/OneDrive/jj/UofC/data_analysis/sylviaCheckley/alyssaButters/eColi_genomics/SQLite_database/hybrid_data/mobile_element_finder'

# Parse and insert data for all mge result files
for filename in os.listdir(directory):
    if filename.endswith('.csv'):
        filepath = os.path.join(directory, filename)
        print(f"Processing file: {filepath}")
        mge_data = parse_mge_file(filepath)
        if mge_data:
            insert_mge_data('Ecoli.db', mge_data)
            print(f"Inserted data into Ecoli.db from file: {filename}")
        else:
            print(f"No data to insert for file: {filename}")
