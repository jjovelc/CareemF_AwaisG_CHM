#!/usr/bin/python
import sqlite3
import pandas as pd

def create_prediction_table(conn):
    c = conn.cursor()
    c.execute(
        """
        CREATE TABLE IF NOT EXISTS rfplasmid_results (
            sample_id TEXT,
            prediction TEXT,
            votes_chromosomal REAL,
            votes_plasmid REAL,
            contigID TEXT,
            UNIQUE (sample_id)
        )
        """
    )
    conn.commit()

def insert_prediction_data(db_name, prediction_data):
    conn = sqlite3.connect(db_name)
    c = conn.cursor()

    sql_insert = """
        INSERT OR IGNORE INTO rfplasmid_results (sample_id, prediction, votes_chromosomal, votes_plasmid, contigID) 
        VALUES (?, ?, ?, ?, ?)
    """

    for index, row in prediction_data.iterrows():
        # Now use the index as sample_id
        record = (index, row['prediction'], row['votes chromosomal'], row['votes plasmid'], row['contigID'])
        print(f"Inserting record: {record}")
        c.execute(sql_insert, record)
    
    conn.commit()
    conn.close()

# Create the database and table
conn = sqlite3.connect("Ecoli.db")
create_prediction_table(conn)
conn.close()


# Load the prediction.csv data
csv_file = '/Users/juanjovel/OneDrive/jj/UofC/data_analysis/sylviaCheckley/alyssaButters/eColi_genomics/SQLite_database/hybrid_data/rfplasmid/prediction.csv'
prediction_data = pd.read_csv(csv_file)
prediction_data = pd.read_csv(csv_file, index_col=0)

# Insert the prediction data into the database
insert_prediction_data("Ecoli.db", prediction_data)

print("Inserted prediction data into Ecoli.db")
