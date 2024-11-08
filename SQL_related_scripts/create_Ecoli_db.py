import sqlite3
import csv
import os

def create_tables(conn):
    c = conn.cursor()
    # Create tables
    c.execute('''CREATE TABLE IF NOT EXISTS samples (
            sample_id TEXT PRIMARY KEY,
            isolate TEXT,
            mash_group TEXT,
            pmrB TEXT,
            source TEXT,
            serotype TEXT,
            sequence_type TEXT 
            )
    ''')

    # Commit the changes
    conn.commit()

def fill_table_from_tsv(db_name, table_name, tsv_file):
    # Connect to the SQLite database
    conn = sqlite3.connect(db_name)
    c = conn.cursor()

    # Read in the data from the csv file
    with open(tsv_file, 'r') as f:
        reader = csv.reader(f, delimiter='\t')
        # Get the header row
        header = next(reader)

        # Prepare the SQL insert statement with placeholders
        placeholders = ', '.join(['?'] * len(header))
        sql_insert = f'INSERT INTO {table_name} ({", ".join(header)}) VALUES ({placeholders})'

        # Insert the data
        for row in reader:
            c.execute(sql_insert, row)
        
    # Commit the changes
    conn.commit()
    conn.close()

# Define paths
current_directory = os.getcwd()
db_path = os.path.join(current_directory, 'Ecoli.db')
tsv_path = os.path.join(current_directory, 'metadata.tsv')

print(f"Current working directory: {current_directory}")
print(f"Database path: {db_path}")
print(f"TSV file path: {tsv_path}")

# Check if the script can create the database file
print(f"Trying to create the database at: {db_path}")

# Create the database and tables
conn = sqlite3.connect(db_path)
create_tables(conn)
conn.close()

# Check if the database file was created
if os.path.exists(db_path):
    print(f"Database created successfully: {db_path}")
else:
    print(f"Failed to create the database: {db_path}")

# Fill the table with data from the TSV file
print(f"Trying to fill the table from the TSV file at: {tsv_path}")

fill_table_from_tsv(db_path, 'samples', tsv_path)

print("Finished processing.")