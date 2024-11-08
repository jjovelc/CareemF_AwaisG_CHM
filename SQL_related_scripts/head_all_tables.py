import sqlite3
import os

os.chdir('/Users/juanjovel/OneDrive/jj/UofC/git_repos/EcoliGenomics/SQLite_database/scripts')

def head_tables(db_name, limit=5):
    conn = sqlite3.connect(db_name)
    c = conn.cursor()

    # Get the list of all tables in the database
    c.execute("SELECT name FROM sqlite_master WHERE type='table';")
    tables = c.fetchall()

    # Iterate over all tables and print the first few rows
    for table_name in tables:
        table_name = table_name[0]  # Extract table name from the tuple
        print(f"\nTable: {table_name}")
        
        try:
            # Select the first `limit` rows from the table
            c.execute(f"SELECT * FROM {table_name} LIMIT {limit};")
            rows = c.fetchall()

            # Fetch the column names
            column_names = [description[0] for description in c.description]
            print(f"Columns: {column_names}")

            # Print the first few rows
            for row in rows:
                print(row)
        except sqlite3.OperationalError as e:
            print(f"Error querying table {table_name}: {e}")

    conn.close()

# Call the function with your SQLite database
db_name = 'Ecoli.db'  # Replace with your actual database name
head_tables(db_name)