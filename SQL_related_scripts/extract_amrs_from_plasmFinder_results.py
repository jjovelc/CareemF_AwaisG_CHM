import sqlite3
import pandas as pd

# Connect to the SQLite database
conn = sqlite3.connect('Ecoli.db')  # Replace with your database name

# SQL query to select necessary columns
query = """
SELECT sample_id, gene
FROM plasmFinder_results
"""

# Read data into a pandas DataFrame
df = pd.read_sql_query(query, conn)

# Close the connection
conn.close()

# Display the first few rows of the dataframe
print("Initial DataFrame:")
print(df.head())

# Create a pivot table with gene_symbol as rows and sample_id as columns
pivot_df = df.pivot_table(index='gene', columns='sample_id', aggfunc='size', fill_value=0)

# Display the pivot table
print("\nPivot Table:")
print(pivot_df)

# Optionally, save the pivot table to a CSV file
pivot_df.to_csv('plasmFinder_results_report.tsv', sep='\t', index=True)

print("\nPivot table saved to 'plasmFinder_results_report.tsv'.")
