import sqlite3

# Connect to your SQLite database
conn = sqlite3.connect('Ecoli.db')
cursor = conn.cursor()

# Define and execute the SQL query
query = """
SELECT
    s.source,
    COUNT(DISTINCT i.int_id) AS num_unique_integrons,
    COUNT(i.int_id) AS total_integron_entries,
    GROUP_CONCAT(DISTINCT i.type_elt) AS unique_integron_types
FROM
    samples s
JOIN
    intFinder_results i ON s.isolate = i.sample_id
GROUP BY
    s.source
ORDER BY
    s.source;
"""

cursor.execute(query)
results = cursor.fetchall()

# Print the results
for row in results:
    print(row)

# Close the connection
conn.close()
