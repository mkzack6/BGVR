import csv
import random

# Define the number of genes and rows
num_genes = 10000
num_rows = 100

# Create header
header = [f"gene{i}" for i in range(1, num_genes + 1)] + ["label"]

# Generate data
data = []
for _ in range(num_rows):
    # Random gene expression values (0.0 to 1.0)
    row = [round(random.uniform(0.0, 1.0), 4) for _ in range(num_genes)]
    # Random binary label
    row.append(random.randint(0, 1))
    data.append(row)

# Write to CSV
with open("gene_expression.csv", "w", newline="") as f:
    writer = csv.writer(f)
    writer.writerow(header)
    writer.writerows(data)