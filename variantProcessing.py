# Process variant data

import os
import pandas as pd
from itertools import combinations

# -------------------------------------------------------------------
# Set working directory
os.chdir("C:/Users/acale/OneDrive/Documents/Waterloo BME/4A/BME 499/BME499/Data/")

# Read in the data
cDifficile_variants = pd.read_csv('mutation_analysis.csv', sep=',', header=0)

# -------------------------------------------------------------------
# Define the mutation columns
mutation_columns = ['gyrA_mutation', 'rpoB_mutation', 'nimB_mutation', 'vanR_mutation']

# Function to check if a cell has a mutation
def has_mutation(cell):
    return cell != 'N'

# Function to check if a row has no mutations
def no_mutations(row):
    return all(row[col] == 'N' for col in mutation_columns)

# Count individual mutations
individual_counts = {col: cDifficile_variants[col].apply(has_mutation).sum() for col in mutation_columns}

# Count combinations of mutations
pair_counts = {}
triplet_counts = {}

for pair in combinations(mutation_columns, 2):
    pair_counts[pair] = cDifficile_variants.apply(lambda row: all(has_mutation(row[col]) for col in pair), axis=1).sum()

for triplet in combinations(mutation_columns, 3):
    triplet_counts[triplet] = cDifficile_variants.apply(lambda row: all(has_mutation(row[col]) for col in triplet), axis=1).sum()

# Count samples with no mutations
no_mutation_count = cDifficile_variants.apply(no_mutations, axis=1).sum()

# -------------------------------------------------------------------
# Output the results
print("Individual Mutation Counts:")
for col, count in individual_counts.items():
    print(f"{col}: {count}")

print("\nPairwise Mutation Counts:")
for pair, count in pair_counts.items():
    print(f"{pair}: {count}")

print("\nTriplet Mutation Counts:")
for triplet, count in triplet_counts.items():
    print(f"{triplet}: {count}")

print("\nNumber of Samples with No Mutations:", no_mutation_count)

# Print sample rows where no_mutations returns True
sample_no_mutations = cDifficile_variants[cDifficile_variants.apply(no_mutations, axis=1)]
print(sample_no_mutations.head())