from Bio import SeqIO, Entrez
import pandas as pd
import numpy as np
from sklearn.cluster import AgglomerativeClustering
import time

# NCBI email

Entrez.email = "samramustafa175@gmail.com"

# Input file

fasta_file = "protein.fasta"

# LOAD SEQUENCES

records = list(SeqIO.parse(fasta_file, "fasta"))

accessions = [r.id for r in records]
seqs = [str(r.seq) for r in records]

n = len(records)

# PROTEIN NAME FUNCTION

def clean_name(description, organism):

    desc = description.lower()

    if "groel" in desc:
        return f"GroEL ({organism})"
    elif "hsp60" in desc or "heat shock protein 60" in desc:
        return f"HSP60 ({organism})"
    elif "cpn60" in desc or "chaperonin" in desc:
        return f"Cpn60 ({organism})"
    else:
        return f"Protein ({organism})"


# Fetch annotation from ncbi

def get_ncbi_info(acc):

    try:
        handle = Entrez.efetch(
            db="protein",
            id=acc,
            rettype="gb",
            retmode="text"
        )

        record = SeqIO.read(handle, "genbank")

        description = record.description
        organism = record.annotations.get("organism", "Unknown")

        return description, organism

    except:
        return "Unknown", "Unknown"

print("Fetching annotation...")

protein_names = []
organisms = []

for acc in accessions:

    desc, org = get_ncbi_info(acc)

    name = clean_name(desc, org)

    protein_names.append(name)
    organisms.append(org)

    time.sleep(0.2)


# Pairwise identity

def sequence_identity(seq1, seq2):

    matches = 0
    length = len(seq1)

    for a, b in zip(seq1, seq2):
        if a == b:
            matches += 1

    return (matches / length) * 100


# Build matrix
matrix = np.zeros((n, n))

for i in range(n):
    for j in range(n):

        if i == j:
            matrix[i][j] = 100.0
        else:
            matrix[i][j] = sequence_identity(seqs[i], seqs[j])



matrix_df = pd.DataFrame(matrix, index=protein_names, columns=protein_names)

# SAVE MATRIX

matrix_df.to_csv("protein_similarity_matrix_named.csv")
matrix_df.to_excel("protein_similarity_matrix_named.xlsx", index=True)

print("Named similarity matrix saved")

# CLUSTERING

distance_matrix = 100 - matrix

clustering = AgglomerativeClustering(
    n_clusters=None,
    metric="precomputed",
    linkage="average",
    distance_threshold=25
)

labels = clustering.fit_predict(distance_matrix)

# FINAL TABLE

final_df = pd.DataFrame({
    "Accession": accessions,
    "Protein_Name": protein_names,
    "Organism": organisms,
    "Cluster_Group": labels
})


# SAVE OUTPUT

final_df.to_csv("ortholog_annotated_groups.csv", index=False)
final_df.to_excel("ortholog_annotated_groups.xlsx", index=False)

print("Ortholog groups saved")


# Print
print("\n=== SIMILARITY MATRIX ===")
print(matrix_df)

print("\n=== ORTHOLOG GROUPS ===")
print(final_df)
