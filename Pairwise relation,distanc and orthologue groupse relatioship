from Bio import SeqIO, Entrez
import pandas as pd
import numpy as np
from sklearn.cluster import AgglomerativeClustering
import time

# Ncbi email

Entrez.email = "samramustafa175@gmail.com"

# Input file

fasta_file = "protein.fasta"


# Load sequences

records = list(SeqIO.parse(fasta_file, "fasta"))

accessions = [r.id for r in records]
seqs = [str(r.seq) for r in records]

n = len(records)


# CLEAN NAME FUNCTION

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


# FETCH NCBI INFO

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


# PAIRWISE IDENTITY

def sequence_identity(seq1, seq2):

    matches = sum(a == b for a, b in zip(seq1, seq2))
    return (matches / len(seq1)) * 100

# SIMILARITY MATRIX

matrix = np.zeros((n, n))

for i in range(n):
    for j in range(n):

        if i == j:
            matrix[i][j] = 100.0
        else:
            matrix[i][j] = sequence_identity(seqs[i], seqs[j])


# DISTANCE MATRIX (MAIN FIX)

distance_matrix = 100 - matrix

distance_df = pd.DataFrame(
    distance_matrix,
    index=protein_names,
    columns=protein_names
)


# SAVE DISTANCE MATRIX

distance_df.to_csv("protein_distance_matrix.csv")
distance_df.to_excel("protein_distance_matrix.xlsx")

print("Distance matrix saved")


# CLUSTERING (ORTHOLOG + ALLELE SEPARATION)

clustering = AgglomerativeClustering(
    n_clusters=None,
    metric="precomputed",
    linkage="average",
    distance_threshold=25   # adjust if needed
)

labels = clustering.fit_predict(distance_matrix)

# CLASSIFICATION FUNCTION

def classify_relationship(sim):

    if sim >= 98:
        return "Allele / very close variant"
    elif sim >= 70:
        return "Ortholog (same gene, different species)"
    else:
        return "Distant homolog"

relationship = [
    classify_relationship(100 - distance_matrix[i][i])
    for i in range(n)
]


# FINAL TABLE

final_df = pd.DataFrame({
    "Accession": accessions,
    "Protein_Name": protein_names,
    "Organism": organisms,
    "Cluster_Group": labels
})

# Add relationship classification
final_df["Relationship_Type"] = [
    classify_relationship(100) for _ in range(n)
]


pairwise_relationships = []

for i in range(n):
    for j in range(i+1, n):

        sim = matrix[i][j]

        pairwise_relationships.append({
            "Protein_A": protein_names[i],
            "Protein_B": protein_names[j],
            "Similarity_%": sim,
            "Distance_%": 100 - sim,
            "Relationship": classify_relationship(sim)
        })

pair_df = pd.DataFrame(pairwise_relationships)


# SAVE ALL OUTPUTS

distance_df.to_csv("protein_distance_matrix.csv")
distance_df.to_excel("protein_distance_matrix.xlsx")

final_df.to_csv("ortholog_annotated_groups.csv", index=False)
final_df.to_excel("ortholog_annotated_groups.xlsx", index=False)

pair_df.to_csv("pairwise_relationships.csv", index=False)
pair_df.to_excel("pairwise_relationships.xlsx", index=False)


# PRINT

print("\n=== DISTANCE MATRIX ===")
print(distance_df)

print("\n=== CLUSTER TABLE ===")
print(final_df)

print("\n=== PAIRWISE RELATIONSHIPS ===")
print(pair_df)
