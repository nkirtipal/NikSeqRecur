import os
import shutil
import time

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from matplotlib_venn import venn2
from Bio import SeqIO
from Bio.Seq import Seq


# -----------------------------
# Create Output Folder
# -----------------------------
OUTPUT_DIR = "example_outputs/sarscov2"

if os.path.exists(OUTPUT_DIR):
    shutil.rmtree(OUTPUT_DIR)

os.makedirs(OUTPUT_DIR)

print("Fresh output folder created.")


# -----------------------------
# Load FASTA
# -----------------------------
fasta_path = "example_data/sarscov2_spike.fasta"

record = next(
    SeqIO.parse(
        fasta_path,
        "fasta"
    )
)

sequence = str(record.seq)

# Spike gene region
sequence = sequence[21563:25384]

# Smaller test region
sequence = sequence[:1000]

print("Spike Gene Loaded")
print("Sequence Length:", len(sequence))
# -----------------------------
# Edit Distance Function
# -----------------------------
def edit_distance(s1, s2):

    m = len(s1)
    n = len(s2)

    dp = [[0 for _ in range(n + 1)] for _ in range(m + 1)]

    for i in range(m + 1):
        dp[i][0] = i

    for j in range(n + 1):
        dp[0][j] = j

    for i in range(1, m + 1):
        for j in range(1, n + 1):

            if s1[i - 1] == s2[j - 1]:
                cost = 0
            else:
                cost = 1

            dp[i][j] = min(
                dp[i - 1][j] + 1,
                dp[i][j - 1] + 1,
                dp[i - 1][j - 1] + cost
            )

    return dp[m][n]


# -----------------------------
# Similarity Function
# -----------------------------
def normalized_similarity(a, b):

    dist = edit_distance(a, b)

    similarity = 1 - (dist / max(len(a), len(b)))

    return similarity

# -----------------------------
# Sequence Entropy
# -----------------------------
def sequence_entropy(seq):

    bases = ['A', 'T', 'C', 'G']

    entropy = 0

    for b in bases:

        p = seq.count(b) / len(seq)

        if p > 0:
            entropy -= p * np.log2(p)

    return entropy

# -----------------------------
# GC Content
# -----------------------------
def gc_content(seq):

    gc = seq.count("G") + seq.count("C")

    return round((gc / len(seq)) * 100, 2)

    # -----------------------------
# Reverse Complement
# -----------------------------
def reverse_complement(seq):

    complement = {
        "A": "T",
        "T": "A",
        "C": "G",
        "G": "C"
    }

    revcomp = "".join(
        complement[b]
        for b in reversed(seq)
    )

    return revcomp

# -----------------------------
# Remove Duplicates
# -----------------------------
def remove_duplicates(repeats):

    unique = []

    seen = set()

    for r in repeats:

        key = (
            r['start1'],
            r['end1'],
            r['start2'],
            r['end2']
        )

        if key not in seen:

            unique.append(r)

            seen.add(key)

    return unique

# -----------------------------
# Similarity Filtering
# -----------------------------
def filter_similarity(repeats, threshold=0.95):

    filtered = []

    for r in repeats:

        if r['similarity'] >= threshold:
            filtered.append(r)

    return filtered

# -----------------------------
# Merge Overlapping Regions
# -----------------------------
def merge_overlaps(repeats):

    if len(repeats) == 0:
        return []

    repeats = sorted(repeats, key=lambda x: x['start1'])

    merged = [repeats[0]]

    for current in repeats[1:]:

        previous = merged[-1]

        if current['start1'] <= previous['end1']:

            previous['end1'] = max(
                previous['end1'],
                current['end1']
            )

            previous['similarity'] = max(
                previous['similarity'],
                current['similarity']
            )

        else:
            merged.append(current)

    return merged

    # -----------------------------
# Mutation Classification
# -----------------------------
def classify_mutation(prot1, prot2):

    if prot1 == "NA" or prot2 == "NA":
        return "noncoding"

    elif prot1 == prot2:
        return "synonymous"

    elif "*" in prot2 and "*" not in prot1:
        return "stop_gain"

    elif "*" in prot1 and "*" not in prot2:
        return "stop_loss"

    else:
        return "nonsynonymous"

# -----------------------------
# Transition / Transversion
# -----------------------------
def classify_substitution(seq1, seq2):

    transitions = {
        ("A", "G"),
        ("G", "A"),
        ("C", "T"),
        ("T", "C")
    }

    mutation_types = []

    for a, b in zip(seq1, seq2):

        if a != b:

            if (a, b) in transitions:
                mutation_types.append("transition")

            else:
                mutation_types.append("transversion")

    if len(mutation_types) == 0:
        return "none"

    return ",".join(sorted(set(mutation_types)))

# -----------------------------
# Codon Position Effect
# -----------------------------
def codon_position_effect(seq1, seq2):

    positions = []

    for i, (a, b) in enumerate(zip(seq1, seq2)):

        if a != b:

            codon_pos = (i % 3) + 1

            positions.append(str(codon_pos))

    if len(positions) == 0:
        return "none"

    return ",".join(sorted(set(positions)))

# -----------------------------
# Repeat Detection
# -----------------------------
def find_repeats(sequence,
                 min_len=6,
                 max_len=12,
                 similarity_threshold=0.8):

    repeats = []

    seq_len = len(sequence)

    for k in range(min_len, max_len + 1, 3):

        for i in range(seq_len - 2 * k):

            seq1 = sequence[i:i+k]

            for j in range(i + k, seq_len - k):

                seq2 = sequence[j:j+k]

                # Reverse complement
                rev_seq2 = reverse_complement(seq2)

                # Forward similarity
                forward_similarity = normalized_similarity(
                    seq1,
                    seq2
                )

                # Reverse-complement similarity
                reverse_similarity = normalized_similarity(
                    seq1,
                    rev_seq2
                )

                # Best similarity
                similarity = max(
                    forward_similarity,
                    reverse_similarity
                )

                # Match type
                match_type = (
                    "reverse_complement"
                    if reverse_similarity > forward_similarity
                    else "forward"
                )

                # Entropy
                entropy = sequence_entropy(seq1)

                if similarity >= similarity_threshold and entropy > 1.0:

                    # Protein Translation
                    protein1 = (
                        "-".join(str(Seq(seq1).translate()))
                        if len(seq1) % 3 == 0
                        else "NA"
                    )

                    protein2 = (
                        "-".join(str(Seq(seq2).translate()))
                        if len(seq2) % 3 == 0
                        else "NA"
                    )

                    # Store Repeat
                    repeats.append({

                        "start1": i,
                        "end1": i+k,

                        "start2": j,
                        "end2": j+k,

                        "length": k,

                        "similarity": round(similarity, 3),

                        "match_type": match_type,

                        "repeat1": seq1,
                        "repeat2": seq2,

                        "gc_content1": gc_content(seq1),
                        "gc_content2": gc_content(seq2),

                        "protein1": protein1,
                        "protein2": protein2,

                        "mutation_type": classify_mutation(
                            protein1,
                            protein2
                        ),
                        "substitution_type": classify_substitution(
                            seq1,
                            seq2
                        ),
                        "codon_position": codon_position_effect(
                           seq1,
                           seq2
                       )

                    })

    return repeats


# -----------------------------
# Run Tool
# -----------------------------
start = time.perf_counter()

repeats = find_repeats(
    sequence,
    min_len= 9,    #   6,
    max_len= 18,  # 12,
    similarity_threshold=0.8
)
# Save RAW results
raw_df = pd.DataFrame(repeats)

# Stop timer AFTER raw detection

end = time.perf_counter()

print("\nRaw Repeats Found:", len(repeats))

filtered_repeats = remove_duplicates(repeats)

filtered_repeats = filter_similarity(
    filtered_repeats,
    threshold=0.95
)

print("After Duplicate Filtering:", len(filtered_repeats))

filtered_df = pd.DataFrame(filtered_repeats)

pd.set_option('display.max_columns', None)
pd.set_option('display.width', None)

print(raw_df.head())

# -----------------------------
# Amino Acid Decoder
# -----------------------------
aa_dict = {
    "A": "Alanine",
    "R": "Arginine",
    "N": "Asparagine",
    "D": "Aspartic Acid",
    "C": "Cysteine",
    "Q": "Glutamine",
    "E": "Glutamic Acid",
    "G": "Glycine",
    "H": "Histidine",
    "I": "Isoleucine",
    "L": "Leucine",
    "K": "Lysine",
    "M": "Methionine",
    "F": "Phenylalanine",
    "P": "Proline",
    "S": "Serine",
    "T": "Threonine",
    "W": "Tryptophan",
    "Y": "Tyrosine",
    "V": "Valine"
}

# -----------------------------
# Amino Acid Legend
# -----------------------------
print("\nAmino Acid Legend:")

print("A = Alanine")
print("R = Arginine")
print("N = Asparagine")
print("D = Aspartic Acid")
print("C = Cysteine")
print("Q = Glutamine")
print("E = Glutamic Acid")
print("G = Glycine")
print("H = Histidine")
print("I = Isoleucine")
print("L = Leucine")
print("K = Lysine")
print("M = Methionine")
print("F = Phenylalanine")
print("P = Proline")
print("S = Serine")
print("T = Threonine")
print("W = Tryptophan")
print("Y = Tyrosine")
print("V = Valine")


# -----------------------------
# Save CSV
# -----------------------------
raw_df.to_csv(
    f"{OUTPUT_DIR}/raw_results.csv",
    index=False
)

filtered_df.to_csv(
    f"{OUTPUT_DIR}/filtered_results.csv",
    index=False
)

print("\nBoth CSV files saved.")

# -----------------------------
# BED Export
# -----------------------------
bed_df = filtered_df[['start1', 'end1']].copy()

bed_df.insert(0, 'chrom', 'sequence1')

bed_df.to_csv(
    f"{OUTPUT_DIR}/nseqrecur_regions.bed",
    sep='\t',
    header=False,
    index=False
)

print("BED file exported.")

# -----------------------------
# Plot Repeat Lengths
# -----------------------------
plt.figure(figsize=(8,5))

# plt.hist(df["length"], bins=10)
counts, bins, patches = plt.hist(
    raw_df["length"],
    bins=10
)

# Color bars using normalized counts
norm = plt.Normalize(counts.min(), counts.max())

for count, patch in zip(counts, patches):

    patch.set_facecolor(
        plt.cm.viridis(norm(count))
    )
plt.xlabel("Repeat Length (bp)")
plt.ylabel("Count")
plt.title(
    f"NSeqRecur Raw Recurrence Distribution\n"
    f"Sequence Length: {len(sequence)} bp | "
    f"Raw Repeats: {len(raw_df)}"
)
plt.subplots_adjust(bottom=0.25)
plt.figtext(
    0.15,
    0.05,
    f"Min Length: 9\n"    # 6\n"
    f"Max Length: 18\n"   # 12\n"
    f"Similarity Threshold: 0.8\n"
    f"Runtime: {round(end-start, 3)} sec",
    fontsize=10
)
plt.savefig(
    f"{OUTPUT_DIR}/Raw_Recurrence_distribution.png"
)
plt.show()

# -----------------------------
# FILTERED Distribution
# -----------------------------
plt.figure(figsize=(8,5))

# plt.hist(filtered_df["length"], bins=10)
counts, bins, patches = plt.hist(
    filtered_df["length"],
    bins=10
)

# Apply gradient colors
norm = plt.Normalize(counts.min(), counts.max())

for count, patch in zip(counts, patches):

    patch.set_facecolor(
        plt.cm.plasma(norm(count))
    )
plt.xlabel("Repeat Length (bp)")
plt.ylabel("Count")

plt.title(
    f"NSeqRecur Filtered Recurrence Distribution\n"
    f"Sequence Length: {len(sequence)} bp | "
    f"Filtered Repeats: {len(filtered_df)}"
)
plt.subplots_adjust(bottom=0.25)

plt.figtext(
    0.15,
    0.05,
    f"Similarity Filter: ≥ 0.95\n"
    f"Runtime: {round(end-start, 3)} sec",
    fontsize=10
)

plt.savefig(
    f"{OUTPUT_DIR}/filtered_distribution.png"
)

plt.show()
# -----------------------------
# Dotplot Visualization
# -----------------------------
#
# -----------------------------
# RAW Dotplot
# -----------------------------
plt.figure(figsize=(6,6))

#plt.scatter(
#    raw_df["start1"],
#    raw_df["start2"],
#    alpha=0.6
#)
plt.scatter(
    raw_df["start1"],
    raw_df["start2"],
    c=raw_df["similarity"],
    alpha=0.7
)

plt.colorbar(label="Similarity")

plt.xlabel("Region 1 Position")
plt.ylabel("Region 2 Position")

plt.title("NSeqRecur Raw Dotplot")

plt.savefig(
    f"{OUTPUT_DIR}/raw_dotplot.png"
)

plt.show()

# -----------------------------
# FILTERED Dotplot
# -----------------------------

plt.figure(figsize=(6,6))

plt.scatter(
    filtered_df["start1"],
    filtered_df["start2"],
    c=filtered_df["similarity"],
    cmap='plasma',
    vmin=0.95,
    vmax=1.0,
    alpha=0.7
)
plt.colorbar(label="Similarity")
plt.xlabel("Region 1 Position")
plt.ylabel("Region 2 Position")

plt.title("NSeqRecur Filtered Dotplot")

plt.savefig(
    f"{OUTPUT_DIR}/NSeqRecur_Filtered_dotplot.png"
)

plt.show()

# -----------------------------
# Venn Diagram
# -----------------------------

synonymous_set = set(
    filtered_df[
        filtered_df["mutation_type"] == "synonymous"
    ].index
)

transition_set = set(
    filtered_df[
        filtered_df["substitution_type"].str.contains(
            "transition",
            na=False
        )
    ].index
)

plt.figure(figsize=(6,6))

venn2(
    [synonymous_set, transition_set],
    set_labels=("Synonymous", "Transition")
)

plt.title("Mutation Overlap Venn Diagram")

plt.savefig(
    f"{OUTPUT_DIR}/mutation_venn.png"
)

plt.show()

# -----------------------------
# Runtime
# -----------------------------
print(f"\nExecution Time: {round(end-start, 3)} seconds")
