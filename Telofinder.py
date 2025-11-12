#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Telomere detection script

- Reads a plain text or .gz FASTQ/sequence file (case-insensitive, converted to upper-case).
- Scans for degenrate telomere sequences of S. cerevisiae
- Writes a TSV with columns:
  Index, Telomere length, Position, Context
- Context handling:
  * We slice from (position - (length + 50)) to (position + 50).
  * We append "-END OF READ" ONLY if the raw slice ENDS with a newline.
  * Tabs are replaced by spaces to avoid TSV column shifts.
"""

import argparse
import gzip
import os

# -------------------------------
# CLI arguments
# -------------------------------
parser = argparse.ArgumentParser(description="Telomere detection script")
parser.add_argument(
    "--sequence", "-s",
    default="./Barcode.fastq",
    help="Path to the input sequence file (plain text or .gz). Default: ./Barcode.fastq"
)
parser.add_argument(
    "--min", "-m",
    type=int,
    default=50,
    help="Minimum telomere length required (default: 50, minimum 30 to exclude TG microsatellites)"
)
args = parser.parse_args()

sequence_file = args.sequence
min_length = args.min

# -------------------------------
# Initialize counters (will be reused per hit)
# -------------------------------
index = 0
length = 0
position = 0
trash = 0
ac = 0
gt = 0
tt = 0
gg = 0
context = ""
clean_context = ""
raw_context = ""

# -------------------------------
# Read input file as ONE string (upper-case to ignore case)
# Works for plain text and .gz
# -------------------------------
def read_all_upper(path: str) -> str:
    if path.endswith(".gz"):
        with gzip.open(path, "rt", newline="") as f:
            return f.read().upper()
    else:
        with open(path, "r", newline="") as f:
            return f.read().upper()

text = read_all_upper(sequence_file)

results = []

# -------------------------------
# Scan for telomeres
# -------------------------------
pos = 0
L = len(text)

while pos < L:
    # reset per-hit counters
    length = 0
    trash = 0
    ac = 0
    gt = 0
    tt = 0
    gg = 0
    position = pos
    while pos<L:
            ch = text[pos]

            if ch == "G":
                length += 1
                if text[pos-1] == "T":
                    gt += 1
                else:
                    gt = 0
                tt = 0
                gg += 1
            elif ch == "T":
                length += 1
                tt += 1
                gg = 0
            elif ch in ("A", "C"):
                ac = text.count("A",pos,pos+15) + text.count("C",pos,pos+15)
                if ac > 2:
                    break
                else:
                    ac=0
                    length +=1
            else:
                trash = 1
                break

            if tt > 3:
                break
            if gg > 5:
                break
            if gt > 10:
                trash = 1
                break
                
            pos += 1

    # Valid telomere?
    if length > min_length and length < 1000 and trash == 0 and text[position-1] != "\n" and text[position+length+1] != "\n":
        index += 1

        # Context window:
        start_context = max(0, position - 50)
        end_context = min(L, position + (length + 50))
        raw_context = text[start_context:end_context]

        # Replace tabs to keep TSV structure intact
        clean_context = raw_context.replace("\n", "END OF READ")
        context = clean_context.replace("\t", " ")

        results.append([index, length, position, context])
    pos += 1
    
pos = 0

while pos < L:
    # reset per-hit counters
    length = 0
    trash = 0
    ac = 0
    gt = 0
    tt = 0
    gg = 0
    position = pos
    while pos<L:
            ch = text[pos]

            if ch == "C":
                length += 1
                if text[pos-1] == "A":
                    gt += 1
                else:
                    gt = 0
                tt = 0
                gg += 1
            elif ch == "A":
                length += 1
                tt += 1
                gg = 0
            elif ch in ("T", "G"):
                ac = text.count("T",pos,pos+15) + text.count("G",pos,pos+15)
                if ac > 2:
                    break
                else:
                    ac=0
                    length +=1
            else:
                trash = 1
                break

            if tt > 3:
                break
            if gg > 5:
                break
            if gt > 10:
                trash = 1
                break
                
            pos += 1

    # Valid telomere?
    if length > min_length and length < 1000 and trash == 0 and text[position-1] != "\n" and text[position+length+1] != "\n":
        index += 1

        # Context window: from (position - (length + 50)) to (position + 50)
        start_context = max(0, position - 50)
        end_context = min(L, position + (length + 50))
        raw_context = text[start_context:end_context]

        # Replace tabs to keep TSV structure intact
        clean_context = raw_context.replace("\n", "END OF READ")
        context = clean_context.replace("\t", " ")

        results.append([index, length, position, context])
    pos += 1

# -------------------------------
# Prepare output path (handle .gz safely)
# -------------------------------
base = sequence_file
if base.endswith(".gz"):
    base = base[:-3]  # remove trailing '.gz' only
outfile = base + "-results.tsv"

# -------------------------------
# Write output TSV
# -------------------------------
with open(outfile, "w", newline="") as out:
    out.write("Index\tTelomere length\tPosition\tContext\n")
    for row in results:
        out.write("\t".join(map(str, row)) + "\n")

print(f"Analysis complete. {index} telomers found. Results written to: {outfile}")
