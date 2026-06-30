# Intron Junction Finder

A Python tool that identifies intron splice junctions from RNA-seq alignment data (SAM format) and maps them onto annotated gene coordinates.

## Background

When mRNA is sequenced and aligned to a genomic reference, reads spanning an exon-exon boundary get split across the intron they were spliced out of. These "split reads" show up in the alignment's CIGAR string as a skipped region (`N`), and they're useful for confirming where introns actually sit — including catching mis-annotated genes or alternative splice isoforms.

This script parses CIGAR strings to find those split positions, counts how many reads support each junction, and reports which junctions fall within each annotated gene.

## What it does

1. Reads a SAM alignment file and, for each uniquely-aligned read (`NH:i:1`), parses the CIGAR string to find any skipped regions (`N`).
2. Converts each skipped region into a genomic junction (start/end coordinates), tallying how many reads support each unique junction.
3. Reads a tab-separated gene location file and parses out each gene's chromosome, start, and end coordinates.
4. Matches junctions to genes on the same chromosome where the junction falls fully within the gene's boundaries.
5. Writes the results to `output.txt`, with one row per junction (gene ID, junction start, junction end, supporting read count) and a blank line separating genes.

## Usage

```bash
python intron_junction_finder.py alignment.sam genes.txt
```

- `alignment.sam` — a SAM-format alignment file (header lines starting with `@` are skipped automatically)
- `genes.txt` — a tab-separated file with a header row and three columns: gene ID, transcript ID, and genomic location in the format `chrName:1,234,567..1,234,789(+)`

Output is written to `output.txt` in the working directory.

## Notes on the CIGAR parsing

- Only `M` (match) and `D` (deletion) operations advance the reference position normally.
- `N` (skipped region) marks an intron — its start/end becomes a junction.
- `I` (insertion) and `S` (soft clip) consume read sequence but **not** reference sequence, so they're excluded from junction position calculations.
- A single read can span multiple introns; each `N` region is extracted as a separate junction.

## Requirements

Python 3, standard library only — no external dependencies.

## Context

Originally written for a bioinformatics coursework assignment (BIOL4292, University of Glasgow) covering string parsing, file I/O, functions, and exception handling. Cleaned up afterward for general readability and to remove assignment-specific details (e.g. hardcoded filenames).
