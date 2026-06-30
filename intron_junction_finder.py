#!/usr/bin/env python3
"""
Intron Junction Finder
-----------------------
Identifies intron splice junctions from RNA-seq alignment data (SAM format)
by parsing CIGAR strings, then maps those junctions onto annotated gene
coordinates from a tab-separated gene location file.

Originally developed as a university coursework assignment (BIOL4292,
University of Glasgow), using only the Python standard library.

Usage:
    python intron_junction_finder.py alignment.sam genes.txt
"""

import sys


# ---------------------------------------------------------------------------
# FUNCTIONS
# ---------------------------------------------------------------------------

def parsing_cigar(cigar):
    """
    Parses a CIGAR string into a list of (length, operation) tuples.
    e.g. '50M40N2M' -> [(50, 'M'), (40, 'N'), (2, 'M')]
    """
    num = ''
    cigar = cigar.upper()  # normalise case
    parsed_cigars = []
    for char in cigar:
        if char.isdigit():
            num += char
        else:
            parsed_cigars.append((int(num), char))
            num = ''
    return parsed_cigars


assert parsing_cigar('100M') == [(100, 'M')]
assert parsing_cigar('50M40N2m') == [(50, 'M'), (40, 'N'), (2, 'M')]
assert parsing_cigar('50M2m') == [(50, 'M'), (2, 'M')]


def extract_introns_sam(sam_file):
    """
    Reads a SAM file and extracts intron junctions from uniquely-aligned,
    split reads (reads whose CIGAR string contains an 'N' operation).

    Returns a dict of {(chromosome, junction_start, junction_end): read_count}.
    """
    junction_counts = {}
    try:
        with open(sam_file) as sam:
            for line in sam:
                if line.startswith('@'):
                    continue  # skip header lines

                contents_sam = line.rstrip('\n').split('\t')
                chromosome = contents_sam[2]

                try:
                    align_start = int(contents_sam[3])
                except ValueError:
                    continue  # skip rows with a non-numeric position

                cigar = contents_sam[5]

                # NH:i:x is always in the last column for this assignment's
                # input files, so no extra parsing/validation is needed here.
                nread = contents_sam[-1]
                lastno = int(nread.split(':')[-1])
                if lastno != 1:
                    continue  # skip reads that aligned more than once

                if 'N' not in cigar:
                    continue  # skip reads with no splice junction

                try:
                    parsed_cigar = parsing_cigar(cigar)
                    ref_pos = align_start
                except ValueError:
                    continue  # skip CIGAR strings that fail to parse

                for bp, letter in parsed_cigar:
                    if letter == 'N':
                        start = ref_pos
                        end = ref_pos + bp
                        key = (chromosome, start, end)
                        # tally read support for this junction
                        junction_counts[key] = junction_counts.get(key, 0) + 1
                        ref_pos += bp
                    elif letter in ('M', 'D'):
                        # match and deletion both consume reference sequence
                        ref_pos += bp
                    else:
                        # insertions (I) and soft clips (S) do not consume
                        # reference sequence, so they are ignored here
                        pass

        return junction_counts  # {(chromosome, start, end): read_count}

    except FileNotFoundError:
        print("SAM file not found:", sam_file)
        sys.exit(1)


def parse_gene_locations(gene_file):
    """
    Reads a tab-separated gene location file and returns a dict of
    {gene_id: (chromosome, start, end)}.

    Expected location format: 'chrName:1,234..5,678(+)'
    """
    gene_data = {}

    with open(gene_file) as locations:
        next(locations)  # skip header row
        for line in locations:
            geneid, source_id, genomic_loc = line.rstrip().split('\t')

            try:
                chrom, location = genomic_loc.split(':')
            except ValueError:
                continue  # skip rows with an unexpected location format

            starts, ends = location.split('..')

            try:
                ends = ends[0:-3].replace(',', '')   # strip strand suffix and commas
                starts = starts.replace(',', '')
            except ValueError:
                continue

            gene_data[geneid] = (chrom, int(starts), int(ends))

    return gene_data


def match_introns(genes, junction_counts, output_name):
    """
    Matches junctions to genes that share a chromosome and fall fully
    within the gene's boundaries, then writes the results to a tab-
    separated output file (one blank line between genes).
    """
    output_filename = output_name + '.txt'
    with open(output_filename, 'w') as out:
        for gene_id, gchrom, gstart, gend in genes:
            for (jchrom, jstart, jend), count in junction_counts.items():
                if gchrom != jchrom:
                    continue  # chromosomes must match
                if gstart <= jstart and jend <= gend:
                    # junction lies within the gene boundaries
                    out.write(f"{gene_id}\t{jstart}\t{jend}\t{count}\n")
            out.write("\n")
    print("Wrote output to", output_filename)


# ---------------------------------------------------------------------------
# MAIN
# ---------------------------------------------------------------------------

def main():
    if len(sys.argv) != 3:
        sys.exit(
            "Please provide two input files at the command line\n"
            "The first file should be in SAM format\n"
            "The second file must be in TXT format"
        )

    sam_file = sys.argv[1]
    gene_file = sys.argv[2]

    if not sam_file.lower().endswith(".sam"):
        sys.exit("Error: input alignment file must end with .sam")

    gene_data = parse_gene_locations(gene_file)
    genes = [
        (geneid, chrom, start, end)
        for geneid, (chrom, start, end) in gene_data.items()
    ]

    junction_counts = extract_introns_sam(sam_file)

    output_name = "output"  # change to a fixed identifier if required
    match_introns(genes, junction_counts, output_name)


if __name__ == "__main__":
    main()
