#!/usr/bin/env python3
# Original script by K. Kortright; refactored with Claude (Anthropic, 2026).
# Changes: PEP-8 formatting, docstrings, glob-based file discovery replacing
# 42 hardcoded calls, and bug fixes for misassigned output variables (files 14–28).
"""
vcf_parse.py

Parses GATK-produced VCF files from an P. aeruginosa evolution experiment,
extracting per-site allele frequencies and annotating variants with
gene information from a reference genome dictionary.

Output columns (tab-separated):
    chrom, pos, ref_bp, freq_ref, alt_bp, freq_alt, gene_number, gene_name, qual

Usage:
    python vcf_parse.py
    (Expects genomedict.txt and all target VCF files in the working directory.)
"""

import glob
import os
import sys


# ---------------------------------------------------------------------------
# Load genome annotation dictionary
# ---------------------------------------------------------------------------

def load_genome_dict(dict_path: str) -> dict:
    """
    Load a tab-delimited genome annotation file mapping genomic positions
    to gene annotations.

    Expected file format (one entry per line):
        position <tab> field1 <tab> gene_number <tab> gene_name

    Args:
        dict_path: Path to the genome dictionary file.

    Returns:
        A dict mapping position strings to [field1, gene_number, gene_name].
    """
    genome_dict = {}
    with open(dict_path, "r") as fh:
        for line in fh:
            fields = line.strip("\n").split("\t")
            # Key: genomic position; value: [field1, gene_number, gene_name]
            genome_dict[fields[0]] = [fields[1], fields[2], fields[3]]
    return genome_dict


# ---------------------------------------------------------------------------
# VCF parser
# ---------------------------------------------------------------------------

def parse_vcf(input_path: str, output_path: str, genome_dict: dict) -> None:
    """
    Parse a single GATK VCF file and write a tab-delimited summary of
    high-quality variant sites with allele frequencies.

    For each non-header VCF line, the function:
      - Skips LowQual sites.
      - Requires the FORMAT field to follow the AD:DP layout used by GATK
        HaplotypeCaller in haploid/pooled mode (index 1 = AD, index 2 = DP).
      - Computes reference and alternate allele frequencies from AD counts.
      - Annotates the site with gene information from genome_dict, or marks
        it as intergenic if the position is absent from the dictionary.

    Args:
        input_path:   Path to the input VCF file.
        output_path:  Path for the parsed output text file.
        genome_dict:  Annotation dictionary from load_genome_dict().
    """
    with open(input_path, "r") as infile, open(output_path, "w") as outfile:
        for line in infile:

            # Skip VCF header lines
            if line.startswith("#"):
                continue

            fields = line.strip("\n").split("\t")
            chrom   = fields[0]
            pos     = fields[1]
            ref_bp  = fields[3]
            alt_bp  = fields[4]
            qual    = fields[5]
            filter_ = fields[6]   # FILTER column ("PASS", "LowQual", etc.)

            # Drop low-quality calls before any further processing
            if filter_ == "LowQual":
                continue

            # Parse the FORMAT field to locate AD and DP subfields.
            # GATK emits FORMAT as e.g. "GT:AD:DP:..." — we require
            # index 1 == "AD" and index 2 == "DP" to match expected layout.
            format_keys = fields[8].split(":")
            if not (len(format_keys) >= 3
                    and format_keys[1] == "AD"
                    and format_keys[2] == "DP"):
                continue  # Unexpected FORMAT layout; skip site

            # Extract AD (allelic depth) and DP (total read depth) from
            # the sample column (index 9), matching FORMAT order above.
            sample_values = fields[9].split(":")
            ad_field = sample_values[1]   # "ref_count,alt_count"
            # dp_field = sample_values[2]  # total depth (unused downstream)

            ref_count = int(ad_field.split(",")[0])
            alt_count = int(ad_field.split(",")[1])
            total     = ref_count + alt_count

            # Compute allele frequencies (guard against zero-depth sites)
            freq_ref = ref_count / total if total > 0 else 0.0
            freq_alt = alt_count / total if total > 0 else 0.0

            # Look up gene annotation; fall back to "intergenic" if absent
            if pos in genome_dict:
                gene_number = genome_dict[pos][1]
                gene_name   = genome_dict[pos][2]
            else:
                gene_number = "intergenic"
                gene_name   = "intergenic"

            # Write parsed record to output file
            outfile.write(
                "\t".join([
                    chrom,
                    pos,
                    ref_bp,
                    str(freq_ref),
                    alt_bp,
                    str(freq_alt),
                    gene_number,
                    gene_name,
                    qual,
                ]) + "\n"
            )

    print(f"Done with {input_path}")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    """
    Discover all VCF files in the working directory and parse each one.

    VCF files are expected to be named:
        bacteria_<sample>_raw_variants.vcf

    Each produces a corresponding:
        bacteria_<sample>_parsed.txt
    """
    # Load genome position-to-gene annotation dictionary
    genome_dict = load_genome_dict("genomedict.txt")

    # Discover all target VCF files rather than hardcoding each name
    vcf_files = sorted(glob.glob("bacteria_*_raw_variants.vcf"))

    if not vcf_files:
        print("No VCF files matching 'bacteria_*_raw_variants.vcf' found.",
              file=sys.stderr)
        sys.exit(1)

    for vcf_path in vcf_files:
        # Derive output path by replacing the suffix
        out_path = vcf_path.replace("_raw_variants.vcf", "_parsed.txt")
        parse_vcf(vcf_path, out_path, genome_dict)


if __name__ == "__main__":
    main()
