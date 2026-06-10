#!/bin/bash

# ==============================================================================
# PIPELINE REQUIREMENTS:
# 1. The script requires two input directories passed as arguments:
#    - Ancestor directory containing the ancestor pileup files.
#    - Endpoint directory containing the endpoint pileup files.
# 2. CRITICAL SORTING REQUIREMENT: Input pileup files MUST be strictly sorted
#    by chromosome, and within each chromosome, sorted by base position.
#    Otherwise, the intersection logic in the 'common' step will completely fail.
#
# FILE NAMING CONVENTION:
# Input files must already be split by Copy Number (CN) and named exactly as:
# - Ancestor: XXX_CN_ancestor_pileup.gz
# - Endpoint: XXX_CN_endpoint_pileup.gz
# (XXXX is the sample name/ID, CN is the Copy Number, e.g., SampleA_2_ancestor_pileup.gz)
#
# FOLDER STRUCTURE:
# For each XXX_CN pair, a dedicated directory named 'XXX_CN' will be created.
# All uncompressed output files (.dat) and subsequent pipeline steps will
# be saved inside this specific pair folder.
# ==============================================================================

ANCESTOR_DIR=$1
ENDPOINT_DIR=$2

# Parameters
MAX_COVERAGE_THRESHOLD=1000  # <-- Change this threshold as needed

# Executables paths
PARSER_EXEC="./parsing/exec/parsing"
COMMON_EXEC="./common/exec/common"
MEAN_COV_EXEC="./mean_coverage/exec/mean_coverage"
MODAL_COV_EXEC="./modal_coverage/exec/modal_coverage"
LIST_COORD_EXEC="./signatures/exec/list_coord_mut"

# Basic check for input directories
if [ -z "$ANCESTOR_DIR" ] || [ -z "$ENDPOINT_DIR" ]; then
    echo "Usage: $0 <ancestor_directory> <endpoint_directory>"
    exit 1
fi

# --- 1. PROCESS ANCESTOR FILES ---
echo "Starting parsing for ANCESTOR files..."

for anc_file in "$ANCESTOR_DIR"/*_ancestor_pileup.gz; do
    # Extract just the filename from the path
    filename=$(basename "$anc_file")

    # Extract the XXX_CN prefix by removing the suffix
    prefix="${filename%_ancestor_pileup.gz}"

    # Create the dedicated folder for this pair
    mkdir -p "$prefix"

    # Define output file inside the new folder
    out_dat="$prefix/${prefix}_ancestor_pseudopileup.dat"

    echo "Processing $filename -> $out_dat"
    $PARSER_EXEC "$anc_file" "$out_dat"
done

# --- 2. PROCESS ENDPOINT FILES ---
echo "Starting parsing for ENDPOINT files..."

for end_file in "$ENDPOINT_DIR"/*_endpoint_pileup.gz; do
    filename=$(basename "$end_file")
    prefix="${filename%_endpoint_pileup.gz}"

    # Ensure folder exists (it should already be created by the ancestor loop)
    mkdir -p "$prefix"

    out_dat="$prefix/${prefix}_endpoint_pseudopileup.dat"

    echo "Processing $filename -> $out_dat"
    $PARSER_EXEC "$end_file" "$out_dat"
done

# --- 3. EXTRACT COMMON LINES, STATS AND MUTATION FREQUENCIES ---
echo "Starting common lines extraction, statistics, and mutations..."

for anc_dat in */*_ancestor_pseudopileup.dat; do

    # Extract the pair folder name (e.g., XXX_CN)
    pair_dir=$(dirname "$anc_dat")
    prefix=$(basename "$pair_dir")

    # Extract Copy Number (CN) by keeping only what is after the last underscore
    CN="${prefix##*_}"

    # Inputs
    end_dat="$pair_dir/${prefix}_endpoint_pseudopileup.dat"

    # Subdirectories
    mkdir -p "$pair_dir/common"
    mkdir -p "$pair_dir/stats"
    mkdir -p "$pair_dir/freq"

    out_common_anc="$pair_dir/common/common_lines_ancestor.dat"
    out_common_end="$pair_dir/common/common_lines_endpoint.dat"

    out_stats_mean_anc="$pair_dir/stats/ancestor_mean_coverage.dat"
    out_stats_mean_end="$pair_dir/stats/endpoint_mean_coverage.dat"

    out_stats_modal_anc="$pair_dir/stats/ancestor_modal_coverage.dat"
    out_stats_modal_end="$pair_dir/stats/endpoint_modal_coverage.dat"

    # Execute Common Lines
    echo "Extracting common lines for $prefix..."
    $COMMON_EXEC "$anc_dat" "$end_dat" "$out_common_anc" "$out_common_end"

    # Execute Mean Coverage
    echo "Calculating mean coverage for $prefix..."
    $MEAN_COV_EXEC "$out_common_anc" "$out_stats_mean_anc" "$MAX_COVERAGE_THRESHOLD"
    $MEAN_COV_EXEC "$out_common_end" "$out_stats_mean_end" "$MAX_COVERAGE_THRESHOLD"

    # Execute Modal Coverage
    echo "Calculating modal coverage for $prefix..."
    $MODAL_COV_EXEC "$out_common_anc" "$out_stats_modal_anc"
    $MODAL_COV_EXEC "$out_common_end" "$out_stats_modal_end"

    # --- 4. LIST COORDINATE MUTATIONS ---
    echo "Extracting mutation frequencies for $prefix (CN: $CN)..."

    { read anc_mean; read anc_std_dev; } < "$out_stats_mean_anc"
    { read end_mean; read end_std_dev; } < "$out_stats_mean_end"

    { read anc_modal; } < "$out_stats_modal_anc"
    { read end_modal; } < "$out_stats_modal_end"

    $LIST_COORD_EXEC "$out_common_anc" \
                     "$out_common_end" \
                     "$pair_dir/freq/" \
                     "$anc_modal" \
                     "$anc_std_dev" \
                     "$end_modal" \
                     "$end_std_dev" \
                     "$CN"

done
