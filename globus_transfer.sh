#!/usr/bin/env bash
#
# globus_transfer.sh
#
# Downloads ISMIP7 Greenland SMB, thermal forcing, and subglacial discharge
# data for the core experiments from Globus:
#
#   Models:    CESM2-WACCM, MRI-ESM2-0 (and OCX)
#   Scenarios: historical (filtered to 2007-2014), ssp370, ssp126, ssp585
#
# Each model/scenario pulls four variables: tf, sgd, acabf-anomaly, dacabfdz.
# Directory naming and version numbers differ slightly by model (see
# get_relpaths() below):
#
#   CESM2-WACCM:
#     ocean-1000m/tf/v2
#     SDBN1-1000m/sgd/v2
#     SDBN1-1000m/acabf-anomaly/v3
#     SDBN1-1000m/dacabfdz/v3
#
#   MRI-ESM2-0:
#     ocean-1000m/tf/v1
#     GEMB-SDBN1-1000m/sgd/v1
#     GEMB-SDBN1-1000m/acabf-anomaly/v2
#     GEMB-SDBN1-1000m/dacabfdz/v2
#
#   OCX products (no scenario level, filtered to 2007-2025 by default;
#   RACMO2.3p2-ERA acabf (SMB) instead pulls 1958-2025):
#     OCX/EN4/ocean-1000m/tf/v1
#     OCX/RACMO2.3p2-ERA/SDBN1-1000m/acabf/v1
#     OCX/RACMO2.3p2-ERA/SDBN1-1000m/sgd/v1
#
# Also grabs one standalone file (see SINGLE_FILES below):
#     ISMIP7/tools/ismip7-gris-ocean-forcing/subglacial_discharge_basins_ismip.nc
#
# ---- Before running (future users) ----
#   1. Start Globus Connect Personal on this machine so it can act as the
#      transfer destination endpoint:
#         ./globusconnectpersonal -start &
#   2. Update DST_EP below to your own Globus Connect Personal endpoint
#      UUID (`globusconnectpersonal -show-config` or the web app will show
#      it), and update LOCAL_BASE to wherever you want the data to land.

set -euo pipefail

# ---- Endpoints ----
SRC_EP="ccc9bbd2-4091-4e35-addd-eeb639cf5332"
DST_EP="735d3da2-a79b-11f1-81e8-0ee7ef9370d9"  # <-- CHANGE ME: your Globus Connect Personal endpoint UUID

# ---- Base paths (model name is appended in the loop below) ----
REMOTE_BASE="ISMIP7/GrIS"
LOCAL_BASE="/data2/issm/shields/ismip7greenland/ModelData/ISMIP7/GrIS"  # <-- CHANGE ME: your local destination path

# ---- ISMIP7 root paths (for items outside GrIS/, e.g. shared tools) ----
ISMIP7_ROOT_REMOTE="ISMIP7"
ISMIP7_ROOT_LOCAL="/data2/issm/shields/ismip7greenland/ModelData/ISMIP7"

# ---- Standalone single files (no directory recursion/filtering) ----
SINGLE_FILES=(
    "tools/ismip7-gris-ocean-forcing/subglacial_discharge_basins_ismip.nc"
)

# ---- Models + scenarios to process ----
MODELS=("CESM2-WACCM" "MRI-ESM2-0")
SSP_SCENARIOS=("ssp370" "ssp126" "ssp585" "ctrl")

# ---- Year filter used for the "historical" pulls (2007-2014 only) ----
YEAR_REGEX='_(200[7-9]|201[0-4])\.nc$'

# ---- OCX products: no scenario level ----
# Default filter: 2007-2025. The RACMO SMB (acabf) product instead uses
# 1958-2025, via the per-entry override below.
OCX_YEAR_REGEX='_(200[7-9]|201[0-9]|202[0-5])\.nc$'
OCX_SMB_YEAR_REGEX='_(195[89]|19[6-9][0-9]|200[0-9]|201[0-9]|202[0-5])\.nc$'

# Each entry: product|rel_path|year_regex_override
# (leave the third field blank to fall back to OCX_YEAR_REGEX)
OCX_PATHS=(
    "EN4|ocean-1000m/tf/v1|"
    "RACMO2.3p2-ERA|SDBN1-1000m/acabf/v1|${OCX_SMB_YEAR_REGEX}"
    "RACMO2.3p2-ERA|SDBN1-1000m/sgd/v1|"
)

# Scratch files for batch transfers (reused/overwritten each iteration)
MATCHED_FILE="$(mktemp)"
BATCH_FILE="$(mktemp)"
trap 'rm -f "$MATCHED_FILE" "$BATCH_FILE"' EXIT

# ---------------------------------------------------------------
# Per-model variable/version layout: tf, sgd, acabf-anomaly
# ---------------------------------------------------------------
get_relpaths() {
    local model="$1"
    case "$model" in
        CESM2-WACCM)
            printf '%s\n' \
                "ocean-1000m/tf/v2" \
                "SDBN1-1000m/sgd/v2" \
                "SDBN1-1000m/acabf-anomaly/v3" \
                "SDBN1-1000m/dacabfdz/v3"
            ;;
        MRI-ESM2-0)
            printf '%s\n' \
                "ocean-1000m/tf/v1" \
                "GEMB-SDBN1-1000m/sgd/v1" \
                "GEMB-SDBN1-1000m/acabf-anomaly/v2" \
                "GEMB-SDBN1-1000m/dacabfdz/v2"
            ;;
        *)
            echo "Unknown model: $model" >&2
            return 1
            ;;
    esac
}

# ---------------------------------------------------------------
# Filtered (year-restricted) transfer, used for "historical"
# ---------------------------------------------------------------
transfer_filtered() {
    local model="$1"
    local rel_path="$2"   # e.g. ocean-1000m/tf/v2
    local remote_dir="${REMOTE_BASE}/${model}/historical/${rel_path}/"
    local local_dir="${LOCAL_BASE}/${model}/historical/${rel_path}/"

    echo "==> [${model}] [historical] ${rel_path}"
    mkdir -p "$local_dir"

    globus ls "${SRC_EP}:${remote_dir}" \
        | grep -E "$YEAR_REGEX" > "$MATCHED_FILE" || true

    if [[ ! -s "$MATCHED_FILE" ]]; then
        echo "    no matching files found for ${model}/${rel_path}, skipping"
        return
    fi

    while IFS= read -r f; do
        echo "$f $f"
    done < "$MATCHED_FILE" > "$BATCH_FILE"

    globus transfer \
        "${SRC_EP}:${remote_dir}" \
        "${DST_EP}:${local_dir}" \
        --batch "$BATCH_FILE"
}

# ---------------------------------------------------------------
# Full recursive transfer, used for the ssp scenarios
# ---------------------------------------------------------------
transfer_recursive() {
    local model="$1"
    local scenario="$2"   # e.g. ssp370, ssp126, ssp585
    local rel_path="$3"   # e.g. ocean-1000m/tf/v2
    local remote_dir="${REMOTE_BASE}/${model}/${scenario}/${rel_path}/"
    local local_dir="${LOCAL_BASE}/${model}/${scenario}/${rel_path}/"

    echo "==> [${model}] [${scenario}] ${rel_path}"
    mkdir -p "$local_dir"

    globus transfer -r \
        "${SRC_EP}:${remote_dir}" \
        "${DST_EP}:${local_dir}"
}

# ---------------------------------------------------------------
# Filtered transfer for OCX products (EN4, RACMO2.3p2-ERA).
# These have no scenario level and use a wider year range (2007-2025).
# ---------------------------------------------------------------
transfer_ocx() {
    local product="$1"    # e.g. EN4, RACMO2.3p2-ERA
    local rel_path="$2"   # e.g. ocean-1000m/tf/v1
    local year_regex="$3" # per-entry year filter (falls back to OCX_YEAR_REGEX)
    local remote_dir="${REMOTE_BASE}/OCX/${product}/${rel_path}/"
    local local_dir="${LOCAL_BASE}/OCX/${product}/${rel_path}/"

    echo "==> [OCX] [${product}] ${rel_path} (years: ${year_regex})"
    mkdir -p "$local_dir"

    globus ls "${SRC_EP}:${remote_dir}" \
        | grep -E "$year_regex" > "$MATCHED_FILE" || true

    if [[ ! -s "$MATCHED_FILE" ]]; then
        echo "    no matching files found for OCX/${product}/${rel_path}, skipping"
        return
    fi

    while IFS= read -r f; do
        echo "$f $f"
    done < "$MATCHED_FILE" > "$BATCH_FILE"

    globus transfer \
        "${SRC_EP}:${remote_dir}" \
        "${DST_EP}:${local_dir}" \
        --batch "$BATCH_FILE"
}

# ---------------------------------------------------------------
# Single-file transfer (no listing/filtering, just one file).
# Destination directory is created first since globus transfer
# won't create intermediate dirs on its own.
# ---------------------------------------------------------------
transfer_single_file() {
    local rel_file="$1"   # path relative to ISMIP7_ROOT_REMOTE / ISMIP7_ROOT_LOCAL
    local remote_file="${ISMIP7_ROOT_REMOTE}/${rel_file}"
    local local_file="${ISMIP7_ROOT_LOCAL}/${rel_file}"
    local local_dir
    local_dir="$(dirname "$local_file")"

    echo "==> [single file] ${rel_file}"
    mkdir -p "$local_dir"

    globus transfer \
        "${SRC_EP}:${remote_file}" \
        "${DST_EP}:${local_file}"
}

# ---------------------------------------------------------------
# Main: for each model, run historical (filtered) + all ssp
# scenarios (recursive) across all three variables
# ---------------------------------------------------------------
for model in "${MODELS[@]}"; do
    echo "===================================================="
    echo " Model: ${model}"
    echo "===================================================="

    # Make sure the top-level model directory exists even before
    # any sub-transfers run.
    mkdir -p "${LOCAL_BASE}/${model}"

    # historical: filtered to 2007-2014
    while IFS= read -r rel_path; do
        transfer_filtered "$model" "$rel_path"
    done < <(get_relpaths "$model")

    # ssp scenarios: full recursive pull
    for scenario in "${SSP_SCENARIOS[@]}"; do
        while IFS= read -r rel_path; do
            transfer_recursive "$model" "$scenario" "$rel_path"
        done < <(get_relpaths "$model")
    done
done

echo "===================================================="
echo " OCX products"
echo "===================================================="
mkdir -p "${LOCAL_BASE}/OCX"

for entry in "${OCX_PATHS[@]}"; do
    product="${entry%%|*}"
    rest="${entry#*|}"
    rel_path="${rest%%|*}"
    year_override="${rest#*|}"
    transfer_ocx "$product" "$rel_path" "${year_override:-$OCX_YEAR_REGEX}"
done

echo "===================================================="
echo " Standalone files"
echo "===================================================="
for rel_file in "${SINGLE_FILES[@]}"; do
    transfer_single_file "$rel_file"
done

echo "All transfers submitted. Check 'globus task list' for status."
