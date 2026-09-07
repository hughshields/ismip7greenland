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
# Each model/scenario pulls three variables: tf, sgd, acabf-anomaly.
# Directory naming and version numbers differ slightly by model (see
# get_relpaths() below):
#
#   CESM2-WACCM:
#     ocean-1000m/tf/v2
#     SDBN1-1000m/sgd/v2
#     SDBN1-1000m/acabf-anomaly/v3 
#
#   MRI-ESM2-0:
#     ocean-1000m/tf/v1
#     GEMB-SDBN1-1000m/sgd/v1
#     GEMB-SDBN1-1000m/acabf-anomaly/v2
#
#   OCX products (no scenario level, filtered to 2007-2025 instead):
#     OCX/EN4/ocean-1000m/tf/v1
#     OCX/RACMO2.3p2-ERA/SDBN1-1000m/acabf/v1
#     OCX/RACMO2.3p2-ERA/SDBN1-1000m/sgd/v1

set -euo pipefail

# ---- Endpoints ----
SRC_EP="ccc9bbd2-4091-4e35-addd-eeb639cf5332"
DST_EP="735d3da2-a79b-11f1-81e8-0ee7ef9370d9"

# ---- Base paths (model name is appended in the loop below) ----
REMOTE_BASE="ISMIP7/GrIS"
LOCAL_BASE="/data2/issm/shields/ismip7greenland/ModelData/ISMIP7/GrIS"

# ---- Models + scenarios to process ----
MODELS=("CESM2-WACCM" "MRI-ESM2-0")
SSP_SCENARIOS=("ssp370" "ssp126" "ssp585" "ctrl")

# ---- Year filter used for the "historical" pulls (2007-2014 only) ----
YEAR_REGEX='_(200[7-9]|201[0-4])\.nc$'

# ---- OCX products: no scenario level, filtered to 2007-2025 ----
OCX_YEAR_REGEX='_(200[7-9]|201[0-9]|202[0-5])\.nc$'
OCX_PATHS=(
    "EN4|ocean-1000m/tf/v1"
    "RACMO2.3p2-ERA|SDBN1-1000m/acabf/v1"
    "RACMO2.3p2-ERA|SDBN1-1000m/sgd/v1"
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
                "SDBN1-1000m/acabf-anomaly/v3"
            ;;
        MRI-ESM2-0)
            printf '%s\n' \
                "ocean-1000m/tf/v1" \
                "GEMB-SDBN1-1000m/sgd/v1" \
                "GEMB-SDBN1-1000m/acabf-anomaly/v2"
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
    local remote_dir="${REMOTE_BASE}/OCX/${product}/${rel_path}/"
    local local_dir="${LOCAL_BASE}/OCX/${product}/${rel_path}/"

    echo "==> [OCX] [${product}] ${rel_path}"
    mkdir -p "$local_dir"

    globus ls "${SRC_EP}:${remote_dir}" \
        | grep -E "$OCX_YEAR_REGEX" > "$MATCHED_FILE" || true

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
    rel_path="${entry#*|}"
    transfer_ocx "$product" "$rel_path"
done

echo "All transfers submitted. Check 'globus task list' for status."
