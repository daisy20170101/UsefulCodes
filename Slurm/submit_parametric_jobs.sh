#!/bin/bash
#
# Submit parametric SLURM jobs with different mu_s values
# Creates model directories, modifies parameters, and submits jobs
#

# ============================================================================
# Configuration
# ============================================================================

BASE_DIR="/Users/DuoL/Documents/SeisSol/slurm/tibet"
OUTPUT_BASE="/nesi/nobackup/gns04005/daisy/output/tibet/t16mod7PL"

# Define model numbers and corresponding mu_s values
# Add or modify these arrays as needed
MODEL_NUMBERS=(1 2 3 4 5)
MU_S_VALUES=(0.5 0.6 0.7 0.8 0.9)

# Whether to submit jobs (1=yes, 0=no, just create directories)
SUBMIT_JOBS=1

# ============================================================================
# Functions
# ============================================================================

create_model() {
    local model_num=$1
    local mu_s_value=$2

    local model_name=$(printf "model_%02d" $model_num)
    local model_dir="${BASE_DIR}/${model_name}"
    local output_path="${OUTPUT_BASE}_${model_name}"

    echo "========================================================================"
    echo "Creating ${model_name}"
    echo "========================================================================"

    # Create model directory
    if [ -d "$model_dir" ]; then
        echo "  WARNING: ${model_dir} already exists!"
        read -p "  Overwrite? (y/n): " -n 1 -r
        echo
        if [[ ! $REPLY =~ ^[Yy]$ ]]; then
            echo "  Skipping ${model_name}"
            return 1
        fi
        rm -rf "$model_dir"
    fi

    mkdir -p "$model_dir"
    echo "  Created directory: ${model_dir}"

    # Copy required files
    echo "  Copying files..."
    cp "${BASE_DIR}/parameters.par" "$model_dir/" || { echo "  ERROR: Failed to copy parameters.par"; return 1; }
    cp "${BASE_DIR}/submit_att.sh" "$model_dir/" || { echo "  ERROR: Failed to copy submit_att.sh"; return 1; }
    cp "${BASE_DIR}/fault2.yaml" "$model_dir/" || { echo "  ERROR: Failed to copy fault2.yaml"; return 1; }
    echo "  ✓ Copied: parameters.par, submit_att.sh, fault2.yaml"

    # Modify OutputFile in parameters.par
    echo "  Modifying parameters.par..."
    sed -i.bak "s|OutputFile = '.*'|OutputFile = '${output_path}'|" "${model_dir}/parameters.par"
    rm -f "${model_dir}/parameters.par.bak"
    echo "    OutputFile: ${output_path}"

    # Modify mu_s in fault2.yaml
    echo "  Modifying fault2.yaml..."
    # This modifies the first occurrence of mu_s=<number>
    sed -i.bak "s/^\(\s*\)mu_s\s*=\s*[0-9.]\+/\1mu_s=${mu_s_value}/" "${model_dir}/fault2.yaml"
    rm -f "${model_dir}/fault2.yaml.bak"
    echo "    mu_s: ${mu_s_value}"

    # Submit job
    if [ $SUBMIT_JOBS -eq 1 ]; then
        echo "  Submitting job..."
        cd "$model_dir" || return 1

        job_output=$(sbatch submit_att.sh 2>&1)

        if [ $? -eq 0 ]; then
            echo "  ✓ Job submitted: ${job_output}"
        else
            echo "  ✗ Job submission failed: ${job_output}"
        fi

        cd "$BASE_DIR"
    else
        echo "  (Skipping job submission)"
    fi

    echo ""
    return 0
}

# ============================================================================
# Main script
# ============================================================================

echo "========================================================================"
echo " PARAMETRIC JOB SUBMISSION"
echo "========================================================================"
echo "Base directory: ${BASE_DIR}"
echo "Number of jobs: ${#MODEL_NUMBERS[@]}"
echo "Submit jobs: ${SUBMIT_JOBS}"
echo ""

# Check if arrays have same length
if [ ${#MODEL_NUMBERS[@]} -ne ${#MU_S_VALUES[@]} ]; then
    echo "ERROR: MODEL_NUMBERS and MU_S_VALUES arrays must have same length!"
    exit 1
fi

# Create each model
created_count=0
failed_count=0

for i in "${!MODEL_NUMBERS[@]}"; do
    model_num=${MODEL_NUMBERS[$i]}
    mu_s=${MU_S_VALUES[$i]}

    if create_model $model_num $mu_s; then
        ((created_count++))
    else
        ((failed_count++))
    fi
done

# Summary
echo "========================================================================"
echo " SUMMARY"
echo "========================================================================"
echo "Successfully created: ${created_count} model directories"
if [ $failed_count -gt 0 ]; then
    echo "Failed: ${failed_count} model directories"
fi

# List created directories
echo ""
echo "Model directories:"
for model_num in "${MODEL_NUMBERS[@]}"; do
    model_name=$(printf "model_%02d" $model_num)
    model_dir="${BASE_DIR}/${model_name}"
    if [ -d "$model_dir" ]; then
        echo "  - ${model_dir}"
    fi
done

echo ""
echo "Done!"
