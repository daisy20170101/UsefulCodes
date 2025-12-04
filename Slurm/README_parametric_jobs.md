# Parametric SLURM Job Submission

Automated scripts to create model directories and submit parametric SLURM jobs with different parameters.

## Files

- **`submit_parametric_jobs.py`** - Python version (more flexible, recommended)
- **`submit_parametric_jobs.sh`** - Bash version (simpler, standalone)

## What These Scripts Do

For each model:
1. Create a new directory `model_XX` (where XX is the model number)
2. Copy `parameters.par`, `submit_att.sh`, and `fault2.yaml` into the new directory
3. Modify `OutputFile` in `parameters.par`
4. Modify `mu_s` value in `fault2.yaml`
5. Optionally submit the SLURM job using `sbatch`

## Python Version (Recommended)

### Quick Start

```bash
# Default example: Create models 1-5 with mu_s from 0.5 to 0.9
./submit_parametric_jobs.py

# Dry run (see what would be done without executing)
./submit_parametric_jobs.py --dry-run
```

### Basic Usage

```bash
# Specify model numbers and mu_s values explicitly
./submit_parametric_jobs.py \
    --models 1 2 3 4 5 \
    --mu-s-values 0.5 0.6 0.7 0.8 0.9

# Or use a range for mu_s (start, stop, step)
./submit_parametric_jobs.py \
    --models 1 2 3 4 5 6 7 8 9 10 \
    --mu-s-range 0.5 0.95 0.05
```

### Advanced Options

```bash
# Custom output path template
./submit_parametric_jobs.py \
    --models 10 11 12 \
    --mu-s-values 0.75 0.80 0.85 \
    --output-template '/nesi/nobackup/gns04005/daisy/output/tibet/newrun_model{model_num:02d}'

# Create directories but don't submit jobs
./submit_parametric_jobs.py \
    --models 1 2 3 \
    --mu-s-values 0.6 0.7 0.8 \
    --no-submit

# Custom base directory
./submit_parametric_jobs.py \
    --base-dir /path/to/templates \
    --models 1 2 3 \
    --mu-s-values 0.5 0.6 0.7
```

### All Options

```
--base-dir PATH           Base directory with template files
                         (default: /Users/DuoL/Documents/SeisSol/slurm/tibet)

--models N1 N2 N3 ...    Model numbers to create

--mu-s-values V1 V2 ...  mu_s values for each model

--mu-s-range START STOP STEP
                         Generate mu_s values as range

--output-template PATH   Template for output paths
                         (use {model_num} placeholder)
                         (default: /nesi/nobackup/.../t16mod7PL_model{model_num:02d})

--no-submit             Create directories but don't submit jobs

--dry-run               Show what would be done without executing
```

## Bash Version

### Quick Start

Edit the configuration section in `submit_parametric_jobs.sh`:

```bash
# Edit these arrays
MODEL_NUMBERS=(1 2 3 4 5)
MU_S_VALUES=(0.5 0.6 0.7 0.8 0.9)

# Set to 0 to skip job submission
SUBMIT_JOBS=1
```

Then run:

```bash
./submit_parametric_jobs.sh
```

### Customization

Open `submit_parametric_jobs.sh` and modify:

```bash
# Base directory containing template files
BASE_DIR="/Users/DuoL/Documents/SeisSol/slurm/tibet"

# Output path prefix
OUTPUT_BASE="/nesi/nobackup/gns04005/daisy/output/tibet/t16mod7PL"

# Model numbers
MODEL_NUMBERS=(1 2 3 4 5 6 7 8 9 10)

# Corresponding mu_s values (must match length of MODEL_NUMBERS)
MU_S_VALUES=(0.50 0.55 0.60 0.65 0.70 0.75 0.80 0.85 0.90 0.95)

# Set to 0 to create directories without submitting jobs
SUBMIT_JOBS=1
```

## Examples

### Example 1: Test Run (No Submission)

```bash
# Python version
./submit_parametric_jobs.py \
    --models 99 \
    --mu-s-values 0.75 \
    --no-submit

# Bash version (edit script to set SUBMIT_JOBS=0)
```

### Example 2: Fine Parameter Sweep

```bash
# Create 20 models with mu_s from 0.50 to 0.95 (step 0.025)
./submit_parametric_jobs.py \
    --models $(seq 1 20) \
    --mu-s-range 0.50 0.975 0.025
```

### Example 3: Specific Values

```bash
# Test specific mu_s values of interest
./submit_parametric_jobs.py \
    --models 1 2 3 4 \
    --mu-s-values 0.38 0.50 0.70 0.85
```

### Example 4: Dry Run First

```bash
# Check what will be created before executing
./submit_parametric_jobs.py \
    --models 1 2 3 4 5 \
    --mu-s-range 0.5 0.9 0.1 \
    --dry-run

# If looks good, run without --dry-run
./submit_parametric_jobs.py \
    --models 1 2 3 4 5 \
    --mu-s-range 0.5 0.9 0.1
```

## Output Structure

After running, your directory structure will be:

```
/Users/DuoL/Documents/SeisSol/slurm/tibet/
├── parameters.par           # Original template
├── submit_att.sh           # Original template
├── fault2.yaml             # Original template
├── submit_parametric_jobs.py
├── submit_parametric_jobs.sh
├── model_01/               # Created by script
│   ├── parameters.par      # Modified (OutputFile changed)
│   ├── submit_att.sh       # Copy
│   └── fault2.yaml         # Modified (mu_s changed)
├── model_02/
│   ├── parameters.par
│   ├── submit_att.sh
│   └── fault2.yaml
└── ...
```

## What Gets Modified

### In `parameters.par`:

```
# Before:
OutputFile = '/nesi/nobackup/gns04005/daisy/output/tibet/t16mod7PL'

# After (for model_01):
OutputFile = '/nesi/nobackup/gns04005/daisy/output/tibet/t16mod7PL_model01'
```

### In `fault2.yaml`:

```
# Before (line 34):
mu_s=0.8

# After (for model with mu_s=0.6):
mu_s=0.6
```

## Checking Job Status

After submission, check your jobs:

```bash
# View all your jobs
squeue -u $USER

# View specific job details
scontrol show job <JOB_ID>

# Cancel a job
scancel <JOB_ID>

# Cancel all your jobs
scancel -u $USER
```

## Monitoring Output

Check output files in the specified output directories:

```bash
# Check if output directory was created
ls -lh /nesi/nobackup/gns04005/daisy/output/tibet/

# Monitor a specific model's output
ls -lh /nesi/nobackup/gns04005/daisy/output/tibet/t16mod7PL_model01/
```

## Troubleshooting

### Script won't run

```bash
# Make executable
chmod +x submit_parametric_jobs.py
chmod +x submit_parametric_jobs.sh
```

### Python version issues

```bash
# Use python3 explicitly
python3 submit_parametric_jobs.py --models 1 2 3 --mu-s-values 0.5 0.6 0.7
```

### Directory already exists

The script will prompt you to overwrite. Answer 'y' to continue or 'n' to skip.

### Job submission fails

- Check that you're on a SLURM system
- Verify `sbatch` command is available: `which sbatch`
- Check your SLURM account and allocations: `sacctmgr show user $USER`

### Wrong mu_s modified in fault2.yaml

The script modifies the first occurrence of `mu_s=<number>` in the file. If your YAML structure is different, you may need to adjust the regex pattern in the script.

## Notes

- **Template files**: The original `parameters.par`, `submit_att.sh`, and `fault2.yaml` in the base directory are not modified
- **Overwriting**: If a model directory already exists, you'll be prompted to overwrite
- **Job names**: All jobs will have the name "seissol" (from `#SBATCH -J seissol` in submit_att.sh)
- **Output paths**: Ensure the output directory path exists on the compute system

## Advanced: Generating Parameter Sets

You can generate model numbers and mu_s values programmatically:

```bash
# Using seq for model numbers
./submit_parametric_jobs.py \
    --models $(seq 1 10) \
    --mu-s-range 0.5 0.95 0.05

# Or in Python:
python3 << EOF
import numpy as np
models = list(range(1, 11))
mu_s_vals = list(np.linspace(0.5, 0.95, 10))
print(f"--models {' '.join(map(str, models))}")
print(f"--mu-s-values {' '.join(map(str, mu_s_vals))}")
EOF
```

## Cleaning Up

To remove all model directories:

```bash
# Be careful! This removes all model_* directories
rm -rf /Users/DuoL/Documents/SeisSol/slurm/tibet/model_*

# Or remove specific models
rm -rf /Users/DuoL/Documents/SeisSol/slurm/tibet/model_{01..05}
```
