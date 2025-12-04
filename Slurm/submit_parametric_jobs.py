#!/usr/bin/env python3
"""
Iteratively submit SLURM jobs with different parameters
Creates model directories, modifies parameters, and submits jobs
"""

import os
import shutil
import subprocess
import re
from pathlib import Path


def modify_parameters_par(file_path, output_file_value):
    """
    Modify OutputFile parameter in parameters.par

    Parameters:
    -----------
    file_path : str
        Path to parameters.par file
    output_file_value : str
        New value for OutputFile
    """
    with open(file_path, 'r') as f:
        content = f.read()

    # Replace OutputFile line
    # Match: OutputFile = '/path/to/output'
    pattern = r"(OutputFile\s*=\s*)'[^']*'"
    replacement = rf"\1'{output_file_value}'"

    new_content = re.sub(pattern, replacement, content)

    with open(file_path, 'w') as f:
        f.write(new_content)

    print(f"  Modified OutputFile to: {output_file_value}")


def modify_fault2_yaml(file_path, mu_s_value):
    """
    Modify mu_s parameter in fault2.yaml

    Parameters:
    -----------
    file_path : str
        Path to fault2.yaml file
    mu_s_value : float
        New value for mu_s
    """
    with open(file_path, 'r') as f:
        lines = f.readlines()

    # Find and replace mu_s in the LuaMap function (line 34)
    # Pattern: mu_s=0.8
    modified = False
    new_lines = []

    for line in lines:
        if re.search(r'^\s*mu_s\s*=\s*[\d.]+\s*$', line):
            # Extract indentation
            indent = re.match(r'^(\s*)', line).group(1)
            new_line = f'{indent}mu_s={mu_s_value}\n'
            new_lines.append(new_line)
            modified = True
            print(f"  Modified mu_s to: {mu_s_value} (line: {line.strip()} -> {new_line.strip()})")
        else:
            new_lines.append(line)

    with open(file_path, 'w') as f:
        f.writelines(new_lines)

    if not modified:
        print(f"  WARNING: Could not find mu_s line to modify!")

    return modified


def create_model_directory(base_dir, model_number,
                          output_path_template, mu_s_value,
                          submit_job=True):
    """
    Create a model directory with modified parameters

    Parameters:
    -----------
    base_dir : str
        Base directory containing template files
    model_number : int
        Model number (for folder name model_XX)
    output_path_template : str
        Template for output path (can include {model_num})
    mu_s_value : float
        Value for mu_s parameter
    submit_job : bool
        Whether to submit the job after setup

    Returns:
    --------
    model_dir : str
        Path to created model directory
    """

    base_dir = Path(base_dir)
    model_name = f"model_{model_number:02d}"
    model_dir = base_dir / model_name

    print(f"\n{'='*70}")
    print(f"Creating {model_name}")
    print(f"{'='*70}")

    # Create model directory
    if model_dir.exists():
        print(f"  WARNING: {model_dir} already exists!")
        response = input(f"  Overwrite? (y/n): ").lower()
        if response != 'y':
            print(f"  Skipping {model_name}")
            return None
        shutil.rmtree(model_dir)

    model_dir.mkdir(parents=True)
    print(f"  Created directory: {model_dir}")

    # Files to copy
    files_to_copy = ['parameters.par', 'submit_att.sh', 'fault2.yaml']

    # Copy files
    for filename in files_to_copy:
        src = base_dir / filename
        dst = model_dir / filename

        if not src.exists():
            print(f"  ERROR: Source file not found: {src}")
            return None

        shutil.copy2(src, dst)
        print(f"  Copied: {filename}")

    # Modify parameters.par
    output_path = output_path_template.format(model_num=model_number)
    params_file = model_dir / 'parameters.par'
    modify_parameters_par(params_file, output_path)

    # Modify fault2.yaml
    fault_file = model_dir / 'fault2.yaml'
    modify_fault2_yaml(fault_file, mu_s_value)

    print(f"\n  Summary:")
    print(f"    Model directory: {model_dir}")
    print(f"    OutputFile: {output_path}")
    print(f"    mu_s: {mu_s_value}")

    # Submit job
    if submit_job:
        print(f"\n  Submitting job...")
        try:
            submit_script = model_dir / 'submit_att.sh'
            result = subprocess.run(
                ['sbatch', str(submit_script)],
                cwd=str(model_dir),
                capture_output=True,
                text=True
            )

            if result.returncode == 0:
                job_id = result.stdout.strip()
                print(f"  ✓ Job submitted: {job_id}")
            else:
                print(f"  ✗ Job submission failed!")
                print(f"    Error: {result.stderr}")
        except Exception as e:
            print(f"  ✗ Error submitting job: {e}")

    return str(model_dir)


def batch_submit_jobs(base_dir,
                     model_numbers,
                     mu_s_values,
                     output_path_template='/nesi/nobackup/gns04005/daisy/output/tibet/t16mod7PL_model{model_num:02d}',
                     submit_jobs=True,
                     dry_run=False):
    """
    Batch create model directories and submit jobs

    Parameters:
    -----------
    base_dir : str
        Base directory containing template files
    model_numbers : list
        List of model numbers to create
    mu_s_values : list
        List of mu_s values (same length as model_numbers)
    output_path_template : str
        Template for output paths
    submit_jobs : bool
        Whether to submit jobs
    dry_run : bool
        If True, only show what would be done without executing
    """

    if len(model_numbers) != len(mu_s_values):
        raise ValueError("model_numbers and mu_s_values must have the same length!")

    print("="*70)
    print(" PARAMETRIC JOB SUBMISSION")
    print("="*70)
    print(f"Base directory: {base_dir}")
    print(f"Number of jobs: {len(model_numbers)}")
    print(f"Submit jobs: {submit_jobs}")
    print(f"Dry run: {dry_run}")
    print()

    if dry_run:
        print("DRY RUN MODE - No files will be created or jobs submitted")
        print()
        for model_num, mu_s in zip(model_numbers, mu_s_values):
            output_path = output_path_template.format(model_num=model_num)
            print(f"  model_{model_num:02d}:")
            print(f"    OutputFile: {output_path}")
            print(f"    mu_s: {mu_s}")
        return

    created_dirs = []

    for model_num, mu_s in zip(model_numbers, mu_s_values):
        try:
            model_dir = create_model_directory(
                base_dir=base_dir,
                model_number=model_num,
                output_path_template=output_path_template,
                mu_s_value=mu_s,
                submit_job=submit_jobs
            )

            if model_dir:
                created_dirs.append(model_dir)

        except Exception as e:
            print(f"\n  ERROR creating model_{model_num:02d}: {e}")
            continue

    print("\n" + "="*70)
    print(" SUMMARY")
    print("="*70)
    print(f"Successfully created {len(created_dirs)} model directories:")
    for d in created_dirs:
        print(f"  - {d}")


# ===========================================================================
# Main script
# ===========================================================================

if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(
        description="Create model directories and submit parametric SLURM jobs"
    )

    parser.add_argument(
        '--base-dir',
        default='/Users/DuoL/Documents/SeisSol/slurm/tibet',
        help='Base directory containing template files'
    )

    parser.add_argument(
        '--models',
        type=int,
        nargs='+',
        help='Model numbers to create (e.g., --models 1 2 3 4 5)'
    )

    parser.add_argument(
        '--mu-s-values',
        type=float,
        nargs='+',
        help='mu_s values for each model (e.g., --mu-s-values 0.5 0.6 0.7 0.8 0.9)'
    )

    parser.add_argument(
        '--mu-s-range',
        type=float,
        nargs=3,
        metavar=('START', 'STOP', 'STEP'),
        help='Generate mu_s values as range (start, stop, step)'
    )

    parser.add_argument(
        '--output-template',
        default='/nesi/nobackup/gns04005/daisy/output/tibet/t16mod7PL_model{model_num:02d}',
        help='Template for output paths (use {model_num} placeholder)'
    )

    parser.add_argument(
        '--no-submit',
        action='store_true',
        help='Create directories but do not submit jobs'
    )

    parser.add_argument(
        '--dry-run',
        action='store_true',
        help='Show what would be done without executing'
    )

    args = parser.parse_args()

    # Determine model numbers and mu_s values
    if args.models and args.mu_s_values:
        model_numbers = args.models
        mu_s_values = args.mu_s_values
    elif args.models and args.mu_s_range:
        import numpy as np
        model_numbers = args.models
        start, stop, step = args.mu_s_range
        mu_s_values = list(np.arange(start, stop + step/2, step))
        if len(mu_s_values) != len(model_numbers):
            print(f"ERROR: Generated {len(mu_s_values)} mu_s values but {len(model_numbers)} models specified")
            print(f"       Adjust range or model list")
            exit(1)
    else:
        # Default example: 5 models with mu_s from 0.5 to 0.9
        print("No models specified, using default example:")
        print("  Models: 1-5")
        print("  mu_s: 0.5, 0.6, 0.7, 0.8, 0.9")
        print()
        model_numbers = [1, 2, 3, 4, 5]
        mu_s_values = [0.5, 0.6, 0.7, 0.8, 0.9]

    # Run batch submission
    batch_submit_jobs(
        base_dir=args.base_dir,
        model_numbers=model_numbers,
        mu_s_values=mu_s_values,
        output_path_template=args.output_template,
        submit_jobs=not args.no_submit,
        dry_run=args.dry_run
    )
