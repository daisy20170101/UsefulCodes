# Troubleshooting Guide

## Error: "No data source found in state file"

### Problem
You're trying to use the state file without providing input data.

### Solution
Use the **XDMF input method** instead:

```bash
# ✓ Correct - provide XDMF file as input
pvpython generate_fault_snapshots.py -i your_fault_data.xdmf

# ✗ Wrong - state file alone (causes error)
pvpython generate_fault_snapshots.py -s rakeNew.pvsm
```

### Understanding the Error

The `rakeNew.pvsm` state file contains:
- Camera angles and view settings
- Color map configurations
- Rendering preferences

It does **NOT** contain:
- The actual fault data
- Mesh geometry
- Field values (SRs, SRd)

### Workflow Options

#### Option 1: XDMF File Only (Recommended)
```bash
pvpython generate_fault_snapshots.py -i fault_output.xdmf
```
This loads data from XDMF and uses default visualization settings.

#### Option 2: XDMF File + State File
```bash
pvpython generate_fault_snapshots.py -i fault_output.xdmf -s rakeNew.pvsm
```
This loads data from XDMF and applies camera/view settings from the state file.

#### Option 3: State File with Embedded Data (Legacy)
If your state file was saved with data included:
```bash
pvpython generate_fault_snapshots.py -s my_state_with_data.pvsm
```
**Note:** The `rakeNew.pvsm` file likely doesn't have embedded data.

---

## Diagnostic Tool

Check what's in your state file:

```bash
pvpython check_state_file.py rakeNew.pvsm
```

This will show:
- Whether the state file loads successfully
- What data sources are present (if any)
- Available data fields
- Camera and view settings

---

## Common Issues

### Issue 1: Missing XDMF File

**Error:**
```
FileNotFoundError: Data file not found: fault_output.xdmf
```

**Solution:**
- Check the file path is correct
- Use absolute path or correct relative path
- Verify the file exists: `ls -l path/to/fault_output.xdmf`

**Example:**
```bash
# Check file exists
ls -l /path/to/simulation/fault_output.xdmf

# Use absolute path
pvpython generate_fault_snapshots.py -i /full/path/to/fault_output.xdmf

# Or relative path from current directory
pvpython generate_fault_snapshots.py -i ../data/fault_output.xdmf
```

---

### Issue 2: Missing SRs or SRd Fields

**Error:**
```
Error: Field 'SRs' not found in data
```

**Solution:**
Check what fields are available in your XDMF file:

```bash
# Use the diagnostic tool
pvpython check_state_file.py your_data.xdmf
```

Or check in ParaView GUI:
1. Open the XDMF file in ParaView
2. Look at the "Information" panel
3. Check "Point Arrays" or "Cell Arrays"

If SRs/SRd fields have different names in your data, you'll need to modify the script.

---

### Issue 3: State File Not Found

**Warning:**
```
Warning: State file not found: ./visualization/rakeNew.pvsm
```

**Solution:**
This is just a warning. The script will continue without camera settings.

To fix:
```bash
# Copy state file to visualization folder
cp /Users/DuoL/Documents/NSHM/Central/Paraviews/rakeNew.pvsm ./visualization/

# Or specify full path
pvpython generate_fault_snapshots.py -i data.xdmf -s /full/path/to/rakeNew.pvsm
```

---

### Issue 4: ParaView/pvpython Not Found

**Error:**
```
pvpython: command not found
```

**Solution:**

**On macOS:**
```bash
# Add ParaView to PATH
export PATH="/Applications/ParaView-5.11.0.app/Contents/bin:$PATH"

# Or use full path
/Applications/ParaView-5.11.0.app/Contents/bin/pvpython generate_fault_snapshots.py -i data.xdmf
```

**On Linux:**
```bash
# Add to PATH
export PATH="/opt/paraview/bin:$PATH"

# Or use full path
/opt/paraview/bin/pvpython generate_fault_snapshots.py -i data.xdmf
```

**Verify installation:**
```bash
which pvpython
pvpython --version
```

---

## Quick Reference

### Correct Usage Patterns

```bash
# Single snapshot from XDMF
pvpython generate_fault_snapshots.py -i fault.xdmf

# With state file for camera settings
pvpython generate_fault_snapshots.py -i fault.xdmf -s rakeNew.pvsm

# High resolution
pvpython generate_fault_snapshots.py -i fault.xdmf -r 3840 2160

# Only SRs field
pvpython generate_fault_snapshots.py -i fault.xdmf -f SRs

# Batch process multiple files
pvpython batch_fault_snapshots.py -d /path/to/data -p "*.xdmf"

# Batch process single file (all timesteps)
pvpython batch_fault_snapshots.py -i fault.xdmf
```

### Using the Wrapper Script

```bash
# Single snapshot
./run_snapshots.sh -i fault.xdmf

# Batch mode
./run_snapshots.sh -m batch -d /path/to/data

# With custom resolution
./run_snapshots.sh -i fault.xdmf -r 3840x2160
```

---

## Getting Help

1. **Check script help:**
   ```bash
   pvpython generate_fault_snapshots.py --help
   pvpython batch_fault_snapshots.py --help
   ```

2. **Use diagnostic tool:**
   ```bash
   pvpython check_state_file.py your_file.pvsm
   ```

3. **Check ParaView installation:**
   ```bash
   pvpython --version
   which pvpython
   ```

4. **Read the README:**
   ```bash
   cat visualization/README.md
   ```

---

## Still Having Issues?

If none of these solutions work, provide:
1. The exact command you're running
2. The complete error message
3. Output of `pvpython check_state_file.py rakeNew.pvsm`
4. Your ParaView version: `pvpython --version`
5. Operating system

This will help diagnose the specific issue.
