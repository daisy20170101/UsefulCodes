# ParaView State File Placeholder

This folder should contain the ParaView state file: `rakeNew.pvsm`

## To add the state file:

Copy the state file from your local machine to this folder:

```bash
cp /Users/DuoL/Documents/NSHM/Central/Paraviews/rakeNew.pvsm ./visualization/
```

Or if you're setting up on a different machine:

```bash
# Copy to this directory
cp /path/to/your/rakeNew.pvsm ./visualization/rakeNew.pvsm
```

## State File Requirements

The state file should contain:
- Camera view settings
- Color map configurations for SRs and SRd fields
- Any filters or transformations
- Rendering preferences

## Note

Once you place `rakeNew.pvsm` in this folder, the scripts will use it as the default state file. You can still override it with the `-s` or `--state-file` option.
