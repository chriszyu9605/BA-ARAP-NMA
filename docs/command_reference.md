# Command reference

Use `ba-arap-nma <command> --help` for complete built-in help.

## Global version command

```bash
ba-arap-nma --version
# or
ba-arap-nma version
```

## `ba-arap-nma generate`

Required:

- `--pdb INPUT.pdb`
- `--output DIRECTORY`

Options:

- `--chain A`
- `--modes 1-10` or `--modes 1,3,5`
- `--orientations -1 1`
- `--max-amplitude 20`
- `--amplitude-step 0.25`
- `--cutoff 15`
- `--overwrite`: replace a non-empty output directory
- `--quiet`: suppress progress messages

## `ba-arap-nma two-state`

Required:

- `--input-pdb INPUT.pdb`
- `--target-pdb SECOND_STATE.pdb`
- `--input-chain A`
- `--target-chain A`
- `--output DIRECTORY`

Additional options:

- `--analysis-modes 50`
- `--summary-k 10`
- `--min-identity 0.85`
- `--min-coverage 0.80`
- all common mode, orientation, amplitude, cutoff, output-safety, and progress options listed above.
