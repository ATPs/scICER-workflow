# scICER Workflow Scripts

This directory follows the structure of the original `ATPs/scICER-workflow` repository and contains simplified example workflows for scICER.

## Core scripts

1. `run_scICER.R`
   General single-sample scICER example for RNA, SCT, or scLENS-style inputs.

2. `Harmony integration + scICER`
   Harmony-based integration followed by scICER.

3. `scVI integration + scICER`
   scVI-based integration followed by scICER.

The scripts now focus on the analysis steps themselves:

- generic input arguments instead of local paths;
- no manuscript-specific file names;
- no figure-export steps by default.

## Figure scripts

The manuscript plotting scripts are stored separately in `codes/figures/` so that the original workflow examples remain easy to find.
