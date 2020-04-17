# Maximum Achievable Potency of a Protein Surface

## Step 1: Determine Target Residue
For example: DNA binding site on mdm2. Select 1 residue (non-buried) on protein surface as target residue.

Or iterate over the all the protein surface residues with SASA > 10 Å

## Step 2: Open pockets at target site (Optional)
Rosetta based application that explores low-energy fluctuations of the protein surface to reveal cryptic pockets.

Sample Command: relax.linuxgccrelease -s input_pdb -relax:fast -pocket_max_spacing 12 -pocket_zero_derivatives -pocket_psp false -pocket_sps -pocket_num_angles 2 -ex1  ex1aro -ex2 -score:patch pocket.wts.patch -nstruct 1 -cst_fa_file constraints

Dependency: Rosetta Software Suite

If there is a good pocket at the target site then skip to Step 3

## Step 3: Build Exemplars
“Exemplar”: a perfect, but nonphysical, pseudoligand that would optimally match the shape and chemical features of the pocket.

Dependency: Rosetta Software Suite

Sample Command: make_exemplar.linuxgccrelease -database ~/Rosetta/main/database -in:file:s input.pdb -central_relax_pdb_num 99:A -pocket_grid_size 12 -pocket_static_grid -pocket_filter_by_exemplar

## Step 4.1: Obtain Pocket
Script: prox4_exemplar.pl

"Pocket": protein surface in 4Å proximity to the exemplar

## Step 4.2: Calculate Pocket Features
Script: prox4_exemplar.pl && features.pl

"Pocket Features": 13 features of protein pocket geometric and  physiochemical properties.

Depedencies: RADI and Naccess

## ------------------- Pocket Features --------------------------
| Property Name  | Description |
| ------------- | ------------- |
| VOLUME_HULL  | volume of convex hull computed using RADI software  |
| hydrophobicity_kyte  | hydrophobicity based properties of residues  |
| SMALLEST_SIZE  | distance separating the two closest slabs enclosing the hull computed using RADI software  |
| INERTIA_3  | smallest eigenvalue of inertia matrix computed using RADI software  |
| p_N_atom  | frequency of N atoms in pocket  |
| hydrophobicity_pocket  | hydrophobicity pocket estimated with solvent accessibility computed using NACCESS software  |
| p_aliphatic_residues  | frequency of positive residues in pocket (I, L, V)   |
| p_aromatic_residues  | frequency of aromatic residues in pocket (F, Y, H, W)   |
| SURFACE_HULL | surface of convex hull computed using RADI software  |
| p_negative_residues | frequency of negative residues in pocket (D, E)  |
| C_RESIDUE   | number of residues in pocket  |
| p_Ooh_atom  | frequency of Ooh atoms in pocket  |
| p_Ccoo_atom  | frequency of Ccoo atoms in pocket  |

## Step 5: Form Set
Script: FormSets.pl 

This script forms the test set. 

Usage: perl FormSets.pl > DataSet


## Step 6: Applying GBM model
Model Name: model_FINAL.rds

Usage: Rscript GBMrunFINAL.R

Predictions adds the predicted activity as the last column of the Set.

# Walk-Through: Example of Mdm2

