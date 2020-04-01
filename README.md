# Maximum Achievable Potency of a Protein Surface

### Pre-processing Step ###

## Step 1: Determine Target Residue
For example: DNA binding site on p53. Select 1 residue (non-buried) on protein surface as target residue.

Or iterate over the all the protein surface residues with SASA > 10 A

## Step 2: Open pockets at target site (Optional)
Rosetta based application that explores low-energy fluctuations of the protein surface to reveal cryptic pockets.

Dependency: Rosetta Software Suite

If there is a good pocket at the target site then skip to Step 3

## Step 3: Build Exemplars
“Exemplar”: a perfect, but nonphysical, pseudoligand that would optimally match the shape and chemical features of the pocket.

Dependency: Rosetta Software Suite

Sample Command: make_exemplar.linuxgccrelease -database ~/Rosetta/main/database -in:file:s input.pdb -central_relax_pdb_num 97:A,143:A -pocket_grid_size 12 -pocket_static_grid -pocket_filter_by_exemplar

## Step 4: Obtain Pocket
Script: prox4_exemplar.pl

"Pocket": protein surface in 4A proximity to the exemplar
