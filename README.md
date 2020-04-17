# Maximum Achievable Potency of a Protein Surface

## Step 1: Determine Target Residue
For example: Protein *Mdm2*. Select 1 residue (non-buried) on the protein surface as target residue.

Or iterate over the all the protein surface residues with SASA > 10 Å

## Step 2: Open pockets at the target site (Optional)
Rosetta based application that explores low-energy fluctuations of the protein surface to reveal cryptic pockets.

**Sample Command**: relax.linuxgccrelease -s input_pdb -relax:fast -pocket_max_spacing 12 -pocket_zero_derivatives -pocket_psp false -pocket_sps -pocket_num_angles 2 -ex1  ex1aro -ex2 -score:patch pocket.wts.patch -nstruct 1 -cst_fa_file constraints

**Dependency**: Rosetta Software Suite

If there is a good pocket at the target site then skip to Step 3

## Step 3: Build Exemplars
**“Exemplar”**: a perfect, but nonphysical, pseudoligand that would optimally match the shape and chemical features of the pocket.

**Dependency**: Rosetta Software Suite

**Sample Command**: make_exemplar.linuxgccrelease -database ~/Rosetta/main/database -in:file:s input.pdb -central_relax_pdb_num 99:A -pocket_grid_size 12 -pocket_static_grid -pocket_filter_by_exemplar

## Step 4.1: Obtain Pocket
**Script**: *prox4_exemplar.pl*

**"Pocket"**: protein surface in 4Å proximity to the exemplar

## Step 4.2: Calculate Pocket Features
**Script**: *prox4_exemplar.pl* && *features.pl*

**"Pocket Features"**: 13 features of protein pocket geometric and  physiochemical properties.

**Dependencies**: RADI and Naccess

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
**Script**: *FormSets.pl* 

This script forms the test set. 

**Usage**: perl FormSets.pl > DataSet


## Step 6: Applying GBM model
**Model Name**: *model_FINAL.rds*

**Usage**: Rscript GBMrunFINAL.R

Predictions are added in the last column of the Set.

# Download Prerequisites 

1. Download RADI
	- Download binary from http://petitjeanmichel.free.fr/itoweb.petitjean.freeware.html#RADI
	- Change path for binary in line 32 script features.pl
2. Download Naccess
	- Download Naccess app from http://wolf.bms.umist.ac.uk/naccess/
	- Change path for the app in line 86 script features.pl

# Walk-Through: Example of Mdm2
**Protein**: *Mdm2*

**Target Site**: Chain *A* Residue *99*

###### \#\# Skip these steps for this particular example \#\#

mkdir ~/PocketDruggability/data/TA99/

cd ~/PocketDruggability/data/TA99/


mkdir complexes/

cd complexes/

\#\# Place Protein-Exemplar Complexes here. Filename: "*_Complex.pdb"

cd ../


mkdir apo/

cd apo/

\#\# Place Protein Apo Structure files here. Filename: "*.pdb"


cd ~/PocketDruggability/

###### \#\# Previous steps have already been done for this particular example. Start Here. \#\#

perl prox4_exemplar.pl

perl FormSets.pl > TA99Set

Rscript GBMrunFINAL.R

**Output file**: *TA99Predictions.txt*

**Column15** : PredictedActivity (Predicted attainable pactivity for a pocket)

## References
1. Borrel A,  Regad L,  Xhaard H,  Petitjean M, Camproux A-C, PockDrug: A Model for Predicting Pocket Druggability That Overcomes Pocket Estimation Uncertainties. *Journal of Chemical Information and Modeling* **2015**, 55 (4), 882-895.
2. Burgoyne NJ, Jackson RM, Predicting protein interaction sites: binding hot-spots in protein–protein and protein–ligand interfaces. *Bioinformatics* **2006**, 22 (11), 1335-1342.
3. Petitjean M, Applications of the radius-diameter diagram to the classification of topological and geometrical shapes of chemical compounds. *Journal of chemical information and computer sciences* **1992**, 32 (4), 331-337.
4. Kyte J, Doolittle RF, A simple method for displaying the hydropathic character of a protein. *Journal of molecular biology* **1982**, 157 (1), 105-132.
5. Milletti F. Vulpetti A, Predicting polypharmacology by binding site similarity: from kinases to the protein universe. *Journal of chemical information and modeling* **2010**, 50 (8), 1418-1431.

