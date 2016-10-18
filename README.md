# Objective

a) Select structure-candidates.

   Select ~10-20 P1 structures from PDB and run SE (mozyeme) refinements without using clustering.
   This requires finding structure-candidates and preparing them for refinement.
   
b) Once line-search is ready (Mark), pick up final or intermediate results from "a)" and
   use it as input to run SE using clustering.
   
c) Analyse results of "b)" and pick one or two most favorable (showing most improvement)
   to run HF based refinements.

Details for each step a-b-c) follow below.

Structure preparation
---------------------

1) Get structures from PDB. 

This is done manually using RCSB site. Selection criteria are: P1, not too many atoms, 
resolution ~3-4 Å, bad validation metrics. The criteria are rather arbitrary: we try
to find low-resolution models that we believe can be improved by quantum refinement.

File were downloaded using 
    phenix.fetch_pdb PDB_CODE --mtz

Some models have unknown to Phenix ligands. Corresponding CIF files were created using 
phenix.ready_set file_name.pdb

Some data file were missing R-free flags. They were added at subsequent (re-refinement).

*.pdb, *.mtz and *.cif files are stored in 00 folder.

2) Refine structures from "1)" using phenix.refine. List of commands per each structure is
   in 01_run_phenix_refine file.

Refinement results are stored in 01 folder. MTZ files from this folder are to be used in all
subsequent refinements.

3) Run structure from "2)" through run_finalise.py using command

   XXX

and store results in 03 folder. Structures in 03 are ready for quantum refinement.

Current issues to be resolved ASAP:

Errors:
ligand or metal
[1ok9 1pag 1u0d 1va7 1y1l 2ghj 2iwe 2oy0 3nm9 4l21]:
run_finalise.py throw error like 
     Sorry: no charge found in the model file for ""HETATM 7669 ZN    ZN A1129 .*.    ZN  ""
or:
  AssertionError: residue "HETATM 8395  N   SAH A 301 .*.     N  " charge 2 is greater than 1

No errors:
[1fh5 2jee 2oeq 3dtj]

No errors for 2jee, while 2jee_complete.pdb has missing atoms.

Structure refinement
--------------------

XXX more results here, in folders 04, 05, etc.
