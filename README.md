## Objective and steps

#### a) Select structure-candidates.

   * Select ~10-20 P1 structures from PDB and prepare them for quantum refinement.
   
#### b) Pre-screen by preliminary SE refinement.
   
   * Run SE (mozyeme) refinements without using clustering. Make sure qrefine outputs 
   intermediate PDB files so that refinement can be stopped any time and best outcome 
   selected.
   
#### c) Further SE refinement using clustering.

   * Once line-search is ready (Mark), pick up final or intermediate results from "b)" and
   use it as input to run SE using clustering.
   
   XXX Using what programs?
   
#### d) Final HF based refinement.

   * Analyse results of "c)" and pick one or two most favorable (showing most improvement)
   to run HF based refinements.

Details for each step a-b-c-d) follow below.

Structure preparation
---------------------

#### 1) Get structures from PDB. 

   * This is done manually using RCSB site. Selection criteria are: P1, not too many atoms, 
   resolution ~3-4 Å, bad validation metrics. The criteria are rather arbitrary: we try
   to find low-resolution models that we believe can be improved by quantum refinement.

   * File were downloaded using 
   phenix.fetch_pdb PDB_CODE --mtz

   * Some models have unknown to Phenix ligands. Corresponding CIF files were 
created using 
  phenix.ready_set file_name.pdb
This also adds hydrogens with are needed for the run_finalise.py script. The
output is file_name.updated.pdb which needs to be renamed? We can use 
phenix.reduce which has it's own problems.

   * Some data file were missing R-free flags. They were added at subsequent (re-refinement).

   * .pdb, .mtz and .cif files are stored in 00 folder.

#### 2) Phenix-Refine 

   * Refine structures from "1)" using phenix.refine. List of commands per each structure is
   in 01_run_phenix_refine file.

   * Refinement results are stored in 01 folder. MTZ files from this folder are to be used in all
   subsequent refinements.

#### 3) Finalize

   * Run structure from "2)" through 02_run_update_pdb.py using command
   phenix.ready_set is used to add hydrogens to amino acids and ligands
   run_finalise.py is used to complete structures, remove altlocs and resset occupancies
   all completed structrues by 02_run_update_pdb.py are in 02

   XXX This needs a careful and clear description what's being done and how. Adding H and missing
   non-H atoms? Remove altlocs? Reset occupancies? XXX

   After checking those completed structures in 02 by visualisaton, all proper completed structures are stored in 03 and  ready for quantum refinement.

Current issues to be resolved ASAP:

  * Errors that have been postponed:
    * 1u0d - multiple models
    * 2jee - has a terminal amino acid with just a N
    * 3nm9 - DNA not currently supported 

  * Errors:

  * No errors: but 100% sure that all is well.
    * 1fh5 2oeq 3dtj 1ok9 1pag 1va7 1y1l 2ghj 2iwe 2oy0 4l21

Structure refinement
--------------------

XXX more results here, in folders 04, 05, etc.
