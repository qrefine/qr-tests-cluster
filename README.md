## Objective and steps

#### a) Select structure-candidates.

   * Select ~10-20 P1 structures from PDB and prepare them for quantum refinement.
   
#### b) Run cheap HF energy and gradients calculation using clustering.

   * pick up results from "a)" and  use it as input to run cheap HF calculation  using clustering.
   
   * write only gradient-based LBFGS (the line-search version)   

#### c) Final HF based refinement.

   * Analyse results of "b)" and pick one or two most favorable (showing most improvement)
   to run HF based refinements.

Details for each step a-b-c) follow below.

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

   All PDBs in folder 01 have no errors in running 02_run_update_pdb.py are listed in folder 02.  

   After checking those completed structures in 02 by visualisaton and comparing their charges with MOZYME calculated charges, all properly completed structures are stored in 03 and  ready for quantum refinement.

Current issues to be resolved ASAP:

  * Errors that have been postponed:
    * 1u0d - multiple models
    * 2jee - has a terminal amino acid with just a N
    * 3nm9 - DNA not currently supported 

  * Errors:
    *  4drw, 3nm9,  2oy0 and 1y1l have different charges from MOZYME

  * No errors: but not 100% sure that all is well.
    *  1va7 2oeq  4xa1 3tz9  3dtj   2ghj

  * Waters
    * Do we make a policy about the inclusion of water below a certain resolution? 

Testing
-------

Testing is very simple. Run tests/run_tests.py and check for Sorry's.

Structure refinement
--------------------

XXX more results here, in folders 04, 05, etc.
