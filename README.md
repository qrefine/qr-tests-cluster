# qr-tests-cluster
## structure prepareation 
0) structure selection: ~15 structures (P1, not too many atoms, resolution ~3 Å, bad validation metrics) were selected from the PDB.

Residues missing: 
[1il5.pdb 2pro.pdb 3kyi.pdb 3tz9.pdb 4ctd.pdb 4drw.pdb 4k2r.pdb 4rnf.pdb 4tql.pdb 4xa1.pdb 5d12.pdb] 

This causes a known problem of run_finalize not working for residue incomplete strucutres. Enhancement issue opened ()

1) run_finalize.py was run

Errors:
ligand or metal
[1ok9.pdb 1pag.pdb 1u0d.pdb 1va7.pdb 1y1l.pdb 2ghj.pdb 2iwe.pdb 2oy0.pdb 3nm9.pdb 4l21.pdb]:
run_finalise.py throw error like 
     Sorry: no charge found in the model file for ""HETATM 7669 ZN    ZN A1129 .*.    ZN  ""   
or:
  AssertionError: residue "HETATM 8395  N   SAH A 301 .*.     N  " charge 2 is greater than 1

No-errors.
[1fh5.pdb 2jee.pdb 2oeq.pdb 3dtj.pdb]: run_finalise.py

2) re-refine using standard cctbx.

1fh5 is not good because data seem to be corrupted, plus published R factors were not repoducible.

3) pre-refine using mozyeme for all structures using full system (i.e. no clustering)
2oeq_refine_001.pdb, 3dtj_refine_001.pdb

4) Run cluster-Q|R refinement -  run clustering at two levels:
 1. semi-empirical  for all structures
 2. HF-2c for best structure (based on metrics) found so far. 

sample command:
python ../qr-core/qrefine.py a87_99_h.pdb.mtz perturbed/1.5/4.pdb restraints=cctbx update_all_scales=False


5) run_analyze.py is to analyze the results.
