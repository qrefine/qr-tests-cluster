## structure preparation 
0) structure download: ~25 structures (P1, not too many atoms, resolution ~3 Å, bad validation metrics) were selected from the PDB. *.pdb files are in 0_pdb, and *.mtz files are in 0_mtz.

1) Residues missing:
[1il5.pdb 2pro.pdb 3kyi.pdb 3tz9.pdb 4ctd.pdb 4drw.pdb 4k2r.pdb 4rnf.pdb 4tql.pdb 4xa1.pdb 5d12.pdb]

It is not available now in run_finalize.py. 
Just keep those structures without missing residues in a chain in 1_pdb

2) Keep completed structures by run_finalise.py

Errors:
ligand or metal
[1ok9.pdb 1pag.pdb 1u0d.pdb 1va7.pdb 1y1l.pdb 2ghj.pdb 2iwe.pdb 2oy0.pdb 3nm9.pdb 4l21.pdb]:
run_finalise.py throw error like 
     Sorry: no charge found in the model file for ""HETATM 7669 ZN    ZN A1129 .*.    ZN  ""
or:
  AssertionError: residue "HETATM 8395  N   SAH A 301 .*.     N  " charge 2 is greater than 1


No errors:
[1fh5.pdb 2jee.pdb 2oeq.pdb 3dtj.pdb]


*Incompleteness:

No errors throw for 2jee, while 2jee_complete.pdb  has missing atoms through visualization.

3) Refine using standard cctbx.

1fh5 is not good because data seem to be corrupted, plus published R factors were not reproducible.

4) Complete those structures from 3) as input structures for Q|R

## Q|R refinement procedure.
a) pre-refine using mozyeme for all structures in 4_pdb using full system (i.e. no clustering) 2oeq_refine_001.pdb, 3dtj_refine_001.pdb


b) Run cluster-Q|R refinement -  run clustering at two levels:
 1. semi-empirical  for all structures
 2. HF-2c for best structure (based on metrics) found so far.

sample command:
phenix.python ../qr-core/qrefine.py a87_99_h.pdb.mtz perturbed/1.5/4.pdb 


_analyze.py is to analyze the results.  

