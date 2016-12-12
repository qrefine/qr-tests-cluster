import os, sys

import iotbx

from run_finalise import get_hetero_charges, default_ion_charges, remove_alt_loc
from run_finalise import write_hierarchy, get_processed_pdb
from run_finalise import get_inter_residue_bonds, complete_pdb_hierarchy
from run_finalise import calculate_pdb_hierarchy_charge

def run(pdb_filename):
  print "run",pdb_filename
  pdb_inp = iotbx.pdb.input(pdb_filename)
  hetero_charges = get_hetero_charges(pdb_inp)
  if not hetero_charges:
    # some defaults
    hetero_charges = default_ion_charges
  hierarchy = pdb_inp.construct_hierarchy()
  #ppf = get_processed_pdb(pdb_filename=pdb_filename)
  #inter_residue_bonds = get_inter_residue_bonds(ppf)
  total_charge = calculate_pdb_hierarchy_charge(
    hierarchy,
    #hetero_charges=hetero_charges,
    #inter_residue_bonds=inter_residue_bonds,
  )
  print "total_charge",total_charge

if __name__=="__main__":
  args = sys.argv[1:]
  del sys.argv[1:]
  run(*tuple(args))
