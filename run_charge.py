import os, sys

import iotbx

from run_finalise import write_hierarchy, get_processed_pdb

def get_total_charge_from_file_name(pdb_filename,
                                    hetero_charges=None,
                                    inter_residue_bonds=None,
                                   ):
  from run_finalise import get_hetero_charges, default_ion_charges
  from run_finalise import get_inter_residue_bonds
  from run_finalise import calculate_pdb_hierarchy_charge
  ppf = get_processed_pdb(pdb_filename=pdb_filename)
  if not hetero_charges:
    hetero_charges = get_hetero_charges(ppf.all_chain_proxies.pdb_inp)
    if not hetero_charges:
      # some defaults
      hetero_charges = default_ion_charges
  if not inter_residue_bonds:
    inter_residue_bonds = get_inter_residue_bonds(ppf)
  total_charge = calculate_pdb_hierarchy_charge(
    ppf.all_chain_proxies.pdb_hierarchy,
    hetero_charges=hetero_charges,
    inter_residue_bonds=inter_residue_bonds,
  )
  return total_charge

def run(pdb_filename):
  print "run",pdb_filename
  total_charge = get_total_charge_from_file_name(pdb_filename)
  print "total_charge",total_charge

if __name__=="__main__":
  args = sys.argv[1:]
  del sys.argv[1:]
  run(*tuple(args))
