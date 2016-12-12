import os, sys

import iotbx

from run_finalise import get_hetero_charges, default_ion_charges, remove_alt_loc
from run_finalise import write_hierarchy, get_processed_pdb
from run_finalise import get_inter_residue_bonds, complete_pdb_hierarchy

def run(pdb_filename):
  print "run",pdb_filename
  pdb_inp = iotbx.pdb.input(pdb_filename)
  hetero_charges = get_hetero_charges(pdb_inp)
  if not hetero_charges:
    # some defaults
    hetero_charges = default_ion_charges
  hierarchy = pdb_inp.construct_hierarchy()
  hierarchy = remove_alt_loc(hierarchy) # should use cctbx
  if 1:
    output = write_hierarchy(pdb_filename,
                             pdb_inp,
                             hierarchy,
                             'no_alt_loc',
                           )
    ppf = get_processed_pdb(pdb_filename=output)
  elif 0:
    raw_records = []
    for atom in hierarchy.atoms():
      raw_records.append(atom.format_atom_record())
    ppf = get_processed_pdb(raw_records=raw_records)
  else:
    ppf = get_processed_pdb(pdb_inp=hierarchy.as_pdb_input())
  inter_residue_bonds = get_inter_residue_bonds(ppf)
  complete_pdb_hierarchy(hierarchy,
                         ppf.geometry_restraints_manager(),
                         use_capping_hydrogens=True,
                        )
  # not required at the moment, no clutering
  if 0:
    write_pdb_hierarchy_qxyz_file(hierarchy,
                                  hetero_charges=hetero_charges,
                                )
  total_charge = calculate_pdb_hierarchy_charge(
    hierarchy,
    hetero_charges=hetero_charges,
    inter_residue_bonds=inter_residue_bonds,
  )
  print "total_charge",total_charge
  write_hierarchy(pdb_filename, pdb_inp, hierarchy, "complete")

if __name__=="__main__":
  args = sys.argv[1:]
  del sys.argv[1:]
  run(*tuple(args))
