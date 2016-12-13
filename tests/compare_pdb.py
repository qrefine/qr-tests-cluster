import os, sys

import iotbx
from iotbx import pdb

from run_finalise import calculate_pdb_hierarchy_charge, get_hetero_charges
from run_finalise import get_processed_pdb, get_inter_residue_bonds
from run_finalise import default_ion_charges

def get_hierarchy(pdb_filename):
  pdb_inp = iotbx.pdb.input(pdb_filename)
  return pdb_inp.construct_hierarchy()

def compare_pdbs(pdb1, pdb2):
  h1 = get_hierarchy(pdb1)
  h2 = get_hierarchy(pdb2)
  for rg1, rg2 in zip(h1.residue_groups(), h2.residue_groups()):
    #print rg1, rg2
    if len(rg1.atoms())!=len(rg2.atoms()):
      print '1'
      for atom in rg1.atoms(): print atom.quote()
      print '2'
      for atom in rg2.atoms(): print atom.quote()

def run(*args):
  print "run",args
  compare_pdbs(args[0], args[1])

if __name__=="__main__":
  args = sys.argv[1:]
  del sys.argv[1:]
  run(*tuple(args))
