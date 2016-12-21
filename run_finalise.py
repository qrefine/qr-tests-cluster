# LIBTBX_SET_DISPATCHER_NAME phenix.development.ready_set
import os, sys
from string import letters

import iotbx
import iotbx.pdb
from libtbx.utils import Sorry
from libtbx import easy_run

from cctbx.array_family import flex

from iotbx.pdb.atom_name_interpretation import \
  interpreters as protein_atom_name_interpreters



#from libtbx.utils import null_out

skip = [ # see README.md for reasons
  '1il5.pdb',
  '2jee.pdb',
  #'1u0d.pdb',
  '3nm9.pdb',
  '3kyi.pdb',
  '4k2r.pdb',
  '4rnf.pdb',
  '5d12.pdb',
  ]


def junk():
  protein_interpreter = protein_atom_name_interpreters.get(ag.resname)
  atom_name_interpretation = None
  if (protein_interpreter is not None):
    atom_name_interpretation = protein_interpreter.match_atom_names(
        atom_names=atom_names)
    #if (atom_name_interpretation is not None):
    #    atom_name_interpretation.d_aa_residue_name = d_aa_rn
    assert not atom_name_interpretation.unexpected
    print atom_names
  else:
    raise Sorry("?????")

  assert 0

def remove_alt_loc(hierarchy):
  # should use cctbx
  for rg in hierarchy.residue_groups():
    if len(rg.atom_groups())==1: continue
    resnames = []
    for ag in rg.atom_groups():
      if ag.resname not in resnames: resnames.append(ag.resname)
    assert len(resnames)==1
    bag = rg.atom_groups()[0]
    bag.altloc = ""
    for i, ag in enumerate(rg.atom_groups()):
      if i==0: continue
      for atom in ag.atoms():
        if bag.get_atom(atom.name.strip()): continue
        da = atom.detached_copy()
        bag.append_atom(da)
      rg.remove_atom_group(ag)
  hierarchy.atoms_reset_serial()
  hierarchy.atoms().reset_i_seq()
  return hierarchy

def run_ready_set(pdb_filename):
  from StringIO import StringIO
  assert pdb_filename.find('.pdb')>-1, 'ReadySet! only works on PDB, not CIF'
  cmd = 'phenix.ready_set %s' % pdb_filename
  print 'Running ReadySet!'
  print cmd
  ero = easy_run.fully_buffered(command=cmd)
  std = StringIO()
  ero.show_stdout(std)
  for line in std.getvalue().splitlines():
    print line
  return os.path.basename(pdb_filename).replace('.pdb', '.updated.pdb') 
  # maybe read from stdout

def run_fetch_pdb(code):
  cmd = 'phenix.fetch_pdb %s' % code
  print 'Fetching files'
  print cmd
  easy_run.call(cmd)

def loop_over_dir(d):
  i=0
  for filename in os.listdir(d):
    if not filename.endswith('.pdb'): continue
    i+=1
    print '%s\n %3d %s\n%s' % ('*'*80, i, os.path.join(d, filename), '*'*80)
    if filename in skip:
      print 'skipping'
      continue
    if os.path.exists(filename.replace('.pdb', '.updated.pdb')):
      run(filename.replace('.pdb', '.updated.pdb'))
    else:
      run(os.path.join(d, filename))

def run(pdb_filename):
  print "run",pdb_filename
  if os.path.isdir(pdb_filename):
    loop_over_dir(pdb_filename)
    return
  if pdb_filename.find(".pdb")==-1:
    if not os.path.exists('%s.pdb' % pdb_filename):
      run_fetch_pdb(pdb_filename)
    pdb_filename = '%s.pdb' % pdb_filename
  if pdb_filename.endswith('.updated.pdb'):
    pdb_filename_h = pdb_filename
  else:
    pdb_filename_h = pdb_filename.replace('.pdb', '.updated.pdb')
  print 'pdb_filename_h',pdb_filename_h
  if not os.path.exists(pdb_filename_h):
    pdb_filename = run_ready_set(pdb_filename)
  else: pdb_filename = pdb_filename_h
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
                         #use_capping_hydrogens=True,
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
