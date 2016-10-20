# LIBTBX_SET_DISPATCHER_NAME phenix.development.ready_set
import os, sys
import math

import iotbx
import iotbx.pdb
from libtbx.utils import Sorry
from scitbx import matrix

from cctbx.array_family import flex

from iotbx.pdb.atom_name_interpretation import \
  interpreters as protein_atom_name_interpreters
get_class = iotbx.pdb.common_residue_names_get_class

from mmtbx.chemical_components import get_cif_dictionary

from mmtbx.monomer_library import server
mon_lib_server = server.server()

from mmtbx import monomer_library
from mmtbx.monomer_library import pdb_interpretation
#import mmtbx.monomer_library.server
#from iotbx import pdb
#from libtbx.utils import null_out

class chemical_component_class(dict):
  def get_total_charge(self):
    total = 0
    for atom in self.get("_chem_comp_atom", []):
      total += getattr(atom, 'charge', 0)
    return total

  def get_hydrogens(self):
    hs = []
    for atom in self.get("_chem_comp_atom", []):
      if getattr(atom, 'type_symbol') in ["H", "D"]:
        hs.append(atom)
    return hs

charge_per_aa_polymer = {}
hydrogens_per_aa_polymer = {}

def d_squared(xyz1, xyz2):
  d2 = 0
  for i in range(3):
    d2 += (xyz2[i]-xyz1[i])**2
  return d2

def get_aa_charge(code):
  # get from cache first
  # then look in the chemical components
  # not sure what to do about novel ligands...
  tmp = charge_per_aa_polymer.get(code, None)
  if tmp: return tmp
  cc = chemical_component_class()
  cc.update(get_cif_dictionary(code))
  tmp = cc.get_total_charge()
  charge_per_aa_polymer[code] = tmp
  return tmp

def get_aa_polymer_hydrogens(code):
  tmp = hydrogens_per_aa_polymer.get(code, None)
  if tmp: return tmp
  cc = chemical_component_class()
  cc.update(get_cif_dictionary(code))
  tmp = cc.get_hydrogens()
  hydrogens_per_aa_polymer[code] = tmp
  return tmp

def get_bond_vector(a1,a2,unit=False):
  vector = []
  l = 0
  for i in range(3):
    vector.append(a1.xyz[i]-a2.xyz[i])
    l+=vector[i]**2
  if unit:
    l=math.sqrt(l)
    for i in range(3):
      vector[i] /= l
  return tuple(vector)

def construct_xyz(ba, bv,
                  aa, av,
                  da, dv,
                  period=3,
                  ):
  assert ba is not None
  assert aa is not None
  assert da is not None
  rn = matrix.col(ba.xyz)
  rca = matrix.col(aa.xyz)
  rc = matrix.col(da.xyz)
  rcca = rc -rca

  e0 = (rn - rca).normalize()
  e1 = (rcca - (rcca.dot(e0))*e0).normalize()
  e2 = e0.cross(e1)

  pi = math.pi
  alpha = math.radians(av)
  phi = math.radians(dv)

  rh_list = []
  for n in range(0, period):
    rh = rn + bv * (math.sin(alpha)*(math.cos(phi + n*2*pi/period)*e1 +
                                     math.sin(phi + n*2*pi/period)*e2) - 
                    math.cos(alpha)*e0)
    rh_list.append(rh)
  return rh_list

def add_n_terminal_hydrogens_to_atom_group(ag):
  n = ag.get_atom("N")
  if n is None: return
  ca = ag.get_atom("CA")
  if ca is None: return
  c = ag.get_atom("C")
  if c is None: return
  if ag.get_atom("H"): # maybe needs to be smarter
    for atom in ag.atoms(): print atom.quote()
    ag.remove_atom(ag.get_atom('H'))
  # add H1
  rh3 = construct_xyz(n, 0.9,
                      ca, 109.5,
                      c, 120.,
                     )
  for i in range(0,3):
    name = " H%d " % (i+1)
    if ag.get_atom(name.strip()): continue
    atom = iotbx.pdb.hierarchy.atom()
    atom.name = name
    atom.element = "H"
    atom.xyz = rh3[i]
    atom.occ = n.occ
    atom.b = n.b
    ag.append_atom(atom)

def add_n_terminal_hydrogens(hierarchy,
                             #residue_selection=None,
                            ):
  # add N terminal hydrogens because Reduce only does it to resseq=1
  # needs to be alt.loc. aware for non-quantum-refine
  for chain_i, chain in enumerate(hierarchy.chains()):
    for res_i, residue_group in enumerate(chain.residue_groups()):
      if len(residue_group.atom_groups())>1: continue
      atom_group = residue_group.atom_groups()[0]
      if get_class(atom_group.resname) not in ["common_amino_acid",
                                               "modified_amino_acid",
                                             ]:
        continue
      if res_i==0: # need better switch
        add_n_terminal_hydrogens_to_atom_group(atom_group)
  hierarchy.atoms_reset_serial()
  hierarchy.atoms().reset_i_seq()
  return hierarchy

def add_c_terminal_oxygens_to_atom_group(ag):
  #
  # do we need ANISOU
  #
  if ag.get_atom("OXT"): return
  c = ag.get_atom("C")
  if c is None: return
  ca = ag.get_atom("CA")
  if ca is None: return
  n = ag.get_atom("N")
  if n is None: return
  ro2 = construct_xyz(c, 1.231,
                      ca, 120.,
                      n, 160.,
                      period=2,
                     )
  oxys = [' O  ', ' OXT'] 
  for i in range(0,2):
    name = oxys[i]
    atom = ag.get_atom(name.strip())
    #if atom:
    #  print name, atom.format_atom_record(), ro2[i]
    #  print 'd2',d_squared(atom.xyz, ro2[i])
    if atom:
      atom.xyz = ro2[i]
    else:
      atom = iotbx.pdb.hierarchy.atom()
      atom.name = name
      atom.element = "O"
      atom.occ = c.occ
      atom.b = c.b
      atom.xyz = ro2[i]
      ag.append_atom(atom)

def add_c_terminal_oxygens(hierarchy,
                             #residue_selection=None,
                            ):
  for chain_i, chain in enumerate(hierarchy.chains()):
    for res_i, residue_group in enumerate(chain.residue_groups()):
      if len(residue_group.atom_groups())>1: continue
      atom_group = residue_group.atom_groups()[0]
      if get_class(atom_group.resname) not in ["common_amino_acid",
                                               "modified_amino_acid",
                                             ]:
        continue
      if res_i==len(chain.residue_groups())-1: # need better switch
        add_c_terminal_oxygens_to_atom_group(atom_group)
  hierarchy.atoms_reset_serial()
  hierarchy.atoms().reset_i_seq()
  return hierarchy

def display_residue_group(rg):
  return '  residue_group: resseq="%s" icode="%s"' % (rg.resseq, rg.icode)

def display_atom_group(ag, verbose=False):
  if verbose:
    outl = display_atom_group(ag)
    for atom in ag.atoms():
      outl += '\n  %s ' % atom.format_atom_record()
    return outl
  return '    atom_group: resname="%s" altloc="%s"' % (ag.resname, ag.altloc)

def generate_residue_groups(hierarchy,
                            assert_no_alt_loc=False,
                            exclude_water=False,
                            verbose=False,
                            ):
  for rg in hierarchy.residue_groups():
    if verbose: display_residue_group(rg)
    if exclude_water:
      ag=rg.atom_groups()[0]
      if get_class(ag.resname)=='common_water': continue
    if assert_no_alt_loc:
      if len(rg.atom_groups())>1:
        raise Sorry("Contains alt. locs.")
    yield rg

def calculate_residue_charge(rg,
                             assert_contains_hydrogens=True,
                             assert_no_alt_loc=True,
                             verbose=False,
                             ):
  def _terminal(names, check):
    for name in check:
      if name not in names:
        break
    else:
      return True
    return False
  def n_terminal(names):
    return _terminal(names, [' H1 ', ' H2 ', ' H3 '])
  def c_terminal(names):
    return _terminal(names, [' OXT'])

  max_charge=1
  if assert_no_alt_loc:
    if len(rg.atom_groups())>1:
      raise Sorry("alt locs in %s" % display_residue_group(rg))
  # ions  
  if get_class(rg.atom_groups()[0].resname)=="common_element":
    atom = rg.atoms()[0]
    if not atom.charge.strip():
      raise Sorry('no charge found in the model file for "%s"' % atom.quote())
    else:
      return atom.charge_as_int()
  # others
  hs=0
  atom_names = []
  for atom in rg.atoms():
    if verbose: print atom.quote()
    if atom.element_is_hydrogen(): hs+=1
    atom_names.append(atom.name)
  if verbose: print get_class(rg.atom_groups()[0].resname)
  if assert_contains_hydrogens:
    if hs==0:
      hydrogens = get_aa_polymer_hydrogens(rg.atom_groups()[0].resname)
      if len(hydrogens)!=0:
        if verbose:
          for atom in rg.atoms(): print atom.quote()
        raise Sorry("no hydrogens: %s" % display_residue_group(rg))
  ag = rg.atom_groups()[0]
  charge = get_aa_charge(ag.resname)
  if ( get_class(ag.resname) in ["common_amino_acid", "modified_amino_acid"] or
       ag.resname in ['DVA', # need complete list of D-amino and non-standard
                      ]
       ):
    poly_hs = len(get_aa_polymer_hydrogens(ag.resname))-2
    diff_hs = hs-poly_hs
    if verbose: print 'charge: %s poly_hs: %s diff_hs: %s' % (charge,
                                                              poly_hs,
                                                              diff_hs,
                                                            )

    if n_terminal(atom_names):
      diff_hs-=1
    elif c_terminal(atom_names):
      diff_hs-=1
      max_charge+=1
    charge+=diff_hs

    assert abs(diff_hs)<=max_charge, 'residue %s charge %s is greater than %s' % (
      rg.atoms()[0].quote(),
      diff_hs,
      max_charge,
    )
  else:
    print "NOT POLY",ag.resname, get_class(ag.resname)
    if ag.resname not in ["NO3",
                          "PXZ",
                          "SRT",
                          "HM7",
                        ]:
      assert 0
  return charge

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

d_amino_acids = {"DVA" : "VAL",
                 }
non_standard_amino_acids = { #"SAR" : None,
                            }

def get_partial_point_charges(rg):
  """
  This function relies only on the residue group and monomer library server
  """
  v2_to_3 = {' HA3':' HA1',
             ' HB3':' HB1',
             ' HG3':' HG1',
             'HG13':'HG11',
             ' HD3':' HD1',
             ' HE3':' HE1',
             }
  misc = {' OXT' : ' O  ',
          }
  tmp = []
  for ag in rg.atom_groups():
    restraints = mon_lib_server.get_comp_comp_id_direct(ag.resname)
    if restraints is None:
      resname = d_amino_acids.get(ag.resname, None)
      if resname is not None:
        restraints = mon_lib_server.get_comp_comp_id_direct(resname)
    atom_dict = restraints.atom_dict()
    for atom in ag.atoms():
      # ions
      if get_class(ag.resname)=="common_element":
        assert len(ag.atoms())==1
        if not atom.charge.strip():
          raise Sorry('no charge found in the model file for "%s"' % atom.quote())
        else:
          tmp.append([atom.charge_as_int()]+list(atom.xyz))
          continue
      # other atoms
      cif = atom_dict.get(atom.name.strip(), None)
      if cif is None:
#        print atom.quote()
        if atom.name in [" H1 ", " H2 ", " H3 "]: # needs calculating...
          tmp.append([0]+list(atom.xyz))
          continue
        if atom.name in v2_to_3:
          cif = atom_dict.get(v2_to_3[atom.name].strip())
        elif atom.name in misc:
          cif = atom_dict.get(misc[atom.name].strip())
        elif atom.name.find("'")>-1:
          name = atom.name.replace("'", "*")
          cif = atom_dict.get(name.strip(), None)
      assert cif, "%s" % atom_dict
      tmp.append([cif.partial_charge]+list(atom.xyz))
  return tmp

def remove_alt_loc(hierarchy):
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

def get_processed_pdb(pdb_filename):
  # need to add params
  mon_lib_srv = monomer_library.server.server()
  ener_lib = monomer_library.server.ener_lib()
  ppf = monomer_library.pdb_interpretation.process(
    mon_lib_srv           = mon_lib_srv,
    ener_lib              = ener_lib,
    file_name             = pdb_filename,
    keep_monomer_mappings = True,
  #  log                   = null_out(),
  )
  return ppf

def remove_acid_side_chain_hydrogens(hierarchy):
  removes = {"GLU" : "HE2",
             "ASP" : "HD2",
             }
  for ag in hierarchy.atom_groups():
    r = removes.get(ag.resname, None)
    if r is None: continue
    atom = ag.get_atom(r)
    if atom:
      ag.remove_atom(atom)
  hierarchy.atoms_reset_serial()
  hierarchy.atoms().reset_i_seq()
  return hierarchy

def write_hierarchy(pdb_filename, pdb_inp, hierarchy, underscore):
  output = "%s_%s.pdb" % (pdb_filename.split(".")[0],
                          underscore,
                          )
  f=file(output, "wb")
  f.write(hierarchy.as_pdb_string(
    crystal_symmetry=pdb_inp.crystal_symmetry()),
          )
  f.close()
  print "\n  Output written to: %s" % output

def complete_pdb_hierarchy(hierarchy):
  hierarchy = remove_alt_loc(hierarchy)
  from mmtbx.building import extend_sidechains
  n_changed = extend_sidechains.extend_protein_model(hierarchy,
                                                     hydrogens=True,
                                                   )
  remove_acid_side_chain_hydrogens(hierarchy)
  add_n_terminal_hydrogens(hierarchy) # in place
  add_c_terminal_oxygens(hierarchy) # in place
  hierarchy.atoms().set_chemical_element_simple_if_necessary()
  hierarchy.sort_atoms_in_place()

def calculate_pdb_hierarchy_charge(hierarchy):
  charge = 0
  for residue in generate_residue_groups(hierarchy,
                                         assert_no_alt_loc=True,
                                         exclude_water=True,
                                        ):
    tmp = calculate_residue_charge(residue)
    charge += tmp
  return charge

def write_pdb_hierarchy_qxyz_file(hierarchy, file_name="qxyz_cctbx.dat",charge_scaling_positions=None,scale=0):
  qxyz_file = open(file_name,"w+")
  qxyz_file.write(str(hierarchy.atoms_size())+ "  \n")
  qxyz_file.write("  \n")
  qxyz = None
  for residue in generate_residue_groups(hierarchy,
                                         assert_no_alt_loc=True,
                                         exclude_water=True,
                                        ):
    if qxyz is None:
      qxyz = get_partial_point_charges(residue) 
    else: qxyz = qxyz + get_partial_point_charges(residue)
  scale_partial_point_charges(qxyz,charge_scaling_positions, scale=0)
  for item in qxyz:
    item_list = item + ["  \n"]
    item_string = "  ".join(str(elm) for elm in item_list)
    qxyz_file.write(item_string)
  qxyz_file.close()

def write_pdb_hierarchy_xyzq_file(hierarchy, file_name="xyzq_cctbx.dat", charge_scaling_positions=None, scale=0):
  qxyz = None
  xyzq_file = open(file_name,"w+")
  for residue in generate_residue_groups(hierarchy,
                                         assert_no_alt_loc=True,
                                         exclude_water=True,
                                        ):
    if qxyz is None:
      qxyz = get_partial_point_charges(residue)
    else: qxyz = qxyz + get_partial_point_charges(residue)
  scale_partial_point_charges(qxyz,charge_scaling_positions, scale=0)
  for item in qxyz:
    item_list = item[1:]+item[0:1] + ["  \n"]
    item_string = "  ".join(str(elm) for elm in item_list)
    xyzq_file.write(item_string)
  xyzq_file.close()

def scale_partial_point_charges(qxyz,charge_scaling_positions=None, scale=0):
  def partial_charge_in_charge_scaling_positions(partial_charge,charge_scaling_positions):
    scaling = False
    for xyz in charge_scaling_positions:
      same_point = abs(xyz[0] - partial_charge[1]) < 1.0E-3 and abs(xyz[1] - partial_charge[2]) < 1.0E-3 and abs(xyz[2] - partial_charge[3]) < 1.0E-3
      if same_point:
        scaling = True
        break
    return scaling 
  if charge_scaling_positions != None:
    for item in qxyz:
      if partial_charge_in_charge_scaling_positions(item,charge_scaling_positions):
        item[0] =  item[0]*scale

def run(pdb_filename):
  print "run",pdb_filename
  pdb_inp = iotbx.pdb.input(pdb_filename)
  hierarchy = pdb_inp.construct_hierarchy()
  complete_pdb_hierarchy(hierarchy)
  write_pdb_hierarchy_qxyz_file(hierarchy)
  total_charge = calculate_pdb_hierarchy_charge(hierarchy)
  print "total_charge",total_charge
  write_hierarchy(pdb_filename, pdb_inp, hierarchy, "complete")

if __name__=="__main__":
  args = sys.argv[1:]
  del sys.argv[1:]
  run(*tuple(args))
