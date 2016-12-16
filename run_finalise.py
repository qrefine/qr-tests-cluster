# LIBTBX_SET_DISPATCHER_NAME phenix.development.ready_set
import os, sys
import math
from string import letters

import iotbx
import iotbx.pdb
from libtbx.utils import Sorry
from scitbx import matrix
from libtbx import easy_run

from cctbx.array_family import flex

from iotbx.pdb.atom_name_interpretation import \
  interpreters as protein_atom_name_interpreters
get_class = iotbx.pdb.common_residue_names_get_class

from mmtbx.chemical_components import get_cif_dictionary

from mmtbx.monomer_library import server
mon_lib_server = server.server()

from mmtbx import monomer_library
from mmtbx.monomer_library import pdb_interpretation
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

default_ion_charges = {
  "PT" : 2,
  "CL" : -1,
  "CD" : 2,
  "ZN" : 2,
  }
allowable_amino_acid_charges = {
  "ARG" : 1,
  "CYS" : 0,
  }

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

def _add_atom_to_residue_group(atom, ag):
  tag = iotbx.pdb.hierarchy.atom_group()
  tag.resname = ag.resname
  tag.append_atom(atom)
  rg = iotbx.pdb.hierarchy.residue_group()
  rg.resseq = ag.parent().resseq
  rg.append_atom_group(tag)
  for i, c in enumerate(letters):
    if c==ag.parent().parent().id:
      break
  atom.tmp = i
  for atom in rg.atoms(): print atom.format_atom_record()
  return rg

def add_n_terminal_hydrogens_to_atom_group(ag,
                                           use_capping_hydrogens=False,
                                           append_to_end_of_model=False,
                                          ):
  rc=[]
  n = ag.get_atom("N")
  if n is None: return
  ca = ag.get_atom("CA")
  if ca is None: return
  c = ag.get_atom("C")
  if c is None: return
  if ag.get_atom("H"): # maybe needs to be smarter or actually work
    #for atom in ag.atoms(): print atom.quote()
    ag.remove_atom(ag.get_atom('H'))
  # add H1
  rh3 = construct_xyz(n, 0.9,
                      ca, 109.5,
                      c, 120.,
                     )
  possible = ['H', 'H1', 'H2', 'H3']
  h_count = 0
  for h in possible:
    if ag.get_atom(h): h_count+=1
  number_of_hydrogens=3
  if use_capping_hydrogens:
    number_of_hydrogens-=1
  print 'h_count',h_count,number_of_hydrogens
  if h_count>=number_of_hydrogens: return []
  for i in range(0, number_of_hydrogens):
    name = " H%d " % (i+1)
    if ag.get_atom(name.strip()): continue
    if ag.resname=='PRO':
      if i==0:
        continue
    atom = iotbx.pdb.hierarchy.atom()
    atom.name = name
    atom.element = "H"
    atom.xyz = rh3[i]
    atom.occ = n.occ
    atom.b = n.b
    if append_to_end_of_model and i+1==number_of_hydrogens:
      rg = _add_atom_to_residue_group(atom, ag)
      rc.append(rg)
    else:
      ag.append_atom(atom)
  #for atom in rc: print atom.quote()
  return rc

def add_n_terminal_hydrogens_to_residue_group(rg,
                                              use_capping_hydrogens=False,
                                              append_to_end_of_model=False,
                                             ):
  rc=[]
  for ag in rg.atom_groups():
    rc += add_n_terminal_hydrogens_to_atom_group(
      ag,
      use_capping_hydrogens=use_capping_hydrogens,
      append_to_end_of_model=append_to_end_of_model,
    )
  return rc

def add_n_terminal_hydrogens(hierarchy,
                             #residue_selection=None,
                             add_to_chain_breaks=False,
                            ):
  assert 0
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

def add_c_terminal_oxygens_to_atom_group(ag,
                                         use_capping_hydrogens=False,
                                         append_to_end_of_model=False,
                                        ):
  #
  # do we need ANISOU
  #
  rc = []
  atom_name=' OXT'
  atom_element = 'O'
  bond_length=1.231
  if use_capping_hydrogens:
    if ag.get_atom(atom_name.strip()): return []
    atom_name=" HC "
    atom_element="H"
    bond_length=1.
  if ag.get_atom(atom_name.strip()): return []
  c = ag.get_atom("C")
  if c is None: return
  ca = ag.get_atom("CA")
  if ca is None: return
  n = ag.get_atom("N")
  if n is None: return
  ro2 = construct_xyz(c, bond_length,
                      ca, 120.,
                      n, 160.,
                      period=2,
                     )
  oxys = [' O  ', atom_name] 
  for i in range(0,2):
    name = oxys[i]
    atom = ag.get_atom(name.strip())
    if atom:
      atom.xyz = ro2[i]
    else:
      atom = iotbx.pdb.hierarchy.atom()
      atom.name = name
      atom.element = atom_element
      atom.occ = c.occ
      atom.b = c.b
      atom.xyz = ro2[i]
      if append_to_end_of_model:
        rg = _add_atom_to_residue_group(atom, ag)
        rc.append(rg)
      else:
        ag.append_atom(atom)
  return rc

def add_c_terminal_oxygens_to_residue_group(rg,
                                            use_capping_hydrogens=False,
                                            append_to_end_of_model=False,
                                          ):
  rc=[]
  for ag in rg.atom_groups():
    rc += add_c_terminal_oxygens_to_atom_group(
      ag,
      use_capping_hydrogens=use_capping_hydrogens,
      append_to_end_of_model=append_to_end_of_model,
    )
  return rc

def add_c_terminal_oxygens(hierarchy,
                          ):
  assert 0
  for chain_i, chain in enumerate(hierarchy.chains()):
    for res_i, residue_group in enumerate(chain.residue_groups()):
      if len(residue_group.atom_groups())>1: continue
      atom_group = residue_group.atom_groups()[0]
      if get_class(atom_group.resname) not in ["common_amino_acid",
                                               "modified_amino_acid",
                                             ]:
        continue
      if capping_hydrogens:
        assert 0
      if res_i==len(chain.residue_groups())-1: # need better switch
        add_c_terminal_oxygens_to_atom_group(atom_group)
  hierarchy.atoms_reset_serial()
  hierarchy.atoms().reset_i_seq()
  return hierarchy

def add_cys_hg_to_atom_group(ag,
                             append_to_end_of_model=False,
                            ):
  #
  # do we need ANISOU
  #
  rc = []
  atom_name=' HG '
  atom_element = 'H'
  bond_length=1.
  if ag.get_atom(atom_name.strip()): return []
  sg = ag.get_atom("SG")
  if sg is None: return
  cb = ag.get_atom("CB")
  if cb is None: return
  ca = ag.get_atom("CA")
  if ca is None: return
  ro2 = construct_xyz(sg, bond_length,
                      cb, 120.,
                      ca, 160.,
                      period=1,
                     )
  atom = iotbx.pdb.hierarchy.atom()
  atom.name = atom_name
  atom.element = atom_element
  atom.occ = sg.occ
  atom.b = sg.b
  atom.xyz = ro2[0]
  if append_to_end_of_model:
    rg = _add_atom_to_residue_group(atom, ag)
    rc.append(rg)
  else:
    ag.append_atom(atom)
  return rc

def add_cys_hg_to_residue_group(rg,
                                append_to_end_of_model=False,
                               ):
  rc=[]
  for ag in rg.atom_groups():
    if ag.resname not in ['CYS']: continue
    rc += add_cys_hg_to_atom_group(
      ag,
      append_to_end_of_model=append_to_end_of_model,
    )
  return rc

def generate_residues_via_conformer(hierarchy,
                                    backbone_only=False,
                                    verbose=False,
                                    ):
  backbone_asc = hierarchy.atom_selection_cache()
  backbone_sel = backbone_asc.selection("name ca or name c or name n or name o or name cb")
  backbone_hierarchy = hierarchy.select(backbone_sel)
  get_class = iotbx.pdb.common_residue_names_get_class
  loop_hierarchy=hierarchy
  if backbone_only: loop_hierarchy=backbone_hierarchy
  for model in loop_hierarchy.models():
    if verbose: print 'model: "%s"' % model.id
    for chain in model.chains():
      if verbose: print 'chain: "%s"' % chain.id
      for conformer in chain.conformers():
        if verbose: print '  conformer: altloc="%s"' % (
          conformer.altloc)
#        while threes: del threes[0]
#        threes.start=None
#        threes.end=None
#        list_of_threes = []
        for residue in conformer.residues():
          if verbose:
            if residue.resname not in ["HOH"]:
              print '    residue: resname="%s" resid="%s"' % (
                residue.resname, residue.resid())
          if verbose: print '      residue class : %s' % get_class(residue.resname)
          if get_class(residue.resname) not in ["common_amino_acid"]:
            continue
          yield residue

def generate_protein_fragments(hierarchy,
                               geometry,
                               backbone_only=False,
                               verbose=False,
                               ):
  from mmtbx.conformation_dependent_library.multi_residue_class import \
    ThreeProteinResidues, RestraintsRegistry
  registry = RestraintsRegistry()
  threes = ThreeProteinResidues(geometry, registry=registry)
  for residue in generate_residues_via_conformer(hierarchy,
                                                 backbone_only=backbone_only,
                                                 verbose=verbose,
                                                 ):
    list.append(threes, residue)
    sub_unit = threes.provide_second_sub_unit_if_unlinked()
    if sub_unit:
      threes.start = True
      threes.end = True
      yield threes
      threes = sub_unit
  threes.start = True
  threes.end = True
  yield threes

def add_terminal_hydrogens(
    hierarchy,
    geometry_restraints_manager,
    add_to_chain_breaks=False,
    use_capping_hydrogens=False,  # instead of terminal H
    append_to_end_of_model=False, # useful for Q|R
    #use_capping_only_on_chain_breaks=False,
    ):
  # add N terminal hydrogens because Reduce only does it to resseq=1
  # needs to be alt.loc. aware for non-quantum-refine
  atoms = hierarchy.atoms()
  def get_residue_group(residue):
    for atom in residue.atoms():
      atom = atoms[atom.i_seq]
      break
    return atom.parent().parent()
  ###
  additional_hydrogens=[]
  for three in generate_protein_fragments(hierarchy,
                                          geometry_restraints_manager,
                                          backbone_only=False,
                                         ):
    assert three.are_linked()
    for i in range(len(three)):
      rg = get_residue_group(three[i])
      add_cys_hg_to_residue_group(rg)
    if three.start:
      print 'start'*10,three
      rg = get_residue_group(three[0])
      rc = add_n_terminal_hydrogens_to_residue_group(
        rg,
        use_capping_hydrogens=use_capping_hydrogens,
        append_to_end_of_model=append_to_end_of_model,
      )
      if rc: additional_hydrogens.append(rc)
      #hierarchy.reset_i_seq_if_necessary()
    if three.end:
      print 'end'*10,three
      rg = get_residue_group(three[-1])
      rc = add_c_terminal_oxygens_to_residue_group(
        rg,
        use_capping_hydrogens=use_capping_hydrogens,
        append_to_end_of_model=append_to_end_of_model,
      )
      if rc: additional_hydrogens.append(rc)
      #hierarchy.reset_i_seq_if_necessary()
    else:
      pass

  if append_to_end_of_model and additional_hydrogens:
    tmp = []
    for group in additional_hydrogens:
      for atom in group:
        tmp.append(atom)
    _add_atoms_to_end_of_hierarchy(hierarchy, tmp)

def _add_atoms_to_end_of_hierarchy(hierarchy, rgs):
  chains = {}
  for rg in rgs:
    for atom in rg.atoms():
      cid = atom.tmp
      if cid not in chains:
        chains[cid] = iotbx.pdb.hierarchy.chain()
        chains[cid].id = letters[cid]
      chains[cid].append_residue_group(rg)
  model = hierarchy.models()[0]
  for i, chain in sorted(chains.items()):
    model.append_chain(chain)

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
                             hetero_charges=None,
                             inter_residue_bonds=None,
                             verbose=False,
                             ):
  def _terminal(names, check):
    for name in check:
      if name not in names:
        break
    else:
      return True
    return False
  def n_terminal(residue_name, atom_names):
    if residue_name in ["PRO"]:
      check_names = [' H2 ',' H3 ']
    else:
      check_names = [' H1 ',' H2 ',' H3 ']
    rc = _terminal(atom_names, check_names)
    if not rc and residue_name in ['PRO']: # CHARMM...
      check_names = [' H 1', ' H 2']
      rc = _terminal(atom_names, check_names)
    if not rc and residue_name in ['PRO']: # BABEL...
      check_names = [' HN1', ' HN2']
      rc = _terminal(atom_names, check_names)
    return rc
  def nh2_terminal(atom_names):
    return _terminal(atom_names, [' HT1', ' HT2'])
  def nh3_terminal(atom_names):
    return _terminal(atom_names, [' HT1', ' HT2', ' HT3'])
  def c_terminal(atom_names):
    rc = _terminal(atom_names, [' OXT'])
    if not rc: rc = _terminal(atom_names, [' OT1', ' OT2'])
    return rc
  def covalent_bond(i_seqs, inter_residue_bonds):
    for i_seq in i_seqs:
      if i_seq in inter_residue_bonds:
        return True
    return False

  max_charge=1
  if assert_no_alt_loc:
    if len(rg.atom_groups())>1:
      raise Sorry("alt locs in %s" % display_residue_group(rg))
  # ions
  # needs to be centralised!!!
  resname = rg.atom_groups()[0].resname
  if get_class(resname)=="common_element":
    atom = rg.atoms()[0]
    if not atom.charge.strip():
      if hetero_charges:
        charge = hetero_charges.get( atom.parent().resname.strip(), None)
        if charge is None:
          raise Sorry('no charge found in the model file or hetero_charges for "%s"' % atom.quote())
        else:
          return charge
      else:
        raise Sorry('no charge found in the model file for "%s"' % atom.quote())
    else:
      return atom.charge_as_int()
  # others
  hs=0
  atom_names = []
  atom_i_seqs = []
  for atom in rg.atoms():
    if verbose: print atom.quote()
    if atom.element_is_hydrogen(): hs+=1
    atom_names.append(atom.name)
    atom_i_seqs.append(atom.i_seq)
  if verbose: print get_class(resname)
  if assert_contains_hydrogens:
    if hs==0:
      hydrogens = get_aa_polymer_hydrogens(resname)
      if len(hydrogens)!=0:
        if verbose:
          for atom in rg.atoms(): print 'H',atom.quote()
        raise Sorry("no hydrogens: %s" % display_residue_group(rg))
  ag = rg.atom_groups()[0]
  charge = get_aa_charge(ag.resname)
  rc = get_aa_charge(ag.resname)
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
    if verbose: print atom_names
    if n_terminal(ag.resname, atom_names):
      diff_hs-=1
      if verbose:
        print 'n_terminal'
        print 'charge: %s poly_hs: %s diff_hs: %s' % (charge,
                                                              poly_hs,
                                                              diff_hs,
                                                            )
    elif nh3_terminal(atom_names):
      diff_hs-=1
      if verbose: print 'nh3_terminal'
    elif nh2_terminal(atom_names):
      diff_hs-=1
      if verbose: print 'nh2_terminal'
    elif c_terminal(atom_names):
      diff_hs-=1
      max_charge+=1
      if verbose:
        print 'c_terminal'
        print 'charge: %s poly_hs: %s diff_hs: %s' % (charge,
                                                    poly_hs,
                                                    diff_hs,
      )
    if covalent_bond(atom_i_seqs, inter_residue_bonds):
      diff_hs+=1
    charge+=diff_hs
    if verbose:
      print 'charge: %s poly_hs: %s diff_hs: %s' % (charge,
                                                    poly_hs,
                                                    diff_hs,
      )

    assert abs(diff_hs)<=max_charge, 'residue %s charge %s is greater than %s' % (
      rg.atoms()[0].quote(),
      diff_hs,
      max_charge,
    )
    if resname in allowable_amino_acid_charges:
      assert allowable_amino_acid_charges[resname]-1 <= charge <= allowable_amino_acid_charges[resname]+1
  else:
    print "NOT POLY",ag.resname, get_class(ag.resname),
    print charge
    
    if 0 and ag.resname not in ["NO3",
                          "PXZ",
                          "SRT",
                          "HM7",
                        ]:
      assert 0
  return charge, rc

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

def get_partial_point_charges(rg, hetero_charges=None):
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
          if hetero_charges:
            key = atom.parent().resname
            print 'using hetero_charges for :%s' % key
            charge = hetero_charges.get(key.strip(), None)
            if charge:
              tmp.append([charge]+list(atom.xyz))
            else:
              raise Sorry('no charge found in the model file or hetero_charges for "%s"' % atom.quote())
          else:
            raise Sorry('no charge found in the model file for "%s"' % atom.quote())
        else:
          tmp.append([atom.charge_as_int()]+list(atom.xyz))
          continue
      # other atoms
      cif = atom_dict.get(atom.name.strip(), None)
      if cif is None:
        if atom.name in [" H1 ", " H2 ", " H3 "]: # needs calculating...
          tmp.append([0.26]+list(atom.xyz))
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

def get_processed_pdb(pdb_filename=None,
                      raw_records=None,
                      pdb_inp=None,
                    ):
  # need to add params
  mon_lib_srv = monomer_library.server.server()
  ener_lib = monomer_library.server.ener_lib()
  ppf = monomer_library.pdb_interpretation.process(
    mon_lib_srv           = mon_lib_srv,
    ener_lib              = ener_lib,
    file_name             = pdb_filename,
    pdb_inp               = pdb_inp,
    raw_records           = raw_records,
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
  return output

def complete_pdb_hierarchy(hierarchy,
                           geometry_restraints_manager,
                           use_capping_hydrogens=False,
                           append_to_end_of_model=False,
                          ):
  for ag in hierarchy.atom_groups():
    if get_class(ag.resname) in ['common_rna_dna']:
      raise Sorry('')
  from mmtbx.building import extend_sidechains
  n_changed = extend_sidechains.extend_protein_model(hierarchy,
                                                     hydrogens=True,
                                                   )
  remove_acid_side_chain_hydrogens(hierarchy)
  add_terminal_hydrogens(hierarchy,
                         geometry_restraints_manager,
                         use_capping_hydrogens=use_capping_hydrogens,
                         append_to_end_of_model=append_to_end_of_model,
                        ) # in place
  hierarchy.atoms().set_chemical_element_simple_if_necessary()
  hierarchy.sort_atoms_in_place()

def calculate_pdb_hierarchy_charge(hierarchy,
                                   hetero_charges=None,
                                   inter_residue_bonds=None,
                                   check=None,
                                   verbose=False,
                                   ):
  charge = 0
  if inter_residue_bonds is None: inter_residue_bonds=[]
  for residue in generate_residue_groups(hierarchy,
                                         assert_no_alt_loc=True,
                                         exclude_water=True,
                                        ):
    tmp, rc = calculate_residue_charge(residue,
                                       hetero_charges=hetero_charges,
                                       inter_residue_bonds=inter_residue_bonds,
                                       verbose=verbose,
                                     )
    if check:
      print residue
      key = 'PRO%s' % residue.parent().id
      key += '.%s' % residue.resseq.strip()
      key += '.%s' % residue.atom_groups()[0].resname
      key = key.replace('HIS', 'HSD')
      print key
      print ' CHARMM %f Phenix %f' % (check[key], tmp)
      if 1:
        print inter_residue_bonds
      assert abs(check[key]-tmp)<0.001
    charge += tmp
    if verbose:
      if tmp:
        outl = '-'*80
        outl += '\nNON-ZERO CHARGE current %2d base %2d total %2d' % (tmp,rc,charge)
        for i, atom in enumerate(residue.atoms()):
          if i==0:
            outl += ' "%s"' % (atom.parent().id_str())
            if tmp!=rc: outl += ' DIFF %2d' % (tmp-rc)
          assert abs(tmp-rc)<=2, outl
          outl += "\n%s" % atom.quote()
        outl += '\n%s' % ('-'*80)
        print outl
  return charge

def get_total_charge_from_pdb(pdb_filename=None,
		              raw_records=None,
                              pdb_inp=None,
                              hetero_charges=None,
                              inter_residue_bonds=None,
                              verbose=False,
                              check = None,
                                   ):
  from run_finalise import get_processed_pdb
  from run_finalise import get_hetero_charges, default_ion_charges
  from run_finalise import get_inter_residue_bonds
  from run_finalise import calculate_pdb_hierarchy_charge
  ppf = get_processed_pdb(pdb_filename=pdb_filename,
                          raw_records=raw_records,
                          pdb_inp=pdb_inp,
                        )
  pdb_inp = ppf.all_chain_proxies.pdb_inp
  cs = pdb_inp.crystal_symmetry_from_cryst1()
  assert cs, 'There is no CRYST1 record in the input file'
  if not hetero_charges:
    hetero_charges = get_hetero_charges(ppf.all_chain_proxies.pdb_inp)
    if not hetero_charges:
      # some defaults
      hetero_charges = default_ion_charges
  if not inter_residue_bonds:
    inter_residue_bonds = get_inter_residue_bonds(ppf)
  if verbose:
    for key, item in inter_residue_bonds.items():
      if type(key)!=type(0) and len(key)==2: print key, item
  total_charge = calculate_pdb_hierarchy_charge(
    ppf.all_chain_proxies.pdb_hierarchy,
    hetero_charges=hetero_charges,
    inter_residue_bonds=inter_residue_bonds,
    check=check,
    verbose=verbose,
  )
  return total_charge

def write_charge_and_coordinates_from_hierarchy(hierarchy,
                                                file_name,
                                                qxyz_order='qxyz',
                                                hetero_charges=None,
                                                exclude_water=True,
                                                charge_scaling_positions=None,
                                                scale=0
                                                ):
  qxyz = None
  for residue in generate_residue_groups(hierarchy,
                                         assert_no_alt_loc=True,
                                         exclude_water=exclude_water,
                                        ):
    if qxyz is None:
      qxyz = get_partial_point_charges(residue, hetero_charges=hetero_charges) 
    else: 
      qxyz = qxyz + get_partial_point_charges(residue,
                                              hetero_charges=hetero_charges)
  if qxyz is None: return
  scale_partial_point_charges(qxyz,charge_scaling_positions, scale=0)
  qxyz_file = open(file_name,"w+")
  if qxyz_order=='qxyz': # tetrachem?
    qxyz_file.write(str(hierarchy.atoms_size())+ "  \n")
    qxyz_file.write("  \n")
  elif qxyz_order=='xyzq':
    pass
  else:
    raise Sorry('invalid qxyz_order parameter "%s"' % qxyz_order)
  for item in qxyz:
    if qxyz_order=='qxyz':
      item_list = item + ["  \n"]
    elif qxyz_order=='xyzq':
      item_list = item[1:]+item[0:1] + ["  \n"]
    else:
      raise Sorry('invalid qxyz_order parameter "%s"' % qxyz_order)
    item_string = "  ".join(str(elm) for elm in item_list)
    qxyz_file.write(item_string)
  qxyz_file.close()

def write_pdb_hierarchy_qxyz_file(hierarchy,
                                  file_name="qxyz_cctbx.dat",
                                  hetero_charges=None,
                                  exclude_water=True,
                                  charge_scaling_positions=None,
                                  scale=0,
                                 ):
  write_charge_and_coordinates_from_hierarchy(hierarchy,
                                              file_name=file_name,
                                              qxyz_order='qxyz',
                                              hetero_charges=hetero_charges,
                                              exclude_water=exclude_water,
                                              charge_scaling_positions=None,
	                                      scale=0,
                                              )

def write_pdb_hierarchy_xyzq_file(hierarchy,
                                  file_name="xyzq_cctbx.dat",
                                  hetero_charges=None,
                                  exclude_water=True,
                                  charge_scaling_positions=None,
                                  scale=0,
                                  ):
  write_charge_and_coordinates_from_hierarchy(hierarchy,
                                              file_name=file_name,
                                              qxyz_order='xyzq',
                                              hetero_charges=hetero_charges,
                                              exclude_water=exclude_water,
                                              charge_scaling_positions=None,
                                              scale=0,
                                              )
def scale_partial_point_charges(qxyz,
                                charge_scaling_positions=None,
                                scale=0):
  def partial_charge_in_charge_scaling_positions(partial_charge,
                                                 charge_scaling_positions):
    scaling = False
    for xyz in charge_scaling_positions:
      same_point = abs(xyz[0] - partial_charge[1]) < 1.0E-3 and abs(xyz[1] - partial_charge[2]) < 1.0E-3 and abs(xyz[2] - partial_charge[3]) < 1.0E-3
      if same_point:
        scaling = True
        break
    return scaling
  #####
  if charge_scaling_positions != None:
    for item in qxyz:
      if partial_charge_in_charge_scaling_positions(item,charge_scaling_positions):
        item[0] =  item[0]*scale

def get_hetero_charges(pdb_inp):
  # get the hetero charges from the FORMUL record
  hetero_charges = {}
  for line in pdb_inp.heterogen_section():
    if line.find("FORMUL")==-1: continue
    tmp = line.split()
    hetero_charges.setdefault(tmp[2], 0)
    charge = tmp[-1].replace(")", '')
    sign = None
    if charge.find('-')>-1: sign=-1
    if charge.find('+')>-1: sign=1
    if sign:
      charge=charge.replace("+", '').replace('-','')
      charge = int(charge)
      hetero_charges[tmp[2]]=charge*sign
  return hetero_charges

def get_inter_residue_bonds(ppf):
  # must use this before changing the hierarchy
  # used for S-S bonds etc...
  inter_residue_bonds = {}
  grm = ppf.geometry_restraints_manager()
  if not hasattr(grm, 'get_all_bond_proxies'): return inter_residue_bonds
  atoms = ppf.all_chain_proxies.pdb_hierarchy.atoms()
  for bond in grm.get_all_bond_proxies():
    if not hasattr(bond, 'get_proxies_with_origin_id'): continue
    for p in bond.get_proxies_with_origin_id():
      assert p.origin_id==0
      r1 = atoms[p.i_seqs[0]]
      r2 = atoms[p.i_seqs[1]]
      # exclude peptide links
      # but maybe should include all for completeness
      if r1.name.strip()=='C' and r2.name.strip()=='N': continue
      #outl = 'bonding %s %s' % ( r1.quote(),r2.quote())
      r1=r1.parent().parent()
      r2=r2.parent().parent()
      if r1.id_str()!=r2.id_str():
        inter_residue_bonds[p.i_seqs] = True
        for i in range(2):
          inter_residue_bonds.setdefault(p.i_seqs[i], [])
          inter_residue_bonds[p.i_seqs[i]].append(inter_residue_bonds[p.i_seqs])
  return inter_residue_bonds

def run_ready_set(pdb_filename):
  from StringIO import StringIO
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
