import os, sys
from scitbx.array_family import flex
import iotbx.pdb
from mmtbx import monomer_library
import mmtbx.monomer_library.server
import mmtbx.monomer_library.pdb_interpretation
from mmtbx import model_statistics
from cctbx import uctbx
from libtbx.utils import null_out
from mmtbx.validation.ramalyze import ramalyze
from mmtbx.validation.rotalyze import rotalyze
from mmtbx.validation.cbetadev import cbetadev
from mmtbx.validation.clashscore import clashscore
from mmtbx.validation.utils import molprobity_score
from mmtbx.validation import omegalyze
from libtbx import group_args
from iotbx import reflection_file_reader
from iotbx import reflection_file_utils
import iotbx.pdb
import mmtbx.f_model

mon_lib_srv = monomer_library.server.server()
ener_lib = monomer_library.server.ener_lib()
  
def all_single_atom_residues(ph):
  sizes = flex.int()
  for r in ph.residue_groups():
    sizes.append(r.atoms().size())
  fr = sizes.count(sizes[0])*100./sizes.size()
  if(fr>90.): return True
  else: return False
  
def get_model_stat(file_name):
  pdb_inp = iotbx.pdb.input(file_name=file_name)  
  atoms = pdb_inp.atoms()
  box = uctbx.non_crystallographic_unit_cell_with_the_sites_in_its_center(
    sites_cart   = atoms.extract_xyz(),
    buffer_layer = 5)
  atoms.set_xyz(new_xyz=box.sites_cart)
  ph = pdb_inp.construct_hierarchy()
  if(all_single_atom_residues(ph=ph)): return None
  raw_recs = ph.as_pdb_string(crystal_symmetry=box.crystal_symmetry()).splitlines()
  #  
  params = monomer_library.pdb_interpretation.master_params.extract()
  params.clash_guard.nonbonded_distance_threshold=None
  params.disable_uc_volume_vs_n_atoms_check=False
  processed_pdb_file = monomer_library.pdb_interpretation.process(
    mon_lib_srv        = mon_lib_srv,
    ener_lib           = ener_lib,
    raw_records        = raw_recs,
    params             = params,
    log                = null_out())
  restraints_manager   = processed_pdb_file.geometry_restraints_manager(
    show_energies      = False,
    plain_pairs_radius = 5.0)
  #
  sites_cart = processed_pdb_file.all_chain_proxies.pdb_hierarchy.atoms(
    ).extract_xyz()
  energies_sites = \
    restraints_manager.energies_sites(
      sites_cart        = sites_cart,
      compute_gradients = False)
  a_mean = energies_sites.angle_deviations()[2]
  b_mean = energies_sites.bond_deviations()[2]
  nonbonded_distances = energies_sites.nonbonded_distances()
  number_of_worst_clashes = (nonbonded_distances<0.5).count(True)
  #
  ramalyze_obj = ramalyze(pdb_hierarchy=ph, outliers_only=False)
  ramachandran_outliers = ramalyze_obj.percent_outliers
  rotamer_outliers = rotalyze(
    pdb_hierarchy=ph, outliers_only=False).percent_outliers
  c_beta_dev = cbetadev(
    pdb_hierarchy = ph,
    outliers_only = True,
    out           = null_out()).get_outlier_count()
  omglz = omegalyze.omegalyze(pdb_hierarchy=ph, quiet=True)
  n_cis_proline     = omglz.n_cis_proline()
  n_cis_general     = omglz.n_cis_general()
  n_twisted_proline = omglz.n_twisted_proline()
  n_twisted_general = omglz.n_twisted_general()
  #
  clsc = clashscore(pdb_hierarchy=ph).get_clashscore()
  mpscore = molprobity_score(
    clashscore = clsc,
    rota_out   = rotamer_outliers,
    rama_fav   = ramalyze_obj.percent_favored)
  #
  occ = atoms.extract_occ()
  bs = atoms.extract_b()
  #
  return group_args(
    b_mean                  = b_mean, 
    a_mean                  = a_mean, 
    number_of_worst_clashes = number_of_worst_clashes, 
    ramachandran_outliers   = ramachandran_outliers, 
    rotamer_outliers        = rotamer_outliers, 
    c_beta_dev              = c_beta_dev, 
    n_cis_proline           = n_cis_proline,
    n_cis_general           = n_cis_general,
    n_twisted_proline       = n_twisted_proline,
    n_twisted_general       = n_twisted_general,
    o                       = occ.min_max_mean().as_tuple(),
    b                       = bs.min_max_mean().as_tuple(),
    mpscore                 = mpscore,
    clsc                    = clsc,
    n_atoms                 = atoms.size())
            
def get_r(pdb_file, mtz_file):
  miller_arrays = reflection_file_reader.any_reflection_file(file_name =
    mtz_file).as_miller_arrays()
  f_obs, r_free_flags = None,None
  for ma in miller_arrays:
    if(ma.info().label_string().count("F-obs-filtered")>0):
      f_obs = ma.deep_copy()
      merged = f_obs.as_non_anomalous_array().merge_equivalents()
      f_obs = merged.array().set_observation_type(f_obs)
    if(ma.info().label_string().count("R-free-flags")):
      r_free_flags = ma.deep_copy()
      merged = r_free_flags.as_non_anomalous_array().merge_equivalents()
      r_free_flags = merged.array().set_observation_type(r_free_flags)
  [f_obs, r_free_flags].count(None) == 0
  fr = r_free_flags.data().count(1)*100./r_free_flags.data().size()
  if(fr<15.):
    r_free_flags = r_free_flags.array(data = r_free_flags.data()==1)
  else:
    r_free_flags = r_free_flags.array(data = r_free_flags.data()==0)
  f_obs, r_free_flags = f_obs.common_sets(r_free_flags)
  xrs = iotbx.pdb.input(file_name = pdb_file).xray_structure_simple()
  fmodel = mmtbx.f_model.manager(
    f_obs          = f_obs,
    r_free_flags   = r_free_flags,
    xray_structure = xrs)
  fmodel.update_all_scales()
  return fmodel.r_work(), fmodel.r_free()
          
def run():
  unique_codes = []
  for f in os.listdir("01/"):
    if(f.endswith(".pdb")):
      if(not f[:4] in unique_codes):
        unique_codes.append(f[:4])
  assert len(unique_codes)==26
  #
  print "PDB        | ORIGINAL                                                   | PHENIX.REFINE"
  print "           | RW     RF      BONDS   ANGL   CLSC  RAMA   ROTA   CB  MPSC | RW     RF      BONDS   ANGL   CLSC  RAMA   ROTA   CB  MPSC"
  for code in unique_codes:
    pdb_0 = "00/%s.pdb"%code
    pdb_1 = "01/%s_refine_001.pdb"%code
    hkl   = "01/%s_refine_001.mtz"%code
    #
    rw_0, rf_0 = get_r(pdb_file=pdb_0, mtz_file=hkl)
    rw_1, rf_1 = get_r(pdb_file=pdb_1, mtz_file=hkl)
    #
    ms_0 = get_model_stat(file_name = pdb_0)
    ms_1 = get_model_stat(file_name = pdb_1)
    #
    fmt="%6.4f %6.4f %6.3f %6.2f %6.2f %5.2f %6.2f %4d %5.2f"
    s0 = fmt%(
      rw_0, rf_0,
      ms_0.b_mean, 
      ms_0.a_mean, 
      ms_0.clsc, 
      ms_0.ramachandran_outliers, 
      ms_0.rotamer_outliers, 
      ms_0.c_beta_dev,
      ms_0.mpscore)
    s1 = fmt%(
      rw_1, rf_1,
      ms_1.b_mean, 
      ms_1.a_mean, 
      ms_1.clsc, 
      ms_1.ramachandran_outliers, 
      ms_1.rotamer_outliers, 
      ms_1.c_beta_dev,
      ms_1.mpscore)
    print "%s %5d"%(code, ms_0.n_atoms), "|", s0, "|", s1
    sys.stdout.flush() 
    
   
if (__name__ == "__main__"):
  run()