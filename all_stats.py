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

def get_bonds_angles_rmsd(restraints_manager, xrs):
  hd_sel = xrs.hd_selection()
  energies_sites = \
    restraints_manager.select(~hd_sel).energies_sites(
      sites_cart        = xrs.sites_cart().select(~hd_sel),
      compute_gradients = False)
  a_mean = energies_sites.angle_deviations()[2]
  b_mean = energies_sites.bond_deviations()[2]
  return a_mean, b_mean

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
  params.use_neutron_distances = True
  params.restraints_library.cdl = False
  processed_pdb_file = monomer_library.pdb_interpretation.process(
    mon_lib_srv        = mon_lib_srv,
    ener_lib           = ener_lib,
    raw_records        = raw_recs,
    params             = params,
    log                = null_out())
  xrs = processed_pdb_file.xray_structure()
  sctr_keys = xrs.scattering_type_registry().type_count_dict().keys()
  has_hd = "H" in sctr_keys or "D" in sctr_keys
  restraints_manager   = processed_pdb_file.geometry_restraints_manager(
    show_energies      = False,
    assume_hydrogens_all_missing = not has_hd,
    plain_pairs_radius = 5.0)
  a_mean, b_mean = get_bonds_angles_rmsd(
    restraints_manager = restraints_manager, xrs = xrs)
  energies_sites = \
    restraints_manager.energies_sites(
      sites_cart        = xrs.sites_cart(),
      compute_gradients = False)
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
    if(ma.info().label_string().count("F-obs")>0):
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
  print "PDB        | ORIGINAL                                                   | PHENIX.REFINE                                              \
                    | COMPLETED(no completed PDB structure found, all values 0)  | COMPLETED REFINED in Q|R                                              "
  print "           | RW     RF      BONDS   ANGL   CLSC  RAMA   ROTA   CB  MPSC | RW     RF      BONDS   ANGL   CLSC  RAMA   ROTA   CB  MPSC \
                    | RW     RF      BONDS   ANGL   CLSC  RAMA   ROTA   CB  MPSC | RW     RF      BONDS   ANGL   CLSC  RAMA   ROTA   CB  MPSC "
  for code in unique_codes:
    pdb_0 = "00/%s.pdb"%code
    pdb_1 = "01/%s_refine_001.pdb"%code
    pdb_2 = "02/%s_complete.pdb"%code
    pdb_3 = "03/%s_complete_refined.pdb"%code
    hkl   = "01/%s.mtz"%code
    #
    rw_0, rf_0 = get_r(pdb_file=pdb_0, mtz_file=hkl)
    rw_1, rf_1 = get_r(pdb_file=pdb_1, mtz_file=hkl)
    if(os.path.isfile(pdb_2)):
      rw_2, rf_2 = get_r(pdb_file=pdb_2, mtz_file=hkl)
      rw_3, rf_3 = get_r(pdb_file=pdb_3, mtz_file=hkl)
    else:
      rw_2, rf_2 = 0, 0
      rw_3, rf_3 = 0, 0
    #
    ms_0 = get_model_stat(file_name = pdb_0)
    ms_1 = get_model_stat(file_name = pdb_1)
    if (os.path.isfile(pdb_2)):
      ms_2 = get_model_stat(file_name = pdb_2)
      ms_3 = get_model_stat(file_name = pdb_3)
    else:
      ms_2 = group_args(
    b_mean                  = 0,
    a_mean                  = 0,
    number_of_worst_clashes = 0,
    ramachandran_outliers   = 0,
    rotamer_outliers        = 0,
    c_beta_dev              = 0,
    n_cis_proline           = 0,
    n_cis_general           = 0,
    n_twisted_proline       = 0,
    n_twisted_general       = 0,
    o                       = 0,
    b                       = 0,
    mpscore                 = 0,
    clsc                    = 0,
    n_atoms                 = 0)
      ms_3 = ms_2
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
    s2 = fmt%(
      rw_2, rf_2,
      ms_2.b_mean,
      ms_2.a_mean,
      ms_2.clsc,
      ms_2.ramachandran_outliers,
      ms_2.rotamer_outliers,
      ms_2.c_beta_dev,
      ms_2.mpscore)
    s3 = fmt%(
      rw_3, rf_3,
      ms_3.b_mean,
      ms_3.a_mean,
      ms_3.clsc,
      ms_3.ramachandran_outliers,
      ms_3.rotamer_outliers,
      ms_3.c_beta_dev,
      ms_3.mpscore)
    print "%s %5d"%(code, ms_0.n_atoms), "|", s0, "|", s1, "|", s2, "|", s3
    sys.stdout.flush()


if (__name__ == "__main__"):
  run()
