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


            
def run_one(mtz_file):
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
  f_obs, r_free_flags = f_obs.common_sets(r_free_flags)
  mtz_dataset = f_obs.as_mtz_dataset(column_root_label="F-obs")
  mtz_dataset.add_miller_array(
    miller_array      = r_free_flags,
    column_root_label = "R-free-flags")
  mtz_object = mtz_dataset.mtz_object()
  mtz_object.write(file_name = "%s.mtz"%mtz_file[3:7])
          
def run():
  unique_codes = []
  for f in os.listdir("01/"):
    if(f.endswith(".mtz")):
      print "01/"+f
      run_one(mtz_file="01/"+f)

    
   
if (__name__ == "__main__"):
  run()
