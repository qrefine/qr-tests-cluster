import os, sys

import iotbx

from run_finalise import get_total_charge_from_file_name

def run(pdb_filename):
  print "run",pdb_filename
  total_charge = get_total_charge_from_file_name(pdb_filename,
                                                 verbose=True,
  )
  print "total_charge",total_charge

if __name__=="__main__":
  args = sys.argv[1:]
  del sys.argv[1:]
  run(*tuple(args))
