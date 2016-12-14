import os, sys

import iotbx

from run_finalise import get_total_charge_from_file_name

def run(pdb_filename):
  print "run",pdb_filename
  data = {}
  if os.path.exists(pdb_filename.replace('.pdb', '.psf')):
    f=file(pdb_filename.replace('.pdb', '.psf'), 'rb')
    lines = f.readlines()
    f.close()
    for line in lines:
      tmp = line.split()
      if not tmp: continue
      if tmp[1].find('PRO')==-1: continue
      print line
      key = '.'.join(tmp[1:4])
      print key
      data.setdefault(key, 0)
      data[key]+=float(tmp[6])
    print data

  total_charge = get_total_charge_from_file_name(pdb_filename,
                                                 check=data,
                                                 verbose=True,
  )
  print "total_charge",total_charge

if __name__=="__main__":
  args = sys.argv[1:]
  del sys.argv[1:]
  run(*tuple(args))
