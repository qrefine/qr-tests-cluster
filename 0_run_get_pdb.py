from __future__ import division
from libtbx import easy_run
import os

pdb_codes=['2jee','1fh5','3dtj','2oeq','1il5', '2pro', '3kyi', '3tz9', '4ctd', '4drw', '4k2r', '4rnf', '4tql', '4xa1', '5d12','1ok9', '1pag', '1u0d', '1va7', '1y1l', '2ghj', '2iwe', '2oy0', '3nm9', '4l21']
def run():
  for pc in pdb_codes:
    if not os.path.exists(pc+'.pdb'):
      print pc+'.pdb'
      cmd = "phenix.fetch_pdb %s --mtz"%pc.lower()
      print cmd
      try: easy_run.call(cmd)
      except Exception, e:
        print "  FAILED:", str(e)
   
  easy_run.call("rm -rf *.fa *.cif")
  easy_run.call("mv *.mtz ../0_mtz")  

if __name__ == "__main__":
  run()
