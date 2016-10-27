import os, sys
from libtbx import easy_run
from iotbx import pdb

import run_finalise

pdbs = {"PRO_terminal" : """
CRYST1   42.664   49.718   66.065 109.25  94.96  99.24 P 1           4
ATOM    874  N   PRO A 115      -9.167  -7.159   4.783  1.00 30.39           N
ATOM    875  CA  PRO A 115      -9.350  -8.630   5.081  1.00 30.17           C
ATOM    876  C   PRO A 115      -9.382  -9.587   3.862  1.00 32.77           C
ATOM    877  O   PRO A 115     -10.069  -9.304   2.880  1.00 33.78           O
ATOM    878  CB  PRO A 115     -10.669  -8.662   5.838  1.00 28.12           C
ATOM    879  CG  PRO A 115     -10.665  -7.351   6.592  1.00 29.06           C
ATOM    880  CD  PRO A 115     -10.013  -6.329   5.660  1.00 30.28           C
ATOM    881  N  AGLY A 116      -8.745 -10.760   4.020  0.50 31.15           N
ATOM    882  CA AGLY A 116      -8.770 -11.833   3.012  0.50 27.99           C
ATOM    883  C  AGLY A 116      -7.959 -13.081   3.370  0.50 26.73           C
ATOM    884  O  AGLY A 116      -7.545 -13.840   2.487  0.50 26.23           O
ATOM    885  N  BGLY A 116      -8.621 -10.686   3.937  0.50 35.25           N
ATOM    886  CA BGLY A 116      -8.079 -11.393   2.744  0.50 36.56           C
ATOM    887  C  BGLY A 116      -8.989 -12.091   1.734  0.50 36.99           C
ATOM    888  O  BGLY A 116      -9.994 -11.547   1.293  0.50 37.57           O
ATOM    889  N  AGLY A 117      -7.748 -13.313   4.661  0.50 25.75           N
ATOM    890  CA AGLY A 117      -6.947 -14.457   5.099  0.50 24.75           C
ATOM    891  C  AGLY A 117      -5.660 -14.520   4.307  0.50 24.36           C
ATOM    892  O  AGLY A 117      -5.284 -13.541   3.660  0.50 23.17           O
ATOM    893  N  BGLY A 117      -8.582 -13.285   1.315  0.50 36.89           N
ATOM    894  CA BGLY A 117      -9.277 -14.011   0.264  0.50 35.40           C
ATOM    895  C  BGLY A 117      -8.318 -14.993  -0.358  0.50 37.53           C
ATOM    896  O  BGLY A 117      -8.109 -14.988  -1.568  0.50 39.46           O
ATOM    897  N  ASER A 118      -4.970 -15.660   4.366  0.50 25.06           N
ATOM    898  CA ASER A 118      -3.753 -15.855   3.578  0.50 25.99           C
ATOM    899  C  ASER A 118      -3.903 -15.046   2.311  0.50 28.77           C
ATOM    900  O  ASER A 118      -4.893 -15.203   1.584  0.50 30.82           O
ATOM    901  CB ASER A 118      -3.569 -17.326   3.226  0.50 24.41           C
ATOM    902  OG ASER A 118      -4.540 -17.742   2.278  0.50 25.05           O
ATOM    903  N  BSER A 118      -7.724 -15.839   0.481  0.50 38.64           N
ATOM    904  CA BSER A 118      -6.635 -16.713   0.055  0.50 38.40           C
ATOM    905  C  BSER A 118      -5.353 -15.909  -0.093  0.50 39.49           C
ATOM    906  O  BSER A 118      -4.263 -16.472  -0.176  0.50 39.26           O
ATOM    907  CB BSER A 118      -6.980 -17.441  -1.256  0.50 37.71           C
ATOM    908  OG BSER A 118      -8.020 -16.788  -1.972  0.50 34.61           O
""",
        }

def test_PRO_terminal_and_alt_loc():
  tf = 'PRO_terminal.pdb'
  f=file(tf, "wb")
  f.write(pdbs["PRO_terminal"])
  f.close()
  cmd = 'iotbx.python ../run_finalise.py %s' % tf
  print cmd
  easy_run.call(cmd)
  pdb_inp = pdb.input(tf.replace('.pdb', '_complete.pdb'))
  hierarchy = pdb_inp.construct_hierarchy()
  must_find = ['H2', 'H3', 'OXT']
  for atom in hierarchy.atoms():
    assert atom.name.strip()!='H1'
    if atom.name.strip() in must_find:
      must_find.remove(atom.name.strip())
  assert not must_find
  print 'OK'

def test_1yjp_charge():
  pdb_inp = pdb.input('1yjp.pdb')
  hierarchy = pdb_inp.construct_hierarchy()
  try:
    charge = run_finalise.calculate_pdb_hierarchy_charge(hierarchy)
    assert 0
  except Exception, e:
    assert e.message.find('no hydrogens')>-1
  print 'OK'
  tf='1yjp.pdb'
  cmd = 'iotbx.python ../run_finalise.py %s' % tf
  print cmd
  easy_run.call(cmd)
  pdb_inp = pdb.input(tf.replace('.pdb', '_complete.pdb'))
  hierarchy = pdb_inp.construct_hierarchy()
  charge = run_finalise.calculate_pdb_hierarchy_charge(hierarchy, verbose=1)
  print 'charge',charge
  assert charge==0, 'change of 1yjp should be zero not %s %s' % charge
  print 'OK'

def run():
  test_PRO_terminal_and_alt_loc()
  test_1yjp_charge()

if __name__=="__main__":
  args = sys.argv[1:]
  del sys.argv[1:]
  run(*tuple(args))
