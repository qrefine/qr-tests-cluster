import os, sys, shutil
from libtbx import easy_run

from multiprocessing import Pool
pdb_dir = './'

def callback(args):
  print args
  return args
  
def _process_pdb_filename(pdb_file):
  updated_file=pdb_file[:-4]+".updated.pdb"
  complete_file=pdb_file[:-4]+"_complete.pdb"
  if ( pdb_file.endswith("pdb")  and not os.path.exists(complete_file) ):
   # pdb_file=pdb_dir+pdb_file
    print pdb_file
    # removing alt loc done in run_finalise.py
    #cmd = "phenix.pdbtools %(pdb_file_path)s" % locals()
    #cmd +=' keep="%s"' % "altloc ' ' or altloc 'A'"
    #cmd +=' stop_for_unknowns=False occupancies.set=1'
    # maybe check for non-unit occ.
    #print cmd
    #easy_run.call(cmd)
    modified_file="%(pdb_file)s_modified.pdb" % locals()
    shutil.copyfile(pdb_file, modified_file)
    cmd ="phenix.reduce -trim  "+modified_file +" >  "+pdb_file
    easy_run.call(cmd)
    if os.path.exists(modified_file): os.remove(modified_file)
    cmd ="phenix.ready_set %(pdb_file)s  add_h_to_water=true" % locals()
 #   cmd += ' ligand_cache_directory=%s' % os.path.abspath(os.path.join(
 #             "..",
 #             "00_ligands_cif",
 #             ))
    print cmd
    easy_run.call(cmd)
    shutil.copyfile(updated_file, pdb_file)

    cmd = "phenix.python ../run_finalise.py %s" % pdb_file + "> "+ pdb_file[:-4]+".log"
    print '\n\t~> %s\n' % cmd
    easy_run.call(cmd)
    if not os.path.exists(complete_file):
      print 'run_finalise failed'
      cmd = "rm   "+pdb_file[:-4]+"*pdb"
      print cmd
      os.system("rm   "+pdb_file[:-4]+"*pdb")
    return pdb_file
  return None

def run(
    nproc=8,
    only_code=None, 
    ):
  os.system("rm *")
  cmd = "cp ../01/*.pdb ./"
  print cmd 
  os.system(cmd)
  try: nproc=int(nproc)
  except: nproc=1  
  pool = None
  if nproc>1:
    pool = Pool(processes=nproc)

  pdb_files = os.listdir(pdb_dir)
  for pdb_file in pdb_files:
    cmd='mv  '+ pdb_file +'  ' + pdb_file[:4]+'.pdb'
    os.system(cmd)
  pdb_files = os.listdir(pdb_dir)
  for pdb_file in pdb_files:
    if not pdb_file.endswith(".pdb"): continue
    if len(pdb_file.split('.'))!=2: continue
    if len(pdb_file.split('.')[0])!=4: continue
    only_code = pdb_file.split('.')[0]
#    if only_code is not None and pdb_file.find(only_code)==-1: continue
    if nproc==1:
      _process_pdb_filename(pdb_file)
    else:
      rc = pool.apply_async(
        _process_pdb_filename,
        [pdb_file],
        callback=callback,
        )
  if pool:
    pool.close()
    pool.join()
  os.system("rm *.pickle")
  os.system("rm *.eff")
  os.system("rm *.dat")

if __name__=="__main__":
  args = sys.argv[1:]
  del sys.argv[1:]
  run(*tuple(args))
