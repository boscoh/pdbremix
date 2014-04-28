import os
import urllib
import gzip

from docopt import docopt
import util


def expand_pdbs(pdbs):
  results = []
  for pdb in pdbs:
    if os.path.isfile(pdb):
      print "Reading PDB codes from", pdb
      new_pdbs = expand_pdbs(util.words_in_file(pdb))
      results.extend(new_pdbs)
    else:
      pdb = pdb.lower()
      if pdb.endswith('.pdb'):
        pdb = pdb[:-4]
      if len(pdb) != 4:
        raise Exception("PDB code is not 4 letters long")
      results.append(pdb)
  return results


def get_pdbs_with_http(pdbs):
  for pdb in expand_pdbs(pdbs):
    fname = pdb if pdb.endswith('pdb') else '%s.pdb' % pdb 
    if os.path.isfile(fname):
      print "Skip: %s exists" % fname
    else:
      try:
        site = 'http://www.rcsb.org/pdb/files'
        url = '%s/%s' % (site, fname)
        sock = urllib.urlopen(url)
        text = sock.read()
        sock.close()
        if '<html>' in text:
          raise Exception
        open(fname, 'w').write(text)
        print fname
      except:
        print "Failed: %s" % fname


template = """#!/bin/bash
ftp -n %(host)s <<END_SCRIPT
quote user guest
quote 
ls
cd pub/pdb/data/structures/all/pdb
binary
%(mget_str)s
quit
END_SCRIPT
exit 0
"""


def get_pdbs_with_ftp(pdbs):
  pdbs = expand_pdbs(pdbs)
  entries = ['pdb%s.ent.gz' % pdb for pdb in pdbs]
  mget_scripts = ["get %s" % entry for entry in entries]
  substitutions = {
    'host': 'ftp.wwpdb.org',
    'mget_str': '\n'.join(mget_scripts),
  }
  ftp_script = template % substitutions
  open('ftp.sh', 'w').write(ftp_script)
  util.run_with_output('chmod +x ftp.sh')
  util.run_with_output('./ftp.sh')
  for pdb, entry in zip(pdbs, entries):
    fname = pdb + '.pdb'
    if os.path.isfile(entry):
      out_f = open(fname, 'w')
      for line in gzip.open(entry):
        out_f.write(line)
      util.clean_fname(entry)
      print fname
    else:
      print "Failed: %s" % fname
  util.clean_fname('ftp.sh')

