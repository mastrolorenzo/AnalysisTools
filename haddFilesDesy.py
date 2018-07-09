#!/bin/env python

import multiprocessing
import glob
import sys
import os

join = os.path.join
FILES_PER_BLOCK = 100

if len(sys.argv) < 2:
    print 'I need the path to job output directory...'
    exit(-1)

base_dir = sys.argv[1]

out_dir = join(base_dir, 'haddjobs')
if len(sys.argv) > 2:
    out_dir = join(sys.argv[2],'haddjobs')
print '>    output goes to', out_dir
os.system('mkdir -p ' + out_dir)

_, sub_dirs, __ = next(os.walk(base_dir))
if len(sys.argv) < 3:
    sub_dirs.remove('haddjobs')

def do_block(args):
    sub_dir, block_number, files = args
    outfile = '%s/sum_%s_%i.root' % (out_dir, sub_dir, block_number)
    if os.path.isfile(outfile):
        import ROOT
        f = ROOT.TFile(outfile)
        if not isinstance(f.Get('Events'), ROOT.TTree):
            print 'REDOING cuz no tree named "Events":', outfile
        else:
            print 'SKIPPING cuz existing:', outfile
            return outfile
    cmd = 'hadd -f %s ' % outfile
    cmd += ' '.join(files)
    print '>    issuing hadd for ', outfile
    os.system(cmd)
    return outfile

def prep_block_args(sub_dir):
    sample_files = glob.glob(join(base_dir, sub_dir, '*.root'))
    blocks = list(
        sample_files[i:i+FILES_PER_BLOCK]
        for i in xrange(0, len(sample_files), FILES_PER_BLOCK)
    )
    return list(
        (sub_dir, block_num, block_files)
        for block_num, block_files in enumerate(blocks)
    )

n_cpu = multiprocessing.cpu_count()/4*3 + 1
print 'running on %s cpus' % n_cpu
pool = multiprocessing.Pool(n_cpu)
block_args = list(
    args
    for sub_dir in sub_dirs
    for args in prep_block_args(sub_dir)
)
for outfile in pool.imap_unordered(do_block, block_args):
    print 'DONE:', os.path.basename(outfile)

# cmd = 'hadd -f {od}/sum_Data.root {od}/sum_Run*.root'.format(od=out_dir)
# print '>    issuing', cmd
# os.system(cmd)
