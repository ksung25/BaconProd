#!/usr/bin/env python

# baconBatch.py #############################################################################
# Python driver for Bacon Analyzer executable
# Original Author N.Wardle (CERN) 

# TODO : Provide output support to EOS
# For now assume output is small enough to store locally. 
# ------------------------------------------------------------------------------------

import sys, commands, os, fnmatch
from optparse import OptionParser
from optparse import OptionGroup
from BaconAna.Utils.makeFilelist import *

# Ok this is dangerous since we pretty much have to assume some arguments for the exec
# Take WAnalysis as the standard command line style
# 'maxevents, input, isGen
default_args = ['10000000','nothing.root','1'] #,output.root -> could add to analyzer

# Options
parser = OptionParser()
parser = OptionParser(usage="usage: %prog analyzer [options] \nrun with --help to get list of options")
parser.add_option("-d","--directory",default='./',help="Pick up files from a particular directory. can also pass from /eos/")
parser.add_option("-o","--outdir",default='bacon',help="output for analyzer. This will always be the output for job scripts.")
parser.add_option("-a","--arg",dest="args",default=[],action="append",help="Pass analyzer args n:arg")

# Make batch submission scripts options
parser.add_option("-n","--njobs",dest="njobs",type='int',default=1,help="Split into n jobs, will automatically produce submission scripts")
parser.add_option("-q","--queue",default='1nh',help="submission queue")

parser.add_option("--dryRun",default=False,action="store_true",help="Do nothing, just create jobs if requested")

# Monitor options (submit,check,resubmit failed)  -- just pass outodir as usual but this time pass --monitor sub --monitor check or --monitor resub
parser.add_option("--monitor",default='',help="Monitor mode (sub/resub/check directory of jobs)")

cwd = os.getcwd()
(options,args) = parser.parse_args()
if len(args)<1 and not options.monitor: sys.exit('Error -- must specify ANALYZER' )
njobs = options.njobs

def write_job(exec_line, out, analyzer, i, n):

	sub_file = open('%s/sub_%s_job%d.sh'%(out,analyzer,i),'w')
	sub_file.write('#!/bin/bash\n')
	sub_file.write('# Job Number %d, running over %d files \n'%(i,n))
	sub_file.write('touch %s.run\n'%os.path.abspath(sub_file.name))
	sub_file.write('cd %s\n'%os.getcwd())
	sub_file.write('eval `scramv1 runtime -sh`\n')
	sub_file.write('cd -\n')
	sub_file.write('mkdir -p scratch\n')
	sub_file.write('cd scratch\n')
	sub_file.write('cp -p $CMSSW_BASE/bin/$SCRAM_ARCH/%s .\n'%analyzer)
	sub_file.write('mkdir -p %s\n'%(out))

	sub_file.write('if ( %s ) then\n'%exec_line)
	sub_file.write('\t hadd -f Output_job%d.root %s/*.root \n'%(i,(out)))
	sub_file.write('\t mv Output_job*.root %s\n'%os.path.abspath(out))
	sub_file.write('\t rm -rf ./bacon ./Output_job* \n')
	sub_file.write('\t touch %s.done\n'%os.path.abspath(sub_file.name))
	sub_file.write('else\n')
	sub_file.write('\t touch %s.fail\n'%os.path.abspath(sub_file.name))
	sub_file.write('fi\n')
	sub_file.write('rm -f %s.run\n'%os.path.abspath(sub_file.name))
	sub_file.close()
	os.system('chmod +x %s'%os.path.abspath(sub_file.name))
  
def submit_jobs(lofjobs):
   for sub_file in lofjobs:
    os.system('rm -f %s.done'%os.path.abspath(sub_file))
    os.system('rm -f %s.fail'%os.path.abspath(sub_file))
    os.system('rm -f %s.log'%os.path.abspath(sub_file))
    os.system('bsub -q %s -o %s.log %s'%(options.queue,os.path.abspath(sub_file),os.path.abspath(sub_file)))
  
if options.monitor: 
  if options.monitor not in ['sub','check','resub']: sys.exit('Error -- Unknown monitor mode %s'%options.monitor)
  dir = options.outdir

  if options.monitor == 'sub' or options.monitor == 'resub': 
    # pick up job scripts in output directory (ends in .sh)
    lofjobs = []
    for root,dirs,files in os.walk(dir):
     for file in fnmatch.filter(files,'*.sh'):
       if options.monitor == 'resub' and not os.path.isfile('%s/%s.fail'%(root,file)): continue
       lofjobs.append('%s/%s'%(os.path.abspath(root),file))
    print 'Submitting %d jobs from directory %s'%(len(lofjobs),dir)
    submit_jobs(lofjobs) 

  if options.monitor == 'check': 
    failjobs = []
    runjobs  = []
    donejobs = []

    for root,dirs,files in os.walk(dir):
     for file in fnmatch.filter(files,'*.sh'):
       if os.path.isfile('%s/%s.fail'%(root,file)): failjobs.append('%s'%file)
       if os.path.isfile('%s/%s.done'%(root,file)):
       		if not '%s.sh'%file in failjobs : donejobs.append('%s'%file)
       if os.path.isfile('%s/%s.run'%(root,file)): runjobs.append('%s'%file)
    print 'Status of jobs directory ', dir
    print '  %d in status Fail -> (resub them with --monitor resub)'%len(failjobs)
    for job in failjobs : print '\t %s'%job
    print '  %d in status Running -> '%len(runjobs)
    for job in runjobs : print '\t %s'%job
    print '  %d in status Done -> '%len(donejobs)
    for job in donejobs : print '\t %s'%job

  sys.exit('Finished Monitor -- %s'%options.monitor)

def parse_to_dict(l_list):
  if len(l_list)<1: return {}
  ret = {}
  for item in l_list:     
    ni,varg = item.split(':') # should put a try here
    ret[int(ni)]=varg
  return ret

def getFilesJob(dir,job,njobs):
  if njobs == 1 : 
  	njobs = -1
	job = 0
  infiles = []
  if '/store/' in options.directory : infiles = makeCaFiles(options.directory,njobs,job)
  else : infiles = makeFiles(options.directory,njobs,job)
  return infiles

# -- MAIN
os.system('mkdir -p %s/%s'%(cwd,options.outdir)) 

analyzer = args[0]
analyzer_args = parse_to_dict(options.args)
exec_line = '%s'%analyzer


for arg_i,arg in enumerate(default_args):
  if arg_i in analyzer_args.keys(): 
  	arg = analyzer_args[arg_i]
  exec_line+=' %s '%arg

print 'running Analyser -- (default call) \n\t%s'%exec_line

if not options.dryRun and njobs > 1:
	print 'Writing %d Submission Scripts to %s (submit after with --monitor sub)'%(njobs,options.outdir)

for job_i in range(njobs):
 
 files = getFilesJob(options.directory,job_i,njobs)
 job_exec = ''

 nfiles_i = 0
 for fil_i,fil in enumerate(files):
   if not fil[1]: continue
   exec_line_i = exec_line.replace('nothing.root',fil[0]) 
   job_exec+=exec_line_i+'; mv Output.root %s/Ouput_job%d_file%d.root; '%(options.outdir,job_i,fil_i) 
   nfiles_i += 1

 if options.dryRun : 
 	print 'job %d/%d -> '%(job_i+1,njobs), job_exec
 elif njobs > 1: 
   write_job(job_exec, options.outdir, analyzer, job_i, nfiles_i)
 else: os.system(job_exec) 
