#! /usr/bin/env python
import commands,sys,os,subprocess,ROOT,numpy
from optparse import OptionParser
import argparse

eos='/afs/cern.ch/project/eos/installation/cms/bin/eos.select'

aparser = argparse.ArgumentParser()
aparser.add_argument('-dir'   ,'--dir'      ,nargs='+',type=str,default='crap')
aparser.add_argument('-base'  ,'--base'     ,nargs='+',type=str,default=['08'])
aparser.add_argument('-mod'   ,'--mod'       ,nargs='+',type=str,default=[50])

args = aparser.parse_args()
count=0
basecount=0

def clear():
    global count
    global basecount
    os.system('hadd %s_%s.root *.root ' % (args.dir[0],basecount))
    #os.system('mv *.root /tmp/pharris/tmp')
    os.system('cmsStage %s_%s.root /store/cmst3/group/monojet/production/%s/%s/' % (args.dir[0],basecount,args.base[0],args.dir[0]))
    os.system('rm *.root')
    count=0
    basecount=basecount+1

def search(dirname):
    global count
    global basecount
    print 'dir',dirname
    dirSearch = '%s ls %s' %(eos,dirname)
    exists = commands.getoutput(dirSearch)
    for label in exists.splitlines():
        if label.find('log') > 0 or label == 'failed':
            os.system('%s rm -r %s/%s' % (eos,dirname,label))
            continue
        if label.find('.root') > 0 and dirname != 'eos/cms/store/cmst3/group/monojet/production/%s/%s' % (args.base[0],args.dir[0]):
            count=count+1
            shortdirname=dirname.replace('eos/cms','')
            os.system('cmsStage %s/%s .' % (shortdirname,label))
            if count % args.mod[0] == 0:
                clear()
            continue
        search('%s/%s' % (dirname,label))

def cleardir(dirname):
    dirSearch = '%s ls %s' %(eos,dirname)
    exists = commands.getoutput(dirSearch)
    for label in exists.splitlines():
        if label.find('.root') > 0:
            continue
        os.system('echo %s rm -r %s/%s >> clearAll.sh' % (eos,dirname,label))
            

os.system('mkdir /tmp/pharris/tmp/')
search('eos/cms/store/cmst3/group/monojet/production/%s/%s' % (args.base[0],args.dir[0]) )
clear()
cleardir('eos/cms/store/cmst3/group/monojet/production/%s/%s' % (args.base[0],args.dir[0]) )
