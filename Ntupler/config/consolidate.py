#! /usr/bin/env python
import commands,sys,os,subprocess,ROOT,numpy
from optparse import OptionParser
import argparse

eos='eos root://cmseos.fnal.gov'

aparser = argparse.ArgumentParser()
aparser.add_argument('-dir'   ,'--dir'      ,nargs='+',type=str,default='crap')
aparser.add_argument('-base'  ,'--base'     ,nargs='+',type=str,default=['13'])
aparser.add_argument('-mod'   ,'--mod'       ,nargs='+',type=str,default=[5])

args = aparser.parse_args()
count=0
basecount=0

def clear():
    global count
    global basecount
    #os.system('hadd %s_%s.root *.root ' % (args.dir[0],basecount))
    print 'hadd %s_%s.root *.root ' % (args.dir[0],basecount)
    print 'ls -lrth %s_%s.root' % (args.dir[0],basecount)
    #os.system('xrdcp %s_%s.root /store/group/lpcbacon/%s/%s/' % (args.dir[0],basecount,args.base[0],args.dir[0]))
    print 'xrdcp %s_%s.root root://cmseos.fnal.gov//store/group/lpcbacon/%s/%s/' % (args.dir[0],basecount,args.base[0],args.dir[0])
    #os.system('rm *.root')
    print 'rm *.root'
    count=0
    basecount=basecount+1

def search(dirname):
    global count
    global basecount
    #print 'dir',dirname
    dirSearch = '%s ls %s' %(eos,dirname)
    exists = commands.getoutput(dirSearch)
    for label in exists.splitlines():
        if label.find('log') > 0 or label == 'failed':
            os.system('%s rm -r %s/%s' % (eos,dirname,label))
            #print '%s rm -r %s/%s' % (eos,dirname,label)
            continue
        if label.find('.root') > 0 and dirname != '/store/group/lpcbacon/%s/%s' % (args.base[0],args.dir[0]):
            count=count+1
            shortdirname=dirname.replace('','')
            #os.system('xrdcp /eos/cms/%s/%s .' % (shortdirname,label))
            print 'xrdcp root://cmseos.fnal.gov/%s/%s .' % (shortdirname,label)
            if count % args.mod[0] == 0:
                clear()
            continue
        search('%s/%s' % (dirname,label))

def cleardir(dirname,iFile=True):
    dirSearch = '%s ls %s' %(eos,dirname)
    exists = commands.getoutput(dirSearch)
    for label in exists.splitlines():
        #print label
        if label.find('.root') > 0 and iFile:
            continue
        if label.find('.root') < 0 and not iFile:
            continue
        if iFile:
            os.system('echo %s rm -r %s/%s >> clearAll.sh' % (eos,dirname,label))
            #print 'echo %s rm -r %s/%s >> clearAll.sh'% (eos,dirname,label)
        else:
            os.system('%s rm %s/%s' % (eos,dirname,label))
            #print 'trying to erase'
            #print '%s rm %s/%s' % (eos,dirname,label)

cleardir('/store/group/lpcbacon/%s/%s' % (args.base[0],args.dir[0]),False )
print '\n'
print 'rm *.root'
print 'ls -lrth'
search('/store/group/lpcbacon/%s/%s' % (args.base[0],args.dir[0]) )
clear()
cleardir('/store/group/lpcbacon/%s/%s' % (args.base[0],args.dir[0]) )
