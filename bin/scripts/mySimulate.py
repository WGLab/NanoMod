
import os;
import sys;

import time;
#import numpy as np;
import random
import subprocess
import getpass

import h5py

import myDetect
import myFast5
from myCom import *

from collections import defaultdict

import rpy2.robjects as robjects
from rpy2.robjects.packages import importr

from pkg_resources import resource_string



target_pos = 3072
target_strand = '-'
target_chr = 'spel'
data_labels = ['simulate_case', 'simulate_case', 'folder_control']
foldersep = 3;

def myBinDefault(seqsize=(6184/3)):
   mpercertiles = [0.001, 0.0025, 0.005, 0.01, 0.02, 0.03, 0.04, 0.05]
   mybin = defaultdict()
   #myRankStr = [('P<=%.4f' % mpercertiles[0])]
   myRankStr = [('(, %.2f%%]' % (mpercertiles[0]*100))]
   for rp in range(1, int(mpercertiles[0]*seqsize)+1):
      mybin[rp] = myRankStr[0]
   split_points = []
   for mp in mpercertiles:
      split_points.append(int(mp*seqsize))
   for i in range(len(split_points)):
      if i==len(split_points)-1:
         #myRankStr.append('P>%.4f' % (mpercertiles[i]))
         #mybin[split_points[i]+1] = ('P>%.4f' % (mpercertiles[i]))
         myRankStr.append('(%.2f%%, )' % ((mpercertiles[i]*100)))
         mybin[split_points[i]+1] = ('(%.2f%%, )' % ((mpercertiles[i]*100)))
      else:
         a, b = split_points[i]+1, split_points[i+1]+1
         #myRankStr.append('%.4f<P<=%.4f' % (mpercertiles[i], mpercertiles[i+1]))
         myRankStr.append('(%.2f%%, %.2f%%]' % (mpercertiles[i]*100, mpercertiles[i+1]*100))
         for j in range(a, b):
            #mybin[j] = ('%.4f<P<=%.4f' % (mpercertiles[i], mpercertiles[i+1]))
            mybin[j] = ('(%.2f%%, %.2f%%]' % (mpercertiles[i]*100, mpercertiles[i+1]*100))
   return (mybin, split_points, myRankStr)
def myBinDefault_old():
#def myBinDefault():
   mybin = defaultdict()
   myRankStr = ['Rank 1']
   mybin[1] = 'Rank 1'
   split_points = [2, 4, 10, 20, 40, 60, 100]
   for i in range(len(split_points)):
      if i==len(split_points)-1: 
         myRankStr.append('Rank %d+' % (split_points[i]))
         mybin[split_points[i]] = ('Rank %d+' % (split_points[i]))
      else:
         a, b = split_points[i], split_points[i+1]
	 myRankStr.append('Rank %d-%d' % (a, b-1))
         for j in range(a, b):
            mybin[j] = ('Rank %d-%d' % (a, b-1))
   split_points.insert(0, 1);
   return (mybin, split_points, myRankStr)

def getSubFolders(moptions):
   ds3 = [moptions["wrkBase1"], moptions["wrkBase2"], moptions["wrkBase1"]]
   subfolders = [[], [], []]
   subf_max_int = []
   for ds_ind in range(len(ds3)):
      cur_ds = ds3[ds_ind]
      cur_max_folder_int = -1;
      cur_subfolder_list = os.listdir(cur_ds)
      for cursbf in cur_subfolder_list:
          if os.path.isdir(cur_ds + cursbf):
             try:
                cur_sbf_int = int(cursbf)
             except:
                if moptions['outLevel']<=OUTPUT_WARNING: print ("\t Warning!!! The folder name (%s) is not from nanopore software for the base folder %s" % (cursbf, cur_ds))
		continue
	     if cur_sbf_int>cur_max_folder_int: cur_max_folder_int = cur_sbf_int
	     subfolders[ds_ind].append(cursbf)
      if cur_max_folder_int == -1:
         if moptions['outLevel']<=OUTPUT_ERROR:
            print ("Error!!! no expected folder for %s (%d)" % (cur_ds, cur_max_folder_int))
            print subfolders[ds_ind]
	 sys.exit()
      subf_max_int.append(cur_max_folder_int)

   moptions['ds3'] = ds3
   moptions['subfolders'], moptions['subf_max_int'] = subfolders, subf_max_int

def readEvents(mfolder, mpkey, moptions):
   f5suf = moptions['.fast5']

   #moptions[mpkey] = defaultdict(set)
   noevent = 0; notfast5 = 0;

   fast5files = os.listdir(mfolder);
   for f5f in fast5files:
      if f5f[-len(f5suf):]==f5suf:
         with h5py.File(mfolder+f5f, 'r') as mf5:
             if mf5.__contains__(myFast5.rawAlignment_full):
                #moptions[mpkey][f5f] = [myFast5.ReadNanoraw_events(mf5), myFast5.ReadMapInfoInRef(mf5)]
		moptions[mpkey][f5f] = (myFast5.ReadNanoraw_events(mf5), myFast5.ReadMapInfoInRef(mf5))
             else:
                #moptions[mpkey][f5f] = [[], ()]
		moptions[mpkey][f5f] = ([], ())
		noevent += 1
      else:
         if moptions['outLevel']<=OUTPUT_WARNING:
            print "Warning!!! not fast5 file", mfolder+f5f
	 notfast5 += 1
   return (noevent, notfast5)

def getGenomeEvents(data_dict, dat_labels, moptions):
   start_time = time.time();
   for curd_ind in range(len(data_dict)):
       curd = data_dict[curd_ind]
       curl = dat_labels[curd_ind]
       
       curdkeys = curd.keys(); 
       for curdk in curdkeys:
	   for i in range(len(curd[curdk][0])):
               if curd[curdk][1][2]=='+':
                   curpos = i + curd[curdk][1][1]
               else: curpos = curd[curdk][1][1] + len(curd[curdk][0])-1 - i
	       moptions[curl]['base'][(curd[curdk][1][0], curd[curdk][1][2])][curpos] = curd[curdk][0]['base'][i]
	       moptions[curl]['norm_mean'][(curd[curdk][1][0], curd[curdk][1][2])][curpos].append(curd[curdk][0]['norm_mean'][i])
   end_time = time.time();
   if moptions['outLevel']<=OUTPUT_INFO: print ("Combine events: consuming time=%d" % (end_time-start_time))

#def mSimulate1(mj, mi, mk, moptions):
#   casefolder = moptions["wrkBase2"] + ('%d/' % mj)
#   controlsample = moptions["wrkBase1"] + ('%d/' % mi)
#   controltest = moptions["wrkBase1"] + ('%d/' % mk)

def getNextFolder(cf, nextadd=1):
   while cf[-1]=='/' or cf[-1]=='\\': cf = cf[:-1]
   posf = cf.rfind('/')
   basef = cf[:posf]; curfolderid = int(cf[posf+1:]);
   allf = []; exsf = []; fid = []
   for i in range(0, nextadd + 2):
      curf = basef+'/'+str(curfolderid + i) + '/'
      allf.append(curf);
      fid.append(curfolderid + i)
      exsf.append(os.path.isdir(curf));
   for i in range(len(exsf)):
      if exsf[i]: maxexist = i
      #print allf[i], exsf[i], fid[i]
   #print maxexist, fid[maxexist], cf, nextadd
   if maxexist>nextadd: return allf[nextadd]
   else: return basef+'/'+str((curfolderid+nextadd)%fid[maxexist]) + '/'


def mSimulate1(casefolder, controlsample, controltest, moptions):
   if moptions['outLevel']<=OUTPUT_INFO: print 'case foler', casefolder
   start_time = time.time();
   mpkey='casefolder'; moptions[mpkey] = defaultdict(set)
   #noevent, notfast5 = readEvents(casefolder, 'casefolder', moptions)
   noevent, notfast5 = 0, 0
   for nadd in range(foldersep):
      curf = getNextFolder(casefolder,nadd)
      no1event, not1fast5 = readEvents(curf, 'casefolder', moptions)
      noevent += no1event
      notfast5 += not1fast5
      if moptions['outLevel']<=OUTPUT_WARNING: print ('Read case %d fast5 files(available=%d, inavail=%d), %d non-fast5 files for %s(%s)' % (len(moptions['casefolder']), len(moptions['casefolder'])-noevent, noevent, notfast5, casefolder, curf))
   end_time = time.time();
   if moptions['outLevel']<=OUTPUT_INFO: print ("consuming time=%d" % (end_time-start_time))

   if moptions['outLevel']<=OUTPUT_INFO: print 'control sample folder', controlsample
   start_time = time.time();
   mpkey='controlsample'; moptions[mpkey] = defaultdict(set)
   #noevent, notfast5 = readEvents(controlsample, 'controlsample', moptions)
   noevent, notfast5 = 0, 0
   for nadd in range(foldersep):
      curf = getNextFolder(controlsample,nadd)
      no1event, not1fast5 = readEvents(curf, 'controlsample', moptions)
      noevent += no1event
      notfast5 += not1fast5
      if moptions['outLevel']<=OUTPUT_WARNING: print ('Read mixed-c %d fast5 files(available=%d, inavail=%d), %d non-fast5 files for %s(%s)' % (len(moptions['controlsample']), len(moptions['controlsample'])-noevent, noevent, notfast5, controlsample, curf))
   end_time = time.time();
   if moptions['outLevel']<=OUTPUT_INFO: print ("consuming time=%d" % (end_time-start_time))

   if moptions['outLevel']<=OUTPUT_INFO: print 'control test folder', controltest
   start_time = time.time();
   mpkey='controltest'; moptions[mpkey] = defaultdict(set)
   noevent, notfast5 = readEvents(controltest, 'controltest', moptions)
   noevent, notfast5 = 0, 0
   for nadd in range(foldersep):
      curf = getNextFolder(controltest,nadd)
      no1event, not1fast5 = readEvents(curf, 'controltest', moptions)
      noevent += no1event
      notfast5 += not1fast5
      if moptions['outLevel']<=OUTPUT_WARNING: print ('Read control %d fast5 files(available=%d, inavail=%d), %d non-fast5 files for %s(%s)' % (len(moptions['controltest']), len(moptions['controltest'])-noevent, noevent, notfast5, controltest, curf))
   end_time = time.time();
   if moptions['outLevel']<=OUTPUT_INFO: print ("consuming time=%d" % (end_time-start_time))

   moptions['ds2'] = ['simulate_case', 'folder_control']

   for mperc in moptions['Percentages']:
      start_time = time.time();
      if not moptions.has_key('PercDis'):
         moptions['PercDis'] = defaultdict(list)
      for rt in range(moptions['random_times']):
         if abs( len(moptions['casefolder']) - len(moptions['controlsample']) ) > 10:
            if moptions['outLevel']<=OUTPUT_WARNING: print 'Warning!!! larger difference', len(moptions['casefolder']), len(moptions['controlsample']), abs( len(moptions['casefolder']) - len(moptions['controlsample']) ), casefolder, controlsample

         #case_rand = [random.uniform(0,1) for _ in range(len(moptions['casefolder']))]
	 #mcase1 = {x:moptions['casefolder'][x] for x,y in zip(moptions['casefolder'], case_rand) if y<=mperc}
	 mcase1 = {x:moptions['casefolder'][x] for x,y in zip(moptions['casefolder'], [random.uniform(0,1) for _ in range(len(moptions['casefolder']))]) if y<=mperc}

	 #con1_rand = [random.uniform(0,1) for _ in range(len(moptions['controlsample']))]
	 #mcon1 = {x:moptions['controlsample'][x] for x,y in zip(moptions['controlsample'], con1_rand) if y<1-mperc}
	 mcon1 = {x:moptions['controlsample'][x] for x,y in zip(moptions['controlsample'], [random.uniform(0,1) for _ in range(len(moptions['controlsample']))]) if y<1-mperc}

	 #con2_rand = [random.uniform(0,1) for _ in range(len(moptions['controltest']))]
	 #mcon2 = {x:moptions['controltest'][x] for x,y in zip(moptions['controltest'], con2_rand) if y<=mperc}

	 if moptions['outLevel']<=OUTPUT_DEBUG: print "randtime, len(case), len(con), len(test) ", rt, len(mcase1), len(mcon1), len(moptions['controltest'])
         moptions['sign_test'] = []
	 moptions['sorted_sign_test'] = []

         moptions[moptions['ds2'][0]] = {}
	 moptions[moptions['ds2'][1]] = {}

	 moptions[moptions['ds2'][0]]['base'] = defaultdict(lambda: defaultdict(str))
	 moptions[moptions['ds2'][1]]['base'] = defaultdict(lambda: defaultdict(str))

	 moptions[moptions['ds2'][0]]['norm_mean'] = defaultdict(lambda: defaultdict(list))
	 moptions[moptions['ds2'][1]]['norm_mean'] = defaultdict(lambda: defaultdict(list))

         getGenomeEvents([mcase1, mcon1, moptions['controltest']], data_labels, moptions)

         myDetect.mfilter_coverage(moptions)

	 myDetect.mtest2(moptions)

	 moptions['PercDis'][mperc].append(getTopRank(moptions))
	 if moptions['outLevel']<=OUTPUT_INFO: print 'Rank ', mperc, rt, moptions['PercDis'][mperc][-1]
	 end_time = time.time();
	 if moptions['outLevel']<=OUTPUT_INFO: print ("consuming time=***%d*** for perc=%.2f_randtime=%d" % (end_time-start_time, mperc, rt))
	 sys.stdout.flush()
   
   if moptions['outLevel']<=OUTPUT_INFO: print ''
   for mperc in moptions['Percentages']:
       if moptions['outLevel']<=OUTPUT_WARNING: print mperc, moptions['PercDis'][mperc]
   sys.stdout.flush()

   cur_file_base = moptions['outFolder'] + '/' + moptions["FileID"] 
   opfile = cur_file_base + '.output'
   mySaveRes(opfile, moptions) 
   os.system('touch '+cur_file_base+'.done')

def mySaveRes(opfile, moptions, format1='%.5f', format2=' %d'):
   fh = open(opfile, 'w')
   #print 'mrank', moptions['PercDis']; #moptions['Percentages']
   for mperc in moptions['Percentages']:
      #print 'mrank', mperc, moptions['PercDis'][mperc]
      fh.write(format1 % mperc)
      for currank in moptions['PercDis'][mperc]:
         #print 'mrank', mperc, currank
         if int(currank)<0: 
            if moptions['outLevel']<=OUTPUT_INFO: print 'Negative rank', currank, int(currank)
            continue;
         
         fh.write(format2 % currank)
      fh.write('\n')
   fh.close();

def printTestInfo1(significant_pos, moptions):
    mstr = [("At pos %4d of %s strand in %s " % (significant_pos[0][2]+1, significant_pos[0][1], significant_pos[0][0]))]
    if moptions["neighborPvalues"]>0 and (not moptions["testMethod"]=="ks"):
       mstr.append('u=%.1E(%.1E) t=%.1E(%.1E) ks=%.1E(%.1E) c=%.1E(%.1E)' % (significant_pos[1][0][1],significant_pos[1][0][0], significant_pos[1][1][1],significant_pos[1][1][0], significant_pos[1][2][1],significant_pos[1][2][0], significant_pos[1][3][1],significant_pos[1][3][0]))
    else:
       mstr.append('u=%.1E(%.1E) t=%.1E(%.1E) ks=%.1E(%.1E)' % (significant_pos[1][0][1],significant_pos[1][0][0], significant_pos[1][1][1],significant_pos[1][1][0], significant_pos[1][2][1],significant_pos[1][2][0]))
    return ' '.join(mstr)

def getTopRank(moptions):
    curn = 0; printn = 0;
    output_pos = []
    rank_default = -1

    closesize = moptions["neighborPvalues"]*2
    if moptions['RegionRankbyST']==1:
       closesize = moptions["window"] #*2-1
       if closesize<1: closesize = 1

    rank_of_interest = [rank_default, -target_pos]
    if rank_of_interest[1]>-1000000: rank_of_interest[1] = -1000000

    for significant_pos in moptions['sorted_sign_test']:
       too_close_to_previous = False;
       if printn<3:
          if moptions['outLevel']<=OUTPUT_DEBUG: #INFO:
             print curn, printTestInfo1(significant_pos, moptions)
	  printn = printn + 1
       for pre_pos in output_pos:
          if pre_pos[0]==significant_pos[0][0] and pre_pos[1]==significant_pos[0][1] and \
             abs(pre_pos[2]-significant_pos[0][2])<closesize: #<=closesize: #moptions["neighborPvalues"]*2:
             too_close_to_previous = True; break;
       if too_close_to_previous: continue;

       cur_ind = moptions['sign_test'].index(significant_pos)
       noenough = False;
       #sk = (significant_pos[0][0], significant_pos[0][1])
       for mind in range(cur_ind-moptions["window"], cur_ind+moptions["window"]+1):
          if not myDetect.pos_check(moptions['sign_test'], cur_ind, mind):
             noenough = True;
          if noenough: break;
       if not noenough:
          output_pos.append((significant_pos[0][0], significant_pos[0][1], significant_pos[0][2]))
          curn = curn + 1
          if significant_pos[0][0]==target_chr and significant_pos[0][1]==target_strand and \
             abs(rank_of_interest[1]-target_pos)>abs(significant_pos[0][2]-target_pos) and abs(significant_pos[0][2]-target_pos)<closesize:  #<=closesize: #moptions["neighborPvalues"]*2:
             rank_of_interest = [curn, significant_pos[0][2], moptions[moptions['ds2'][0]]['base'][(significant_pos[0][0], significant_pos[0][1])][significant_pos[0][2]]]

       #if curn>=moptions["topN"]: break;  
       if rank_of_interest[-1]>-1: break;
    return rank_of_interest[0]

def mSimulate(moptions):
   moptions['random_times'] = 10; #3; #10
   #moptions['outLevel']=OUTPUT_INFO

   if moptions['outLevel']<=OUTPUT_INFO: print 'Percentages', moptions['Percentages']
   random.seed(1)

   if moptions["wrkBase3"]==None:
      current_username = getpass.getuser();
      qscmd = 'qstat -u "'+current_username+'" | awk \'NR>2 && ($5=="r" || $5=="qw") && substr($3, 1, 4)=="sim_" {print $4}\' | sort | wc -l'

      getSubFolders(moptions)
      if moptions['outLevel']<=OUTPUT_INFO: print moptions['subfolders'], '\n', moptions['subf_max_int']

      total_jobs = defaultdict();
      for mj in range(moptions['subf_max_int'][1]):
         for mi in range(moptions['subf_max_int'][0]):
            #mk = (mi + 1)%(moptions['subf_max_int'][0]);
	    mk = (mi + foldersep)%(moptions['subf_max_int'][0]);
            casefolder = moptions["wrkBase2"] + ('%d/' % mj)
            controlsample = moptions["wrkBase1"] + ('%d/' % mi)
            controltest = moptions["wrkBase1"] + ('%d/' % mk)

            for cur_perc in moptions['Percentages']:
               cur_file_id = ('%s_%d_%d_%d_%.5f' % (moptions["FileID"], mi, mj, mk, cur_perc))
               #cmd = ("echo 'python modDetection.py simulate --wrkBase1 %s --wrkBase2 %s --wrkBase3 %s --Percentages %.5f --FileID %s --outFolder %s' | qsub -V -cwd -l h_vmem=8G -N sim_%s -e %s/sim_%s.e -o %s/sim_%s.o -q all.q,bigmem" % (controlsample, casefolder, controltest, cur_perc, cur_file_id, moptions['outFolder'], cur_file_id, moptions['outFolder'], cur_file_id, moptions['outFolder'], cur_file_id))
	       cmd = ("echo 'python modDetection.py simulate --testMethod %s --wrkBase1 %s --wrkBase2 %s --wrkBase3 %s --Percentages %.5f --FileID %s --outFolder %s' | qsub -V -cwd -l h_vmem=12G -N sim_%s -e %s/sim_%s.e -o %s/sim_%s.o -q all.q" % (moptions['testMethod'], controlsample, casefolder, controltest, cur_perc, cur_file_id, moptions['outFolder'], cur_file_id, moptions['outFolder'], cur_file_id, moptions['outFolder'], cur_file_id))
	       total_jobs[cur_file_id] = [cmd, ('%s/sim_%s.e' % (moptions['outFolder'], cur_file_id)), ('%s/sim_%s.o' % (moptions['outFolder'], cur_file_id)), ('%s/%s.output' % (moptions['outFolder'], cur_file_id)), ('%s/%s.done' % (moptions['outFolder'], cur_file_id))]
               
	       #rmcmd = ' '.join(['rm', total_jobs[cur_file_id][1], total_jobs[cur_file_id][2], total_jobs[cur_file_id][3],total_jobs[cur_file_id][4]])
	       #os.system(rmcmd)
               for rmf in total_jobs[cur_file_id][1:]:
                  if os.path.isfile(rmf):
                     os.system('rm '+rmf)

	       print total_jobs[cur_file_id]; sys.stdout.flush()
	       os.system(cmd)
	       time.sleep(0.5)
	 #   if mi>1: break;
         #if mj>1: break;	  

      if not moptions.has_key('PercDis'):
          moptions['PercDis'] = defaultdict(list)
      
      efile = ('%s/sim_%s_merge.e' % (moptions['outFolder'], moptions["FileID"]))
      ofile = ('%s/sim_%s_merge.o' % (moptions['outFolder'], moptions["FileID"]))
      opfile = ('%s/%s_merge.output' % (moptions['outFolder'], moptions["FileID"]))
      rmfiles = [efile, ofile, opfile]
      for rmf in rmfiles:
         if os.path.isfile(rmf):
            os.system('rm '+rmf)

      default_wait_time = 100;
      wait_time = default_wait_time; 
      max_wait_time = 300; min_wait_time = 30;
      zero_qstat_jobs_times = 0;
      jobkeys = total_jobs.keys();
      while len(jobkeys)>0:
         finished_jobs = 0;
         for jk in jobkeys:
            if os.path.isfile(total_jobs[jk][-1]):
               time.sleep(1)
               finished_jobs = finished_jobs + 1

               fr = open(total_jobs[jk][-2], 'r')
	       lines = fr.readlines();
               for l in lines:
                  lsp = l.split()
                  curp = float(lsp[0])
		  if moptions['outLevel']<=OUTPUT_DEBUG: print jk, lsp; sys.stdout.flush()
		  moptions['PercDis'][curp].extend(lsp[1:])
               
	       if os.path.getsize(total_jobs[jk][1])>0:
                  if moptions['outLevel']<=OUTPUT_DEBUG: print ('tail -vn +1 %s >> %s' % (total_jobs[jk][1], efile)); sys.stdout.flush()
                  os.system('tail -vn +1 %s >> %s' % (total_jobs[jk][1], efile))
               if os.path.getsize(total_jobs[jk][2])>0:
                  if moptions['outLevel']<=OUTPUT_DEBUG: print ('tail -vn +1 %s >> %s' % (total_jobs[jk][2], ofile)); sys.stdout.flush()
		  os.system('tail -vn +3 %s >> %s' % (total_jobs[jk][2], ofile))
               if os.path.getsize(total_jobs[jk][3])>0:
                  if moptions['outLevel']<=OUTPUT_DEBUG: print ('tail -vn +1 %s >> %s' % (total_jobs[jk][3], opfile)); sys.stdout.flush()
	          os.system('tail -vn +1 %s >> %s' % (total_jobs[jk][3], opfile))
               rmcmd = ' '.join(['rm', total_jobs[jk][1], total_jobs[jk][2], total_jobs[jk][3], total_jobs[jk][4]])
	       if moptions['outLevel']<=OUTPUT_DEBUG: 
                  rmcmd = ' '.join(['rm', total_jobs[jk][1], total_jobs[jk][4]])
                  print 'RM', rmcmd
	       os.system(rmcmd)
	       del total_jobs[jk]
	 if finished_jobs>0:
            wait_time = reduce_wait_time(wait_time, max_wait_time, min_wait_time)
         else:
            wait_time = double_wait_time(wait_time, max_wait_time, min_wait_time)

	 runjobs = int(subprocess.check_output(qscmd, shell=True))
	 print ('Total running jobs=%d(qstat=%d)' % (len(jobkeys), runjobs))+' ---'+moptions["wrkBase1"].split('/')[0], moptions["wrkBase2"].split('/')[0],'---Waiting time', wait_time, '3s/per dot'
	 sys.stdout.flush()

	 jobkeys = total_jobs.keys();

	 if len(jobkeys)>0 and runjobs==0:
            if moptions['outLevel']<=OUTPUT_WARNING: 
               print ('Warning!!! Expected %d jobs running. But %d jobs running using qstat' % (len(jobkeys), runjobs))
	       print 'Please check whether error occur. The program will exit after this message appears in three consecutive times.'
            zero_qstat_jobs_times = zero_qstat_jobs_times + 1
            if zero_qstat_jobs_times>=3:
               print 'Unfished jobs:'
	       jobkeys.sort()
	       print '\t', jobkeys
	       if wait_time>default_wait_time: 
                  wait_time = default_wait_time
               break;    
         else: zero_qstat_jobs_times = 0

         if len(jobkeys)>0:
            for i in range(wait_time/3+1):
               time.sleep(3);
	       sys.stdout.write('>'); sys.stdout.flush()

	       #cur_hint = ('> %ds left' % (wait_time-i*3))
	       #sys.stdout.write(cur_hint); #print '.',; 
	       #sys.stdout.flush()
	       #sys.stdout.write("\b"*(len(cur_hint)-1))
            print ''
            #time.sleep(wait_time)

      jobkeys = total_jobs.keys();
      runjobs = int(subprocess.check_output(qscmd, shell=True))
      #print 'Info: Total running jobs', str(len(jobkeys))+' ---'+moptions["wrkBase1"], moptions["wrkBase2"]+'--- '+'qstat running jobs', runjobs, '---Waiting time', wait_time
      print ('Total running jobs=%d(qstat=%d)' % (len(jobkeys), runjobs))+' ---'+moptions["wrkBase1"].split('/')[0], moptions["wrkBase2"].split('/')[0],'---Waiting time', wait_time

      cur_file_base = moptions['outFolder'] + '/' + moptions["FileID"]
      opfile = cur_file_base + '_all.output'
      if moptions['outLevel']<=OUTPUT_DEBUG: print moptions['PercDis']
      #print moptions['PercDis']
      mySaveRes(opfile, moptions, format2=' %s')
      mplotHis(moptions);
   else:
      mSimulate1(moptions["wrkBase2"], moptions["wrkBase1"], moptions["wrkBase3"], moptions)
      print '\n'

def double_wait_time(wait_time, max_wait_time, min_wait_time):
   wait_time = int(wait_time * 1.5)
   if wait_time>max_wait_time: wait_time=max_wait_time
   return wait_time
def reduce_wait_time(wait_time, max_wait_time, min_wait_time):
   wait_time = wait_time/3
   if wait_time<min_wait_time: wait_time=min_wait_time
   return wait_time

def group_rank(moptions):
   #print moptions['PercDis']
   mgroupRanks = defaultdict(lambda: defaultdict(int))

   mybin, split_points, myRankStr = myBinDefault()

   percKeys = moptions['PercDis'].keys(); 
   for pk in percKeys:
      for curr in moptions['PercDis'][pk]:
         curr = int(curr)
         if 0<curr<=split_points[-1]:
            cur_r_str = mybin[curr]
         else: 
            if curr==0: print "Warning!!! Rank at pos 0"
	    if curr<0: print "Warning!!! Negative Rank", curr
            cur_r_str = myRankStr[-1]
         mgroupRanks[pk][cur_r_str] += 1

   for pk in percKeys:
      total = 0;
      for k in mgroupRanks[pk]:
         total += mgroupRanks[pk][k]
      if total==0:
         print 'Warning!!! total=0 for ', pk;
      else:
         for k in mgroupRanks[pk]:
             mgroupRanks[pk][k] = mgroupRanks[pk][k]/float(total)

   perclist = []
   rankgrouplist = []
   rankperclist = []
   for pk in percKeys:
      #print pk
      for k in mgroupRanks[pk]:
         #print '\t', k, ('%.3f' % (mgroupRanks[pk][k]))
         perclist.append(pk)
	 rankgrouplist.append(k)
	 rankperclist.append(mgroupRanks[pk][k])

   return (perclist, rankgrouplist, rankperclist, split_points, myRankStr)

def mplotHis(moptions):
   perclist, rankgrouplist, rankperclist, split_points, myRankStr = group_rank(moptions)

   figname = moptions["FileID"]
   mresfolder = moptions['outFolder']

   ggplot = importr('ggplot2')
   importr('gridExtra')

   spvector = robjects.IntVector(split_points)
   rankstrvector = robjects.StrVector(myRankStr)

   moptions["Percentages"].sort()
   percvector = robjects.FloatVector(moptions["Percentages"])

   mdfperc = robjects.DataFrame({"MixedPerc":robjects.FactorVector(robjects.FloatVector(perclist), levels=percvector, labels=percvector), "Percentile":robjects.FactorVector(robjects.StrVector(rankgrouplist), levels=rankstrvector, labels=rankstrvector), "Fraction":robjects.FloatVector(rankperclist)})

   robjects.r(resource_string(__name__, 'Rscript/Hist_sim_plot.R'))
   robjects.r('pdf("'+mresfolder+'/hist_'+figname+'.pdf", width='+("%.0f" % (len(moptions["Percentages"])*1.5))+', height=4, onefile = TRUE)')

   robjects.globalenv['Hist_sim_plot'](mdfperc, spvector, rankstrvector)

   robjects.r('dev.off()')



if __name__=='__main__':
   moptions = {}
   if len(sys.argv)<2:
      print 'Usages: python mySimulate.py 1/(more 9perpage/27perpage)'
      sys.exit()
   if sys.argv[1]=='1':
      moptions["FileID"] = 'mod'
      moptions['outFolder'] = '../mRes/'
      moptions["Percentages"] = [0.3, 0.35, 0.5]

      mresfile = ('%s/%s_all.output' % (moptions['outFolder'], moptions["FileID"]))
      moptions['PercDis'] = {}
      fr = open(mresfile, 'r')
      lines = fr.readlines()
      for l in lines:
         lsp = l.split();
         moptions['PercDis'][float(lsp[0])] = lsp[1:]
         for i in range(len(moptions['PercDis'][float(lsp[0])])):
            moptions['PercDis'][float(lsp[0])][i] = int(moptions['PercDis'][float(lsp[0])][i])
      fr.close()

      mplotHis(moptions)
   elif sys.argv[1]=='more':
      more_usage_str = 'Usages: python mySimulate.py more 9perpage/27perpage'
      moptions['outFolder'] = '../logsim/'
      moptions['outFolder'] = '../logsim_dif_perc/'
      moptions["Percentages"] = [0.15,0.2,0.25,0.3,0.4,0.5,0.6,0.7]

      control = {'vector':'vector_SpeI_cut', 'oligo':'ctrl_oligo_SpeI_cut'}
      #control = {'oligo':'ctrl_oligo_SpeI_cut'}

      cases = {'Glucose':'Glucose_SpeI_cut', \
               'BrdU':'BrdU_SpeI_cut_data', \
               'A647':'A647_Click_SpeI_cut', \
               'Bio':'Bio_pur_SpeI_cut_data', \
               'Edu':'Edu_SpeI_cut_data', \
               'IdU':'IdU', \
               'A488':'A488_Click_SpeI_cut_data', \
               'Nicotine':'Nicotine_SpeI_cut', \
               'Azide':'Azide_SpeI_cut'}

      conkeys = control.keys(); conkeys.sort()
      caskeys = cases.keys(); caskeys.sort(); caskeys = caskeys[::-1]

      xlabs = robjects.StrVector(["MixedPerc"])

      #testmethod3 = ['fisher', 'stouffer', 'ks']
      testmethod3 = ['ks', 'stouffer', 'fisher']

      for conk in conkeys:
         perclist_all = []
         rankgrouplist_all = []
         rankperclist_all = []
	 case_control_all = []
	 testmethod_all = []

         for testmethod in testmethod3:
            for cask in caskeys:
               moptions["FileID"] = ('%s_%s_%s' % (conk, cask, testmethod))
               mresfile = ('%s/%s_all.output' % (moptions['outFolder'], moptions["FileID"]))
               moptions['PercDis'] = {}
               fr = open(mresfile, 'r')
               lines = fr.readlines()
               for l in lines:
                  lsp = l.split();
                  moptions['PercDis'][float(lsp[0])] = lsp[1:]
                  for i in range(len(moptions['PercDis'][float(lsp[0])])):
                     moptions['PercDis'][float(lsp[0])][i] = int(moptions['PercDis'][float(lsp[0])][i])
               fr.close()

               perclist, rankgrouplist, rankperclist, split_points, myRankStr = group_rank(moptions)

	       if True: # False: #True:
                  print '>>>', moptions["FileID"]
                  opdict = defaultdict(lambda: defaultdict(float));
                  max_rankgrp_len = 0
	          for perc, rankgrp, rankperc in zip(perclist, rankgrouplist, rankperclist):
		     opdict[perc][rankgrp] = rankperc
		     if len(rankgrp)>max_rankgrp_len: max_rankgrp_len = len(rankgrp)
                     #print len(perclist), len(rankgrouplist), len(rankperclist), perc, rankgrp, rankperc
                     #print '\t', ('%.5f %s %.3f' % (perc, rankgrp, rankperc))
                  print '   ',
	          for rankgrp in myRankStr:
                     print (('   %'+str(max_rankgrp_len)+'s') % rankgrp),
                  print ''
	          opdictkeys = opdict.keys(); opdictkeys.sort()
	          for perc in opdictkeys:
                     print ('%.2f' % perc),
		     for rankgrp in myRankStr:
                        print (('   %'+str(max_rankgrp_len)+'s') % ('%.3f' % opdict[perc][rankgrp])),
                     print ''

	       perclist_all.extend(perclist)
	       rankgrouplist_all.extend(rankgrouplist)
	       rankperclist_all.extend(rankperclist);
	       #case_control_all.extend(['%s_%s' % (conk, cask)]*len(perclist))
	       if len(sys.argv)>=3 and sys.argv[2]=='9perpage':
                  testmethod_all.extend([testmethod]*len(perclist));
                  case_control_all.extend(['%s_%s_%s' % (conk, cask, testmethod)]*len(perclist))
               elif len(sys.argv)>=3 and sys.argv[2]=='27perpage':
                  testmethod_all.extend(['%s_%s_%s' % (conk, cask, testmethod)]*len(perclist))
		  case_control_all.extend(['%s_%s' % (conk, cask)]*len(perclist))
               else:
                  print more_usage_str; #'Usages: python mySimulate.py more 9perpage/27perpage'
                  sys.exit()


         mresfolder = moptions['outFolder']

         ggplot = importr('ggplot2')
         importr('gridExtra')

         spvector = robjects.IntVector(split_points)
         rankstrvector = robjects.StrVector(myRankStr)

         moptions["Percentages"].sort()
         percvector = robjects.FloatVector(moptions["Percentages"])

	 mdfperc = robjects.DataFrame({"MixedPerc":robjects.FactorVector(robjects.FloatVector(perclist_all), levels=percvector, labels=percvector), "Percentile":robjects.FactorVector(robjects.StrVector(rankgrouplist_all), levels=rankstrvector, labels=rankstrvector), "Fraction":robjects.FloatVector(rankperclist_all), "CaseControl":robjects.StrVector(case_control_all), "Method":robjects.StrVector(testmethod_all)})

         if len(sys.argv)<3:
            print more_usage_str; #'Usages: python mySimulate.py more 9perpage/27perpage' 
            sys.exit()

         figname = conk+'_' + sys.argv[2]
         if sys.argv[2]=='9perpage':
            robjects.r(resource_string(__name__, 'Rscript/Hist_sim_plot9.R'))
            robjects.r('pdf("'+mresfolder+'/hist_'+figname+'.pdf", width='+("%.0f" % (len(moptions["Percentages"])*2))+', height=9, onefile = TRUE)')

            robjects.globalenv['Hist_sim_plot9'](mdfperc, spvector, rankstrvector, xlabs)
         elif sys.argv[2]=='27perpage':
            robjects.r(resource_string(__name__, 'Rscript/Hist_sim_plot27.R'))
	    robjects.r('pdf("'+mresfolder+'/hist_'+figname+'.pdf", width='+("%.0f" % (len(moptions["Percentages"])*1.5))+', height=4, onefile = TRUE)')
	    robjects.globalenv['Hist_sim_plot27'](mdfperc, spvector, rankstrvector, xlabs)
         else:   
            print more_usage_str; #'Usages: python mySimulate.py more 9perpage/27perpage'
	    sys.exit()

         robjects.r('dev.off()')

      

