
import os;
import sys;

import time;
import numpy as np;
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

import mySimulate


target_pos = 3072
target_strand = '-'
target_chr = 'spel'
data_labels = ['simulate_case', 'simulate_case', 'folder_control']
tofile = False; 
tofile = True;
test=False;
#test=True;

def getSubFolders(moptions):
   ds2 = [moptions["wrkBase1"], moptions["wrkBase2"]]
   subfolders = [[], []]
   subf_max_int = []
   for ds_ind in range(len(ds2)):
      cur_ds = ds2[ds_ind]
      cur_max_folder_int = -1;
      cur_subfolder_list = os.listdir(cur_ds)
      for cursbf in cur_subfolder_list:
          if os.path.isdir(cur_ds + cursbf):
             try:
                cur_sbf_int = int(cursbf)
             except:
                if moptions['outLevel']<=OUTPUT_WARNING: print ("\t Warning!!! The folder name is not from nanopore software for the base folder %s" % (cursbf, cur_ds))
		continue
	     if cur_sbf_int>cur_max_folder_int: cur_max_folder_int = cur_sbf_int
	     subfolders[ds_ind].append(cur_ds+cursbf+'/')
      if cur_max_folder_int == -1:
         if moptions['outLevel']<=OUTPUT_ERROR:
            print ("Error!!! no expected folder for %s (%d)" % (cur_ds, cur_max_folder_int))
            print subfolders[ds_ind]
	 sys.exit()
      subf_max_int.append(cur_max_folder_int)

   moptions['subfolders'] = {moptions["wrkBase1"]:subfolders[0], moptions["wrkBase2"]:subfolders[1]}
   moptions['subf_max_int'] = subf_max_int

def readEvents(mfolders, mpkey, moptions):
   f5suf = moptions['.fast5']

   moptions[mpkey] = defaultdict(set)
   noevent = 0; notfast5 = 0;

   start_time = time.time();
   for mfolder_ind in range(len(mfolders)):
      mfolder = mfolders[mfolder_ind]
      if test and mfolder_ind>2: break;
      fast5files = os.listdir(mfolder);
      for f5f in fast5files:
         if f5f[-len(f5suf):]==f5suf:
            #print mfolder+f5f
            with h5py.File(mfolder+f5f, 'r') as mf5:
                if mf5.__contains__(myFast5.rawAlignment_full):
	           moptions[mpkey][mfolder+f5f] = (myFast5.ReadNanoraw_events(mf5), myFast5.ReadMapInfoInRef(mf5))
                else:
                   moptions[mpkey][mfolder+f5f] = ([], ())
                   noevent += 1
         else:
            if moptions['outLevel']<=OUTPUT_WARNING:
               print "Warning!!! not fast5 file", mfolder+f5f
            notfast5 += 1
 
      end_time = time.time();
      if tofile:
         sys.stdout.write('>'); 
      else:
         cur_hint = ('> %d folders left; used time: %d' % (len(mfolders)-mfolder_ind-1, end_time-start_time))
	 sys.stdout.write(cur_hint); 
	 sys.stdout.flush()
         sys.stdout.write("\b"*(len(cur_hint)-1))
      sys.stdout.flush()	 

   return (noevent, notfast5)


def mSimulate1(moptions):
   getSubFolders(moptions)

   if moptions['outLevel']<=OUTPUT_INFO: print 'case foler', moptions["wrkBase2"]
   start_time = time.time();
   noevent, notfast5 = readEvents(moptions['subfolders'][moptions["wrkBase2"]], 'casefolder', moptions)
   if moptions['outLevel']<=OUTPUT_WARNING: print ('Read case %d fast5 files(available=%d, inavail=%d), %d non-fast5 files for %s' % (len(moptions['casefolder']), len(moptions['casefolder'])-noevent, noevent, notfast5, moptions["wrkBase2"]))
   end_time = time.time();
   if moptions['outLevel']<=OUTPUT_INFO: print ("consuming time=%d" % (end_time-start_time))

   if moptions['outLevel']<=OUTPUT_INFO: print 'control sample folder', moptions["wrkBase1"]
   start_time = time.time();
   noevent, notfast5 = readEvents(moptions['subfolders'][moptions["wrkBase1"]], 'controlsample', moptions)
   if moptions['outLevel']<=OUTPUT_WARNING: print ('Read control %d fast5 files(available=%d, inavail=%d), %d non-fast5 files for %s' % (len(moptions['controlsample']), len(moptions['controlsample'])-noevent, noevent, notfast5, moptions["wrkBase1"]))
   end_time = time.time();
   if moptions['outLevel']<=OUTPUT_INFO: print ("consuming time=%d" % (end_time-start_time))

   moptions['ds2'] = ['simulate_case', 'folder_control']

   allcasekeys = moptions['casefolder'].keys(); 
   allcontkeys = moptions['controlsample'].keys();

   if moptions['outLevel']<=OUTPUT_INFO: print ('INFO: Total subreads with modifications (%d) and without modifications (%d)' % (len(allcasekeys), len(allcontkeys)))

   start_time = time.time();
   if not moptions.has_key('PercDis'):
      moptions['PercDis'] = []

   cur_num_mod_subreads = moptions["CaseSize"]
   cur_num_unmod_subreads1 = int(cur_num_mod_subreads*(1-moptions['Percentage'])/moptions['Percentage'])
   cur_num_unmod_subreads2 = int(cur_num_mod_subreads/moptions['Percentage'])
   if moptions['outLevel']<=OUTPUT_WARNING: #OUTPUT_INFO: 
      print ('Will simulate %d subreads with modifications and %d subreads without modifications: total=%d, control subreads without modification is %d' % (cur_num_mod_subreads, cur_num_unmod_subreads1, cur_num_mod_subreads+cur_num_unmod_subreads1, cur_num_unmod_subreads2))

   for rt in range(moptions['random_times']):
      mcase_rand = np.random.choice(len(allcasekeys), cur_num_mod_subreads, replace=False);
      mcase1 = {allcasekeys[x]:moptions['casefolder'][allcasekeys[x]] for x in mcase_rand}

      con_random = np.random.choice(len(allcontkeys), cur_num_unmod_subreads1+cur_num_unmod_subreads2, replace=False);
      mcon1 = {allcontkeys[x]:moptions['controlsample'][allcontkeys[x]] for x in con_random[:cur_num_unmod_subreads1]}
      mcon2 = {allcontkeys[x]:moptions['controlsample'][allcontkeys[x]] for x in con_random[cur_num_unmod_subreads1:]}

      if moptions['outLevel']<=OUTPUT_WARNING: #OUTPUT_DEBUG: 
         print "randtime, len(case), len(con), len(test) ", rt, len(mcase1), len(mcon1), len(mcon2)

      if abs(len(mcase1)+len(mcon1) - len(mcon2)) > 10:
         if moptions['outLevel']<=OUTPUT_WARNING: print 'Warning!!! larger difference', len(mcase1), len(mcon1), len(mcon2), abs(len(mcase1)+len(mcon1) - len(mcon2) )

      moptions['sign_test'] = []
      moptions['sorted_sign_test'] = []

      moptions[moptions['ds2'][0]] = {}
      moptions[moptions['ds2'][1]] = {}

      moptions[moptions['ds2'][0]]['base'] = defaultdict(lambda: defaultdict(str))
      moptions[moptions['ds2'][1]]['base'] = defaultdict(lambda: defaultdict(str))

      moptions[moptions['ds2'][0]]['norm_mean'] = defaultdict(lambda: defaultdict(list))
      moptions[moptions['ds2'][1]]['norm_mean'] = defaultdict(lambda: defaultdict(list))

      mySimulate.getGenomeEvents([mcase1, mcon1, mcon2], data_labels, moptions)

      myDetect.mfilter_coverage(moptions)

      myDetect.mtest2(moptions)

      moptions['PercDis'].append(mySimulate.getTopRank(moptions))
      if moptions['outLevel']<=OUTPUT_INFO: print 'Rank ', moptions['Percentage'], rt, moptions['PercDis'][-1]
      end_time = time.time();
      if moptions['outLevel']<=OUTPUT_INFO: print ("consuming time=***%d*** for perc=%.2f_randtime=%d" % (end_time-start_time, moptions['Percentage'], rt))
      sys.stdout.flush()
   
   if moptions['outLevel']<=OUTPUT_INFO: print ''
   if moptions['outLevel']<=OUTPUT_WARNING: print moptions['Percentage'], moptions['PercDis']
   sys.stdout.flush()

   #cur_file_base = moptions['outFolder'] + '/' + moptions["FileID"] + ('_%.5f' % (moptions['Percentage'])) 
   cur_file_base = moptions['outFolder'] + '/' + moptions["FileID"]
   opfile = cur_file_base + '.output'
   mySaveRes(opfile, moptions) 
   os.system('touch '+cur_file_base+'.done')

def mySaveRes(opfile, moptions, format1='%d', format2=' %d'):
   fh = open(opfile, 'w');
   if isinstance(moptions['PercDis'], list):
      fh.write(format1 % moptions["CaseSize"])
      for currank in moptions['PercDis']:
         if int(currank)<0: 
            if moptions['outLevel']<=OUTPUT_INFO: print 'Negative rank', currank, int(currank)
            continue;
         
         fh.write(format2 % currank)
      fh.write('\n')
   else:
     allcasesize = moptions['PercDis'].keys(); allcasesize.sort()
     for curcs in allcasesize:
        fh.write(format1 % curcs)
	for currank in moptions['PercDis'][curcs]:
           if int(currank)<0:
              if moptions['outLevel']<=OUTPUT_INFO: print 'Negative rank', currank, int(currank)
	      continue;
           fh.write(format2 % currank)
	fh.write('\n')   
   fh.close();

def mSimulat2(moptions):
   if test: moptions['random_times'] = 3; #100; #3; #10
   else: moptions['random_times'] = 100
   moptions['outLevel']=OUTPUT_INFO

   if moptions['outLevel']<=OUTPUT_INFO: print 'Percentage', moptions['Percentage']
   random.seed(1)

   if moptions["runType"]==3:
      moptions['perpage']='9perpage'
      moptions['pdfortif'] = 'pdf'
      moptions['pdfortif'] = 'tif'
      mplotall(moptions)

   elif moptions["runType"]==2:
      mSimulate1(moptions)
      print '\n'
   elif moptions["runType"]==1:
      current_username = getpass.getuser();
      if len(moptions["FileID"])>9: job_name_base = moptions["FileID"][:10]
      else: job_name_base = moptions["FileID"]+'_'
      qscmd = 'qstat -u "'+current_username+'" | awk \'NR>2 && ($5=="r" || $5=="qw") && substr($3, 1, '+str(len(job_name_base))+')=="'+job_name_base+'" {print $4}\' | sort | wc -l'

      getSubFolders(moptions)
      if moptions['outLevel']<=OUTPUT_INFO: print moptions['subfolders'], '\n', moptions['subf_max_int']

      total_control_subreads = moptions['subf_max_int'][0]*4000
      #print moptions['subf_max_int'], total_control_subreads
      max_case_subreads = total_control_subreads*moptions['Percentage']/(2-moptions['Percentage'])
      start_case_subreads = 400; mstep = 300;
      start_case_subreads = 1000; mstep = 1000;
      if moptions['outLevel']<=OUTPUT_INFO: print ('Total control_subreads=%d, max_case_subreads=%d; Percentage=%.5f' % (total_control_subreads, max_case_subreads, moptions['Percentage']))

      max_case_subreads = int(max_case_subreads)

      moptions['CaseSizes'] = []
      total_jobs = defaultdict();
      for cur_case_subreads in range(start_case_subreads, max_case_subreads, mstep):
          if test and cur_case_subreads>=1100: break;
          moptions['CaseSizes'].append(cur_case_subreads)
          cur_file_id = ("%s_%d" % (moptions["FileID"], cur_case_subreads))
	  cmd = ("echo 'python modDetection.py simulat2 --CaseSize %d --testMethod %s --wrkBase1 %s --wrkBase2 %s --Percentage %.5f --FileID %s --outFolder %s' | qsub -V -cwd -l h_vmem=20G -N %s -e %s/sim2_%s.e -o %s/sim2_%s.o -q all.q" % (cur_case_subreads, moptions['testMethod'], moptions["wrkBase1"], moptions["wrkBase2"], moptions['Percentage'], cur_file_id, moptions['outFolder'], cur_file_id, moptions['outFolder'], cur_file_id, moptions['outFolder'], cur_file_id))
	  total_jobs[cur_file_id] = [cmd, ('%s/sim2_%s.e' % (moptions['outFolder'], cur_file_id)), ('%s/sim2_%s.o' % (moptions['outFolder'], cur_file_id)), ('%s/%s.output' % (moptions['outFolder'], cur_file_id)), ('%s/%s.done' % (moptions['outFolder'], cur_file_id))]

          for rmf in total_jobs[cur_file_id][1:]:
             if os.path.isfile(rmf):
                 os.system('rm '+rmf)

	  print total_jobs[cur_file_id]; sys.stdout.flush()
	  os.system(cmd)
	  time.sleep(0.5)

      time.sleep(10)
      if not moptions.has_key('PercDis'):
          moptions['PercDis'] = defaultdict(list)
      
      efile = ('%s/sim2_%s_merge.e' % (moptions['outFolder'], moptions["FileID"]))
      ofile = ('%s/sim2_%s_merge.o' % (moptions['outFolder'], moptions["FileID"]))
      opfile = ('%s/sim2_%s_merge.output' % (moptions['outFolder'], moptions["FileID"]))
      rmfiles = [efile, ofile, opfile]
      for rmf in rmfiles:
         if os.path.isfile(rmf):
            os.system('rm '+rmf)

      cur_file_base = moptions['outFolder'] + '/' + moptions["FileID"]
      opfile = cur_file_base + '_all.output'
      #if moptions['outLevel']<=OUTPUT_DEBUG: print moptions['PercDis']
      #mySaveRes(opfile, moptions, format2=' %s')


      print qscmd
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
                  curp = int(lsp[0])
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
            wait_time = mySimulate.reduce_wait_time(wait_time, max_wait_time, min_wait_time)
         else:
            wait_time = mySimulate.double_wait_time(wait_time, max_wait_time, min_wait_time)

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
               if tofile:
                  sys.stdout.write('>');
               else:
	          cur_hint = ('> %ds left' % (wait_time-i*3))
	          sys.stdout.write(cur_hint); #
	          sys.stdout.flush()
	          sys.stdout.write("\b"*(len(cur_hint)-1))
               sys.stdout.flush()
            print ''

      jobkeys = total_jobs.keys();
      runjobs = int(subprocess.check_output(qscmd, shell=True))
      print ('Total running jobs=%d(qstat=%d)' % (len(jobkeys), runjobs))+' ---'+moptions["wrkBase1"].split('/')[0], moptions["wrkBase2"].split('/')[0],'---Waiting time', wait_time

      #cur_file_base = moptions['outFolder'] + '/' + moptions["FileID"]
      #opfile = cur_file_base + '_all.output'
      if moptions['outLevel']<=OUTPUT_DEBUG: print moptions['PercDis']
      mySaveRes(opfile, moptions, format2=' %s')
      mplotHis(moptions);

def group_rank(moptions):
   #print moptions['PercDis']
   mgroupRanks = defaultdict(lambda: defaultdict(int))

   mybin, split_points, myRankStr = mySimulate.myBinDefault()

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

   #moptions["Percentages"].sort()
   #percvector = robjects.FloatVector(moptions["Percentages"])
   moptions['CaseSizes'].sort();
   csvector = robjects.IntVector(moptions['CaseSizes'])

   #mdfperc = robjects.DataFrame({"MixedPerc":robjects.FactorVector(robjects.FloatVector(perclist), levels=percvector, labels=percvector), "Rank":robjects.FactorVector(robjects.StrVector(rankgrouplist), levels=rankstrvector, labels=rankstrvector), "Fraction":robjects.FloatVector(rankperclist)})
   mdfperc = robjects.DataFrame({"MixedPerc":robjects.FactorVector(robjects.IntVector(perclist), levels=csvector, labels=csvector), "Percentile":robjects.FactorVector(robjects.StrVector(rankgrouplist), levels=rankstrvector, labels=rankstrvector), "Fraction":robjects.FloatVector(rankperclist)})

   robjects.r(resource_string(__name__, 'Rscript/Hist_sim_plot.R'))
   #robjects.r('pdf("'+mresfolder+'/hist_'+figname+'.pdf", width='+("%.0f" % (len(moptions["Percentages"])*1.5))+', height=4, onefile = TRUE)')
   robjects.r('pdf("'+mresfolder+'/hist2_'+figname+'.pdf", width='+("%.0f" % (len(moptions["CaseSizes"])*0.8))+', height=4, onefile = TRUE)')

   robjects.globalenv['Hist_sim_plot'](mdfperc, spvector, rankstrvector)

   robjects.r('dev.off()')


def mplotall(moptions):
   #moptions['outFolder'] = 'logsim_step1000_sm2/'


   #moptions["Percentages"] = [0.15,0.2,0.25,0.3,0.4,0.5,0.6,0.7]
   #control = {'vector':'vector_SpeI_cut', 'oligo':'ctrl_oligo_SpeI_cut'}
   #control = {'oligo':'ctrl_oligo_SpeI_cut'}
   #cases = {'Glucose':'Glucose_SpeI_cut', \
   #         'BrdU':'BrdU_SpeI_cut_data', \
   #         'A647':'A647_Click_SpeI_cut', \
   #         'Bio':'Bio_pur_SpeI_cut_data', \
   #         'Edu':'Edu_SpeI_cut_data', \
   #         'IdU':'IdU', \
   #         'A488':'A488_Click_SpeI_cut_data', \
   #         'Nicotine':'Nicotine_SpeI_cut', \
   #         'Azide':'Azide_SpeI_cut'}
   #conkeys = control.keys(); conkeys.sort()
   #caskeys = cases.keys(); caskeys.sort(); caskeys = caskeys[::-1]
   ##testmethod3 = ['fisher', 'stouffer', 'ks']
   #testmethod3 = ['ks', 'stouffer', 'fisher']

   control = {'ol':'ctrl_oligo_SpeI_cut'}
   cases = {'Id':'IdU', \
            'Ed':'Edu_SpeI_cut_data', \
            'Br':'BrdU_SpeI_cut_data', \
            'Bi':'Bio_pur_SpeI_cut_data', \
            'Ni':'Nicotine_SpeI_cut', \
            'Az':'Azide_SpeI_cut'}
   control = {'ol':'oligo'}
   cases = {'Id':'IdU', \
            'Ed':'Edu', \
            'Br':'BrdU', \
            'Bi':'Bio', \
            'Ni':'Nicotine', \
            'Az':'Azide', \
	    'Gl':'Glucose', \
	    'A6':'A647', \
	    'A4':'A488' }


   #Id Ed Br Bi Ni Az
   conkeys = control.keys(); conkeys.sort()
   caskeys = cases.keys(); caskeys.sort(); caskeys = caskeys[::-1]

   testmethod3 = {'st':'stouffer', 'ks':'ks', 'f':'fisher'}
   testmethod3keys = ['st', 'ks', 'f']

   for conk_ind in conkeys:
      conk = control[conk_ind]
      for perc in [0.2, 0.25, 0.15]:
         perclist_all = []
	 rankgrouplist_all = []
	 rankperclist_all = []
	 case_control_all = []
	 testmethod_all = []

         for testmethod_ind in testmethod3keys:
            testmethod = testmethod3[testmethod_ind]
            for cask_ind in caskeys:
               cask = cases[cask_ind]
               #moptions["FileID"] = ('%s_%s_%s' % (conk, cask, testmethod))
	       #logsim_step1000_sm2/olNist0.25_all.output
	       moptions["FileID"] = ('%s%s%s%.2f' % (conk_ind, cask_ind, testmethod_ind, perc))
               mresfile = ('%s/%s_all.output' % (moptions['outFolder'], moptions["FileID"]))
               moptions['PercDis'] = {}
               fr = open(mresfile, 'r')
               lines = fr.readlines()
               for l in lines:
                  lsp = l.split();
		  lsp[0] = int(lsp[0])/2
                  moptions['PercDis'][int(lsp[0])] = lsp[1:]
                  for i in range(len(moptions['PercDis'][int(lsp[0])])):
                     moptions['PercDis'][int(lsp[0])][i] = int(moptions['PercDis'][int(lsp[0])][i])
               fr.close()

               perclist, rankgrouplist, rankperclist, split_points, myRankStr = group_rank(moptions)

	       if True: # False: #True:
                  print '>>>', moptions["FileID"]
                  opdict = defaultdict(lambda: defaultdict(float));
                  max_rankgrp_len = 0
	          for nperc, rankgrp, rankperc in zip(perclist, rankgrouplist, rankperclist):
		     opdict[nperc][rankgrp] = rankperc
		     if len(rankgrp)>max_rankgrp_len: max_rankgrp_len = len(rankgrp)
                  print '   ',
	          for rankgrp in myRankStr:
                     print (('   %'+str(max_rankgrp_len)+'s') % rankgrp),
                  print ''
	          opdictkeys = opdict.keys(); opdictkeys.sort()
	          for nperc in opdictkeys:
                     print ('%d' % nperc),
		     for rankgrp in myRankStr:
                        print (('   %'+str(max_rankgrp_len)+'s') % ('%.3f' % opdict[nperc][rankgrp])),
                     print ''

	       perclist_all.extend(perclist)
	       rankgrouplist_all.extend(rankgrouplist)
	       rankperclist_all.extend(rankperclist);
	       #case_control_all.extend(['%s_%s' % (conk, cask)]*len(perclist))
	       if moptions['perpage']=='9perpage':
                  testmethod_all.extend([testmethod]*len(perclist));
                  case_control_all.extend(['%s_%s_%s' % (conk, cask, testmethod)]*len(perclist))
               elif moptions['perpage']=='27perpage':
                  testmethod_all.extend(['%s_%s_%s' % (conk, cask, testmethod)]*len(perclist))
		  case_control_all.extend(['%s_%s' % (conk, cask)]*len(perclist))
               else:
                  print 'need perpage parameters', moptions['perpage'] 
                  sys.exit()


         mresfolder = moptions['outFolder']

         ggplot = importr('ggplot2')
         importr('gridExtra')

         spvector = robjects.IntVector(split_points)
         rankstrvector = robjects.StrVector(myRankStr)

         #moptions["Percentages"].sort()
         #percvector = robjects.FloatVector(moptions["Percentages"])

	 #mdfperc = robjects.DataFrame({"MixedPerc":robjects.FactorVector(robjects.FloatVector(perclist_all), levels=percvector, labels=percvector), "Rank":robjects.FactorVector(robjects.StrVector(rankgrouplist_all), levels=rankstrvector, labels=rankstrvector), "Fraction":robjects.FloatVector(rankperclist_all), "CaseControl":robjects.StrVector(case_control_all), "Method":robjects.StrVector(testmethod_all)})
	 mdfperc = robjects.DataFrame({"MixedPerc":robjects.FactorVector(robjects.IntVector(perclist_all)), "Percentile":robjects.FactorVector(robjects.StrVector(rankgrouplist_all), levels=rankstrvector, labels=rankstrvector), "Fraction":robjects.FloatVector(rankperclist_all), "CaseControl":robjects.StrVector(case_control_all), "Method":robjects.StrVector(testmethod_all)})

         figname = conk+'_' + ('%.2f' % perc) + '_' + moptions['perpage']

         if moptions['pdfortif']=='pdf':
            robjects.r(resource_string(__name__, 'Rscript/Hist_sim_plot9.R'))
            xlabstr =  robjects.StrVector(["#subreads with modifications"])
            robjects.r('pdf("'+mresfolder+'/hist_'+figname+'.pdf", width='+("%.0f" % (10)+', height=6, onefile = TRUE)'))
            robjects.globalenv['Hist_sim_plot9'](mdfperc, spvector, rankstrvector, xlabstr)
            robjects.r('dev.off()')
         elif moptions['pdfortif']=='tif': 
            robjects.r(resource_string(__name__, 'Rscript/Hist_sim_plot9tif.R'))
            xlabstr =  robjects.StrVector(["#subreads with modifications"])
            robjects.globalenv['Hist_sim_plot9'](mdfperc, spvector, rankstrvector, xlabstr, mresfolder+'/hist_'+figname)
         else: print 'Error Wrong output file type', moptions['pdfortif']

         ##if len(sys.argv)<3:
         ##   print more_usage_str; #'Usages: python mySimulate.py more 9perpage/27perpage' 
         ##   sys.exit()
         ##figname = conk+'_' + sys.argv[2]
         ##if sys.argv[2]=='9perpage':
         ##   robjects.r(resource_string(__name__, 'Rscript/Hist_sim_plot9.R'))
         ##   robjects.r('pdf("'+mresfolder+'/hist_'+figname+'.pdf", width='+("%.0f" % (len(moptions["Percentages"])*2))+', height=9, onefile = TRUE)')
         ##
         ##   robjects.globalenv['Hist_sim_plot9'](mdfperc, spvector, rankstrvector)
         ##elif sys.argv[2]=='27perpage':
         ##   robjects.r(resource_string(__name__, 'Rscript/Hist_sim_plot27.R'))
	 ##   robjects.r('pdf("'+mresfolder+'/hist_'+figname+'.pdf", width='+("%.0f" % (len(moptions["Percentages"])*1.5))+', height=4, onefile = TRUE)')
	 ##   robjects.globalenv['Hist_sim_plot27'](mdfperc, spvector, rankstrvector)
         ##else:   
         ##   print more_usage_str; #'Usages: python mySimulate.py more 9perpage/27perpage'
	 ##   sys.exit()
         #
         #robjects.r('dev.off()')

      

