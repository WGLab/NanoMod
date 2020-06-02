
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

import myDetect
import mySimulate


target_pos = 3072
target_strand = '-'
target_chr = 'spel'
data_labels = ['simulate_case', 'folder_control']
tofile = False;
tofile = True;
test=False;
#test=True;


def mSimulate1(moptions):
   if moptions['outLevel']<=OUTPUT_INFO: print 'read fast5', moptions["wrkBase1"], moptions["wrkBase2"]
   start_time = time.time();
   moptions['casefolder'] = defaultdict(set)
   noevent, notfast5 = mySimulate.readEvents(moptions["wrkBase1"], 'casefolder', moptions)
   if moptions['outLevel']<=OUTPUT_WARNING: print ('Read case %d fast5 files(available=%d, inavail=%d), %d non-fast5 files for %s' % (len(moptions['casefolder']), len(moptions['casefolder'])-noevent, noevent, notfast5, moptions["wrkBase1"]))

   moptions['controlsample'] = defaultdict(set)
   noevent, notfast5 = mySimulate.readEvents(moptions["wrkBase2"], 'controlsample', moptions)
   if moptions['outLevel']<=OUTPUT_WARNING: print ('Read control %d fast5 files(available=%d, inavail=%d), %d non-fast5 files for %s' % (len(moptions['controlsample']), len(moptions['controlsample'])-noevent, noevent, notfast5, moptions["wrkBase2"]))
   end_time = time.time();
   if moptions['outLevel']<=OUTPUT_INFO: print ("consuming time=%d" % (end_time-start_time))

   moptions['ds2'] = data_labels

   start_time = time.time();
   if not moptions.has_key('CovgDis'):
      moptions['CovgDis'] = []

   allcasekeys = moptions['casefolder'].keys();
   allcontkeys = moptions['controlsample'].keys();

   #for rt in range(moptions['random_times']):
   rt = 0; repeat_time = 0; cur_repeat_time = 0;
   while rt<moptions['random_times']:
      moresample = repeat_time
      if moresample>15: moresample = 15;

      case_sample_num = int(moptions['CaseSize']*(1+moresample*0.02))
      if len(allcasekeys) > case_sample_num:
          mcase_rand = np.random.choice(len(allcasekeys), case_sample_num, replace=False)
          mcase1 = {allcasekeys[x]:moptions['casefolder'][allcasekeys[x]] for x in mcase_rand}
      else:
          mcase1 =  moptions['casefolder']

      cont_sample_num = int(moptions['CaseSize']*(1+moresample*0.02))
      if len(allcontkeys) > cont_sample_num:
          con_random = np.random.choice(len(allcontkeys), cont_sample_num, replace=False)
          mcon1 = {allcontkeys[x]:moptions['controlsample'][allcontkeys[x]] for x in con_random}
      else:
          mcon1 = moptions['controlsample']

      if moptions['outLevel']<=OUTPUT_WARNING: #OUTPUT_DEBUG:
         print "randtime, len(case), len(con)) ", rt, len(mcase1), len(mcon1)

      moptions['sign_test'] = []
      moptions['sorted_sign_test'] = []

      moptions[moptions['ds2'][0]] = {}
      moptions[moptions['ds2'][1]] = {}
      moptions[moptions['ds2'][0]]['base'] = defaultdict(lambda: defaultdict(str))
      moptions[moptions['ds2'][1]]['base'] = defaultdict(lambda: defaultdict(str))
      moptions[moptions['ds2'][0]]['norm_mean'] = defaultdict(lambda: defaultdict(list))
      moptions[moptions['ds2'][1]]['norm_mean'] = defaultdict(lambda: defaultdict(list))
      mySimulate.getGenomeEvents([mcase1, mcon1], data_labels, moptions)

      has_enough_coverage_for_region_of_interest = 0;
      coverage_list = []
      for cur_wrkBase_ind in range(len(moptions['ds2'])):
         cur_wrkBase = moptions['ds2'][cur_wrkBase_ind]
         chrstr = moptions[cur_wrkBase]['norm_mean'].keys();
         for curcs in chrstr:
            if curcs[1] == target_strand and curcs[0]==target_chr:
               for pos in range(target_pos-3, target_pos+4):
                  if len(moptions[cur_wrkBase]['norm_mean'][curcs][pos])<(0.95*moptions['CaseSize']/5):
                     has_enough_coverage_for_region_of_interest += 1;
                  coverage_list.append('%d:%d' % (pos, len(moptions[cur_wrkBase]['norm_mean'][curcs][pos])))
      if moptions['outLevel']<=OUTPUT_WARNING:
         print rt, 'Coverage list', ' '.join(coverage_list), "o", repeat_time, int(0.95*moptions['CaseSize']/5), has_enough_coverage_for_region_of_interest, cur_repeat_time
         sys.stdout.flush()
      if has_enough_coverage_for_region_of_interest > 2:
         print '\tNo enough coverage', rt
         if has_enough_coverage_for_region_of_interest > 3 and cur_repeat_time > 5:
            repeat_time += 1;
         cur_repeat_time += 1
         continue;

      myDetect.mfilter_coverage(moptions)
      myDetect.mtest2(moptions)

      moptions['CovgDis'].append(mySimulate.getTopRank(moptions))
      if moptions['outLevel']<=OUTPUT_INFO: print 'Rank ', moptions['CaseSize'], rt, moptions['CovgDis'][-1]
      end_time = time.time();
      if moptions['outLevel']<=OUTPUT_INFO: print ("consuming time=***%d*** for cov=%d_randtime=%d" % (end_time-start_time, moptions['CaseSize'], rt))
      sys.stdout.flush()
      rt = rt + 1; cur_repeat_time = 0; #repeat_time=0

   if moptions['outLevel']<=OUTPUT_INFO: print ''
   if moptions['outLevel']<=OUTPUT_WARNING: print moptions['CaseSize'], moptions['CovgDis']
   sys.stdout.flush()

   cur_file_base = moptions['outFolder'] + '/' + moptions["FileID"]
   opfile = cur_file_base + '.output'
   mySaveRes(opfile, moptions)
   os.system('touch '+cur_file_base+'.done')

def mySaveRes(opfile, moptions, format1='%d', format2=' %d'):
   fh = open(opfile, 'w');
   if isinstance(moptions['CovgDis'], list):
      fh.write(format1 % moptions["CaseSize"])
      for currank in moptions['CovgDis']:
         if int(currank)<0:
            if moptions['outLevel']<=OUTPUT_INFO: print 'Negative rank', currank, int(currank)
            continue;

         fh.write(format2 % currank)
      fh.write('\n')
   else:
     allcasesize = moptions['CovgDis'].keys(); allcasesize.sort()
     for curcs in allcasesize:
        fh.write(format1 % curcs)
	for currank in moptions['CovgDis'][curcs]:
           if int(currank)<0:
              if moptions['outLevel']<=OUTPUT_INFO: print 'Negative rank', currank, int(currank)
	      continue;
           fh.write(format2 % currank)
	fh.write('\n')
   fh.close();

def mSimulat2(moptions):
   if test: moptions['random_times'] = 3; #100; #3; #10
   else: moptions['random_times'] = 100
   #moptions['outLevel']=OUTPUT_INFO

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
      ##full: "qstat -xml | tr '\n' ' ' | sed 's#<job_list[^>]*>#\n#g' | sed 's#<[^>]*>##g' | grep ' ' | column -t"
      qscmd = 'qstat -u "'+current_username+'" | awk \'NR>2 && ($5=="r" || $5=="qw") && substr($3, 1, '+str(len(job_name_base))+')=="'+job_name_base+'" {print $4}\' | sort | wc -l'

      moptions['CaseSizes'] = []
      moptions['CaseSizes'].append(60);
      moptions['CaseSizes'].append(80);
      moptions['CaseSizes'].append(100);
      moptions['CaseSizes'].append(200);
      moptions['CaseSizes'].append(400);
      moptions['CaseSizes'].append(1000);
      moptions['CaseSizes'].append(2000);
      moptions['CaseSizes'].append(3000);


      total_jobs = defaultdict();
      for cur_case_subreads in moptions['CaseSizes']:
          if test and cur_case_subreads>=1100: break;
          cur_file_id = ("%s_%d" % (moptions["FileID"], cur_case_subreads))
	  cmd = ("echo 'python modDetection.py DownSampling --CaseSize %d --testMethod %s --wrkBase1 %s --wrkBase2 %s --FileID %s --outFolder %s --SaveTest %d --neighborPvalues %d --RegionRankbyST %d --percentile %.3f --WindOvlp %d --window %d' | qsub -V -cwd -l h_vmem=10G -N %s -e %s/ds_%s.e -o %s/ds_%s.o -q all.q" % (cur_case_subreads, moptions['testMethod'], moptions["wrkBase1"], moptions["wrkBase2"], cur_file_id, moptions['outFolder'], moptions['SaveTest'], moptions['neighborPvalues'], moptions['RegionRankbyST'], moptions['percentile'], moptions['WindOvlp'], (moptions['window']+1)*2, cur_file_id, moptions['outFolder'], cur_file_id, moptions['outFolder'], cur_file_id))
	  total_jobs[cur_file_id] = [cmd, ('%s/ds_%s.e' % (moptions['outFolder'], cur_file_id)), ('%s/ds_%s.o' % (moptions['outFolder'], cur_file_id)), ('%s/%s.output' % (moptions['outFolder'], cur_file_id)), ('%s/%s.done' % (moptions['outFolder'], cur_file_id))]

          for rmf in total_jobs[cur_file_id][1:]:
             if os.path.isfile(rmf):
                 os.system('rm '+rmf)

	  print total_jobs[cur_file_id]; sys.stdout.flush()
	  os.system(cmd)
	  time.sleep(0.5)

      time.sleep(10)
      if not moptions.has_key('CovgDis'):
          moptions['CovgDis'] = defaultdict(list)

      efile = ('%s/ds_%s_merge.e' % (moptions['outFolder'], moptions["FileID"]))
      ofile = ('%s/ds_%s_merge.o' % (moptions['outFolder'], moptions["FileID"]))
      opfile = ('%s/ds_%s_merge.output' % (moptions['outFolder'], moptions["FileID"]))
      rmfiles = [efile, ofile, opfile]
      for rmf in rmfiles:
         if os.path.isfile(rmf):
            os.system('rm '+rmf)

      cur_file_base = moptions['outFolder'] + '/' + moptions["FileID"]
      opfile = cur_file_base + '_all.output'
      #if moptions['outLevel']<=OUTPUT_DEBUG: print moptions['CovgDis']
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
		  moptions['CovgDis'][curp].extend(lsp[1:])

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
      if moptions['outLevel']<=OUTPUT_DEBUG: print moptions['CovgDis']
      mySaveRes(opfile, moptions, format2=' %s')
      mplotHis(moptions);

def group_rank(moptions):
   #print moptions['CovgDis']
   mgroupRanks = defaultdict(lambda: defaultdict(int))

   mybin, split_points, myRankStr = mySimulate.myBinDefault()

   percKeys = moptions['CovgDis'].keys();
   for pk in percKeys:
      for curr in moptions['CovgDis'][pk]:
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

   moptions['CaseSizes'].sort();
   csvector = robjects.IntVector(moptions['CaseSizes'])

   #mdfperc = robjects.DataFrame({"MixedPerc":robjects.FactorVector(robjects.FloatVector(perclist), levels=percvector, labels=percvector), "Rank":robjects.FactorVector(robjects.StrVector(rankgrouplist), levels=rankstrvector, labels=rankstrvector), "Fraction":robjects.FloatVector(rankperclist)})
   mdfperc = robjects.DataFrame({"MixedPerc":robjects.FactorVector(robjects.IntVector(perclist), levels=csvector, labels=csvector), "Percentile":robjects.FactorVector(robjects.StrVector(rankgrouplist), levels=rankstrvector, labels=rankstrvector), "Fraction":robjects.FloatVector(rankperclist)})

   robjects.r(resource_string(__name__, 'Rscript/Hist_sim_plot.R'))
   robjects.r('pdf("'+mresfolder+'/hist2_'+figname+'.pdf", width='+("%.0f" % (len(moptions["CaseSizes"])*0.8))+', height=4, onefile = TRUE)')

   robjects.globalenv['Hist_sim_plot'](mdfperc, spvector, rankstrvector)

   robjects.r('dev.off()')


def mplotall(moptions):
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
   testmethod3keys = ['st', 'ks']

   for conk_ind in conkeys:
         conk = control[conk_ind]
         perclist_all = []
	 rankgrouplist_all = []
	 rankperclist_all = []
	 case_control_all = []
	 testmethod_all = []

         for testmethod_ind in testmethod3keys:
            testmethod = testmethod3[testmethod_ind]
            for cask_ind in caskeys:
               cask = cases[cask_ind]
	       moptions["FileID"] = ('%s%s%s%s' % (moptions['mprefix'], conk_ind, cask_ind, testmethod_ind))
               mresfile = ('%s/%s_all.output' % (moptions['outFolder'], moptions["FileID"]))
               moptions['CovgDis'] = {}
               fr = open(mresfile, 'r')
               lines = fr.readlines()
               for l in lines:
                  lsp = l.split();
		  lsp[0] = int(lsp[0])/2
                  moptions['CovgDis'][int(lsp[0])] = lsp[1:]
                  for i in range(len(moptions['CovgDis'][int(lsp[0])])):
                     moptions['CovgDis'][int(lsp[0])][i] = int(moptions['CovgDis'][int(lsp[0])][i])
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

	 mdfperc = robjects.DataFrame({"MixedPerc":robjects.FactorVector(robjects.IntVector(perclist_all)), "Percentile":robjects.FactorVector(robjects.StrVector(rankgrouplist_all), levels=rankstrvector, labels=rankstrvector), "Fraction":robjects.FloatVector(rankperclist_all), "CaseControl":robjects.StrVector(case_control_all), "Method":robjects.StrVector(testmethod_all)})

         figname = conk+ '_' + moptions['perpage']

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
