
import os;
import sys;
import copy
import string;
import math;

import h5py
import numpy as np

import time

import rpy2.robjects as robjects
from rpy2.robjects.packages import importr

from pkg_resources import resource_string

from scipy.stats import mannwhitneyu
from scipy.stats import ttest_ind
from scipy.stats import ks_2samp
from scipy.stats import combine_pvalues


from collections import defaultdict
from collections import deque

import myFast5
from myCom import *

has_boxplot=0;
has_ut = 0;

def mReadSignalBase(moptions):
   fn = moptions['fast5filename']
   if moptions['outLevel']<=OUTPUT_DEBUG: print "Read:", fn
   if not os.path.isfile(fn):
      if moptions['outLevel']<=OUTPUT_ERROR: print "Error!! no such file", fn
      return

   tocon = True;
   with h5py.File(fn, 'r') as mf5:
      #if fn[:len(moptions['ds2'][1])] == moptions['ds2'][1]:
      #   print fn, mf5, myFast5.rawAlignment_full
      #if mf5.__contains__(myFast5.rawAlignment_full):
      
      noerror = True;
      try:
         if not mf5.__contains__(myFast5.rawAlignment_full): noerror = False
      except:
         if moptions['outLevel']<=OUTPUT_ERROR: print ('Cannot find %s in %s' % (myFast5.rawAlignment_full, fn))
         noerror = False;
         #print 'Error'
      if noerror:

         mapped_chrom, mapped_start, mapped_strand = myFast5.ReadMapInfoInRef(mf5)

         #if not moptions[moptions['cur_wrkBase']].has_key((mapped_chrom, "-")):
         #   moptions[moptions['cur_wrkBase']][(mapped_chrom, "-")] = defaultdict(list) #{}
         #if not moptions[moptions['cur_wrkBase']].has_key((mapped_chrom, "+")):
         #   moptions[moptions['cur_wrkBase']][(mapped_chrom, "+")] = defaultdict(list) #{}

         #genome_fa = myFast5.ReadNanoraw_Genome_alignment(mf5)
         #genome_fa_pure = genome_fa.replace('-', '')
        
         #print fn,
         nanoraw_events = myFast5.ReadNanoraw_events(mf5)

         if moptions.has_key('Chr') and (not moptions["Chr"]==mapped_chrom):
            tocon = False;
         if tocon and moptions.has_key("start_pos") and moptions.has_key("end_pos"):
            if moptions.has_key("checkN") and moptions[moptions['cur_wrkBase']]['norm_mean'].has_key((moptions["Chr"], mapped_strand)):
               mposkeys = moptions[moptions['cur_wrkBase']]['norm_mean'][(moptions["Chr"], mapped_strand)].keys();
               if len(mposkeys)>0 and len(moptions[moptions['cur_wrkBase']]['norm_mean'][(moptions["Chr"], mapped_strand)][mposkeys[0]])>=moptions["checkN"]:
                  tocon = False;
            #if tocon and (mapped_start>moptions["start_pos"] or mapped_start+len(genome_fa_pure)<moptions["end_pos"]): 
            if tocon and (mapped_start>moptions["start_pos"] or mapped_start+len(nanoraw_events)<moptions["end_pos"]):
              tocon = False;
    
         #print mapped_chrom, mapped_start, mapped_strand    

         if tocon:
            #nanoraw_events = myFast5.ReadNanoraw_events(mf5)
            #read_al_fa = myFast5.ReadNanoraw_read_alignment(mf5)

            for i in range(len(nanoraw_events)):
               if mapped_strand=='+':
                  curpos = i + mapped_start
               else: curpos = mapped_start + len(nanoraw_events)-1 - i
               if moptions.has_key("start_pos") and moptions.has_key("end_pos"):
                  if curpos<moptions["start_pos"] or curpos>moptions["end_pos"]: 
                     continue;
               #if mapped_strand=='-': cur_n = na_bp[genome_fa_pure[i]]
               #else: cur_n = genome_fa_pure[i]
               #cur_n = genome_fa_pure[i]
               #if not moptions[moptions['cur_wrkBase']][(mapped_chrom, mapped_strand)].has_key(curpos):
               #   moptions[moptions['cur_wrkBase']][(mapped_chrom, mapped_strand)][curpos] = [[], nanoraw_events['base'][i]]
               #moptions[moptions['cur_wrkBase']][(mapped_chrom, mapped_strand)][curpos][0].append(nanoraw_events['norm_mean'][i])

               moptions[moptions['cur_wrkBase']]['base'][(mapped_chrom, mapped_strand)][curpos] = nanoraw_events['base'][i]
               moptions[moptions['cur_wrkBase']]['basedict'][(mapped_chrom, mapped_strand)][curpos][nanoraw_events['base'][i]] += 1
               moptions[moptions['cur_wrkBase']]['norm_mean'][(mapped_chrom, mapped_strand)][curpos].append(nanoraw_events['norm_mean'][i])
      else:
         #print "INFO: no alignment info for", myFast5.rawAlignment_full, myFast5.raw_event_ful
         if moptions['outLevel']<=OUTPUT_INFO: print "INFO: no alignment info for", fn  

def plot1(moptions, significant_pos, curn):
   m_signal = [] #deque() #[]
   m_pos = [] #deque() #[]
   m_ds = [] #deque() #[]

   curchr = significant_pos[0][0];
   curstrand = significant_pos[0][1];
   curpos = significant_pos[0][2];

   if moptions["neighborPvalues"]>0 and (not moptions["testMethod"]=="ks"):
      mtitle = ("1=%s VS\n 2=%s:\n p-value=%.1E (ks test p=%.1E) at pos %d of %s strand in %s. Rank %d " % (moptions['ds2'][0], moptions['ds2'][1], significant_pos[1][3][1], significant_pos[1][2][1], curpos+1, curstrand, curchr, curn+1))
   else:
      mtitle = ("1=%s VS\n 2=%s:\n p-value=%.1E at pos %d of %s strand in %s. Rank %d  " % (moptions['ds2'][0], moptions['ds2'][1], significant_pos[1][2][1], curpos+1, curstrand, curchr, curn+1))

   ds0 = moptions[moptions['ds2'][0]]
   ds1 = moptions[moptions['ds2'][1]]

   ds2 = [ds0, ds1]

   sk = (curchr, curstrand)
   noenough = False;
   pv3 = {}
   cur_ind = moptions['sign_test'].index(significant_pos)
   print significant_pos, cur_ind, curn
   nearybysize = moptions["window"]
   if moptions['RegionRankbyST']==1: nearybysize = int(nearybysize*2)
   #for mind in range(cur_ind-moptions["window"], cur_ind+moptions["window"]+1):
   for mind in range(cur_ind-nearybysize, cur_ind+nearybysize+1):
      if pos_check(moptions['sign_test'], cur_ind, mind): 
         #print len(moptions['sign_test']), cur_ind, mind
         pk = moptions['sign_test'][mind][0][2]
         pv = moptions['sign_test'][mind][1]
         pv3[(pk, ds0['base'][sk][pk])] = pv
      else:
         noenough = True;
      if noenough: break;
      for mds_ind in range(len(ds2)):
         mna = ds2[mds_ind]['base'][sk][pk]
         for sg in ds2[mds_ind]['norm_mean'][sk][pk]:
            m_ds.append("%d" % (mds_ind+1))
            if moptions["neighborPvalues"]>0 and (not moptions["testMethod"]=="ks"):
               if has_ut==1:
                  m_pos.append('%d/%s\n%.1E\n%.1E\n%.1E\n%.1E' % (pk+1, mna, pv[0][1], pv[1][1],pv[2][1],pv[3][1]))
               else:
                  m_pos.append('%d/%s\n%.1E\n%.1E' % (pk+1, mna, pv[2][1],pv[3][1])) 
            else:
               if has_ut==1:
                  m_pos.append('%d/%s\n%.1E\n%.1E\n%.1E' % (pk+1, mna, pv[0][1], pv[1][1],pv[2][1]))
               else:
                  m_pos.append('%d/%s\n%.1E' % (pk+1, mna, pv[2][1]))
            m_signal.append(round(sg,3))
 
   #for pk in range(curpos-moptions["window"], curpos+moptions["window"]+1):
   #   pv = None;
   #   if pk==curpos: pv = significant_pos[1]
   #   else:
   #      if ds1['norm_mean'].has_key(sk) and ds1['norm_mean'][sk].has_key(pk) and ds0['norm_mean'].has_key(sk) and ds0['norm_mean'][sk].has_key(pk):
   #         pv = getUtest(ds0['norm_mean'][sk][pk], ds1['norm_mean'][sk][pk])
   #   if pv==None:
   #      noenough = True;
   #   else:
   #      cur_comb_pv = get_fisher_comb_pvalues(moptions, significant_pos)
   #      if not cur_comb_pv==None:
   #         pv.append(cur_comb_pv)
   #      pv3[(pk, ds0['base'][sk][pk])] = pv
   #   if noenough: break;
   #
   #   for mds_ind in range(len(ds2)):
   #      mna = ds2[mds_ind]['base'][sk][pk]
   #      for sg in ds2[mds_ind]['norm_mean'][sk][pk]:
   #         m_ds.append("%d" % (mds_ind+1))
   #        if moptions["neighborPvalues"]>0:
   #            m_pos.append('%d/%s\n%.1E\n%.1E\n%.1E\n%.1E' % (pk+1, mna, pv[0], pv[1],pv[2],pv[3]))
   #         else:
   #            m_pos.append('%d/%s\n%.1E\n%.1E\n%.1E' % (pk+1, mna, pv[0], pv[1],pv[2]))
   #        m_signal.append(round(sg,3))
      
   if not noenough:
      closesize = moptions["neighborPvalues"]*2
      if moptions['RegionRankbyST']==1:
         closesize = moptions["window"]
         if closesize<1: closesize = 1

      #if significant_pos[0][1]=='-' and 3072-moptions["neighborPvalues"]*3<=significant_pos[0][2]<=3072+moptions["neighborPvalues"]*3:
      if significant_pos[0][1]=='-' and 3072-closesize<significant_pos[0][2]<3072+closesize:
         print 'Rank', curn+1, moptions["testMethod"], moptions["FileID"], significant_pos[0][0], significant_pos[0][1], significant_pos[0][2]+1, significant_pos[0][3]

      #poskeys = deque(); pvsp3 = [deque(), deque(), deque()]
      poskeys = []; pvsp3 = [[], [], [], []]
      #print 'pvsp3', pvsp3
      pv3keys = pv3.keys(); pv3keys.sort()
      for pv3k in pv3keys:
          if moptions["neighborPvalues"]>0 and (not moptions["testMethod"]=="ks"):
             print ('%d/%s' % (pv3k[0]+1, pv3k[1])), ('u=%.3E(%.3E) t=%.3E(%.3E) ks=%.3E(%.3E) pv5=%.3E(%.3E)' % (pv3[pv3k][0][1],pv3[pv3k][0][0], pv3[pv3k][1][1],pv3[pv3k][1][0],  pv3[pv3k][2][1],pv3[pv3k][2][0], pv3[pv3k][3][1],pv3[pv3k][3][0]))
          else:
             print ('%d/%s' % (pv3k[0]+1, pv3k[1])), ('u=%.3E(%.3E) t=%.3E(%.3E) ks=%.3E(%.3E)' % (pv3[pv3k][0][1],pv3[pv3k][0][0], pv3[pv3k][1][1],pv3[pv3k][1][0],  pv3[pv3k][2][1],pv3[pv3k][2][0]))
          poskeys.append('%d/%s' % (pv3k[0]+1, pv3k[1]))
          #pvsp3[0].append(pv3[pv3k][0])
          #pvsp3[1].append(pv3[pv3k][1])
          #pvsp3[2].append(pv3[pv3k][2])
          pvsp3[0].append(round(math.log10(pv3[pv3k][0][1]), 3))
          pvsp3[1].append(round(math.log10(pv3[pv3k][1][1]), 3))
          pvsp3[2].append(round(math.log10(pv3[pv3k][2][1]), 3))
          if moptions["neighborPvalues"]>0 and (not moptions["testMethod"]=="ks"):
             pvsp3[3].append(round(math.log10(pv3[pv3k][3][1]), 3))          
      print ''

      stu = {"Position":robjects.StrVector(poskeys), "Pvalue":robjects.FloatVector(pvsp3[0])}; stru = robjects.DataFrame(stu)
      stt = {"Position":robjects.StrVector(poskeys), "Pvalue":robjects.FloatVector(pvsp3[1])}; strt = robjects.DataFrame(stt)
      stks ={"Position":robjects.StrVector(poskeys), "Pvalue":robjects.FloatVector(pvsp3[2])}; strks= robjects.DataFrame(stks)
      if moptions["neighborPvalues"]>0 and (not moptions["testMethod"]=="ks"):
         stcb ={"Position":robjects.StrVector(poskeys), "Pvalue":robjects.FloatVector(pvsp3[3])}; 
      else:
         stcb ={"Position":robjects.StrVector([]), "Pvalue":robjects.FloatVector(pvsp3[3])};
      strcb= robjects.DataFrame(stcb)

      pydf = {"Signal":robjects.FloatVector(m_signal), "Position":robjects.StrVector(m_pos), "DS":robjects.FactorVector(robjects.StrVector(m_ds))}
      plotDat = robjects.DataFrame(pydf)

      mrtitle = robjects.StrVector([mtitle])
      mhasbox = robjects.IntVector([has_boxplot])

      sys.stdout.flush()
      robjects.globalenv['Base_Most_Significant_Plot'](plotDat, stru, strt, strks, strcb, mrtitle, mhasbox)

   return noenough

def mboxplot(moptions):
   figname = moptions["FileID"]
   mresfolder = moptions['outFolder']

   ggplot = importr('ggplot2')
   importr('gridExtra')
   importr('scales')

   robjects.r(resource_string(__name__, 'Rscript/Base_Most_Significant_Plot.R'))
   wdenlarge = 1.7
   if moptions['RegionRankbyST']==1: wdenlarge = 2.5
   if has_boxplot==1:
      robjects.r('pdf("'+mresfolder+'/rplot_'+figname+'.pdf", width='+("%.0f" % (moptions["window"]*wdenlarge))+', height=10, onefile = TRUE)')
   else:
      robjects.r('pdf("'+mresfolder+'/rplot_'+figname+'.pdf", width='+("%.0f" % (moptions["window"]*wdenlarge))+', height=4.5, onefile = TRUE)')

   closesize = moptions["neighborPvalues"]*2
   if moptions['RegionRankbyST']==1:
      closesize = moptions["window"]
      if closesize<1: closesize = 1

   curn = 0;
   output_pos = []
   for mostp in moptions['sorted_sign_test']:
      too_close_to_previous = False;
      for pre_pos in output_pos:
         if pre_pos[0]==mostp[0][0] and pre_pos[1]==mostp[0][1] and \
            abs(pre_pos[2]-mostp[0][2])<closesize: #moptions["neighborPvalues"]*2: 
            too_close_to_previous = True; break;
      if too_close_to_previous: continue;

      if not plot1(moptions, mostp, curn): 
         output_pos.append((mostp[0][0], mostp[0][1], mostp[0][2]))
         curn = curn + 1
      if curn==moptions["topN"]: break;

   robjects.r('dev.off()')

def mfilter_coverage(moptions):
   for dsn in moptions['ds2']:
      curds = moptions[dsn]['norm_mean']
      strandkeys = curds.keys(); strandkeys.sort()
      for sk in strandkeys:
         if moptions['outLevel']<=OUTPUT_DEBUG: print "Info: ", dsn, sk, len(curds[sk])
         poskeys = curds[sk].keys(); poskeys.sort();
         for pk in poskeys:
             if len(curds[sk][pk])<moptions["MinCoverage"]: 
                del curds[sk][pk]
                del moptions[dsn]['base'][sk][pk]
         if len(curds[sk])==0: 
            del curds[sk]
            del moptions[dsn]['base'][sk]


def m_min_float(fv):
   if fv<sys.float_info.min:
      return sys.float_info.min
   else: return fv;

def m_max_float(fv):
   if fv>sys.float_info.max:
      return sys.float_info.max
   else: return fv

def getUtest(a, b):
  st, pu = mannwhitneyu(a, b)
  pu = m_min_float(pu)
  stu = m_max_float(st);

  st, pt = ttest_ind(a, b, equal_var=False)
  pt = m_min_float(pt)
  stt = m_max_float(st);

  st, pks = ks_2samp(a, b)
  pks = m_min_float(pks)
  stks = m_max_float(st);

  return [(stu, pu), (stt, pt), (stks, pks)]

#@inline
def pos_check(mlist, i, j):
    if j<0 or j>len(mlist)-1: return False
    if i==j or \
       (mlist[i][0][0]==mlist[j][0][0] and mlist[i][0][1]==mlist[j][0][1] and i-j==mlist[i][0][2]-mlist[j][0][2]):
       return True;
    else: return False;

def combin_pvalues(moptions):
   for i in range(len(moptions['sign_test'])):
      comb_pv = get_combin_pvalue(moptions, i)
      if not comb_pv==None:
         moptions['sign_test'][i][1].append(comb_pv)

def get_combin_pvalue(moptions, i):
   if moptions["neighborPvalues"]>0 and len(moptions['sign_test'])>0:
      #enoughNeighbor = True;
      pvalue_neighbors = []
      for j in range(i-moptions["neighborPvalues"], i+moptions["neighborPvalues"]+1):
          if j<0 or j>len(moptions['sign_test'])-1 or (not pos_check(moptions['sign_test'], i, j)):
             #enoughNeighbor = False;
             #continue;
             pvalue_neighbors.append(1.0)
          else:
             pvalue_neighbors.append(moptions['sign_test'][j][1][2][1])

      #default: fisher, no weights
      if moptions["testMethod"]=='fisher':
         comb_p_st, comb_p_p = combine_pvalues(pvalue_neighbors)
      #stouffer, weights
      if moptions["testMethod"]=='stouffer':
         midweight = 100;
         mweights = [midweight]
         for k in range(moptions["neighborPvalues"]):
           mweights.insert(0, mweights[0]/moptions["WeightsDif"])
           mweights.append(mweights[-1]/moptions["WeightsDif"])
         comb_p_st, comb_p_p = combine_pvalues(pvalue_neighbors, method='stouffer', weights=mweights); #[1,2,4,2,1])

      comb_p_p = m_min_float(comb_p_p)
      comb_p_st = m_max_float(comb_p_st)

      return (comb_p_st, comb_p_p)

      #if not enoughNeighbor:
      #   return combine_pvalues(pvalue_neighbors)[1]
      #else:
      #   return 1.0
   else: 
      if moptions["neighborPvalues"] == 0: return moptions['sign_test'][i][1][2]
      else: return None;

def mtest2(moptions):
   ds0 = moptions[moptions['ds2'][0]]
   ds1 = moptions[moptions['ds2'][1]]

   strandkeys = ds0['norm_mean'].keys(); strandkeys.sort()
   #https://docs.python.org/3/library/collections.html#collections.deque
   #Indexed access is O(1) at both ends but slows to O(n) in the middle. For fast random access, use lists instead
   moptions['sign_test'] = [] #deque() #[]
   start_time = time.time();
   for sk in strandkeys:
      if ds1['norm_mean'].has_key(sk):
         poskeys = ds0['norm_mean'][sk].keys(); poskeys.sort();
         for pk in poskeys:
             if ds1['norm_mean'][sk].has_key(pk):
                if not ds1['base'][sk][pk] == ds0['base'][sk][pk]: 
                   if moptions['outLevel']<=OUTPUT_ERROR:
                      print 'Error not equal', sk, pk, ds1['base'][sk][pk], ds0['base'][sk][pk], ds1['basedict'][sk][pk].items(), ds0['basedict'][sk][pk].items()
                #moptions['sign_test'].append([(sk[0], sk[1], pk, ds1['base'][sk][pk]), getUtest(ds0['norm_mean'][sk][pk], ds1['norm_mean'][sk][pk])])
                moptions['sign_test'].append(((sk[0], sk[1], pk, ds1['base'][sk][pk]), getUtest(ds0['norm_mean'][sk][pk], ds1['norm_mean'][sk][pk])))
   end_time = time.time();
   if moptions['outLevel']<=OUTPUT_INFO: print ("Producing pvalues: consuming time %d" % (end_time-start_time))

   start_time = time.time();
   if not moptions["testMethod"]=="ks": combin_pvalues(moptions)
   end_time = time.time();
   if moptions['outLevel']<=OUTPUT_DEBUG: print ("Combining pvalues: consuming time %d" % (end_time-start_time))

   if moptions['rankUse']=='pv':
      rank_use_p_or_st = 'pvalue'; use_pind=1;
   else:
      rank_use_p_or_st = 'st'; use_pind=0; 

   sorted_ind = 2;
   if not moptions["testMethod"]=="ks": #moptions["neighborPvalues"]>0 and len(moptions['sign_test'])>0: 
      sorted_ind = 3
   start_time = time.time();

   save_test(moptions)

   if moptions['RegionRankbyST']==0:
      moptions['sorted_sign_test'] = sorted(moptions['sign_test'], key=lambda mpv: (mpv[1][sorted_ind][use_pind], mpv[1][2][use_pind], mpv[1][0][use_pind]));
      if rank_use_p_or_st=='st':
         moptions['sorted_sign_test'] = moptions['sorted_sign_test'][::-1]
   else:
      windseg = [];
      moptions['window'] = moptions['window'] + 1
      windlist = range(-moptions['window'], moptions['window']+1)
      if moptions['outLevel']<=OUTPUT_INFO: print 'windlist', windlist, len(moptions['sign_test']), moptions['window']
      movesize = (moptions['window'])
      if moptions['WindOvlp'] == 1: movesize = 1

      strand_specific = defaultdict(lambda: defaultdict(set))
      strans_specific_max = defaultdict(int)
      for mp in moptions['sign_test']:
         strand_specific[(mp[0][0], mp[0][1])][mp[0][2]] = (mp[0][3], mp[1])
         strans_specific_max[(mp[0][0], mp[0][1])] = mp[0][2]
      strandskeys = strand_specific.keys(); strandskeys.sort()
      for sk in strandskeys:
         if moptions['outLevel']<=OUTPUT_INFO:
            print sk, strans_specific_max[sk], len(strand_specific[sk])
            print strand_specific[sk][3078][0], strand_specific[sk][3079][0],strand_specific[sk][3080][0]
      for sk in strandskeys:
         for pk in range(0, strans_specific_max[sk], movesize):
            pvlist = []; not50 = True
            for wind in windlist:
               curposk = pk + wind
               if curposk<0 or curposk>=strans_specific_max[sk] or (not strand_specific[sk].has_key(curposk)):
                  not50 = True;  break;
               not50 = False;
               if (not moptions.has_key('NA')) or len(moptions['NA'])==0 or moptions['NA']==strand_specific[sk][curposk][0]:
                  pvlist.append(strand_specific[sk][curposk][1][sorted_ind][use_pind])
            if not not50:
               opvlist = copy.deepcopy(pvlist)
               pvlist.sort()
               if len(pvlist)>5:
                  windseg.append(((sk[0], sk[1], pk, strand_specific[sk][pk][0]), pvlist, opvlist))
      moptions['sorted_sign_test'] = []
      windseg_sort = sorted(windseg, key=lambda mpv: (mpv[1][int(moptions['percentile']*(len(mpv[1])-1)+0.5)], abs(moptions['window']-mpv[2].index(mpv[1][0]))))

      if False: #True:
         for r1 in range(10):
            print windseg_sort[r1][0], windseg_sort[r1][1][int(moptions['percentile']*(len(windseg_sort[r1][1])-1)+0.5)], moptions['percentile'], moptions['percentile']*(len(windseg_sort[r1][1])-1)+0.5, len(windseg_sort[r1][1]), abs(moptions['window']-windseg_sort[r1][2].index(windseg_sort[r1][1][0])), windseg_sort[r1][2].index(windseg_sort[r1][1][0])

      if moptions['WindOvlp'] == 1: ovlap_higher_rank = []
      for rmp_ind in range(len(windseg_sort)):
         if moptions['WindOvlp'] == 1:
            closeneighbor = False;
            for pre_ind in ovlap_higher_rank:
               if windseg_sort[pre_ind][0][0]==windseg_sort[rmp_ind][0][0] and windseg_sort[pre_ind][0][1]==windseg_sort[rmp_ind][0][1] and abs(windseg_sort[pre_ind][0][2]-windseg_sort[rmp_ind][0][2])<(moptions['window']): #*2:
                  closeneighbor = True; break;
            if closeneighbor: continue;
            ovlap_higher_rank.append(rmp_ind)
         moptions['sorted_sign_test'].append((windseg_sort[rmp_ind][0],strand_specific[(windseg_sort[rmp_ind][0][0], windseg_sort[rmp_ind][0][1])][windseg_sort[rmp_ind][0][2]][1]))
 
   end_time = time.time();
   if moptions['outLevel']<=OUTPUT_DEBUG: print ("Sorting according to pvalues: consuming time %d" % (end_time-start_time))

   if moptions['outLevel']<=OUTPUT_DEBUG: print "Info in sign_test", len(moptions['sorted_sign_test'])

def save_test(moptions):
   #print 'SaveTest',  moptions['SaveTest']
   if moptions['SaveTest']==0: return;

   txtfile = moptions['outFolder'] + '/' + moptions["FileID"] + '_sign_test.txt'
   if moptions['outLevel']<=OUTPUT_ERROR: print 'Test data is saved in', txtfile
   txtwriter = open(txtfile, 'w')
   for mostp in moptions['sign_test']:
      txtwriter.write('%s %s %d %s %.3f %.3E %.3f %.3E %.3f %.3E' % (mostp[0][0], mostp[0][1], mostp[0][2], mostp[0][3], mostp[1][0][0], mostp[1][0][1], mostp[1][1][0], mostp[1][1][1], mostp[1][2][0], mostp[1][2][1]))
      if (moptions["neighborPvalues"]>0 and (not moptions["testMethod"]=="ks")):
         txtwriter.write(' %.3f %.3E\n' % (mostp[1][3][0], mostp[1][3][1]))
      else:
         txtwriter.write('\n')


def ReadAllFast5(moptions):
   #neighbors =  moptions["window"]

   if moptions.has_key("Pos"):
      neighbors =  moptions["window"]
      pos_of_interest = moptions["Pos"]
      start_pos = pos_of_interest - neighbors # -1
      if start_pos<0: start_pos = 0
      end_pos = pos_of_interest + neighbors
      moptions['start_pos'] = start_pos
      moptions["end_pos"] = end_pos
      print moptions["Pos"], moptions["window"], moptions['start_pos'], moptions["end_pos"]

   f5suf = moptions['.fast5']

   moptions['ds2'] = [moptions["wrkBase1"], moptions["wrkBase2"]]
   for cur_wrkBase_ind in range(len(moptions['ds2'])):
      start_time = time.time();
      cur_wrkBase = moptions['ds2'][cur_wrkBase_ind]
      moptions['cur_wrkBase'] = cur_wrkBase
      #moptions[moptions['cur_wrkBase']] = {}

      moptions[moptions['cur_wrkBase']] = {}
      moptions[moptions['cur_wrkBase']]['base'] = defaultdict(lambda: defaultdict(str))
      moptions[moptions['cur_wrkBase']]['norm_mean'] = defaultdict(lambda: defaultdict(list))
      moptions[moptions['cur_wrkBase']]['basedict'] = defaultdict(lambda: defaultdict(lambda: defaultdict(int)))
      #moptions[moptions['cur_wrkBase']]['norm_mean'] = defaultdict(lambda: defaultdict(deque))

      f5num = 0;
      f5sub = [cur_wrkBase]
      while len(f5sub)>0:
         f5subnew = [] #deque() #[]
         if moptions['outLevel']<=OUTPUT_WARNING: print '.sub_fast5_folder', f5sub
         for cursub in f5sub:
            f5subadd, f5num= readsubfolder(cursub, moptions, f5num, cur_wrkBase_ind, start_time, f5suf)
            f5subnew.extend(f5subadd)
         f5sub = f5subnew


def readsubfolder(cursub, moptions, f5num, cur_wrkBase_ind, start_time, f5suf='.fast5'):
   f5sub = [] #deque() #[]
   f5list = os.listdir(cursub)
   for f5_ind in range(len(f5list)):
     f5 = f5list[f5_ind]
     if moptions.has_key("checkN") and moptions[moptions['cur_wrkBase']]['norm_mean'].has_key((moptions["Chr"], '-')) and moptions[moptions['cur_wrkBase']]['norm_mean'].has_key((moptions["Chr"], '+')): 
        mposkeys1 = moptions[moptions['cur_wrkBase']]['norm_mean'][(moptions["Chr"], '-')].keys();
        mposkeys2 = moptions[moptions['cur_wrkBase']]['norm_mean'][(moptions["Chr"], '+')].keys();
        if len(mposkeys1)>0 and len(moptions[moptions['cur_wrkBase']]['norm_mean'][(moptions["Chr"], '-')][mposkeys1[0]])>=moptions["checkN"] and len(mposkeys2)>0 and len(moptions[moptions['cur_wrkBase']]['norm_mean'][(moptions["Chr"], '+')][mposkeys2[0]])>=moptions["checkN"]:
           break;

     if f5[-len(f5suf):]==f5suf:
        moptions['fast5filename'] = cursub+'/'+f5
        if f5num>0 and f5num%1000==0 and moptions['outLevel']<=OUTPUT_WARNING:
           chrkeys = moptions[moptions['cur_wrkBase']]['norm_mean'].keys(); chrkeys.sort()
           print 'chrkeys', chrkeys
           mposkeys1 = moptions[moptions['cur_wrkBase']]['norm_mean'][chrkeys[0]].keys();
           mposkeys2 = moptions[moptions['cur_wrkBase']]['norm_mean'][chrkeys[1]].keys();

           print cur_wrkBase_ind, f5num, f5[-50:], 
           if len(mposkeys1)>0: 
              print chrkeys[0], len(mposkeys1), mposkeys1[0], len(moptions[moptions['cur_wrkBase']]['norm_mean'][chrkeys[0]][mposkeys1[0]]), 
           else: print 0, 
           if len(mposkeys2)>0:
              print chrkeys[1], len(mposkeys2), mposkeys2[0], len(moptions[moptions['cur_wrkBase']]['norm_mean'][chrkeys[1]][mposkeys2[0]]),
           else: print 0,
           end_time = time.time();
           print ("consuming time=%d" % (end_time-start_time)), 
           #start_time = end_time
           if moptions.has_key("checkN"): print moptions["checkN"]
           else: print ''
           sys.stdout.flush()
        #print cursub+'/'+f5,
        mReadSignalBase(moptions)
        f5num = f5num + 1
     elif os.path.isdir(cursub+'/'+f5): 
        if not f5=='mall':
           f5sub.append(cursub+'/'+f5)
   if len(f5sub)>0:
     print 'Under', cursub
     print '\t find', f5sub
   return [f5sub, f5num]

def mDetect(moptions):
   print myFast5.rawAlignment_full
   ReadAllFast5(moptions)

   mfilter_coverage(moptions)

   mtest2(moptions)

   mboxplot(moptions)


