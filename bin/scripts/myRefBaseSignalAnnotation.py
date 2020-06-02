
import os;
import sys;
import string;
import glob;
import time
import copy

import h5py
import numpy as np
import multiprocessing

from collections import defaultdict
from distutils.version import LooseVersion

import tempfile
import subprocess

import re;

import myCom

fast5_channel_id= 'UniqueGlobalKey/channel_id'
fast5_analysis = ''.join(['/', myCom.analyses_base]) #
fast5_events = myCom.basecall_events_base #
fast5_move = myCom.basecall_move_base
fast5_rawReads = ''.join(['/', myCom.raw_base, '/', myCom.reads_base]) #
fast5_basecall_fq = myCom.basecall_fastq_base #
fast5_signal = myCom.signal_base #


#for line 1157:
moresignalperc = 0.3 #0.4

#fast5_mo_correction = "NanomoCorrected_000"

def get_channel_info(moptions, sp_param):
   if not sp_param['f5status']=="": return;
   try:
      channel_info = sp_param['f5reader'][fast5_channel_id].attrs
      sp_param["channel_info"] = {'digitisation':channel_info['digitisation'], 'offset':channel_info['offset'], 'range':channel_info['range'], 'sampling_rate':channel_info['sampling_rate'], 'channel_number':channel_info['channel_number']}
   except:
      raiseError("No Channel Info", sp_param, "No Channel Info")

def raiseError(sp_info, sp_param, errk):
   sp_param['f5status'] = errk
   print ('Error!!! %s in %s' % (sp_info, sp_param['mfile_path']))
   #raise RuntimeError, ('No %s in %s' % (sp_info, moptions['mfile_path']))

def getAlbacoreVersion(moptions, sp_param):
   if not sp_param['f5status']=="": return;
   try:
      if 'Guppy' in sp_param['f5reader'][''.join([fast5_analysis,'/',moptions['basecall_1d']])].attrs['name']:
          sp_param['basecaller'] = 'Guppy'
      used_version = LooseVersion(sp_param['f5reader'][''.join([fast5_analysis,'/',moptions['basecall_1d']])].attrs['version'] if 'version' in sp_param['f5reader'][''.join([fast5_analysis,'/',moptions['basecall_1d']])].attrs else "0.0")
      if used_version < LooseVersion("1.0"): #
         #raiseError("Not supported albacore version", sp_param, "Not supported albacore version")
         sp_param['used_albacore_version'] = 1;
      elif used_version < LooseVersion("2.0"): sp_param['used_albacore_version'] = 1;
      elif used_version >= LooseVersion("2.0"): sp_param['used_albacore_version'] = 2;
   except:
      sp_param['used_albacore_version'] = 1;
      #raiseError(''.join(['No version in', fast5_analysis,'/',moptions['basecall_1d']]), sp_param, "No version")

#https://github.com/jts/nanopolish/tree/master/etc/r9-models
def get_kmer_corrected_info(moptions):
   if moptions['kmer_model_file']==None or (not os.path.isfile(moptions['kmer_model_file'])): return;

   fr = open(moptions['kmer_model_file'], 'r')
   moptions['kmer_model_dict'] = defaultdict()
   line = fr.readline();
   while line:
      line = string.strip(line);
      if len(line)>0 and (not line[0]=='#'):
         try:
            c_kmer, c_level_mean, c_level_stdv = line.split()[:3]
            c_level_mean, c_level_stdv = float(c_level_mean), float(c_level_stdv)
            moptions['kmer_model_dict'][c_kmer] = (c_level_mean, 1/(c_level_stdv*c_level_stdv))
         except:
            pass;
      line = fr.readline();
   fr.close();


# https://community.nanoporetech.com/posts/squiggle-plot-for-raw-data
def get_cur_shift_scale(moptions, sp_param):
   if not sp_param['f5status']=="": return;
   if not moptions.has_key("kmer_model_dict"): return;

   event_key = 'm_event'
   event_key = 'events_data'

   try:
      cur_model = np.array([moptions['kmer_model_dict'][c_model_state] for c_model_state in sp_param[event_key]['model_state']], dtype=[('level_mean', np.float), ('level_stdv', np.float)]);
      c_mean_stdv = cur_model['level_mean']*cur_model['level_stdv']
      c_mean_stdv_sum = c_mean_stdv.sum()
      model_coef_matrix = np.array(( (cur_model['level_stdv'].sum(), c_mean_stdv_sum), \
                                     (c_mean_stdv_sum, (c_mean_stdv*cur_model['level_mean']).sum()) \
                                  ))
      c_event_stdv = sp_param[event_key]['mean'] * cur_model['level_stdv']
      c_event_stdv_mean = c_event_stdv * cur_model['level_mean']
      dependent_array = np.array((c_event_stdv.sum(), c_event_stdv_mean.sum()));

      sp_param['shift_scale'] = {}
      sp_param['shift_scale']['cal_shift'], sp_param['shift_scale']['cal_scale'] = np.linalg.solve(model_coef_matrix, dependent_array)
      sp_param['shift_scale']['chn_shift'], sp_param['shift_scale']['chn_scale'] = -sp_param["channel_info"]['offset'], sp_param["channel_info"]['digitisation']/sp_param["channel_info"]['range']

      sp_param['shift_scale']['shift']=sp_param['shift_scale']['chn_shift']+sp_param['shift_scale']['chn_scale']*sp_param['shift_scale']['cal_shift']
      sp_param['shift_scale']['scale']=sp_param['shift_scale']['chn_scale']*sp_param['shift_scale']['cal_scale']

      sp_param['raw_signals'] = np.round(sp_param['raw_signals']/sp_param['shift_scale']['cal_scale'] - sp_param['shift_scale']['cal_shift']/sp_param['shift_scale']['cal_scale'], 6)
   except:
      raiseError('Cannot nanopore correction', sp_param, "Cannot nanopore correction")

def getEvent(moptions, sp_param):
  if not sp_param['f5status']=="": return;

  try:
     if sp_param['basecaller'] == 'Guppy':
        event_path = ''.join([fast5_analysis, '/', moptions['basecall_1d'], '/', moptions['basecall_2strand'], '/', fast5_move])
     else:
        event_path = ''.join([fast5_analysis, '/', moptions['basecall_1d'], '/', moptions['basecall_2strand'], '/', fast5_events])
     events_data = sp_param['f5reader'][event_path][()]
  except:
     raiseError('No events/move data', sp_param, "No events/move data")
     return;

  convertError = False;

  if sp_param['f5status'] == "":
     sp_param['events_data'] = events_data
     if sp_param['basecaller'] == 'Guppy': #use move table
        fq_path = ''.join([fast5_analysis,'/',moptions['basecall_1d'],'/',moptions['basecall_2strand'],'/',fast5_basecall_fq])
        fq_data = sp_param['f5reader'][fq_path][()]
        fq_data = (fq_data.decode(encoding="utf-8")).split('\n')
        sp_param['fq_seq'] = fq_data[1]

        sp_param['m_event'] = getMove_Info(moptions, sp_param, events_data)
        sp_param['m_event_basecall'] = sp_param['fq_seq']

     elif sp_param['used_albacore_version']==1:
        move0_left = 0; move0_right = len(events_data)-1;
        while move0_left<move0_right:
           if events_data['move'][move0_left]==0: move0_left += 1;
           #elif not (events_data['move'][move0_left] + events_data['move'][move0_left+1] > 0 and \
           #     events_data['stdv'][move0_left]<2 and events_data['stdv'][move0_left+1]<2): move0_left += 1;
           else: break;
        if move0_left>move0_right-20:
           #if moptions['outLevel']<=myCom.OUTPUT_INFO:
           if moptions['outLevel']<=myCom.OUTPUT_DEBUG:
              print 'left:', move0_left, move0_right, len(events_data)
              for mli in range(move0_left+2):
                 print ('%.1f' % events_data['stdv'][mli]),
              print ''
           raiseError(("Too many move0 at 3'(l%d, r%d)" % (move0_left, move0_right)), sp_param, "Remove too many bases on left")
           return;
        while move0_right>move0_left:
           if events_data['move'][move0_right]==0: move0_right -= 1
           #elif not (events_data['move'][move0_right] + events_data['move'][move0_right-1] > 0 and \
           #     events_data['stdv'][move0_right]<2 and events_data['stdv'][move0_right-1]<2): move0_right -= 1
           else: break;
        if move0_right<move0_left+20:
           #if moptions['outLevel']<=myCom.OUTPUT_INFO:
           if moptions['outLevel']<=myCom.OUTPUT_DEBUG:
              print 'right:', move0_left, move0_right, len(events_data)
              for mli in range(move0_right-2, len(events_data['stdv'])):
                  print ('%.1f' % events_data['stdv'][mli]),
              print ''
           raiseError(("Too many move0 at 5'(l%d, r%d)" % (move0_left, move0_right)), sp_param, 'Remove too many bases on right')
           return

        first_base_index_in_raw_signal = np.round(events_data['start'][move0_left].astype(np.float64)*sp_param["channel_info"]["sampling_rate"]).astype(np.int64) - sp_param['raw_attributes']['start_time']

        if first_base_index_in_raw_signal<-2:
           raiseError(('The index of the first base is less than -2(%d=%.6f*%d-%d)' % (first_base_index_in_raw_signal, events_data['start'][move0_left].astype(np.float64), sp_param["channel_info"]["sampling_rate"], sp_param['raw_attributes']['start_time'])), sp_param, "The index of the first base is less than -2")
           return;
        elif first_base_index_in_raw_signal<0:
           first_base_index_in_raw_signal = 0
           if moptions['outLevel']<=myCom.OUTPUT_INFO: print 'Warning!!! first_base_index_in_raw_signal less than 0', sp_param['mfile_path']

        first_base_index_in_raw_signal = np.uint64(first_base_index_in_raw_signal)

        if moptions['outLevel']<=myCom.OUTPUT_DEBUG: print 'first_base_index_in_raw_signal', move0_left, move0_right, events_data['start'][move0_left], events_data['start'][move0_left].astype(np.float64)*sp_param["channel_info"]["sampling_rate"], np.round(events_data['start'][move0_left].astype(np.float64)*sp_param["channel_info"]["sampling_rate"]).astype(np.uint64), sp_param['raw_attributes']['start_time'], first_base_index_in_raw_signal

        m_event = []; pre_i = move0_left;
        cur_length=(events_data['length'][pre_i]*sp_param["channel_info"]["sampling_rate"]).astype('uint64');
        if moptions['outLevel']<=myCom.OUTPUT_DEBUG:
           for i in range(move0_left+1):
              print 'x', i, events_data['move'][i], ('mean=%.3f stdv=%.3f start=%.6f length=%.6f(%d) %s' % (events_data['mean'][i], events_data['stdv'][i], events_data['start'][i], events_data['length'][i], (events_data['length'][i]*sp_param["channel_info"]["sampling_rate"]).astype('uint64'), events_data['model_state'][i]))
        for i in range(move0_left+1, move0_right+1):
           if moptions['outLevel']<=myCom.OUTPUT_DEBUG and (i-move0_left<30 or move0_right-i<10):
              print 'y', i, events_data['move'][i], ('mean=%.3f stdv=%.3f start=%.6f length=%.6f(%d) %s' % (events_data['mean'][i], events_data['stdv'][i], events_data['start'][i], events_data['length'][i], (events_data['length'][i]*sp_param["channel_info"]["sampling_rate"]).astype('uint64'), events_data['model_state'][i]))

           if events_data['move'][i]>0:
              if pre_i==move0_left:
                 m_event.append((round(events_data['mean'][pre_i],3), round(events_data['stdv'][pre_i],3), first_base_index_in_raw_signal, cur_length, events_data['model_state'][pre_i]))
              else:
                 #try:
                 m_event.append((round(events_data['mean'][pre_i],3), round(events_data['stdv'][pre_i],3), m_event[-1][2]+ m_event[-1][3], cur_length, events_data['model_state'][pre_i]))
                 if m_event[-1][2]>np.iinfo(np.int64).max-2 or m_event[-1][2]<0:
                    if not convertError:
                       print ('ex: %.7f*%d=%.0f' % (events_data['start'][move0_left].astype(np.float64), sp_param["channel_info"]["sampling_rate"], events_data['start'][move0_left].astype(np.float64)*sp_param["channel_info"]["sampling_rate"])), sp_param['raw_attributes']['start_time'], sp_param['mfile_path'], m_event[-1][2], m_event[-1][3]
                    convertError = True;
                 #except Exception:
                 #   if not convertError:
                 #      print ('ex: %.7f*%d=%.0f' % (events_data['start'][move0_left].astype(np.float64), sp_param["channel_info"]["sampling_rate"], events_data['start'][move0_left].astype(np.float64)*sp_param["channel_info"]["sampling_rate"])), sp_param['raw_attributes']['start_time'], sp_param['mfile_path'], m_event[-1][2], m_event[-1][3]
                 #   convertError = True;
              pre_i = i;
              cur_length=(events_data['length'][i]*sp_param["channel_info"]["sampling_rate"]).astype('uint64');
           else:
              cur_length += (events_data['length'][i]*sp_param["channel_info"]["sampling_rate"]).astype('uint64')
        if sp_param['f5status'] == "":
           m_event.append((round(events_data['mean'][pre_i],3), round(events_data['stdv'][pre_i],3), m_event[-1][2]+ m_event[-1][3], cur_length, events_data['model_state'][pre_i]))

        m_event = np.array(m_event, dtype=[('mean', '<f4'), ('stdv', '<f4'), ('start', np.uint64), ('length', np.uint64), ('model_state', 'S5')])
        sp_param['m_event'] = m_event
        sp_param['m_event_basecall'] = ''.join([event_model_state[2] for event_model_state in m_event['model_state']]);

        if moptions['outLevel']<=myCom.OUTPUT_DEBUG:
           for i in range(len(m_event)):
              print ('%d mean=%8.3f stdv=%8.3f' % (move0_left+i, m_event['mean'][i], m_event['stdv'][i])),
              cmean = np.mean(sp_param['raw_signals'][m_event['start'][i]:(m_event['start'][i]+m_event['length'][i])])
              cstdv = np.std(sp_param['raw_signals'][m_event['start'][i]:(m_event['start'][i]+m_event['length'][i])])
              print ('cmean=%8.3f cstd=%8.3f for start=%10d length=%10d' % (cmean, cstdv, m_event['start'][i], m_event['length'][i]))

     elif sp_param['used_albacore_version']==2:
        m_event = [];
        pre_i = 0; pre_length = events_data['length'][pre_i].astype('uint64');
        for cur_i in range(1, len(events_data)):
           if events_data['move'][cur_i]>0:
              m_event.append( (round(events_data['mean'][pre_i],3), round(events_data['stdv'][pre_i],3), events_data['start'][pre_i], pre_length, events_data['model_state'][pre_i]) )

              pre_i = cur_i; pre_length = events_data['length'][pre_i].astype('uint64');
           else:
              pre_length += events_data['length'][cur_i].astype('uint64');
        m_event.append( (round(events_data['mean'][pre_i],3), round(events_data['stdv'][pre_i],3), events_data['start'][pre_i], pre_length, events_data['model_state'][pre_i]) )

        m_event = np.array(m_event, dtype=[('mean', '<f4'), ('stdv', '<f4'), ('start', np.uint64), ('length', np.uint64), ('model_state', 'S5')])
        sp_param['m_event'] = m_event
        sp_param['m_event_basecall'] = ''.join([event_model_state[2] for event_model_state in m_event['model_state']]);
        #raise RuntimeError, ("Albacore 2.x is not supported yet will be implemented later")
     else:
        raise RuntimeError, ("This version of Albacore is not supported. Please use the version of Albacore 1.x or 2.x")

def mnormalized(moptions, sp_param):
   get_cur_shift_scale(moptions, sp_param)

   if not sp_param['m_event']['start'][0] < (sp_param['m_event']['start'][-1]+sp_param['m_event']['length'][-1]):
      print 'Fatal error signal start position is less than the end position', sp_param['mfile_path'], sp_param['m_event']['start'][0], sp_param['m_event']['start'][-1], sp_param['m_event']['length'][-1]

   mshift = np.median(sp_param['raw_signals'][sp_param['m_event']['start'][0]:(sp_param['m_event']['start'][-1]+sp_param['m_event']['length'][-1])])
   mscale = np.median(np.abs(sp_param['raw_signals'][sp_param['m_event']['start'][0]:(sp_param['m_event']['start'][-1]+sp_param['m_event']['length'][-1])]-mshift));
   sp_param['raw_signals'] = (sp_param['raw_signals'] - mshift)/mscale
   read_med = np.median(sp_param['raw_signals'][sp_param['m_event']['start'][0]:(sp_param['m_event']['start'][-1]+sp_param['m_event']['length'][-1])])
   read_mad = np.median(np.abs(sp_param['raw_signals'][sp_param['m_event']['start'][0]:(sp_param['m_event']['start'][-1]+sp_param['m_event']['length'][-1])] - read_med))
   print sp_param['mfile_path'], mshift, mscale, read_med, read_mad; sys.stdout.flush()
   lower_lim = read_med - (read_mad * 5)
   upper_lim = read_med + (read_mad * 5)
   sp_param['raw_signals'] = np.round(np.array([upper_lim if sp_param['raw_signals'][i]>upper_lim else (lower_lim if sp_param['raw_signals'][i]<lower_lim  else sp_param['raw_signals'][i]) for i in range(np.size(sp_param['raw_signals']))]), 3)

def getMove_Info(moptions, sp_param, move_data):
   '''
    sp_param.keys: fq_seq, raw_signals, first_sample_template, duration_template
   '''

    #sp_param['first_sample_template'] = sp_param['f5reader']['/Analyses/Segmentation_001/Summary/segmentation'].attrs['first_sample_template']
    #sp_param['duration_template'] = sp_param['f5reader']['/Analyses/Segmentation_001/Summary/segmentation'].attrs['duration_template']

   seg = "Segmentation_" + moptions['basecall_1d'].split('_')[-1]


   attr_path = '/'.join(['', 'Analyses', seg, 'Summary', 'segmentation'])
    #mv_str = '/'.join(['', 'Analyses', moptions['basecall_1d'], moptions['basecall_2strand'], 'Move'])
   sp_param['first_sample_template'] = sp_param['f5reader'][attr_path].attrs['first_sample_template']
   sp_param['duration_template'] = sp_param['f5reader'][attr_path].attrs['duration_template']
    #move_data = sp_param['f5reader'][mv_str][()]
   nrow = len(sp_param['fq_seq']) # row number of event_info; equals to the base number
   nsig = len(sp_param['raw_signals'])
   first = int(sp_param['first_sample_template'])
   duration = int(sp_param['duration_template'])
   move_info = np.empty([nrow], dtype=[('mean', '<f4'), ('stdv', '<f4'), ('start', np.uint64), ('length', np.uint64), ('model_state', 'U5')])
   effect_sig_index = list(range(first, nsig))
   pivot = first
   seg_count = 0 #which segmentation
   for i in range(1, len(move_data)):
      if move_data[i] == 1:
         move_info[seg_count]['mean'] = np.mean(sp_param['raw_signals'][pivot:(2*i + first)])
         move_info[seg_count]['length'] = 2*i + first - pivot
         move_info[seg_count]['stdv'] = np.std(sp_param['raw_signals'][pivot:(2*i + first)])
         move_info[seg_count]['start'] = pivot
         if seg_count == 0:
            move_info[seg_count]['model_state'] = 'N'*2 + sp_param['fq_seq'][seg_count:seg_count+3]
         elif seg_count == 1:
            move_info[seg_count]['model_state'] = 'N' + sp_param['fq_seq'][seg_count-1:seg_count+3]
         elif seg_count == nrow-2:
            move_info[seg_count]['model_state'] = sp_param['fq_seq'][seg_count-2:seg_count+2] + 'N'
         else:
            move_info[seg_count]['model_state'] = sp_param['fq_seq'][seg_count-2 : seg_count+3]
         pivot = 2*i + first
         seg_count += 1
   move_info[seg_count]['mean'] = np.mean(sp_param['raw_signals'][pivot:nsig])
   move_info[seg_count]['length'] = nsig - pivot
   move_info[seg_count]['stdv'] = np.std(sp_param['raw_signals'][pivot:nsig])
   move_info[seg_count]['start'] = pivot
   move_info[seg_count]['model_state'] = sp_param['fq_seq'][seg_count-2:seg_count+1] + 'N'*2

   #############################################################


   return move_info

def getRawInfo(moptions, sp_param):
   if not sp_param['f5status']=="": return;

   try:
      raw_data = sp_param['f5reader'][fast5_rawReads].values()[0]
      sp_param['raw_attributes'] = dict(raw_data.attrs.items())
      #print 'Test', raw_data['Signal'].value[:5]

      sp_param['raw_signals'] = raw_data['Signal'].value
      #sp_param['raw_signals'] = np.round((raw_data['Signal'].value + sp_param["channel_info"]['offset']) * sp_param["channel_info"]['range']/ sp_param["channel_info"]['digitisation'], 6)
      #print 'Test', sp_param['raw_signals'][:5], sp_param["channel_info"]['offset'], sp_param["channel_info"]['range'], sp_param["channel_info"]['digitisation'], sp_param["channel_info"]['range']/ sp_param["channel_info"]['digitisation']

      '''
      mshift = np.median(sp_param['raw_signals'])
      mscale = np.median(np.abs(sp_param['raw_signals']-mshift));
      sp_param['raw_signals'] = (sp_param['raw_signals'] - mshift)/mscale
      read_med = np.median(sp_param['raw_signals'])
      read_mad = np.median(np.abs(sp_param['raw_signals'] - read_med))
      print sp_param['mfile_path'], mshift, mscale, read_med, read_mad; sys.stdout.flush()
      lower_lim = read_med - (read_mad * 5)
      upper_lim = read_med + (read_mad * 5)
      sp_param['raw_signals'] = np.round(np.array([upper_lim if sp_param['raw_signals'][i]>upper_lim else (lower_lim if sp_param['raw_signals'][i]<lower_lim  else sp_param['raw_signals'][i]) for i in range(np.size(sp_param['raw_signals']))]), 3)
      '''

      #sp_param["channel_info"] = {'digitisation':channel_info['digitisation'], 'offset':channel_info['offset'], 'range':channel_info['range'], 'sampling_rate':channel_info['sampling_rate'], 'channel_number':channel_info['channel_number']}
   except:
      raiseError(("No Raw_reads/Signal data %s for %s" % (fast5_rawReads, sp_param['mfile_path'])), sp_param, "No Raw_reads/Signal")

def getFast5Info(moptions, sp_param):
   get_channel_info(moptions, sp_param)
   if not sp_param.has_key("channel_info"):
      raiseError(("Channel information could not be found in %s " % fast5_channel_id), sp_param, "Channel information could not be found")
      return;
   getAlbacoreVersion(moptions, sp_param)
   if not sp_param.has_key('used_albacore_version'):
      return

   try:
      fq_path = ''.join([fast5_analysis, '/', moptions['basecall_1d'], '/', moptions['basecall_2strand'], '/', fast5_basecall_fq])
      fq_data =sp_param['f5reader'][fq_path][()].split('\n')
      sp_param['read_id'] = (fq_data[0][1:] if fq_data[0][0]=='@' else fq_data[0]).replace(" ", ":::").replace("\t", "|||")
   except:
      raiseError('No Fastq data', sp_param, "No Fastq data")

   getRawInfo(moptions, sp_param)
   if sp_param['f5status']=="": getEvent(moptions, sp_param)
   if sp_param['f5status']=="": mnormalized(moptions, sp_param)


def getEventAnnotation(moptions, f5files):
   if moptions['outLevel']<=myCom.OUTPUT_DEBUG:
      start_time = time.time(); runnum = 0;

   f5data = {}
   moptions["Error"] = defaultdict(list)
   for f5f in f5files:
      with h5py.File(f5f, 'r') as mf5:
         sp_param = {}
         sp_param['mfile_path'] = f5f
         sp_param['f5reader'] = mf5
         sp_param['f5status'] = "";
         getFast5Info(moptions, sp_param)

         if sp_param['f5status'] == "":
            if sp_param['read_id'] in f5data:
               print 'Duplicate id', sp_param['read_id'], f5f
            f5data[sp_param['read_id']] = (sp_param['m_event_basecall'], sp_param['m_event'], sp_param['raw_signals'], f5f)
         else:
            moptions["Error"][sp_param['f5status']].append(f5f)

         if moptions['outLevel']<=myCom.OUTPUT_DEBUG:
            runnum += 1;
            if runnum%500==0:
               end_time = time.time();
               print ("%d consuming time %d" % (runnum, end_time-start_time))
   if moptions['outLevel']<=myCom.OUTPUT_DEBUG:
      end_time = time.time();
      print ("All consuming time %d" % (end_time-start_time))

   return f5data;

def correctAndAnnotate(moptions, f5files):
   f5data = getEventAnnotation(moptions, f5files)

   if moptions['outLevel']<=myCom.OUTPUT_DEBUG: start_time = time.time();
   temp_fa = tempfile.NamedTemporaryFile(suffix='.fa')
   f5keys = f5data.keys(); f5keys.sort()
   for f5k in f5keys:
      temp_fa.write(''.join(['>', f5k, '\n', f5data[f5k][0], '\n']))
   temp_fa.flush();
   if moptions['outLevel']<=myCom.OUTPUT_DEBUG:
      end_time = time.time();
      print ("Write consuming time %d" % (end_time-start_time))

   temp_sam = tempfile.NamedTemporaryFile()
   if moptions['alignStr']=='bwa':
      cmd_opt = ['mem', '-x', 'ont2d', '-v', '1', '-t', '1', moptions['Ref'], temp_fa.name]
   else:
      cmd_opt = ['-ax', 'map-ont', moptions['Ref'], temp_fa.name]
   returncode = subprocess.call([moptions['alignStr'],]+cmd_opt, stdout=temp_sam)
   if not returncode==0:
      print ('Fatal Error!!! returncode is non-zero(%d) for "%s"' % (returncode, curcmd))
      errkey = "Cannot running aligment"
      for f5k in f5keys:
         moptions["Error"][errkey].append(f5data[f5k][3])
      return;

   temp_fa.close();
   temp_sam.seek(0);
   align_info = temp_sam.readlines()
   align_info = map(string.strip, align_info);
   temp_sam.close();

   sp_param = {};
   sp_param['f5data'] = f5data

   f5align = defaultdict()
   f5keydict = defaultdict();
   sp_param['ref_info'] = defaultdict()

   if moptions['outLevel']<=myCom.OUTPUT_DEBUG:start_time = time.time();
   ilid = 0;
   while ilid < len(align_info):
      if len(align_info[ilid])==0 or align_info[ilid][0]=='@':
         ilid += 1
         continue;

      sp_param['f5status'] = "";
      sp_param['line'] = align_info[ilid]
      qname = handle_line(moptions, sp_param, f5align)
      if sp_param['f5status'] == "":
         f5keydict[qname] = True;
      ilid += 1

   for f5k in f5keys:
      if f5k not in f5keydict:
         moptions["Error"]["Not in alignment sam"].append(f5data[f5k][3])

   if moptions['outLevel']<=myCom.OUTPUT_DEBUG:
      end_time = time.time();
      print ("Get BAM consuming time %d" % (end_time-start_time))

   sp_param['f5status']= ""
   sp_param['line'] = ""
   if moptions['outLevel']<=myCom.OUTPUT_DEBUG:start_time = time.time();
   handle_record(moptions, sp_param, f5align, f5data)
   if moptions['outLevel']<=myCom.OUTPUT_DEBUG:
      end_time = time.time();
      print ("Analyze & annotate & save consuming time %d" % (end_time-start_time))

def getRefSeq(moptions, sp_param, rname):
   temp_seq = tempfile.NamedTemporaryFile()
   cmd_opt = ['faidx', moptions['Ref'], rname]
   returncode = subprocess.call(['samtools',]+cmd_opt, stdout=temp_seq)
   if not returncode==0:
      print ('Fatal Error!!! cannot find the chrosome sequence %s' % rname)
   else:
      temp_seq.seek(0);
      seqinfo = map(string.strip, temp_seq.readlines());
      temp_seq.close();

      sp_param['ref_info'][rname] = string.strip(''.join(seqinfo[1:])).upper()

def handle_record(moptions, sp_param, f5align, f5data):
   alignkeys = f5align.keys();
   numreg = re.compile('\d+')
   mdireg = re.compile('[MIDNSHPX=]{1}')

   signalneighbor = defaultdict(int);

   if moptions['outLevel']<=myCom.OUTPUT_DEBUG:
      start_time = time.time();
      handlenum = 0;
   for readk in alignkeys:
     if moptions['outLevel']<=myCom.OUTPUT_DEBUG:
        handlenum += 1;
        if handlenum%500==0:
           end_time = time.time();
           print ("%d consuming time %d" % (handlenum, end_time-start_time))

     _, flag, rname, pos, cigar, readseq = f5align[readk]
     if rname not in sp_param['ref_info']:
        getRefSeq(moptions, sp_param, rname)
     refseq = sp_param['ref_info'][rname]

     pos = pos - 1
     numinfo = numreg.findall(cigar);      mdiinfo = mdireg.findall(cigar)

     forward_reverse = '-' if flag&0x10 else '+'
     numinfo = map(int, numinfo)

     if moptions['outLevel']<=myCom.OUTPUT_DEBUG:
        firstn = 7;
        print zip(numinfo[:firstn], mdiinfo[:firstn]); print 'b ref:', refseq[pos:pos+50]; print 'bread:', readseq[:50]; print 'bcall:', f5data[readk][0][:50];
        print zip(numinfo[-firstn:], mdiinfo[-firstn:]);

     leftclip = 0; rightclip = 0;
     while mdiinfo[0] in ['I', 'D', 'N', 'S', 'H', 'P', 'X']:
         if mdiinfo[0] in ['I', 'S', 'X']:
            leftclip += numinfo[0];  readseq = readseq[numinfo[0]:]
         if mdiinfo[0] in ['H']: leftclip += numinfo[0]
         if mdiinfo[0] in ['D', 'N', 'X']:
            pos += numinfo[0]
         numinfo = numinfo[1:];  mdiinfo = mdiinfo[1:]
     while mdiinfo[-1] in ['I', 'D', 'N', 'S', 'H', 'P', 'X']:
         if mdiinfo[-1] in ['I', 'S', 'X']:
            rightclip += numinfo[-1]; readseq = readseq[:-numinfo[-1]]
         if mdiinfo[-1] in ['H']: rightclip += numinfo[-1]
         numinfo = numinfo[:-1]; mdiinfo = mdiinfo[:-1]
     if forward_reverse=='+':
        if rightclip>0: m_event = f5data[readk][1][leftclip:-rightclip]
        else: m_event = f5data[readk][1][leftclip:]
     else:
        if leftclip>0: m_event = f5data[readk][1][rightclip:-leftclip]
        else: m_event = f5data[readk][1][rightclip:]

     if moptions['outLevel']<=myCom.OUTPUT_DEBUG:
        print 'clip', leftclip, rightclip, pos #, zip(numinfo[:firstn], mdiinfo[:firstn])
        print ' ref:', refseq[pos:pos+50];  print 'read:', readseq[:50]

     lastmatch = None; firstmatch = None; first_match_pos = None;
     last_al_match = None; lasmtind = 0;
     base_map_info = []; #indel_groups = defaultdict()
     nummismatch = 0; numinsert = 0; numdel = 0;
     read_ind = 0;
     for n1ind in range(len(numinfo)):
        mdi = mdiinfo[n1ind];
        for n1i in range(numinfo[n1ind]):
           if mdi=='M':
              base_map_info.append((refseq[pos], readseq[read_ind]))
              if refseq[pos]==readseq[read_ind]:
                 if firstmatch==None: firstmatch = read_ind
                 if lastmatch==None or lastmatch<read_ind: lastmatch = read_ind; lasmtind=n1ind
                 if last_al_match==None or last_al_match<len(base_map_info): last_al_match=len(base_map_info)-1
                 if first_match_pos==None: first_match_pos  = pos
              else: nummismatch += 1
              pos += 1; read_ind += 1;
           elif mdi =='I':
              base_map_info.append(('-', readseq[read_ind]))
              #indel_groups[len(base_map_info)-1] = 'I'
              read_ind += 1;
              numinsert += 1
           elif mdi == 'D':
              base_map_info.append((refseq[pos], '-'))
              #indel_groups[len(base_map_info)-1] = 'D'
              pos += 1;
              numdel += 1
           elif mdi == 'N':
              base_map_info.append((refseq[pos], '-'))
              pos += 1;
              if moptions['outLevel']<=myCom.OUTPUT_WARNING:
                 print 'CIGAR-Error N exist', f5data[readk][3]
           elif mdi == 'S':
              #################### !!!!!!!!!!missing
              read_ind += 1;
              if moptions['outLevel']<=myCom.OUTPUT_WARNING:
                 print 'CIGAR-Error!!! S in the middle of the sequence', f5data[readk][3]
           elif mdi == 'H':
              if moptions['outLevel']<=myCom.OUTPUT_WARNING:
                 print 'CIGAR-Error!!! H in the middle of the sequence', f5data[readk][3]
           elif mdi == 'P':
              if moptions['outLevel']<=myCom.OUTPUT_WARNING:
                 print 'CIGAR-Error!!! P exist', f5data[readk][3]
           elif mdi == '=':
             base_map_info.append((refseq[pos], readseq[read_ind]))
             if first_match_pos==None: first_match_pos  = pos
             pos += 1; read_ind += 1;
             if firstmatch==None: firstmatch = read_ind - 1
             if lastmatch==None or lastmatch<read_ind-1: lastmatch = read_ind - 1; lasmtind=n1ind
             if last_al_match==None or last_al_match<len(base_map_info): last_al_match=len(base_map_info)-1
           elif mdi == 'X':
             base_map_info.append((refseq[pos], readseq[read_ind]))
             pos += 1; read_ind += 1;
             nummismatch += 1
           else:
             if moptions['outLevel']<=myCom.OUTPUT_WARNING:
                print 'CIGAR-Error!!!', 'Warning unknow CIGAR element ' + str(numinfo[n1ind]) + ' ' + mdi, f5data[readk][3]
     #print ('Statistics: Mism=%.3f, ins=%.3f, del=%.3f' % (nummismatch/float(len(base_map_info))*100, numinsert/float(len(base_map_info))*100, numdel/float(len(base_map_info))*100))
     if firstmatch==None or lastmatch==None or firstmatch<0 or lastmatch<0:
        if moptions['outLevel']<=myCom.OUTPUT_WARNING:
           print 'match-Error!!! no first and/or last match', f5data[readk][3],
           if firstmatch==None: print('firstmatch=None'),
           elif firstmatch<0: print ('firstmatch=%d' % firstmatch),
           if lastmatch==None: print ('lastmatch=None'),
           elif lastmatch<0: print ('lastmatch=%d' % lastmatch),
           print ''
        moptions["Error"]['Incorrect Alignment'].append(f5data[readk][3])
        continue

     #if moptions['outLevel']<=myCom.OUTPUT_WARNING:
     if moptions['outLevel']<=myCom.OUTPUT_DEBUG:
        print 'pre info', leftclip, rightclip, 'new', firstmatch, lastmatch, len(m_event), len(m_event)-lastmatch-1
     if not firstmatch==None: leftclip += firstmatch
     if (not lastmatch==None) and len(m_event)-lastmatch>1: rightclip += len(m_event)-lastmatch-1

     if forward_reverse=='+':
        if len(m_event)-lastmatch>1:
           m_event = m_event[firstmatch:(lastmatch+1-len(m_event))]
        elif firstmatch>0: m_event = m_event[firstmatch:]
     else:
        if firstmatch>0: m_event = m_event[(len(m_event)-1-lastmatch):-firstmatch]
        elif len(m_event)-lastmatch>1: m_event = m_event[(len(m_event)-1-lastmatch):]
     if firstmatch>0 or len(base_map_info)-last_al_match>1:
        if moptions['outLevel']<=myCom.OUTPUT_WARNING:
           print 'Warning!!! first not match', firstmatch, last_al_match, len(base_map_info), numinfo[lasmtind-2:(lasmtind+5)], mdiinfo[lasmtind-2:(lasmtind+5)], lasmtind, len(numinfo)
        reflast= []; readlast = [];
        for i in range(last_al_match-10, (len(base_map_info) if last_al_match+30>len(base_map_info) else last_al_match+30)):
            reflast.append(base_map_info[i][0]); readlast.append(base_map_info[i][1]);
            if i==last_al_match: reflast.append('*'); readlast.append('*')
        if moptions['outLevel']<=myCom.OUTPUT_DEBUG: #myCom.OUTPUT_INFO:
           print '\t ref', ''.join(reflast)
           print '\tread', ''.join(readlast)

        if len(base_map_info)-last_al_match>1: base_map_info = base_map_info[firstmatch:(last_al_match+1-len(base_map_info))]
        elif firstmatch>0: base_map_info = base_map_info[firstmatch:]

     base_map_info = np.array(base_map_info, dtype=[('refbase', 'S1'), ('readbase', 'S1')])
     #if moptions['outLevel']<=myCom.OUTPUT_INFO: # myCom.OUTPUT_DEBUG: #myCom.OUTPUT_INFO:
     if moptions['outLevel']<=myCom.OUTPUT_DEBUG:
        print '  ref:', ''.join(base_map_info['refbase'][:(50 if len(base_map_info)>50 else len(base_map_info))]), leftclip, rightclip
        print forward_reverse+'read:', ''.join(base_map_info['readbase'][:(50 if len(base_map_info)>50 else len(base_map_info))]), len(base_map_info), len(m_event)
        if forward_reverse=='-':
           if leftclip==0: origianl_bc = f5data[readk][0][(-leftclip-(50 if len(f5data[readk][0])>50+leftclip else len(f5data[readk][0])-leftclip)):][::-1]
           else: origianl_bc = f5data[readk][0][(-leftclip-(50 if len(f5data[readk][0])>50+leftclip else len(f5data[readk][0])-leftclip)):-leftclip][::-1]
        else: origianl_bc = f5data[readk][0][leftclip:(leftclip+(50 if len(f5data[readk][0])>50+leftclip else len(f5data[readk][0])-leftclip))]
        redblist = []; rind = 0
        for i in range(50 if len(base_map_info)>50 else len(base_map_info)):
           if base_map_info['readbase'][i]=='-': redblist.append('-');
           else:
              redblist.append(origianl_bc[rind]);
              rind += 1
        print forward_reverse+'redb:', ''.join(redblist), len(m_event)
        if forward_reverse=='-': origianl_me = m_event[-(50 if len(m_event)>50 else len(m_event)):][::-1]
        else: origianl_me = m_event[:(50 if len(m_event)>50 else len(m_event))]
        origianl_bc = ''.join([event_model_state[2] for event_model_state in origianl_me['model_state']])
        redblist = []; rind = 0
        for i in range(50 if len(base_map_info)>50 else len(base_map_info)):
           if base_map_info['readbase'][i]=='-': redblist.append('-');
           else:
              redblist.append(origianl_bc[rind]);
              rind += 1
        print forward_reverse+'mevt:', ''.join(redblist); #''.join([event_model_state[2] for event_model_state in origianl_me['model_state']]);
        print 'clip', (leftclip, rightclip, nummismatch, numinsert, numdel), len(base_map_info), '<3 5>'
        print '  ref:', ''.join(base_map_info['refbase'][-(50 if len(base_map_info)>50 else len(base_map_info)):]), leftclip, rightclip
        print forward_reverse+'read:', ''.join(base_map_info['readbase'][-(50 if len(base_map_info)>50 else len(base_map_info)):]), len(base_map_info), len(m_event)
     indel_pos = fix_repeat_del(base_map_info, moptions, sp_param, f5data[readk][3])
     if moptions['outLevel']<=myCom.OUTPUT_DEBUG: #myCom.OUTPUT_INFO:
        print len(base_map_info), len(m_event), len(indel_pos), len(m_event)+len(indel_pos)
        print ''

     print 'handling', forward_reverse, f5data[readk][3]
     group_indel_pos = group_indel(indel_pos, m_event, base_map_info, forward_reverse, moptions, sp_param)

     annotated_info, signalnum = annotate1(group_indel_pos, m_event, base_map_info, forward_reverse, f5data[readk][2], f5data[readk][3], moptions, sp_param); #sys.exit(1)
     snkeys = signalnum.keys();
     for snk in snkeys:
        signalneighbor[snk] += signalnum[snk]

     #myCom.fast5_mo_correction
     save_annotation(annotated_info, m_event, base_map_info, f5data[readk][3], (leftclip, rightclip, nummismatch, numinsert, numdel, len(base_map_info)-nummismatch-numinsert-numdel), (forward_reverse, rname, first_match_pos), moptions, sp_param);
   #if moptions['outLevel']<=myCom.OUTPUT_INFO:
   if moptions['outLevel']<=myCom.OUTPUT_DEBUG:
      print signalneighbor

   if moptions['outLevel']<=myCom.OUTPUT_DEBUG:
      end_time = time.time();
      print ("%d consuming time %d" % (handlenum, end_time-start_time))
   moptions['signalneighbor'] = signalneighbor

#
#
#
def get_complement(na):
   #if na in ['-', '*', '+']: return '-'
   #el
   if na in myCom.acgt: return myCom.na_bp[na]
   else: return na;
def save_annotation(annotated_info, m_event, base_map_info, f5f, match_info, map_info, moptions, sp_param):
   m_event_data = [];
   annotated_info_keys = annotated_info.keys(); annotated_info_keys.sort()
   #m_ref_nogap = ''.join(base_map_info['refbase'] if map_info[0]=='+' else map(get_complement, base_map_info['refbase'])[::-1])
   for aim in annotated_info_keys:
      m_event_data.append( ( annotated_info[aim][1], annotated_info[aim][2], \
      annotated_info[aim][3], \
# m_event['start'][annotated_info[aim][0]] if map_info[0]=='+' else m_event['start'][-annotated_info[aim][0]-1], \
      annotated_info[aim][4], \
# m_event['length'][annotated_info[aim][0]] if map_info[0]=='+' else m_event['length'][-annotated_info[aim][0]-1], \
       base_map_info['refbase'][aim] if map_info[0]=='+' else get_complement(base_map_info['refbase'][aim]) \
       #(m_event['model_state'][annotated_info[aim][0]])[2] if map_info[0]=='+' else (m_event['model_state'][-annotated_info[aim][0]-1])[2] \
                          ) )
   m_event_data = np.array((m_event_data if map_info[0]=='+' else m_event_data[::-1]), dtype=[('norm_mean', '<f8'), ('norm_stdev', '<f8'), ('start', '<u4'), ('length', '<u4'), ('base', 'S1')])

   with h5py.File(f5f, 'r+') as save_data:
      base_group = save_data[fast5_analysis]
      if myCom.fast5_mo_correction in base_group:
         try:
            del base_group[myCom.fast5_mo_correction]
            save_data.flush();
         except:
            save_data.close();
            moptions["Error"]['Cannot del old data'].append(f5f)
            print ('Error!!! %s in %s' % ("Cannot del old data", f5f))
            return
      nanomo_correct = base_group.create_group(myCom.fast5_mo_correction)

      #nanomo_correct_bc_group = nanomo_correct.create_group(myCom.rawAlignment_base)
      nanomo_correct_bc_group = nanomo_correct.create_group(myCom.basecall_template_base)
      #rawBaseCalled_template_base = "BaseCalled_template"
      #rawAlignment_base = "Alignment"
      #rawRead_segments_base = "read_segments"

      nanomo_correct_alignment_group = nanomo_correct_bc_group.create_group(myCom.rawAlignment_base)
      nanomo_correct_alignment_group.attrs[myCom.map_start_str] = map_info[2]
      nanomo_correct_alignment_group.attrs[myCom.map_strand_str] = map_info[0]
      nanomo_correct_alignment_group.attrs[myCom.map_chr_str] = map_info[1]
      if map_info[0]=='+':
         nanomo_correct_alignment_group.attrs['clipped_bases_start'] = match_info[0]
         nanomo_correct_alignment_group.attrs['clipped_bases_end'] = match_info[1]
      else:
         nanomo_correct_alignment_group.attrs['clipped_bases_start'] = match_info[1]
         nanomo_correct_alignment_group.attrs['clipped_bases_end'] = match_info[0]
      nanomo_correct_alignment_group.attrs['num_insertions'] = match_info[3]
      nanomo_correct_alignment_group.attrs['num_deletions'] = match_info[4]
      nanomo_correct_alignment_group.attrs['num_matches'] = match_info[5]
      nanomo_correct_alignment_group.attrs['num_mismatches'] = match_info[2]
      nanomo_correct_alignment_group.attrs['Bcinfo'] = moptions['basecall_1d']

      nanomo_correct_alignment_group.create_dataset(myCom.rawRead_alignment_base, data=(base_map_info['readbase'] if map_info[0]=='+' else map(get_complement, base_map_info['readbase'])[::-1]), compression="gzip")
      nanomo_correct_alignment_group.create_dataset(myCom.rawGenome_alignment_base, data=(base_map_info['refbase'] if map_info[0]=='+' else map(get_complement, base_map_info['refbase'])[::-1]), compression="gzip")
      #nanomo_correct_alignment_group.create_dataset(myCom.raw_event_base, data=m_event_data, compression="gzip");
      nanomo_correct_bc_group.create_dataset(myCom.raw_event_base, data=m_event_data, compression="gzip");

      try:
         save_data.flush();
         save_data.close();
      except:
         moptions["Error"]['Cannot save data'].append(f5f)
         print ('Error!!! %s in %s' % ("Cannot save data", f5f))



#
#
#
def annotate1(group_indel_pos, m_event, base_map_info, forward_reverse, raw_pv, f5f, moptions, sp_param):
   cur_level = myCom.OUTPUT_DEBUG; #myCom.OUTPUT_INFO
   annotate_info = defaultdict();

   #if forward_reverse=='+':
   #   print f5f
   #   for i in range(5):
   #      prlst = []
   #      for j in range(m_event[i][2], m_event[i][2]+m_event[i][3]):
   #          prlst.append(('%.5f' % raw_pv[j]))
   #      print i, m_event[i], ' '.join(prlst)

   #for bmi in range(10):
   #   print bmi, base_map_info['refbase'][bmi], base_map_info['readbase'][bmi], forward_reverse
   #print forward_reverse, f5f
   group_indel_pos_keys = group_indel_pos.keys(); group_indel_pos_keys.sort()
   bmi = 0; event_ind = -1;
   if moptions['outLevel']<=cur_level:
      dif2 = [0, 0]
   for gipk in group_indel_pos_keys:
      #print gipk, group_indel_pos[gipk]
      if moptions['outLevel']<=cur_level and forward_reverse=='-':
      #if moptions['outLevel']<=myCom.OUTPUT_INFO:
         print bmi, gipk, group_indel_pos[gipk][2], 'range', gipk-group_indel_pos[gipk][3][0], group_indel_pos[gipk][2]+group_indel_pos[gipk][3][1]+1, group_indel_pos[gipk], event_ind
      #while bmi < (gipk-moptions['Resegment_wind'] if gipk-moptions['Resegment_wind']>-1 else 0):
      while bmi < (gipk-group_indel_pos[gipk][3][0] if gipk-group_indel_pos[gipk][3][0]>-1 else 0):
         if event_ind>-1 and (not (base_map_info['readbase'][bmi]==(m_event['model_state'][event_ind+1][2] if forward_reverse=='+' else get_complement(m_event['model_state'][-event_ind-2][2])))): print 'b', gipk, group_indel_pos[gipk], bmi, base_map_info['refbase'][bmi], base_map_info['readbase'][bmi], (m_event['model_state'][event_ind+1][2] if forward_reverse=='+' else get_complement(m_event['model_state'][-event_ind-2][2])), event_ind+1, ("********" if not (base_map_info['readbase'][bmi]==(m_event['model_state'][event_ind+1][2] if forward_reverse=='+' else get_complement(m_event['model_state'][-event_ind-2][2]))) else "")  # ((m_event['model_state'][event_ind][2], m_event['model_state'][event_ind+1][2]) if forward_reverse=='+' else (get_complement(m_event['model_state'][-event_ind-1][2]), get_complement(m_event['model_state'][-event_ind-2][2])))
         event_ind += 1
         cmean, cstd, cdif = calculate_mean_std(m_event, event_ind, forward_reverse, raw_pv, moptions, sp_param)
         annotate_info[bmi] = (event_ind, cmean, cstd, m_event['start'][event_ind] if forward_reverse=='+' else m_event['start'][-event_ind-1], m_event['length'][event_ind] if forward_reverse=='+' else m_event['length'][-event_ind-1]);
# m_event['start'][annotated_info[aim][0]] if map_info[0]=='+' else m_event['start'][-annotated_info[aim][0]-1]
# m_event['length'][annotated_info[aim][0]] if map_info[0]=='+' else m_event['length'][-annotated_info[aim][0]-1]
         if moptions['outLevel']<=cur_level:
           if cdif:
              dif2[0] += 1
           else: dif2[1] += 1
         bmi += 1
      #while bmi < group_indel_pos[gipk][2]+moptions['Resegment_wind']+1 and bmi < len(base_map_info):
      while bmi < group_indel_pos[gipk][2]+group_indel_pos[gipk][3][1]+1 and bmi < len(base_map_info):
         #print 'i', gipk, group_indel_pos[gipk], bmi, base_map_info['refbase'][bmi], base_map_info['readbase'][bmi], (m_event['model_state'][event_ind+1][2] if forward_reverse=='+' else get_complement(m_event['model_state'][-event_ind-2][2])), (event_ind+1 if (base_map_info['readbase'][bmi] in myCom.acgt) else event_ind) #, ((m_event['model_state'][event_ind][2], m_event['model_state'][event_ind+1][2]) if forward_reverse=='+' else (get_complement(m_event['model_state'][-event_ind-1][2]), get_complement(m_event['model_state'][-event_ind-2][2])))
         if base_map_info['readbase'][bmi] in myCom.acgt:
            event_ind += 1
         if base_map_info['refbase'][bmi] in myCom.acgt:
            annotate_info[bmi] = (event_ind, False)
         bmi += 1
   while bmi < len(base_map_info):
      if not (base_map_info['readbase'][bmi]==(m_event['model_state'][event_ind+1][2] if forward_reverse=='+' else get_complement(m_event['model_state'][-event_ind-2][2]))): print 'e', bmi, base_map_info['refbase'][bmi], base_map_info['readbase'][bmi], (m_event['model_state'][event_ind+1][2] if forward_reverse=='+' else get_complement(m_event['model_state'][-event_ind-2][2])), event_ind+1, ("********" if not (base_map_info['readbase'][bmi]==(m_event['model_state'][event_ind+1][2] if forward_reverse=='+' else get_complement(m_event['model_state'][-event_ind-2][2]))) else "")  #, ((m_event['model_state'][event_ind][2], m_event['model_state'][event_ind+1][2]) if forward_reverse=='+' else (get_complement(m_event['model_state'][-event_ind-1][2]), get_complement(m_event['model_state'][-event_ind-2][2])))
      event_ind += 1
      cmean,cstd,cdif=calculate_mean_std(m_event,event_ind,forward_reverse,raw_pv,moptions,sp_param,isprint=False)#True)
      annotate_info[bmi] = (event_ind, cmean, cstd, m_event['start'][event_ind] if forward_reverse=='+' else m_event['start'][-event_ind-1], m_event['length'][event_ind] if forward_reverse=='+' else m_event['length'][-event_ind-1]);
      if moptions['outLevel']<=cur_level:
         if cdif:
            dif2[0] += 1
         else: dif2[1] += 1
      bmi += 1
   if moptions['outLevel']<=cur_level: print dif2

   #
   signalnum = defaultdict(int);
   for gipk in group_indel_pos_keys:
      if forward_reverse=='+':
         mstart1, mlen1 = m_event[group_indel_pos[gipk][0]][2], m_event[group_indel_pos[gipk][0]][3]
         mstart2, mlen2 = m_event[group_indel_pos[gipk][1]][2], m_event[group_indel_pos[gipk][1]][3]
      else:
         mstart1, mlen1 = m_event[-group_indel_pos[gipk][1]-1][2], m_event[-group_indel_pos[gipk][1]-1][3]
         mstart2, mlen2 = m_event[-group_indel_pos[gipk][0]-1][2], m_event[-group_indel_pos[gipk][0]-1][3]
      if moptions['outLevel']<=myCom.OUTPUT_ERROR:
         if mstart1>=mstart2: print 'Error left start is larger than right start', mstart1, mstart2
         if mstart1+mlen1>mstart2+mlen2: print 'Error left end is larger than right end', mstart1+mlen1, mstart2+mlen2

      expectna = 0;
      refstr = []; readstr = []
      #for bmi in range(gipk-moptions['Resegment_wind'], group_indel_pos[gipk][2]+moptions['Resegment_wind']+1):
      for bmi in range(gipk-group_indel_pos[gipk][3][0], group_indel_pos[gipk][2]+group_indel_pos[gipk][3][1]+1):
         if bmi < 0: continue;
         if not bmi<len(base_map_info): break;
         refstr.append(base_map_info['refbase'][bmi]); readstr.append(base_map_info['readbase'][bmi])
         if base_map_info['refbase'][bmi]=='-':
            continue;
         #if base_map_info['readbase'][bmi] in ['+', '*']: continue;
         #if gipk<10 and forward_reverse=='+':
         #   print '\t', base_map_info['refbase'][bmi], base_map_info['readbase'][bmi]
         if base_map_info['readbase'][bmi]=='~':
            if (bmi>0 and base_map_info['readbase'][bmi-1]=='~'): continue;
         expectna += 1;
      pvsignals = raw_pv[mstart1:mstart2+mlen2]

      for currsw in range(moptions['Resegment_signal_wind'], 1, -1):
         split_pos = find_sp(pvsignals, gipk, refstr, readstr, expectna, currsw, group_indel_pos, mstart1, mstart2, mlen2, f5f, cur_level, moptions)
         if not split_pos==None: break;
      if currsw<moptions['Resegment_signal_wind'] and (not split_pos==None):
         if moptions['outLevel']<=cur_level: #myCom.OUTPUT_INFO:
            print 'Successful!!!', currsw, moptions['Resegment_signal_wind']

      #if gipk<10 and forward_reverse=='+':
      #   print f5f
      #   print pvsignals
      #   print gipk, expectna, mstart1, mstart2, mlen2, moptions['Resegment_wind'], gipk-moptions['Resegment_wind'], group_indel_pos[gipk][2]+moptions['Resegment_wind']+1
      #   print split_pos

      if not split_pos==None:
         signalnum[currsw] += 1
      if split_pos==None:
         signalnum[1] += 1
         all_mean = np.mean(pvsignals);
         all_std = np.std(pvsignals)
      else:
         if moptions['outLevel']<=cur_level:
            spind = 0;
            for pvprind in range(len(pvsignals)):
               if spind<len(split_pos) and split_pos[spind][0]==pvprind:
                  print ('|%.3f|' % split_pos[spind][1]),
                  spind += 1
               print ('%.3f ' % pvsignals[pvprind]),
            print ''

      #if moptions['outLevel']<=myCom.OUTPUT_INFO:
      if moptions['outLevel']<=cur_level:
         #print ' ref', ''.join(refstr), gipk-moptions['Resegment_wind'], gipk, group_indel_pos[gipk][2]+moptions['Resegment_wind']
         print ' ref', ''.join(refstr), gipk-group_indel_pos[gipk][3][0], gipk, group_indel_pos[gipk][2]+group_indel_pos[gipk][3][1]
         print 'read', ''.join(readstr), mstart1, mstart2+mlen2, len(pvsignals), (len(split_pos) if not split_pos==None else ''), expectna

      #bmi = gipk-moptions['Resegment_wind']
      bmi = gipk-group_indel_pos[gipk][3][0]
      if bmi<0: bmi = 0;
      if forward_reverse=='-' and (not split_pos==None):
         spind = len(split_pos)-1
      else: spind = -1;
      #while bmi < group_indel_pos[gipk][2]+moptions['Resegment_wind']+1:
      while bmi < group_indel_pos[gipk][2]+group_indel_pos[gipk][3][1]+1:
         if not bmi<len(base_map_info): break;
         if base_map_info['refbase'][bmi]=='-':
            bmi += 1
            continue;

         if not split_pos==None:
            if spind==-1: start_in_pv = 0;
            else: start_in_pv = split_pos[spind][0]
            if spind==len(split_pos)-1: cur_pv_signal = pvsignals[start_in_pv:]
            else: cur_pv_signal = pvsignals[start_in_pv:split_pos[spind+1][0]]

         if base_map_info['readbase'][bmi]=='~':
            if bmi>0 and base_map_info['readbase'][bmi-1]=='~':
               annotate_info[bmi] = annotate_info[bmi-1]
            else:
               annotate_info[bmi] = (annotate_info[bmi][0], round(np.mean(cur_pv_signal),3), round(np.std(cur_pv_signal),3), mstart1+start_in_pv, split_pos[spind+1][0]-start_in_pv if spind+1<len(split_pos) else mstart2+mlen2-mstart1-start_in_pv) if not split_pos==None else (annotate_info[bmi][0], all_mean, all_std, mstart1, mstart2+mlen2-mstart1)
            if bmi<len(base_map_info['readbase'])-1 and (not base_map_info['readbase'][bmi+1]=='~'):
               if forward_reverse=='+': spind += 1;
               else: spind -= 1
            bmi += 1
         elif base_map_info['readbase'][bmi]=='*':
            annotate_info[bmi] = (annotate_info[bmi][0], round(np.mean(cur_pv_signal),3), round(np.std(cur_pv_signal),3), mstart1+start_in_pv, split_pos[spind+1][0]-start_in_pv if spind+1<len(split_pos) else mstart2+mlen2-mstart1-start_in_pv) if not split_pos==None else (annotate_info[bmi][0], all_mean, all_std, mstart1, mstart2+mlen2-mstart1)
           #if moptions['outLevel']<=myCom.OUTPUT_INFO:
            if moptions['outLevel']<=cur_level:
               #naind = bmi - (gipk - moptions['Resegment_wind'] if gipk - moptions['Resegment_wind']>-1 else 0)
               naind = bmi - (gipk - group_indel_pos[gipk][3][0] if gipk - group_indel_pos[gipk][3][0]>-1 else 0)
               if naind<0: naind=0;
               #print naind, spind, bmi, gipk - moptions['Resegment_wind'], refstr[naind], readstr[naind], (myCom.na_bp[m_event[-annotate_info[bmi][0]-1][4][2]] if forward_reverse=='-' else m_event[annotate_info[bmi][0]][4][2]), annotate_info[bmi]
               print naind, spind, bmi, gipk - group_indel_pos[gipk][3][0], refstr[naind], readstr[naind], (myCom.na_bp[m_event[-annotate_info[bmi][0]-1][4][2]] if forward_reverse=='-' else m_event[annotate_info[bmi][0]][4][2]), annotate_info[bmi]
            bmi += 1
            while base_map_info['readbase'][bmi]=='*':
               annotate_info[bmi] = (annotate_info[bmi][0], round(np.mean(cur_pv_signal),3), round(np.std(cur_pv_signal),3), mstart1+start_in_pv, split_pos[spind+1][0]-start_in_pv if spind+1<len(split_pos) else mstart2+mlen2-mstart1-start_in_pv) if not split_pos==None else (annotate_info[bmi][0], all_mean, all_std, mstart1, mstart2+mlen2-mstart1)
               #if moptions['outLevel']<=myCom.OUTPUT_INFO:
               if moptions['outLevel']<=cur_level:
                  #naind = bmi - (gipk - moptions['Resegment_wind'] if gipk - moptions['Resegment_wind']>-1 else 0)
                  naind = bmi - (gipk - group_indel_pos[gipk][3][0] if gipk - group_indel_pos[gipk][3][0]>-1 else 0)
                  if naind<0: naind=0;
                  #print naind, spind, bmi, gipk - moptions['Resegment_wind'], refstr[naind], readstr[naind], (myCom.na_bp[m_event[-annotate_info[bmi][0]-1][4][2]] if forward_reverse=='-' else m_event[annotate_info[bmi][0]][4][2]), annotate_info[bmi]
                  print naind, spind, bmi, gipk - group_indel_pos[gipk][3][0], refstr[naind], readstr[naind], (myCom.na_bp[m_event[-annotate_info[bmi][0]-1][4][2]] if forward_reverse=='-' else m_event[annotate_info[bmi][0]][4][2]), annotate_info[bmi]
               bmi += 1
            if base_map_info['readbase'][bmi] in myCom.acgt:
               annotate_info[bmi] = (annotate_info[bmi][0], round(np.mean(cur_pv_signal),3), round(np.std(cur_pv_signal),3), mstart1+start_in_pv, split_pos[spind+1][0]-start_in_pv if spind+1<len(split_pos) else mstart2+mlen2-mstart1-start_in_pv) if not split_pos==None else (annotate_info[bmi][0], all_mean, all_std, mstart1, mstart2+mlen2-mstart1)
               #if moptions['outLevel']<=myCom.OUTPUT_INFO:
               if moptions['outLevel']<=cur_level:
                  #naind = bmi - (gipk - moptions['Resegment_wind'] if gipk - moptions['Resegment_wind']>-1 else 0)
                  naind = bmi - (gipk - group_indel_pos[gipk][3][0] if gipk - group_indel_pos[gipk][3][0]>-1 else 0)
                  if naind<0: naind=0;
                  #print naind, spind, bmi, gipk - moptions['Resegment_wind'], refstr[naind], readstr[naind], (myCom.na_bp[m_event[-annotate_info[bmi][0]-1][4][2]] if forward_reverse=='-' else m_event[annotate_info[bmi][0]][4][2]), annotate_info[bmi],
                  print naind, spind, bmi, gipk - group_indel_pos[gipk][3][0], refstr[naind], readstr[naind], (myCom.na_bp[m_event[-annotate_info[bmi][0]-1][4][2]] if forward_reverse=='-' else m_event[annotate_info[bmi][0]][4][2]), annotate_info[bmi],
                  if forward_reverse=='-': print round(m_event[-annotate_info[bmi][0]-1][0],3), round(m_event[-annotate_info[bmi][0]-1][1],3),m_event[-annotate_info[bmi][0]-1][2],m_event[-annotate_info[bmi][0]-1][3],m_event[-annotate_info[bmi][0]-1][4], m_event[-annotate_info[bmi][0]-1][2]+m_event[-annotate_info[bmi][0]-1][3], start_in_pv+mstart1, mstart1+(len(pvsignals) if spind==len(split_pos)-1 else split_pos[spind+1][0]), (round(split_pos[spind][1],3) if spind>-1 else '')
                  else: print round(m_event[annotate_info[bmi][0]][0],3),round(m_event[annotate_info[bmi][0]][1],3),m_event[annotate_info[bmi][0]][2],m_event[annotate_info[bmi][0]][3],m_event[annotate_info[bmi][0]][4],m_event[annotate_info[bmi][0]][2]+m_event[annotate_info[bmi][0]][3], start_in_pv+mstart1, mstart1+(len(pvsignals) if spind==len(split_pos)-1 else split_pos[spind+1][0]), (round(split_pos[spind][1],3) if spind>-1 else '')
               bmi += 1
            else:
               if moptions['outLevel']<=myCom.OUTPUT_ERROR:
                  print 'Error!!! * followed by none nucleotide', base_map_info['readbase'][bmi], gipk, group_indel_pos[gipk], f5f
                  #print ' ref', ''.join(refstr), gipk-moptions['Resegment_wind'], gipk, group_indel_pos[gipk][2]+moptions['Resegment_wind']
                  print ' ref', ''.join(refstr), gipk-group_indel_pos[gipk][3][0], gipk, group_indel_pos[gipk][2]+group_indel_pos[gipk][3][1]
                  print 'read', ''.join(readstr), mstart1, mstart2+mlen2, len(pvsignals)
            if forward_reverse=='+': spind += 1;
            else: spind -= 1
         elif (base_map_info['readbase'][bmi] in myCom.acgt) or base_map_info['readbase'][bmi]=='-':
            annotate_info[bmi] = (annotate_info[bmi][0], round(np.mean(cur_pv_signal),3), round(np.std(cur_pv_signal),3), mstart1+start_in_pv, split_pos[spind+1][0]-start_in_pv if spind+1<len(split_pos) else mstart2+mlen2-mstart1-start_in_pv) if not split_pos==None else (annotate_info[bmi][0], all_mean, all_std, mstart1, mstart2+mlen2-mstart1)
            #if moptions['outLevel']<=myCom.OUTPUT_INFO:
            if moptions['outLevel']<=cur_level:
               #naind = bmi - (gipk - moptions['Resegment_wind'] if gipk - moptions['Resegment_wind']>-1 else 0)
               naind = bmi - (gipk - group_indel_pos[gipk][3][0] if gipk - group_indel_pos[gipk][3][0]>-1 else 0)
               if naind<0: naind=0;
               #print naind, spind, bmi, gipk - moptions['Resegment_wind'], refstr[naind], readstr[naind], (myCom.na_bp[m_event[-annotate_info[bmi][0]-1][4][2]] if forward_reverse=='-' else m_event[annotate_info[bmi][0]][4][2]), annotate_info[bmi],
               print naind, spind, bmi, gipk - group_indel_pos[gipk][3][0], refstr[naind], readstr[naind], (myCom.na_bp[m_event[-annotate_info[bmi][0]-1][4][2]] if forward_reverse=='-' else m_event[annotate_info[bmi][0]][4][2]), annotate_info[bmi],
               if base_map_info['readbase'][bmi] in myCom.acgt:
                  if forward_reverse=='-': print round(m_event[-annotate_info[bmi][0]-1][0],3), round(m_event[-annotate_info[bmi][0]-1][1],3),m_event[-annotate_info[bmi][0]-1][2],m_event[-annotate_info[bmi][0]-1][3],m_event[-annotate_info[bmi][0]-1][4], m_event[-annotate_info[bmi][0]-1][2]+m_event[-annotate_info[bmi][0]-1][3], start_in_pv+mstart1, mstart1+(len(pvsignals) if spind==len(split_pos)-1 else split_pos[spind+1][0]), (round(split_pos[spind][1],3) if spind>-1 else '')
                  else: print round(m_event[annotate_info[bmi][0]][0],3),round(m_event[annotate_info[bmi][0]][1],3),m_event[annotate_info[bmi][0]][2],m_event[annotate_info[bmi][0]][3],m_event[annotate_info[bmi][0]][4], start_in_pv+mstart1, m_event[annotate_info[bmi][0]][2]+m_event[annotate_info[bmi][0]][3], mstart1+(len(pvsignals) if spind==len(split_pos)-1 else split_pos[spind+1][0]), (round(split_pos[spind][1],3) if spind>-1 else '')
               else: print ''
            bmi += 1
            while bmi<len(base_map_info['readbase']) and base_map_info['readbase'][bmi]=='+':
                annotate_info[bmi] = (annotate_info[bmi][0], round(np.mean(cur_pv_signal),3), round(np.std(cur_pv_signal),3), mstart1+start_in_pv, split_pos[spind+1][0]-start_in_pv if spind+1<len(split_pos) else mstart2+mlen2-mstart1-start_in_pv) if not split_pos==None else (annotate_info[bmi][0], all_mean, all_std, mstart1, mstart2+mlen2-mstart1)
                #if moptions['outLevel']<=myCom.OUTPUT_INFO:
                if moptions['outLevel']<=cur_level:
                   #naind = bmi - (gipk - moptions['Resegment_wind'] if gipk - moptions['Resegment_wind']>-1 else 0)
                   naind = bmi - (gipk - group_indel_pos[gipk][3][0] if gipk - group_indel_pos[gipk][3][0]>-1 else 0)
                   if naind<0: naind=0;
                   #print naind, spind, bmi, gipk - moptions['Resegment_wind'], refstr[naind], readstr[naind], (myCom.na_bp[m_event[-annotate_info[bmi][0]-1][4][2]] if forward_reverse=='-' else m_event[annotate_info[bmi][0]][4][2]), annotate_info[bmi]
                   print naind, spind, bmi, gipk - group_indel_pos[gipk][3][0], refstr[naind], readstr[naind], (myCom.na_bp[m_event[-annotate_info[bmi][0]-1][4][2]] if forward_reverse=='-' else m_event[annotate_info[bmi][0]][4][2]), annotate_info[bmi]
                bmi += 1
            if forward_reverse=='+': spind += 1;
            else: spind -= 1
         else:
            if moptions['outLevel']<=myCom.OUTPUT_ERROR:
                print 'Error!!! Unsupported symbols in reads', base_map_info['readbase'][bmi], gipk, group_indel_pos[gipk], f5f
                #print ' ref', ''.join(refstr), gipk-moptions['Resegment_wind'], gipk, group_indel_pos[gipk][2]+moptions['Resegment_wind']
                print ' ref', ''.join(refstr), gipk-group_indel_pos[gipk][3][0], gipk, group_indel_pos[gipk][2]+group_indel_pos[gipk][3][1]
                print 'read', ''.join(readstr), mstart1, mstart2+mlen2, len(pvsignals)
            break;
   #if moptions['outLevel']<=myCom.OUTPUT_INFO:
   if moptions['outLevel']<=cur_level:
      print signalnum.items(),
   if moptions['outLevel']<=cur_level:
   #if moptions['outLevel']<=myCom.OUTPUT_INFO:
      aikeys = annotate_info.keys();
      for aik in aikeys:
         if len(annotate_info[aik])<3:
            print base_map_info['refbase'][aik-4:aik+4]
            print base_map_info['readbase'][aik-4:aik+4]
            print aik, annotate_info[aik],
            if forward_reverse=='+': print round(m_event[annotate_info[aik][0]][0],3),round(m_event[annotate_info[aik][0]][1],3),m_event[annotate_info[aik][0]][2],m_event[annotate_info[aik][0]][3],m_event[annotate_info[aik][0]][4]
            else: print round(m_event[-annotate_info[aik][0]-1][0],3),round(m_event[-annotate_info[aik][0]-1][1],3),m_event[-annotate_info[aik][0]-1][2],m_event[-annotate_info[aik][0]-1][3],m_event[-annotate_info[aik][0]-1][4]
            for gprink in group_indel_pos_keys:
               print '\t', gprink, group_indel_pos[gprink]

   return (annotate_info, signalnum)

#
#
#
def find_sp(pvsignals, gipk, refstr,readstr, expectna, Resegment_signal_wind, group_indel_pos, mstart1,mstart2,mlen2, f5f,cur_level,moptions):
      #cur_level=myCom.OUTPUT_ERROR

      #pairdif = defaultdict();
      #for pvi in range(len(pvsignals)):
      #   for pvj in range(pvi+1, (pvi+Resegment_signal_wind*2 if pvi+Resegment_signal_wind*2<len(pvsignals) else len(pvsignals))):
      #      pairdif[(pvi, pvj)] = np.absolute(pvsignals[pvi]-pvsignals[pvj])
      dif_rank = []
      if moptions['outLevel']<=cur_level: #myCom.OUTPUT_INFO:
      #if moptions['outLevel']<=myCom.OUTPUT_INFO:
         #print ' ref', ''.join(refstr), gipk-moptions['Resegment_wind'], gipk, group_indel_pos[gipk][2]+moptions['Resegment_wind']
         print ' ref', ''.join(refstr), gipk-group_indel_pos[gipk][3][0], gipk, group_indel_pos[gipk][2]+group_indel_pos[gipk][3][1]+1
         print 'read', ''.join(readstr), mstart1, mstart2+mlen2, len(pvsignals)
      for pvi in range(Resegment_signal_wind, len(pvsignals)-Resegment_signal_wind+1):
         #withdif = 0; btwdif = 0; withinnum = 0; btwnum = 0;
         #for pvj in range(pvi-Resegment_signal_wind, pvi):
         #   for pvk in range(pvj+1, pvi):
         #      withdif += pairdif[(pvj, pvk)]; withinnum+=1;
         #   for pvk in range(pvi, pvi+Resegment_signal_wind):
         #      btwdif += pairdif[(pvj, pvk)]; btwnum += 1;
         #for pvj in range(pvi, pvi+Resegment_signal_wind):
         #   for pvk in range(pvj+1, pvi+Resegment_signal_wind):
         #      withdif += pairdif[(pvj, pvk)]; withinnum+=1;
         #if withdif<0.001: withdif=0.001
         ### consider both within dif
         #dif_rank.append((pvi, round((btwdif/(Resegment_signal_wind*Resegment_signal_wind))/(withdif/(Resegment_signal_wind*(Resegment_signal_wind-1))), 5)))
         ### only considered larger within dif
         #withdif1 = 0; withdif2 = 0; btwdif = 0; withinnum1 = 0; withinnum2 = 0; btwnum = 0;
         #for pvj in range(pvi-Resegment_signal_wind, pvi):
         #   for pvk in range(pvj+1, pvi):
         #      withdif1 += pairdif[(pvj, pvk)]; withinnum1+=1;
         #   for pvk in range(pvi, pvi+Resegment_signal_wind):
         #      btwdif += pairdif[(pvj, pvk)]; btwnum += 1;
         #for pvj in range(pvi, pvi+Resegment_signal_wind):
         #   for pvk in range(pvj+1, pvi+Resegment_signal_wind):
         #      withdif2 += pairdif[(pvj, pvk)]; withinnum2+=1;
         #withdif = max(withdif1/withinnum1, withdif2/withinnum2)
         #btwdif = btwdif/btwnum
         #if withdif<0.001: withdif=0.001
         ######
         ##dif_rank.append((pvi, round(btwdif/withdif, 5)))
         #dif_rank.append((pvi, round(btwdif, 5)))
         dif_rank.append((pvi, round(np.absolute(np.mean(pvsignals[pvi-Resegment_signal_wind:pvi])-np.mean(pvsignals[pvi:pvi+Resegment_signal_wind])), 9)))
         if moptions['outLevel']<=cur_level: #myCom.OUTPUT_INFO:
         #if moptions['outLevel']<=myCom.OUTPUT_INFO:
            #print '\tInfo', pvi, ('%.3f ' % pvsignals[pvi]), len(pvsignals), Resegment_signal_wind, 'within', withinnum, np.round(withdif/(Resegment_signal_wind*(Resegment_signal_wind-1)), 6), 'betw', btwnum, np.round(btwdif/(Resegment_signal_wind*Resegment_signal_wind), 6), 'expected-cur', expectna, dif_rank[-1]
            print '\tInfo', pvi, ('%.3f ' % pvsignals[pvi]), len(pvsignals), Resegment_signal_wind, 'expected-cur', expectna, dif_rank[-1]

      #if moptions['outLevel']<=myCom.OUTPUT_INFO:
      if moptions['outLevel']<=cur_level:
         for pvi in range(Resegment_signal_wind): print ('%d:%.3f ' % (pvi, pvsignals[pvi])),
         for dr in range(len(dif_rank)):
            print ('*%.3f* %d:%.3f ' % (dif_rank[dr][1], pvi, pvsignals[dif_rank[dr][0]])),
            pvi +=1
         for pvi in range(len(pvsignals)-Resegment_signal_wind-1, len(pvsignals)):
            print ('%d:%.3f ' % (pvi, pvsignals[pvi])),
            pvi += 1
         print ''

      sorted_dif_rank = sorted(dif_rank, key=lambda mpv: (-mpv[1]))
      split_pos = []
      for dr in sorted_dif_rank:
         contain_neighbor = False;
         for sp in split_pos:
            #if -Resegment_signal_wind<dr[0]-sp[0]<Resegment_signal_wind-1: contain_neighbor = True;
            if -moptions['MinNumSignal']<dr[0]-sp[0]<moptions['MinNumSignal']: contain_neighbor = True
         if not contain_neighbor:
            split_pos.append(dr)
         if len(split_pos)==expectna-1: break;
      if len(split_pos)<expectna-1:
         if moptions['outLevel']<=myCom.OUTPUT_DEBUG or Resegment_signal_wind==2: #myCom.OUTPUT_ERROR:
            print 'Warning Cannot find large enough split', len(split_pos), expectna, len(pvsignals), f5f, len(dif_rank), Resegment_signal_wind, len(pvsignals)-Resegment_signal_wind+1
            #print pvsignals
            print split_pos
            print sorted(split_pos, key=lambda mpv: (mpv[0]))
            print gipk, group_indel_pos[gipk]
            #print ' ref', ''.join(refstr), gipk-moptions['Resegment_wind'], gipk, group_indel_pos[gipk][2]+moptions['Resegment_wind']
            print ' ref', ''.join(refstr), gipk-group_indel_pos[gipk][3][0], gipk, group_indel_pos[gipk][2]+group_indel_pos[gipk][3][1]+1
            print 'read', ''.join(readstr), mstart1, mstart2+mlen2, len(pvsignals)
            for pvi in range(Resegment_signal_wind): print ('%.1f ' % pvsignals[pvi]),
            for dr in range(len(dif_rank)):
               if dif_rank[dr] in split_pos:
                  print ('***|%.1f|*** %.1f ' % (dif_rank[dr][1], pvsignals[dif_rank[dr][0]])),
               else:
                  print ('<%.1f> %.1f ' % (dif_rank[dr][1], pvsignals[dif_rank[dr][0]])),
            for pvi in range(len(pvsignals)-Resegment_signal_wind+1, len(pvsignals)):
               print ('%.1f ' % pvsignals[pvi]),
            print ''

      split_pos = sorted(split_pos, key=lambda mpv: (mpv[0]))
      if moptions['outLevel']<=cur_level: #or Resegment_signal_wind<moptions['Resegment_signal_wind']: #myCom.OUTPUT_INFO:
         for dr in split_pos:
            print 'sp', dr

      if len(split_pos)==expectna-1: return split_pos

def calculate_mean_std(m_event, event_ind, forward_reverse, raw_pv, moptions, sp_param, isprint=False):
   cur_level = myCom.OUTPUT_DEBUG; #myCom.OUTPUT_INFO
   if forward_reverse=='-':
      if isprint:
         print '\tindextest', len(raw_pv), len(m_event), -event_ind-1, m_event[-event_ind-1], len(raw_pv), m_event[-event_ind-1][2], m_event[-event_ind-1][2]+m_event[-event_ind-1][3]
      pvsignal = raw_pv[m_event[-event_ind-1][2]:(m_event[-event_ind-1][2]+m_event[-event_ind-1][3])]
      if moptions['outLevel']<=cur_level: p_mean = m_event[-event_ind-1][0]; p_std = m_event[-event_ind-1][1]
   else:
      if isprint: print '\tindextest', len(raw_pv), len(m_event), event_ind, m_event[event_ind], len(raw_pv), m_event[event_ind][2], m_event[event_ind][2]+m_event[event_ind][3]
      pvsignal = raw_pv[m_event[event_ind][2]:(m_event[event_ind][2]+m_event[event_ind][3])]
      if moptions['outLevel']<=cur_level: p_mean = m_event[event_ind][0]; p_std = m_event[event_ind][1]

   c_mean = round(np.mean(pvsignal), 3)
   c_std = round(np.std(pvsignal), 3)
   if moptions['outLevel']<=cur_level and forward_reverse=='-': #:
      print forward_reverse, len(raw_pv), event_ind, len(m_event), 'start_len',
      if forward_reverse=='-':
         print m_event[::-1][event_ind][2], m_event[::-1][event_ind][3],
      else:
         print m_event[event_ind][2], m_event[event_ind][3],
      print 'mean-std',c_mean, c_std,
      if forward_reverse=='-':
         print m_event[::-1][event_ind][0],  m_event[::-1][event_ind][1],
      else: print m_event[event_ind][0],  m_event[event_ind][1],
      print 'dif', np.round(np.absolute(c_mean-p_mean), 3), np.round(np.absolute(c_std-p_std),3)

   if moptions['outLevel']<=cur_level:
      if (np.absolute(c_mean-p_mean)>1 or np.absolute(c_std-p_std)>1):
         return (c_mean, c_std, False);
      else: return (c_mean, c_std, True);
   else: return (c_mean, c_std, True)

#
#
#
def fix_repeat_del(base_map_info, moptions, sp_param, f5f):
   if moptions['outLevel']<=myCom.OUTPUT_DEBUG:
      print 'test in fix_repeat_del', len(base_map_info), len(base_map_info['refbase'])

   last_non_indel = 1; last_is_repeat = False;
   for rvind in range(1, len(base_map_info)+1):
      if not base_map_info['readbase'][-rvind]=='-':
         last_non_indel = rvind; last_is_repeat = False;
      else:
         if base_map_info['refbase'][-rvind]==base_map_info['refbase'][-last_non_indel] and base_map_info['refbase'][-rvind] in myCom.acgt:
            if last_non_indel==rvind-1 and (base_map_info['readbase'][-last_non_indel]==base_map_info['refbase'][-last_non_indel]): last_is_repeat = True;
            if last_is_repeat:
               # does not consider now
               #base_map_info['readbase'][-rvind]='*' #'-'+base_map_info['readbase'][-last_non_indel]
               if moptions['outLevel']<=myCom.OUTPUT_DEBUG:
                  if last_non_indel==rvind-1 and (not base_map_info['readbase'][-last_non_indel]==base_map_info['refbase'][-last_non_indel]): print 'Warning!!! not equal', base_map_info['readbase'][-last_non_indel], base_map_info['refbase'][-last_non_indel]
         else:
            if last_is_repeat:
               if moptions['outLevel']<=myCom.OUTPUT_DEBUG:
                  if last_non_indel==1:
                     print ' Refr', ''.join(base_map_info['refbase'][(-rvind-3):]), -rvind, len(base_map_info)
                     print 'readr', ''.join(base_map_info['readbase'][(-rvind-3):]), -rvind, len(base_map_info)
                  else:
                     print ' Refr', ''.join(base_map_info['refbase'][(-rvind-3):(-last_non_indel+3)]), -rvind, -last_non_indel, len(base_map_info)
                     print 'readr', ''.join(base_map_info['readbase'][(-rvind-3):(-last_non_indel+3)]), -rvind, -last_non_indel, len(base_map_info)
            last_is_repeat = False;
      if moptions['outLevel']<=myCom.OUTPUT_DEBUG:
         print rvind, base_map_info['refbase'][-rvind], base_map_info['readbase'][-rvind], last_non_indel, last_is_repeat

   last_non_indel = 0; last_is_repeat = False; event_ind = -1;
   # del_pos = defaultdict(); ins_pos = defaultdict();
   indel_pos = defaultdict();
   if moptions['outLevel']<=myCom.OUTPUT_DEBUG: firstdif = True;
   for bmi in range(len(base_map_info)):
      if base_map_info['readbase'][bmi] in myCom.acgt:
         event_ind += 1
         if base_map_info['refbase'][bmi]=='-':
            indel_pos[bmi] = (event_ind, 1)
            #ins_pos[bmi] = event_ind
            if moptions['outLevel']<=myCom.OUTPUT_DEBUG:
               print '\tNonequal', bmi, event_ind, len(indel_pos), event_ind+len(indel_pos);

      if not base_map_info['readbase'][bmi]=='-':
         last_non_indel = bmi; last_is_repeat = False;
         if base_map_info['readbase'][bmi]=='*':
            indel_pos[bmi] = (event_ind, 0)
            #del_pos[bmi] = event_ind
      else:
         if base_map_info['refbase'][bmi]==base_map_info['refbase'][last_non_indel] and base_map_info['refbase'][bmi] in myCom.acgt:
            if last_non_indel==bmi-1 and (base_map_info['readbase'][last_non_indel]==base_map_info['refbase'][last_non_indel]): last_is_repeat = True;
            if last_is_repeat:
               #### does not considered now
               #base_map_info['readbase'][bmi]='+' #'-'+base_map_info['readbase'][last_non_indel]
               if moptions['outLevel']<=myCom.OUTPUT_DEBUG:
                  if last_non_indel==bmi-1 and (not base_map_info['readbase'][last_non_indel]==base_map_info['refbase'][last_non_indel]): print 'Warning!!! not equal', base_map_info['readbase'][last_non_indel], base_map_info['refbase'][last_non_indel]
         else:
            if last_is_repeat:
               if moptions['outLevel']<=myCom.OUTPUT_DEBUG:
                  print ' Ref', ''.join(base_map_info['refbase'][last_non_indel-3:(bmi+3)]), bmi, last_non_indel, len(base_map_info)
                  print 'read', ''.join(base_map_info['readbase'][last_non_indel-3:(bmi+3)]), bmi, last_non_indel, len(base_map_info), '|', bmi, event_ind, len(indel_pos), event_ind+len(indel_pos)
            last_is_repeat = False;

         if base_map_info['refbase'][bmi] in myCom.acgt:
            #del_pos[bmi] = event_ind
            if last_is_repeat: indel_pos[bmi] = (event_ind, 0)
            else: indel_pos[bmi] = (event_ind, -1)
         else:
            print 'Warning!!! both unknown', base_map_info['refbase'][bmi], base_map_info['readbase'][bmi], f5f
      # for test purpose
      if moptions['outLevel']<=myCom.OUTPUT_DEBUG:
         if not bmi==event_ind+len(indel_pos) and firstdif:
            print '\tNonequal', bmi, event_ind, len(indel_pos), event_ind+len(indel_pos);
            print '\t', ''.join(base_map_info['refbase'][bmi-10:bmi]), ''.join(base_map_info['refbase'][bmi:bmi+10])
            print '\t', ''.join(base_map_info['readbase'][bmi-10:bmi]), ''.join(base_map_info['readbase'][bmi:bmi+10])
            firstdif = False

   for bmi in range(3, len(base_map_info)-2):
      if base_map_info['readbase'][bmi] in ['-','+','*']:
         if ''.join(base_map_info['refbase'][bmi-2:bmi+3]) == ''.join(base_map_info['refbase'][bmi-3:bmi+2]):
            base_map_info['readbase'][bmi] = '~'
            if base_map_info['readbase'][bmi-1] in ['-','+','*']:
               base_map_info['readbase'][bmi-1] = '~'

   # for test purpose
   if moptions['outLevel']<=myCom.OUTPUT_DEBUG: #myCom.OUTPUT_INFO:
      gapinref = 0; gapinread = 0;
      for bmi in range(len(base_map_info)):
         if base_map_info['refbase'][bmi] not in myCom.acgt: gapinref += 1;
         if base_map_info['readbase'][bmi] not in myCom.acgt: gapinread += 1;
      print 'gapinfo', gapinref, gapinread, len(base_map_info)
   return indel_pos

#
#
def group_indel(indel_pos, m_event, base_map_info, forward_reverse, moptions, sp_param):
   indel_pos_keys = indel_pos.keys(); indel_pos_keys.sort();
   group_indel_Ipos = defaultdict();   pre_ipk = None;
   for ipk in indel_pos_keys: #### merge close indel ong bmi of base_map_info
      if pre_ipk==None or (not ipk-group_indel_Ipos[pre_ipk][1]<=2):
         group_indel_Ipos[ipk] = (ipk, ipk)
         pre_ipk = ipk
      else: #ipk-group_indel_Ipos[pre_ipk][1]<3:
         group_indel_Ipos[pre_ipk] = (group_indel_Ipos[pre_ipk][0], ipk)
   group_indel_pos = defaultdict();   pre_ipk = None; lastipk = []
   indel_pos_keys = group_indel_Ipos.keys(); indel_pos_keys.sort();
   for ipk_ind in range(len(indel_pos_keys)): #####
      ipk = indel_pos_keys[ipk_ind]
      i1pk, i2pk = group_indel_Ipos[ipk]
      if not i1pk==ipk: print 'ipk not equal', ipk, group_indel_Ipos[ipk]
      leftnum = 0; rightnum = 0
      if base_map_info['refbase'][i1pk]=='-':
         if indel_pos[i1pk][0]-1>=0: start_ev_ind = indel_pos[i1pk][0]-1; leftnum=1
         else: start_ev_ind = 0
      else:
         start_ev_ind = indel_pos[i1pk][0]; leftnum=1
      if indel_pos[i2pk][0]+1<len(m_event): end_evnt_ind = indel_pos[i2pk][0]+1; rightnum=1
      else: end_evnt_ind = len(m_event)-1
      #####################################################
      while True:
             if forward_reverse=='+':
                mstart1,mend1=m_event[start_ev_ind][2], m_event[start_ev_ind][2]+m_event[start_ev_ind][3]
                mstart2,mend2=m_event[end_evnt_ind][2], m_event[end_evnt_ind][2]+m_event[end_evnt_ind][3]
             else:
                mstart1,mend1=m_event[-end_evnt_ind-1][2],m_event[-end_evnt_ind-1][2]+m_event[-end_evnt_ind-1][3]
                mstart2,mend2=m_event[-start_ev_ind-1][2],m_event[-start_ev_ind-1][2]+m_event[-start_ev_ind-1][3]
             numsignals = mend2-mstart1
             if numsignals<1:print "FATAL error!!! negative number of signals",start_ev_ind,end_evnt_ind,pre_ipk,ipk; sys.stdout.flush()
             else:
                expectna = 0;
                for bmi in range(i1pk-leftnum, i2pk+rightnum+1):
                   if bmi < 0: continue;
                   if not bmi<len(base_map_info): break;
                   if base_map_info['refbase'][bmi]=='-': continue;
                   #if base_map_info['readbase'][bmi] in ['+', '*']:continue;
                   if base_map_info['readbase'][bmi]=='~':
                      if (bmi>0 and base_map_info['readbase'][bmi-1]=='~'): continue;
                   expectna += 1;
                if numsignals>(expectna+(1 if expectna*moresignalperc<1 else int(expectna*moresignalperc+0.5)))*moptions['MinNumSignal']:
                   #print 'Test', ipk, i2pk,(start_ev_ind, end_evnt_ind, i2pk, (leftnum, rightnum)), expectna, numsignals
                   #if 1000<ipk<1120:
                   #   print 'Test', ipk, i2pk,(start_ev_ind, end_evnt_ind, i2pk, (leftnum, rightnum)), expectna, numsignals
                   break;
             #print '\tTest', start_ev_ind, pre_ipk, end_evnt_ind, len(m_event), leftnum,rightnum, (group_indel_pos[pre_ipk] if not pre_ipk==None else "")
             if (start_ev_ind==0 or ((not pre_ipk==None) and start_ev_ind<=group_indel_pos[pre_ipk][1])) and end_evnt_ind==len(m_event)-1: break;
             if (pre_ipk==None and start_ev_ind>0) or ((not pre_ipk==None) and start_ev_ind>group_indel_pos[pre_ipk][1]):
                start_ev_ind -= 1; leftnum += 1
             elif (not pre_ipk==None):
                #print 'mb1', pre_ipk, group_indel_pos[pre_ipk],'<<<', start_ev_ind, i1pk, leftnum, len(group_indel_pos),i2pk, rightnum
                start_ev_ind = group_indel_pos[pre_ipk][0];
                i1pk= pre_ipk;
                leftnum=group_indel_pos[pre_ipk][3][0];
                del group_indel_pos[pre_ipk]
                pre_ipk = lastipk.pop()
                #print 'mb2',pre_ipk, (group_indel_pos[pre_ipk][1] if not pre_ipk==None else []), start_ev_ind, i1pk, leftnum, len(group_indel_pos)
             if end_evnt_ind<len(m_event)-1:
                rightnum+=1
                while True:
                   if (base_map_info['readbase'][i2pk+rightnum] in myCom.acgt) and (base_map_info['refbase'][i2pk+rightnum] in myCom.acgt):
                      end_evnt_ind += 1; break;
                   elif (base_map_info['readbase'][i2pk+rightnum] in myCom.acgt) and (base_map_info['refbase'][i2pk+rightnum] not in myCom.acgt):
                      end_evnt_ind += 1; rightnum+=1
                   elif (base_map_info['readbase'][i2pk+rightnum] not in myCom.acgt) and (base_map_info['refbase'][i2pk+rightnum] in myCom.acgt):
                      rightnum+=1
                   else: print 'Fatal error both unknow', i1pk, start_ev_ind, end_evnt_ind, i2pk, (leftnum, rightnum)
      #print 'Test', (start_ev_ind, end_evnt_ind, i2pk, (leftnum, rightnum)), expectna, numsignals
      #######################################################
      if pre_ipk==None or (start_ev_ind>group_indel_pos[pre_ipk][1]):
         group_indel_pos[i1pk] = (start_ev_ind, end_evnt_ind, i2pk, (leftnum, rightnum))
         lastipk.append(pre_ipk)
         pre_ipk = i1pk
      elif start_ev_ind<=group_indel_pos[pre_ipk][1]:
         if end_evnt_ind>=group_indel_pos[pre_ipk][1]:
            group_indel_pos[pre_ipk] = (group_indel_pos[pre_ipk][0], end_evnt_ind, i2pk, (group_indel_pos[pre_ipk][3][0], rightnum))
         #else:
         #   group_indel_pos[pre_ipk] = (group_indel_pos[pre_ipk][0], end_evnt_ind if end_evnt_ind>group_indel_pos[pre_ipk][1] else group_indel_pos[pre_ipk][1], i2pk, (group_indel_pos[pre_ipk][3][0], rightnum if end_evnt_ind>group_indel_pos[pre_ipk][1] else group_indel_pos[pre_ipk][3][1]))


   ''' new but might not corrected.
   group_indel_pos = defaultdict();   pre_ipk = None;
   for ipk in indel_pos_keys:
      leftnum = 0; rightnum = 0
      if base_map_info['refbase'][ipk]=='-':
          if indel_pos[ipk][0]-1>=0: start_ev_ind = indel_pos[ipk][0]-1; leftnum=1
          else: start_ev_ind = 0;
          if indel_pos[ipk][0]+1<len(m_event): end_evnt_ind = indel_pos[ipk][0]+1; rightnum=1
          else: end_evnt_ind = len(m_event)-1
      else:
         firstrun=True; start_ev_ind=indel_pos[ipk][0]; end_evnt_ind=indel_pos[ipk][0]
         while True:
             if not firstleft:
                if (pre_ipk==None) or start_ev_ind>group_indel_pos[pre_ipk][1]:
                   (start_ev_ind, leftnum) = (start_ev_ind-1, leftnum+1) if start_ev_ind-1>=0 else (0, leftnum);
                #elif (not pre_ipk==None) and start_ev_ind==group_indel_pos[pre_ipk][1]:
                #   ipk = pre_ipk; rightnum = leftnum+rightnum+group_indel_pos[pre_ipk][3][1]
                #   leftnum = group_indel_pos[pre_ipk][3][0]
             (end_evnt_ind,rightnum)=(end_evnt_ind+1,rightnum+1) if end_evnt_ind+1<len(m_event) else (len(m_event)-1, rightnum)
             if forward_reverse=='+':
                mstart1, mend1 = m_event[start_ev_ind][2], m_event[start_ev_ind][2]+m_event[start_ev_ind][3]
                mstart2, mend2 = m_event[end_evnt_ind][2], m_event[end_evnt_ind][2]+m_event[end_evnt_ind][3]
             else:
                mstart1, mend1 = m_event[-end_evnt_ind-1][2],m_event[-end_evnt_ind-1][2]+m_event[-end_evnt_ind-1][3]
                mstart2, mend2 = m_event[-start_ev_ind-1][2],m_event[-start_ev_ind-1][2]+m_event[-start_ev_ind-1
             numsignals = mend2-mstart1
             if numsignals<1:print "FATAL error!!! negative number of signals",start_ev_ind,end_evnt_ind,pre_ipk,ipk; sys.stdout.flush()
             else:
                expectna = 0;
                for bmi in range(ipk-leftnum, ipk+rightnum+1):
                   if bmi < 0: bmi += 1; continue;
                   if not bmi<len(base_map_info): break;
                   if base_map_info['refbase'][bmi]=='-': continue;
                   if base_map_info['readbase'][bmi] in ['+', '*']:continue;
                   expectna += 1;
                if numsignals>(expectna+2)*moptions['Resegment_signal_wind']: break;

             if (start_ev_ind==0 or ((not pre_ipk==None) and start_ev_ind==group_indel_pos[pre_ipk][1])) and end_evnt_ind==len(m_event)-1: break;
             firstrun = False;
      if pre_ipk==None or (start_ev_ind>group_indel_pos[pre_ipk][1]):
         group_indel_pos[ipk] = (start_ev_ind, end_evnt_ind, ipk, (leftnum, rightnum))
         pre_ipk = ipk
      elif start_ev_ind<=group_indel_pos[pre_ipk][1]:
         group_indel_pos[pre_ipk] = (group_indel_pos[pre_ipk][0], end_evnt_ind, ipk, (group_indel_pos[pre_ipk][3][0], rightnum))
      '''
      # old fix length for the number of bases for error corrections
      #if base_map_info['refbase'][ipk]=='-':
      #   start_ev_ind = indel_pos[ipk][0]-moptions['Resegment_wind'] if indel_pos[ipk][0]-moptions['Resegment_wind']>=0 else 0
      #else:
      #   start_ev_ind = indel_pos[ipk][0]-moptions['Resegment_wind']+1 if indel_pos[ipk][0]-moptions['Resegment_wind']+1>=0 else 0
      #end_evnt_ind = indel_pos[ipk][0]+moptions['Resegment_wind'] if indel_pos[ipk][0]+moptions['Resegment_wind']<len(m_event) else len(m_event)-1
      #if pre_ipk==None or (start_ev_ind>=group_indel_pos[pre_ipk][1]):
      #   group_indel_pos[ipk] = (start_ev_ind, end_evnt_ind, ipk)
      #   pre_ipk = ipk
      #elif start_ev_ind<group_indel_pos[pre_ipk][1]:
      #   group_indel_pos[pre_ipk] = (group_indel_pos[pre_ipk][0], end_evnt_ind, ipk)

   #if moptions['outLevel']<=myCom.OUTPUT_INFO: #myCom.OUTPUT_DEBUG: #myCom.OUTPUT_INFO:
   if moptions['outLevel']<=myCom.OUTPUT_DEBUG:
   #if moptions['outLevel']<=myCom.OUTPUT_INFO and forward_reverse=='-':
      group_indel_pos_keys = group_indel_pos.keys(); group_indel_pos_keys.sort();
      #print 'Indel   Ref', ''.join(base_map_info['refbase'][1000:1120]), np.sum(base_map_info['refbase'][1000:1120]=='-')
      #print 'Indel read1', ''.join(base_map_info['readbase'][1000:1120]), np.sum(base_map_info['readbase'][1000:1120]=='-')

      for gipk in group_indel_pos_keys:
         #if not (1000<gipk<1120): continue
         #numn = moptions['Resegment_wind']-1
         #start_bmi = gipk-1-numn if gipk-1-numn>0 else 0;
         #end_bmi = group_indel_pos[gipk][2]+2+numn
         start_bmi = gipk - group_indel_pos[gipk][3][0] if gipk - group_indel_pos[gipk][3][0]>0 else 0
         end_bmi = group_indel_pos[gipk][2] + group_indel_pos[gipk][3][1] + 1

         print 'Indel   Ref', ''.join(base_map_info['refbase'][start_bmi:end_bmi]), gipk, group_indel_pos[gipk][2], group_indel_pos[gipk][3]
         print 'Indel read1', ''.join(base_map_info['readbase'][start_bmi:end_bmi]), start_bmi, end_bmi, group_indel_pos[gipk]

         if forward_reverse=='-':
            mseq = ''.join([myCom.na_bp[event_model_state[2]] for event_model_state in m_event[::-1][group_indel_pos[gipk][0]:(group_indel_pos[gipk][1]+1)]['model_state']])
         else:
            mseq = ''.join([event_model_state[2] for event_model_state in m_event[group_indel_pos[gipk][0]:(group_indel_pos[gipk][1]+1)]['model_state']])
         #print 'Indel read2',''.join([mseq[:moptions['Resegment_wind']], (' '*(group_indel_pos[gipk][2]-gipk+1)), mseq[moptions['Resegment_wind']:]]), '::', group_indel_pos[gipk][0], group_indel_pos[gipk][1]+1
         #print 'Indel read2',''.join([mseq[:group_indel_pos[gipk][3][0]+1], (' '*(group_indel_pos[gipk][2]-gipk+1)), mseq[group_indel_pos[gipk][3][0]:]]), '::', group_indel_pos[gipk][0], group_indel_pos[gipk][1]+1
         #print 'Indel read2',mseq
      #sys.exit(1)
   return group_indel_pos

#
#
def handle_line(moptions, sp_param, f5align):
   lsp = sp_param['line'].split('\t')
   qname, flag, rname, pos, mapq, cigar, _, _, _, seq, _ = lsp[:11]
   if qname=='*': sp_param['f5status'] = "qname is *"
   elif int(mapq)==255: sp_param['f5status'] = "mapq is 255"
   elif int(pos)==0: sp_param['f5status'] = "pos is 0"
   elif cigar=='*': sp_param['f5status'] = "cigar is *"
   elif rname=='*': sp_param['f5status'] = "rname is *"
   if not sp_param['f5status']=="": return qname

   if (qname not in f5align) or f5align[qname][0]<int(mapq):
      f5align[qname] = (int(mapq), int(flag), rname, int(pos), cigar, seq)

   return qname

def correctAndAnnotate_handler(moptions, h5files_Q, failed_Q, resegment_neighobr_Q):
   while not h5files_Q.empty():
      try:
         f5files = h5files_Q.get(block=False)
      except:
         break;

      moptions_dc = copy.deepcopy(moptions)
      correctAndAnnotate(moptions_dc, f5files)

      for errtype, errfiles in moptions_dc["Error"].iteritems():
         failed_Q.put((errtype, errfiles));
         if moptions['outLevel']<=myCom.OUTPUT_DEBUG: print errtype, len(errfiles)

      for snk, snv in moptions_dc['signalneighbor'].iteritems():
         resegment_neighobr_Q.put((snk, snv))


def correctAndAnnotate_manager(moptions):
   start_time = time.time();

   #moptions['Resegment_signal_wind.5'] = int(moptions['Resegment_signal_wind']/2.0+0.5)
   #if moptions['Resegment_signal_wind.5']<3: moptions['Resegment_signal_wind.5'] = 3
   #f5files = glob.glob(os.path.join(moptions['wrkBase1'],"*/*.fast5" ))
   #f5files = glob.glob(os.path.join(moptions['wrkBase1'],"0/*.fast5" ))
   #f5files = glob.glob(os.path.join(moptions['wrkBase1'],"*.fast5" ))
   #f5files = ['tofixerror/VAIODG_20170628_FNFAH09129_MN19183_mux_scan_TrainSeq3_Bio_pur_SpeI_cut_31257_ch166_read31_strand.fast5']

   f5files = glob.glob(os.path.join(moptions['wrkBase1'],"*.fast5" ))
   if moptions['recursive']==1:
      #f5files = glob.glob(os.path.join(moptions['wrkBase1'],"*/*.fast5" ))
      f5files.extend(glob.glob(os.path.join(moptions['wrkBase1'],"*/*.fast5" )))
      f5files.extend(glob.glob(os.path.join(moptions['wrkBase1'],"*/*/*.fast5" )))
      f5files.extend(glob.glob(os.path.join(moptions['wrkBase1'],"*/*/*/*.fast5" )))
   #else: f5files = glob.glob(os.path.join(moptions['wrkBase1'],"*.fast5" ))
   #f5files = ['tofixerror/VAIODG_20170628_FNFAH09129_MN19183_mux_scan_TrainSeq3_Bio_pur_SpeI_cut_31257_ch166_read31_strand.fast5']
   #f5files = ['SingleModDs/IdU/ds400/26/VAIODG_20170809_FNFAH20109_MN19183_mux_scan_TrainSeq3_IdU_SpeI_cut_96151_ch473_read32_strand.fast5']


   get_kmer_corrected_info(moptions)

   pmanager = multiprocessing.Manager();
   h5files_Q = pmanager.Queue();
   failed_Q = pmanager.Queue()
   resegment_neighobr_Q = pmanager.Queue()

   h5_batch = [];
   for f5f in f5files:
      h5_batch.append(f5f);
      if len(h5_batch)==moptions['files_per_thread']:
         h5files_Q.put(h5_batch)
         h5_batch = [];
   if len(h5_batch)>0:
      h5files_Q.put(h5_batch)

   share_var = (moptions, h5files_Q, failed_Q, resegment_neighobr_Q)
   handlers = []
   for hid in xrange(moptions['threads']):
      p = multiprocessing.Process(target=correctAndAnnotate_handler, args=share_var);
      p.start();
      handlers.append(p);

   failed_files = defaultdict(list);
   reseg_n = defaultdict(int);
   while any(p.is_alive() for p in handlers):
      try:
         errk, fns = failed_Q.get(block=False);
         failed_files[errk].extend(fns)
         snk, snv = resegment_neighobr_Q.get(block=False);
         reseg_n[snk] += snv
      except:
         time.sleep(1);
         continue;

   print 'Error information for different fast5 files:'
   #for errtype, errfiles in failed_files.iteritems():
   #   print '\t', errtype, errfiles
   for errtype, errfiles in failed_files.iteritems():
      print '\t', errtype, len(errfiles)

   print 'Resegment information for signal neighbors:',
   for snk, snv in reseg_n.iteritems():
      print ('%d(%d)' % (snk, snv)),
   print ''

   end_time = time.time();
   print ("Total consuming time %d" % (end_time-start_time))

if __name__=='__main__':
   if len(sys.argv)>4:
      moptions = {}
      moptions['basecall_1d'] = 'Basecall_1D_000'
      moptions['basecall_2strand'] = 'BaseCalled_template'

      basefolder = sys.argv[1];
      moptions['wrkBase1'] = basefolder
      moptions['targetfolder'] = sys.argv[2]
      moptions['targetfileid'] = sys.argv[3]
      moptions['Ref'] = sys.argv[4]
      moptions['kmer_model_file'] = 'scripts/kmer_model/r9.4_450bps.nucleotide.5mer.template.model'
      moptions['alignStr'] = 'bwa'
      #moptions['alignStr'] = moptions['alignStr'] + ' %s %s'

      moptions['outLevel'] = myCom.OUTPUT_WARNING
      moptions['outLevel'] = myCom.OUTPUT_INFO
      moptions['Resegment_wind'] = 4;
      moptions['Resegment_signal_wind'] = 3; #4;

      moptions['threads'] = 8
      moptions['files_per_thread'] = 500

      #correctAndAnnotate_flowchart(moptions)
      correctAndAnnotate_manager(moptions)

      #'python myGetSignalToEvents.py ctrl_oligo_SpeI_cut/workspace/ seqchecklog oligo ref/trainseq3.fa'
   elif False:
      control = {'vector':'vector_SpeI_cut', 'oligo':'ctrl_oligo_SpeI_cut'}
      cases = {'oligo':'ctrl_oligo_SpeI_cut', \
               'Glucose':'Glucose_SpeI_cut', \
               'BrdU':'BrdU_SpeI_cut_data', \
               'A647':'A647_Click_SpeI_cut', \
               'Bio':'Bio_pur_SpeI_cut_data', \
               'Edu':'Edu_SpeI_cut_data', \
               'IdU':'IdU', \
               'A488':'A488_Click_SpeI_cut_data', \
               'Nicotine':'Nicotine_SpeI_cut', \
               'Azide':'Azide_SpeI_cut'}
      ckeys = cases.keys();
      for ck in ckeys:
         cmd = 'python myGetSignalToEvents.py %s/workspace/ seqchecklog %s' % (cases[ck], ck)
         cmdqsub = ('echo "%s" | qsub -V -cwd -N %s -e seqchecklog/%s.e -o seqchecklog/%s.o ' % (cmd, ck, ck, ck))
         print cmdqsub
         #os.system(cmdqsub)
   ##
   #python myGetFQAlign.py ctrl_oligo_SpeI_cut/workspace/ seqchecklog oligo ref/trainseq3.fa
