
import os;
import sys;

import string;
import math;

import h5py
import numpy as np

from pkg_resources import resource_string

from myCom import *

#/                        Group
#/Analyses                Group
#/Analyses/Basecall_1D_000 Group
#/Analyses/Basecall_1D_000/BaseCalled_template Group
#/Analyses/Basecall_1D_000/BaseCalled_template/Events Dataset {755}
#/Analyses/Basecall_1D_000/BaseCalled_template/Fastq Dataset {SCALAR}
#/Analyses/Basecall_1D_000/Configuration Group
#/Analyses/Basecall_1D_000/Configuration/basecall_1d Group
#/Analyses/Basecall_1D_000/Summary Group
#/Analyses/Basecall_1D_000/Summary/basecall_1d_template Group
#/Analyses/RawGenomeCorrected_000 Group
#/Analyses/RawGenomeCorrected_000/BaseCalled_template Group
#/Analyses/RawGenomeCorrected_000/BaseCalled_template/Alignment Group
#/Analyses/RawGenomeCorrected_000/BaseCalled_template/Alignment/genome_alignment Dataset {419}
#/Analyses/RawGenomeCorrected_000/BaseCalled_template/Alignment/read_alignment Dataset {419}
#/Analyses/RawGenomeCorrected_000/BaseCalled_template/Alignment/read_segments Dataset {369}
#/Analyses/RawGenomeCorrected_000/BaseCalled_template/Events Dataset {417}
#/Analyses/Segmentation_000 Group
#/Analyses/Segmentation_000/Configuration Group
#/Analyses/Segmentation_000/Configuration/event_detection Group
#/Analyses/Segmentation_000/Configuration/stall_removal Group
#/Analyses/Segmentation_000/Summary Group
#/Analyses/Segmentation_000/Summary/segmentation Group
#/Raw                     Group
#/Raw/Reads               Group
#/Raw/Reads/Read_7        Group
#/Raw/Reads/Read_7/Signal Dataset {4558/Inf}


#analyses_base = "Analyses"
#basecall_1D_base = "Basecall_1D_000"
#basecall_template_base = "BaseCalled_template"
#basecall_events_base = "Events"
#basecall_fastq_base = "Fastq"

basecall_event_full = ("/%s/%s/%s/%s" % (analyses_base,basecall_1D_base,basecall_template_base,basecall_events_base))
basecall_fastq_full = ("/%s/%s/%s/%s" % (analyses_base,basecall_1D_base,basecall_template_base,basecall_fastq_base))
basecall_move_full = ("/%s/%s/%s/%s" % (analyses_base,basecall_1D_base,basecall_template_base,basecall_move_base))

def ReadEvents(mf5):
   if mf5.__contains__(basecall_event_full):
      return mf5[analyses_base][basecall_1D_base][basecall_template_base][basecall_events_base]

def ReadMove(mf5):
   if mf5.__contains__(basecall_move_full):
      return mf5[analyses_base][basecall_1D_base][basecall_template_base][basecall_move_base]

def ReadFastQ(mf5):
   if mf5.__contains__(basecall_fastq_full):
      return mf5[analyses_base][basecall_1D_base][basecall_template_base][basecall_fastq_base][()].split("\n")

#raw_base = 'Raw'
#reads_base = "Reads"
#signal_base = "Signal"

reads_full = ("/%s/%s" % (raw_base, reads_base))

def ReadRawSignal(mf5):
   if mf5.__contains__(reads_full):
      readKeys = mf5.__getitem__(reads_full).keys();
      if len(readKeys)>0:
         if len(readKeys)>1: print "Warning!!! more than 1 read keys", readKeys
         if mf5.__contains__("/%s/%s/%s/%s" % (raw_base, reads_base, readKeys[0], signal_base)):
            return mf5[raw_base][reads_base][readKeys[0]][signal_base]

#rawGenomeCorrected_base = "RawGenomeCorrected_000"
#rawBaseCalled_template_base = "BaseCalled_template"
#rawAlignment_base = "Alignment"
#rawRead_segments_base = "read_segments"

rawReads_segments_full = ("/%s/%s/%s/%s/%s" % (analyses_base, rawGenomeCorrected_base, rawBaseCalled_template_base, rawAlignment_base, rawRead_segments_base))

def ReadNanoraw_SingalIndexForBase(mf5):
   if mf5.__contains__(rawReads_segments_full): #
      return mf5[rawReads_segments_full]

#raw_event_base = "Events"
raw_event_ful  = ("/%s/%s/%s/%s" % (analyses_base, rawGenomeCorrected_base, rawBaseCalled_template_base, raw_event_base))

def ReadNanoraw_events(mf5):
   #print raw_event_ful
   if mf5.__contains__(raw_event_ful):
      return mf5[raw_event_ful].value

#rawGenome_alignment_base = "genome_alignment"
rawGenome_alignment_full = ("/%s/%s/%s/%s/%s" % (analyses_base, rawGenomeCorrected_base, rawBaseCalled_template_base, rawAlignment_base, rawGenome_alignment_base))

def ReadNanoraw_Genome_alignment(mf5):
   if mf5.__contains__(rawGenome_alignment_full):
      return ''.join(mf5[rawGenome_alignment_full])

#rawRead_alignment_base = "read_alignment"
rawRead_alignment_full = ("/%s/%s/%s/%s/%s" % (analyses_base, rawGenomeCorrected_base, rawBaseCalled_template_base, rawAlignment_base, rawRead_alignment_base))

def ReadNanoraw_read_alignment(mf5):
   if mf5.__contains__(rawRead_alignment_full):
      return ''.join(mf5[rawRead_alignment_full])

rawAlignment_full = ("/%s/%s/%s/%s" % (analyses_base, rawGenomeCorrected_base, rawBaseCalled_template_base, rawAlignment_base))

#map_chr_str = "mapped_chrom"
#map_start_str = "mapped_start"
#map_strand_str = "mapped_strand"

def ReadMapInfoInRef(mf5):
   align_info = dict(mf5[rawAlignment_full].attrs.items())

   mapped_chrom = align_info[map_chr_str]
   mapped_start = align_info[map_start_str]
   mapped_strand= align_info[map_strand_str]

   return [mapped_chrom, int(mapped_start), mapped_strand]
