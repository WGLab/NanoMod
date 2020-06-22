#!/usr/bin/env python

import os;
import sys;

import string;

import argparse;
from argparse import RawTextHelpFormatter

from scripts import myDetect
from scripts.myCom import *

from scripts import mySimulate
from scripts import mySimulat2
from scripts import myDownSampling0
from scripts import myRefBaseSignalAnnotation


try:
   from rpy2.robjects.packages import importr
except:
   print("Error occurred to import 'rpy2'. \nPlease install 'rpy2' by using 'pip install rpy2'");
   sys.exit()

#detect simulate
parser = argparse.ArgumentParser(description="Detect nucleotide modification from nanopore signals data.", epilog="For example, \n \
\tpython %(prog)s detect: detect the modificaton of nucleotides\n \
\tpython %(prog)s simulate: simulate different percentage of modifications\n \
\tpython %(prog)s simulat2: simulate different number of subreads with modifications\n \
\tpython %(prog)s DownSampling: simulate downsampling reads from a folder\n \
\tpython %(prog)s Annotate: annotate a known sequence using fast5 \n \
", formatter_class=RawTextHelpFormatter);

def printParameters(moptions):
   mpkeys = moptions.keys(); mpkeys.sort()
   for mpk in mpkeys:
      print ('%30s: %s' % (mpk, str(moptions[mpk])))

def mCommonParam(margs):

   ErrorMessage = ""
   moptions = {}
   moptions['outLevel'] = margs.outLevel
   moptions["wrkBase1"] = margs.wrkBase1
   if moptions["wrkBase1"]==None or (not os.path.isdir(moptions["wrkBase1"])):
      ErrorMessage = ErrorMessage + ("\n\tThe directory of the first dataset (%s) does not exist" % moptions["wrkBase1"])
   #moptions["wrkBase2"] = margs.wrkBase2
   #if moptions["wrkBase2"]==None or ( not os.path.isdir(moptions["wrkBase2"])):
   #   ErrorMessage = ErrorMessage + ("\n\tThe directory of the second dataset (%s) does not exist" % moptions["wrkBase2"])
   moptions["window"] = (margs.window-1)/2
   if moptions["window"]<1:
      ErrorMessage = ErrorMessage + ("\n\tWindow size (%d) is too small" % (margs.window))
   moptions["FileID"] = margs.FileID

   moptions['outFolder'] = margs.outFolder
   if not os.path.isdir(moptions['outFolder']):
      os.system('mkdir -p '+moptions['outFolder'])
   #if moptions['outFolder']==None or (not os.path.isdir(moptions['outFolder'])):
   #   ErrorMessage = ErrorMessage + ("\n\tThe output folder (%s) does not exist" % moptions['outFolder'])
   #else:
   #   #if moptions['outFolder'][-1] not in ['/', '\\']: moptions['outFolder'] = moptions['outFolder'] + '/'
   #   moptions['outFolder'] = format_last_letter_of_folder(moptions['outFolder'])

   moptions["MinCoverage"] = margs.MinCoverage
   if moptions["MinCoverage"]<3:
      ErrorMessage = ErrorMessage + ("\n\tThe coverage (%d) is too small" % moptions["MinCoverage"])

   moptions["topN"] = margs.topN
   if moptions["topN"]<1:
      ErrorMessage = ErrorMessage + ("\n\tThe topN (%d) is too smaller" % moptions["topN"])

   moptions["neighborPvalues"] = margs.neighborPvalues
   if moptions["neighborPvalues"] < 0:
      ErrorMessage = ErrorMessage + ("\n\tThe neighborPvalues (%d) cannot be smaller than 0" % moptions["neighborPvalues"])
   moptions["WeightsDif"] = margs.WeightsDif
   if moptions["WeightsDif"]<1:
      moptions["WeightsDif"] = 1.0;

   moptions["testMethod"] = margs.testMethod

   moptions['rankUse'] = margs.rankUse

   moptions['SaveTest'] = margs.SaveTest

   moptions['.fast5'] = '.fast5'

   #moptions[''] = margs.
   moptions['RegionRankbyST'] = margs.RegionRankbyST
   moptions['percentile'] = margs.percentile
   if moptions['percentile']<0: moptions['percentile'] = 0.0;
   if moptions['percentile']>=1: moptions['percentile'] = 0.99
   moptions['WindOvlp'] = margs.WindOvlp
   moptions['NA'] = margs.NA
   #RegionRankbyST percentile WindOvlp NA

   return [moptions, ErrorMessage]

def testPackage():
   try:
      importr('ggplot2')
   except:
      print("Error occurred to import 'ggplot2'. \nPlease install 'ggplot2' in your R using 'install.packages(\"ggplot2\")'.");
      sys.exit()
   try:
      importr('gridExtra')
   except:
      print ("Error occurred to import 'gridExtra'. \nPlease install 'gridExtra' in your R using 'install.packages(\"gridExtra\")'.");
      sys.exit()


def mDetect(margs):
   testPackage()

   moptions, ErrorMessage = mCommonParam(margs)
   moptions['mstd'] = margs.mstd
   if not margs.Pos=="":
      mpos = margs.Pos.split(":")
      moptions["Chr"] = mpos[0]
      if len(mpos)>1:
         moptions["Pos"] = int(mpos[1])-1
         if moptions["Pos"]<0:
            ErrorMessage = ErrorMessage + ("\n\tThe position (%d) of interest should not be less than 0" % moptions["Pos"])
      if len(mpos)>2: ###########
         moptions["Pos2"] = int(mpos[2])-1
         if moptions["Pos2"]<0:
            ErrorMessage = ErrorMessage + ("\n\tThe position (%d) of interest should not be less than 0" % moptions["Pos2"])
         if moptions["Pos2"] - moptions["Pos"] < 1:
            ErrorMessage = ErrorMessage + ("\n\tThe end position (%d) is not larger than the start position (%d)" % moptions["Pos2"], moptions["Pos"])

   #moptions["MinCoverage"] = margs.MinCoverage
   #if moptions["MinCoverage"]<3:
   #   ErrorMessage = ErrorMessage + ("\n\tThe coverage (%d) is too small" % moptions["MinCoverage"])
   #moptions["topN"] = margs.topN
   #if moptions["topN"]<1:
   #   ErrorMessage = ErrorMessage + ("\n\tThe topN (%d) is too smaller" % moptions["topN"])
   #
   #moptions["neighborPvalues"] = margs.neighborPvalues
   #if moptions["neighborPvalues"] < 0:
   #   ErrorMessage = ErrorMessage + ("\n\tThe neighborPvalues (%d) cannot be smaller than 0" % moptions["neighborPvalues"])
   #
   #moptions["testMethod"] = margs.testMethod

   moptions["wrkBase2"] = margs.wrkBase2
   if moptions["wrkBase2"]==None or ( not os.path.isdir(moptions["wrkBase2"])):
      ErrorMessage = ErrorMessage + ("\n\tThe directory of the second dataset (%s) does not exist" % moptions["wrkBase2"])

   if not ErrorMessage=="":
      ErrorMessage = "Please provide correct parameters" + ErrorMessage
      print ErrorMessage, '\n'
      parser.print_help();
      parser.parse_args(['detect', '--help']);
      sys.exit(1)

   #moptions['checkN'] = 50

   moptions['wrkBase1'] = format_last_letter_of_folder(moptions['wrkBase1'])
   moptions['wrkBase2'] = format_last_letter_of_folder(moptions['wrkBase2'])
   moptions['plotType'] = margs.plotType
   moptions['downsampling'] = margs.downsampling
   moptions['downsampling_quantile'] = margs.downsampling_quantile
   moptions['min_lr'] = margs.min_lr
   moptions['min_lr_nb'] = margs.min_lr_nb

   #if moptions['downsampling']:
   #    moptions['Percentages'] = margs.Percentages.split(',')
   #    for ipind in range(len(moptions['Percentages'])):
   #        moptions['Percentages'][ipind] = float(moptions['Percentages'][ipind])
   #        if moptions['Percentages'][ipind]<10**-5:
   #            ErrorMessage = ErrorMessage + ("\n\tThe Percentage (%.6f) is too small or large than %.6f" % (moptions['Percentages'][ipind], 10**-5))
   #            if int(moptions['Percentages'][ipind]) >= 1:
   #                moptions['Percentages'][ipind] = 1.1

   moptions['coverages'] = map(int, margs.coverages.split('-'))
   if len(moptions['coverages'])==1:
      moptions['coverages'] = [moptions['coverages'][0],moptions['coverages'][0]]
   printParameters(moptions)
   myDetect.mDetect(moptions)


def mSimulate(margs):
   testPackage()

   moptions, ErrorMessage = mCommonParam(margs)

   moptions["wrkBase3"] = margs.wrkBase3
   if (not moptions["wrkBase3"]==None) and (not os.path.isdir(moptions["wrkBase3"])):
      ErrorMessage = ErrorMessage + ("\n\tThe directory of the third dataset (%s) does not exist" % moptions["wrkBase3"])
   moptions['Percentages'] = margs.Percentages.split(',')
   for ipind in range(len(moptions['Percentages'])):
      moptions['Percentages'][ipind] = float(moptions['Percentages'][ipind])
      if moptions['Percentages'][ipind]<10**-5:
         ErrorMessage = ErrorMessage + ("\n\tThe Percentage (%.6f) is too small or large than %.6f" % (moptions['Percentages'][ipind], 10**-5))
   moptions['Percentages'].sort()

   moptions["wrkBase2"] = margs.wrkBase2
   if moptions["wrkBase2"]==None or ( not os.path.isdir(moptions["wrkBase2"])):
      ErrorMessage = ErrorMessage + ("\n\tThe directory of the second dataset (%s) does not exist" % moptions["wrkBase2"])

   #ErrorMessage = ErrorMessage.replace(("\n\tThe directory of the first dataset (%s) does not exist" % moptions["wrkBase1"]), '')

   if not ErrorMessage=="":
      ErrorMessage = "Please provide correct parameters" + ErrorMessage
      print ErrorMessage, '\n'
      parser.print_help();
      parser.parse_args(['simulate', '--help']);
      sys.exit(2)

   moptions['run1'] = True;
   moptions["FileID_base"] = moptions["FileID"]

   moptions['wrkBase1'] = format_last_letter_of_folder(moptions['wrkBase1'])
   moptions['wrkBase2'] = format_last_letter_of_folder(moptions['wrkBase2'])
   if (not moptions["wrkBase3"]==None):
      moptions['wrkBase3'] = format_last_letter_of_folder(moptions['wrkBase3'])

   printParameters(moptions)
   mySimulate.mSimulate(moptions)


def mSimulat2(margs):
   testPackage()

   moptions, ErrorMessage = mCommonParam(margs)

   moptions['Percentage'] = margs.Percentage
   if (not margs.runType==3) and moptions['Percentage']==None:
      ErrorMessage = ErrorMessage + ("\n\tNeed a Percentage parameter")
   elif (not margs.runType==3) and moptions['Percentage']<10**-5:
      ErrorMessage = ErrorMessage + ("\n\tThe Percentage (%.6f) is too small or large than %.6f" % (moptions['Percentage'], 10**-5))

   moptions["wrkBase2"] = margs.wrkBase2
   if (not margs.runType==3) and (moptions["wrkBase2"]==None or ( not os.path.isdir(moptions["wrkBase2"]))):
      ErrorMessage = ErrorMessage + ("\n\tThe directory of the second dataset (%s) does not exist" % moptions["wrkBase2"])

   moptions["CaseSize"] = margs.CaseSize
   if margs.runType==2 and (moptions["CaseSize"]==None):
      ErrorMessage = ErrorMessage + ("\n\tNeed a CaseSize parameter")
   elif margs.runType==2 and moptions["CaseSize"]<50:
      ErrorMessage = ErrorMessage + ("\n\tThe number (%d) of subreads with modifications is too small" % moptions["CaseSize"])
   moptions["runType"] = margs.runType

   ErrorMessage = ErrorMessage.replace(("\n\tThe directory of the first dataset (%s) does not exist" % moptions["wrkBase1"]), '')

   if not ErrorMessage=="":
      ErrorMessage = "Please provide correct parameters" + ErrorMessage
      print ErrorMessage, '\n'
      parser.print_help();
      parser.parse_args(['simulat2', '--help']);
      sys.exit(2)

   moptions["FileID_base"] = moptions["FileID"]

   moptions['wrkBase1'] = format_last_letter_of_folder(moptions['wrkBase1'])
   moptions['wrkBase2'] = format_last_letter_of_folder(moptions['wrkBase2'])

   printParameters(moptions)
   mySimulat2.mSimulat2(moptions)

def mDownSampling(margs):
   testPackage()

   moptions, ErrorMessage = mCommonParam(margs)

   moptions["wrkBase2"] = margs.wrkBase2
   if (not margs.runType==3) and (moptions["wrkBase2"]==None or ( not os.path.isdir(moptions["wrkBase2"]))):
      ErrorMessage = ErrorMessage + ("\n\tThe directory of the second dataset (%s) does not exist" % moptions["wrkBase2"])

   moptions["CaseSize"] = margs.CaseSize
   if margs.runType==2 and (moptions["CaseSize"]==None):
      ErrorMessage = ErrorMessage + ("\n\tNeed a CaseSize parameter")
   elif margs.runType==2 and moptions["CaseSize"]<50:
      ErrorMessage = ErrorMessage + ("\n\tThe number (%d) of subreads with modifications is too small" % moptions["CaseSize"])
   moptions["runType"] = margs.runType
   moptions["mprefix"] = margs.mprefix

   ErrorMessage = ErrorMessage.replace(("\n\tThe directory of the first dataset (%s) does not exist" % moptions["wrkBase1"]), '')

   if not ErrorMessage=="":
      ErrorMessage = "Please provide correct parameters" + ErrorMessage
      print ErrorMessage, '\n'
      parser.print_help();
      parser.parse_args(['DownSampling', '--help']);
      sys.exit(2)

   moptions["FileID_base"] = moptions["FileID"]

   moptions['wrkBase1'] = format_last_letter_of_folder(moptions['wrkBase1'])
   moptions['wrkBase2'] = format_last_letter_of_folder(moptions['wrkBase2'])

   printParameters(moptions)
   myDownSampling0.mSimulat2(moptions)


def mAnnotate(margs):
   testPackage()

   ErrorMessage = ""
   moptions = {}
   moptions['outLevel'] = margs.outLevel
   moptions["wrkBase1"] = margs.wrkBase1
   if moptions["wrkBase1"]==None or (not os.path.isdir(moptions["wrkBase1"])):
      ErrorMessage = ErrorMessage + ("\n\tThe directory of the first dataset (%s) does not exist" % moptions["wrkBase1"])

   #Ref kmer_model_file Resegment_wind Resegment_signal_wind threads files_per_thread
   moptions['Ref'] = margs.Ref
   if moptions['Ref']==None or (not os.path.isfile(moptions['Ref'])):
     ErrorMessage = ErrorMessage + ("\n\tThe Ref file (%s) is NONE or does not exist" % moptions["Ref"])
   moptions['kmer_model_file'] = margs.kmer_model_file
   moptions['Resegment_wind'] = margs.Resegment_wind
   if moptions['Resegment_wind']<2:
      ErrorMessage = ErrorMessage + ("\n\t The window size (%d) of nucleotides for regsegmentation is too small" % moptions['Resegment_wind'])
   moptions['Resegment_signal_wind'] = margs.Resegment_signal_wind
   if moptions['Resegment_signal_wind']<2:
      ErrorMessage = ErrorMessage + ("\n\t The window size (%d) of raw signals for regsegmentation is too small" % moptions['Resegment_signal_wind'])
   moptions['MinNumSignal'] = margs.MinNumSignal
   if moptions['MinNumSignal']<2:
      ErrorMessage = ErrorMessage + ("\n\t The number of raw signals of an event for regsegmentation is too small" % moptions['MinNumSignal'])
   moptions['threads'] = margs.threads
   if moptions['threads']<1:
      ErrorMessage = ErrorMessage + ("\n\t The number (%d) of threads is too small" % moptions['threads'])
   moptions['files_per_thread'] = margs.files_per_thread
   if moptions['files_per_thread']<1:
      ErrorMessage = ErrorMessage + ("\n\t The number (%d) of fast5 files per thread is too small" % moptions['files_per_thread'])

   moptions['basecall_1d'] = margs.basecall_1d
   moptions['basecall_2strand'] = margs.basecall_2strand
   moptions['recursive'] = margs.recursive
   moptions['alignStr'] = margs.alignStr; #'bwa'
   #moptions['continue'] = margs.continue
   if not ErrorMessage=="":
      ErrorMessage = "Please provide correct parameters" + ErrorMessage
      print ErrorMessage, '\n'
      parser.print_help();
      parser.parse_args(['Annotate', '--help']);
      sys.exit(2)

   printParameters(moptions)
   myRefBaseSignalAnnotation.correctAndAnnotate_manager(moptions)


#####################################################################################

subparsers = parser.add_subparsers()
parent_parser = argparse.ArgumentParser(add_help=False)

com_group_for_comparison = parent_parser.add_argument_group('Common options for the comparison between two groups of signals.')
com_group_for_comparison.add_argument("--outLevel", type=int, choices=[OUTPUT_DEBUG, OUTPUT_INFO, OUTPUT_WARNING, OUTPUT_ERROR], default=OUTPUT_WARNING, help=("The level for output: %d for DEBUG, %d for INFO, %d for WARNING, %d for ERROR. Default: %d" % (OUTPUT_DEBUG, OUTPUT_INFO, OUTPUT_WARNING, OUTPUT_ERROR, OUTPUT_WARNING)))
com_group_for_comparison.add_argument("--wrkBase1", help="The working base folder for the first group.")
#com_group_for_comparison.add_argument("--wrkBase2", help="The working base folder for the second group.")
com_group_for_comparison.add_argument("--window", type=int, default=21, help="The window size to plot. It is better to be an odd number, such as 11, 21, 31, ... Default: 21")
com_group_for_comparison.add_argument("--FileID", default="mod", help="The unique string for intermediate files and final output files. Default: 'mod'")
com_group_for_comparison.add_argument("--outFolder", default=mresfolder_def, help="The default folder for outputing the results. Default: '"+mresfolder_def+"'")
com_group_for_comparison.add_argument("--MinCoverage", type=int, default=5, help="The minimum coverage for the base of interest and its neighbors. Default: 5.")
#com_group_for_comparison.add_argument("--MinCoverage", type=int, default=30, help="The minimum coverage for the base of interest and its neighbors. Default: 30.")
com_group_for_comparison.add_argument("--topN", type=int, default=30, help="The topN most significance test. Default: 30")
com_group_for_comparison.add_argument("--neighborPvalues", type=int, default=2, help="The number of neighbor p-values are used for combining p-values using fisher's or stouffer's methods. Default:2.")
com_group_for_comparison.add_argument("--WeightsDif", type=float, default=2.0, help="The difference of two adjacent neighbors when weight stouffer's methods are used. Maximum weight for the center position is 100, and the weights of a position are obtained by dividing the weights of position closest to the center by WeightsDif. Default:2.0.")
com_group_for_comparison.add_argument("--testMethod", default="stouffer", choices=['fisher', 'stouffer', 'ks'], help="Which method is used for test statistics: fisher, stouffer (default) or KS-test. ")
#com_group_for_comparison.add_argument("--rankUse", default='st', choices=['st', 'pv'], help="Which criterion is used for ranking: 'st':statistics(default); 'pv':p-value")
com_group_for_comparison.add_argument("--rankUse", default='pv', choices=['st', 'pv'], help="Which criterion is used for ranking: 'st':statistics; 'pv':p-value(default)")
com_group_for_comparison.add_argument("--SaveTest", type=int, default=1, choices=[0,1], help="Whether significant test would be save. Default: 0 (not save)")
com_group_for_comparison.add_argument("--RegionRankbyST", type=int, default=0, choices=[0,1], help="Rank region of window according to the stastitcs of the region."); #RegionRankbyST percentile WindOvlp NA
com_group_for_comparison.add_argument("--percentile", type=float, default=0.1, help="The smallest percentile of p-values in the region is used for ranking. Only used when RegionRankbyST=1.")
com_group_for_comparison.add_argument("--WindOvlp", type=int, default=0, choices=[0,1], help="Whether two windows are overlapped. Default: 0 (not save). The overlapping size will be half of window. Only used when RegionRankbyST=1.")
com_group_for_comparison.add_argument("--NA", type=str, default="", choices=["", 'A', 'C', 'G', 'T'], help="The nucleotide type of interest. Default: ''(all 4 nucleotides). Only used  when RegionRankbyST=1.")


# detect simulate
parser_detect = subparsers.add_parser('detect', parents=[parent_parser], help="Detect nucleotide modifications from nanopore signal data", description="Detect nucleotide modifications from nanopore signal data", epilog="For example, \n \
python %(prog)s --wrkBase1 ctrl_oligo_SpeI_cut --wrkBase2 Nicotine_SpeI_cut --Pos spel:3073 \n \
python %(prog)s --wrkBase1 ctrl_oligo_SpeI_cut/workspace/0 --wrkBase2 Nicotine_SpeI_cut/workspace/0 --Pos spel:3073 \n \
python %(prog)s --wrkBase1 vector_SpeI_cut/workspace/0 --wrkBase2 A488_Click_SpeI_cut_data/workspace/0 \n \
python %(prog)s --wrkBase1 ctrl_oligo_SpeI_cut --wrkBase2 Nicotine_SpeI_cut \n \
", formatter_class=RawTextHelpFormatter)

parser_detect.add_argument("--Pos", default="", help="The position of interest: chr:pos. Default: ''")
parser_detect.add_argument("--mstd", type=bool, default=False, help="Whether save mean and std. Default: False")
#parser_detect.add_argument("--MinCoverage", type=int, default=30, help="The minimum coverage for the base of interest and its neighbors. Default: 30.")
#parser_detect.add_argument("--topN", type=int, default=30, help="The topN most significance test. Default: 10")
#parser_detect.add_argument("--neighborPvalues", type=int, default=1, help="The number of neighbor p-values are used for combining p-values using fisher. Default:2.")
#parser_detect.add_argument("--testMethod", default="ks", choices=['fisher', 'stouffer', 'ks'], help="Which method is used for test statistics: fisher, stouffer or KS-test. Leave empty for KS-test.")

parser_detect.add_argument("--wrkBase2", help="The working base folder for the second group.")
parser_detect.add_argument("--plotType", type=str, default="Density", choices=["Violin", "Density"], help="The type of the plot.")

parser_detect.add_argument("--min_lr", default=500, type=int, help="The minimum length of long reads.")
parser_detect.add_argument("--min_lr_nb", default=0, type=int, help="The variance of requried length of long reads: [min_lr-min_lr_nb, [min_lr+min_lr_nb]. Default: 0(>min_lr).")
parser_detect.add_argument("--downsampling_quantile", default=0.25, type=float, help="The smallest quantile p-value.")
parser_detect.add_argument("--downsampling", default=100, type=int, help="How many times for  downsampling: >1.")
#parser_detect.add_argument("--Percentages", type=str, default='0.3', help="The percentages (seperated by ',', for exmaple '0.3,0.2,0.4,0.5,0.1') of reads with modifications. Default: '0.3'.")
parser_detect.add_argument("--coverages", type=str, default='0-0', help="The coverage threshold. DownSampling will be preformed if the covrage of a site is over the threshold. Default: '0-0'(no downsampling)")
parser_detect.set_defaults(func=mDetect)

parser_simulate = subparsers.add_parser('simulate', parents=[parent_parser], help="Simulate with different percentage of modifications", description="Simulation with different percentage of modifications", epilog="For example, \n \
python %(prog)s --wrkBase1 ctrl_oligo_SpeI_cut/workspace --wrkBase2 Nicotine_SpeI_cut/workspace --Percentages 0.15,0.2,0.25,0.3,0.4,0.5 \n \
python %(prog)s --wrkBase1 ctrl_oligo_SpeI_cut/workspace/0 --wrkBase2 Nicotine_SpeI_cut/workspace/0 --wrkBase3 ctrl_oligo_SpeI_cut/workspace/1  --Percentages 0.3 \n \
python %(prog)s --wrkBase1 ctrl_oligo_SpeI_cut/workspace --wrkBase2 Nicotine_SpeI_cut/workspace --Percentages 0.5,0.3 \n \
", formatter_class=RawTextHelpFormatter)

parser_simulate.add_argument("--wrkBase3", help="The working base folder for the second control group (--wrkBase1 is assumed to be the first control group).")
parser_simulate.add_argument("--Percentages", type=str, default='0.3', help="The percentages (seperated by ',', for exmaple '0.3,0.2,0.4,0.5,0.1') of reads with modifications. Default: '0.3'.")
parser_simulate.add_argument("--wrkBase2", help="The working base folder for the second group.")

parser_simulate.set_defaults(func=mSimulate)


parser_simulat2 = subparsers.add_parser('simulat2', parents=[parent_parser], help="Simulate with different percentage of modifications", description="Simulation with different percentage of modifications", epilog="For example, \n \
python %(prog)s --wrkBase1 ctrl_oligo_SpeI_cut/workspace --wrkBase2 Nicotine_SpeI_cut/workspace --Percentage 0.2 --runType 1 \n \
python %(prog)s --wrkBase1 ctrl_oligo_SpeI_cut/workspace --wrkBase2 Nicotine_SpeI_cut/workspace --Percentage 0.2 --runType 2 --CaseSize 3000 \n \
python %(prog)s --wrkBase1 ctrl_oligo_SpeI_cut/workspace --wrkBase2 Nicotine_SpeI_cut/workspace --Percentage 0.2 --runType 2 --CaseSize 500 \n \
python %(prog)s --outFolder logsim_step1000_sm2/ --runType 3 \n \
", formatter_class=RawTextHelpFormatter)

parser_simulat2.add_argument("--Percentage", type=float, default=None, help="The percentage (for example, 0.2) of reads with modifications. Default: 0.2.")
parser_simulat2.add_argument("--wrkBase2", help="The working base folder for the second group.")
parser_simulat2.add_argument("--CaseSize", type=int, default=None, help="The number of subreads with modifications. Default: 2000")
parser_simulat2.add_argument("--runType", type=int, default=2, choices=[1,2,3], help="The sub-functions to run simulation. 1:run with automatical determination of CaseSize (start from 1000 with step of 1000); 2:run with a specific CaseSize; 3: summarize the results.")

parser_simulat2.set_defaults(func=mSimulat2)

parser_DownSampling = subparsers.add_parser('DownSampling', parents=[parent_parser], help="Simulate with DownSampling from a folder", description="Simulation with DownSampling from a folder", epilog="For example, \n \
python %(prog)s --wrkBase1 ctrl_oligo_SpeI_cut/workspace --wrkBase2 Nicotine_SpeI_cut/workspace --outFolder logwind --runType 1 \n \
python %(prog)s --wrkBase1 ctrl_oligo_SpeI_cut/workspace --wrkBase2 Nicotine_SpeI_cut/workspace --outFolder logwind --runType 2 --CaseSize 100 \n \
python %(prog)s --outFolder logwind/ --runType 3 \n \
python %(prog)s --outFolder logwind0/ --runType 3 \n \
python %(prog)s --outFolder logwind_nocheckCov --runType 3 \n \
", formatter_class=RawTextHelpFormatter)

parser_DownSampling.add_argument("--wrkBase2", help="The working base folder for the second group.")
parser_DownSampling.add_argument("--CaseSize", type=int, default=None, help="The number of subreads with modifications. Default: 2000")
parser_DownSampling.add_argument("--runType", type=int, default=2, choices=[1,2,3], help="The sub-functions to run simulation. 1:run with a set of CaseSize; 2:run with a specific CaseSize; 3: summarize the results.")
parser_DownSampling.add_argument("--mprefix", type=str, default='',  help="The prefix for runType 3.")

parser_DownSampling.set_defaults(func=mDownSampling)


parser_annotate = subparsers.add_parser('Annotate', help="Annotate a known sequence using each fast5", description="Annotate a known sequence using each fast5", epilog="For example, \n \
python %(prog)s --wrkBase1 A488_Click_SpeI_cut_data/workspace --Ref ref/trainseq3.fa --kmer_model_file scripts/kmer_model/r9.4_450bps.nucleotide.5mer.template.model \n \
python %(prog)s --wrkBase1 A488_Click_SpeI_cut_data/workspace --Ref ref/trainseq3.fa --threads 24\n \
", formatter_class=RawTextHelpFormatter)

parser_annotate.add_argument("--outLevel", type=int, choices=[OUTPUT_DEBUG, OUTPUT_INFO, OUTPUT_WARNING, OUTPUT_ERROR], default=OUTPUT_WARNING, help=("The level for output: %d for DEBUG, %d for INFO, %d for WARNING, %d for ERROR. Default: %d" % (OUTPUT_DEBUG, OUTPUT_INFO, OUTPUT_WARNING, OUTPUT_ERROR, OUTPUT_WARNING)))
parser_annotate.add_argument("--wrkBase1", help="The working base folder for the first group.")
parser_annotate.add_argument("--Ref", help="The known sequence")
parser_annotate.add_argument("--kmer_model_file", help="The kmer_model_file")
parser_annotate.add_argument("--Resegment_wind", type=int, default=4, help="The number of neighbor nucleotides used for resegmentation. Default:4")
#parser_annotate.add_argument("--Resegment_signal_wind", type=int, default=3, help="The window size of raw signals used to calculate the difference of consecutive regions. Default:3")
parser_annotate.add_argument("--Resegment_signal_wind", type=int, default=4, help="The window size of raw signals used to calculate the difference of consecutive regions. Default:4")
parser_annotate.add_argument("--threads", type=int, default=12, help="The number of threads used. Default:12")
parser_annotate.add_argument("--files_per_thread", type=int, default=300, help="The number of fast5 files for each thread. Default:300")
#Ref kmer_model_file Resegment_wind Resegment_signal_wind threads files_per_thread
parser_annotate.add_argument("--basecall_1d", default="Basecall_1D_000", help="Path for basecall_1d. Default: Basecall_1D_000")
parser_annotate.add_argument("--basecall_2strand", default="BaseCalled_template", help="Path for basecall_2strand. Default: BaseCalled_template")
parser_annotate.add_argument("--MinNumSignal", type=int, default=4, help="Mininum number of signals for an event. Default:4")
parser_annotate.add_argument("--recursive", type=int, default=1, choices=[0,1], help="Recurise to find fast5 files. Default:0")
parser_annotate.add_argument("--alignStr", type=str, default='bwa', choices=["bwa","minimap2"], help="Alignment tools (bwa or minimap2 is supported). Default: bwa")
#parser_annotate.add_argument("--continue", default=False, action="store_true", type=str, default='bwa',  help="Continue performing annotation")

parser_annotate.set_defaults(func=mAnnotate);



if len(sys.argv)<2:
   parser.print_help();
else:
   args = parser.parse_args()
   args.func(args);
