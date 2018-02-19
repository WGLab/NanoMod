The inputs of NanoMod are two groups of FAST5 files and a known sequence. One group of FAST5 files are from a tested sample with modifications, while the other group of FAST5 files are from a control sample without modifications. FAST5 files are the output of Nanopore sequencing equipment. 

NanoMod needs `Events` data in FAST5 files. It is no `Event` data in FAST5 files, one need to run Albacore (1.0<version<v2.0 is strongly recommended) on their data before running NanoMod. Please refer to Nanopore community for help about how to download/install/run Albacore.

# The step of NanoMod

## 1. `python NanoMod.py Annotate`

First of all, run `python NanoMod.py Annotate` to obtain signal annotation information for the known sequence. The parameters for this command is described below. One usually need to provide where is the FAST5 files (specified by `--wrkBase1`) and a reference sequence(specified by `--Ref`), and keep other parameters default. If you need less or more threads to run this command, `--threads` should be changed. By default 12 threads will be used. 

```
usage: NanoMod.py [-h] {detect,Annotate} ...

Detect nucleotide modification from nanopore signals data.

optional arguments:
  -h, --help            show this help message and exit

For example,
        python NanoMod.py detect: detect the modificaton of nucleotides
        python NanoMod.py Annotate: annotate a known sequence using fast5

usage: NanoMod.py Annotate [-h] [--outLevel {0,1,2,3}] [--wrkBase1 WRKBASE1]
                           [--Ref REF] [--kmer_model_file KMER_MODEL_FILE]
                           [--Resegment_wind RESEGMENT_WIND]
                           [--Resegment_signal_wind RESEGMENT_SIGNAL_WIND]
                           [--threads THREADS]
                           [--files_per_thread FILES_PER_THREAD]
                           [--basecall_1d BASECALL_1D]
                           [--basecall_2strand BASECALL_2STRAND]
                           [--MinNumSignal MINNUMSIGNAL] [--recursive {0,1}]

Annotate a known sequence using each fast5

optional arguments:
  -h, --help            show this help message and exit
  --outLevel {0,1,2,3}  The level for output: 0 for DEBUG, 1 for INFO, 2 for WARNING, 3 for ERROR. Default: 2
  --wrkBase1 WRKBASE1   The working base folder for the first group.
  --Ref REF             The known sequence
  --kmer_model_file KMER_MODEL_FILE
                        The kmer_model_file
  --Resegment_wind RESEGMENT_WIND
                        The number of neighbor nucleotides used for resegmentation. Default:4
  --Resegment_signal_wind RESEGMENT_SIGNAL_WIND
                        The window size of raw signals used to calculate the difference of consecutive regions. Default:4
  --threads THREADS     The number of threads used. Default:12
  --files_per_thread FILES_PER_THREAD
                        The number of fast5 files for each thread. Default:300
  --basecall_1d BASECALL_1D
                        Path for basecall_1d. Default: Basecall_1D_000
  --basecall_2strand BASECALL_2STRAND
                        Path for basecall_2strand. Default: BaseCalled_template
  --MinNumSignal MINNUMSIGNAL
                        Mininum number of signals for an event. Default:4
  --recursive {0,1}     Recurise to find fast5 files. Default:1

For example,
 python NanoMod.py Annotate --wrkBase1 A488_Click_SpeI_cut_data/workspace --Ref ref/trainseq3.fa
 python NanoMod.py Annotate --wrkBase1 A488_Click_SpeI_cut_data/workspace --Ref ref/trainseq3.fa --threads 24
```
### Example
```
For example,
    python NanoMod.py detect --wrkBase1 ctrl_oligo_SpeI_cut --wrkBase2 Nicotine_SpeI_cut
```

### Output:
The annotation information will be writen back to each FAST5 file. One can use `h5dump FAST5file` or `h5ls -r FAST5file` to see the detail. 


## 2. `python NanoMod.py detect`
The second step of modification detection is to run `python NanoMod.py detect`. Similar to running `python NanoMod.py Annotate`, one can use the default values for majority of parameters. Only two parameters are necessary: ` --wrkBase1` and `--wrkBase2 ` which specify where are the FAST5 files for the tested sample and the control sample. One can also use `--testMethod` to tell NanoMod which statistical method should be used. By deault, `Kolmogorov-Smirnov test` is used. `--FileID` should be given to different values if you run NanoMod multiple times at the same time. By default, it is 'mod'


```
usage: NanoMod.py [-h] {detect,Annotate} ...

Detect nucleotide modification from nanopore signals data.

optional arguments:
  -h, --help            show this help message and exit

For example,
        python NanoMod.py detect: detect the modificaton of nucleotides
        python NanoMod.py Annotate: annotate a known sequence using fast5

usage: NanoMod.py detect [-h] [--outLevel {0,1,2,3}] [--wrkBase1 WRKBASE1]
                         [--window WINDOW] [--FileID FILEID]
                         [--outFolder OUTFOLDER] [--MinCoverage MINCOVERAGE]
                         [--topN TOPN] [--neighborPvalues NEIGHBORPVALUES]
                         [--WeightsDif WEIGHTSDIF]
                         [--testMethod {fisher,stouffer,ks}]
                         [--rankUse {st,pv}] [--SaveTest {0,1}]
                         [--RegionRankbyST {0,1}] [--percentile PERCENTILE]
                         [--WindOvlp {0,1}] [--NA {,A,C,G,T}] [--Pos POS]
                         [--wrkBase2 WRKBASE2]

Detect nucleotide modifications from nanopore signal data

optional arguments:
  -h, --help            show this help message and exit
  --Pos POS             The position of interest: chr:pos. Default: ''
  --wrkBase2 WRKBASE2   The working base folder for the second group.

Common options for the comparison between two groups of signals.:
  --outLevel {0,1,2,3}  The level for output: 0 for DEBUG, 1 for INFO, 2 for WARNING, 3 for ERROR. Default: 2
  --wrkBase1 WRKBASE1   The working base folder for the first group.
  --window WINDOW       The window size to plot. It is better to be an odd number, such as 11, 21, 31, ... Default: 21
  --FileID FILEID       The unique string for intermediate files and final output files. Default: 'mod'
  --outFolder OUTFOLDER
                        The default folder for outputing the results. Default: 'mRes/'
  --MinCoverage MINCOVERAGE
                        The minimum coverage for the base of interest and its neighbors. Default: 5.
  --topN TOPN           The topN most significance test. Default: 30
  --neighborPvalues NEIGHBORPVALUES
                        The number of neighbor p-values are used for combining p-values using fisher's or stouffer's methods. Default:2.
  --WeightsDif WEIGHTSDIF
                        The difference of two adjacent neighbors when weight stouffer's methods are used. Maximum weight for the center position is 100, and the weights of a position are obtained by dividing the weights of position closest to the center by WeightsDif. Default:2.0.
  --testMethod {fisher,stouffer,ks}
                        Which method is used for test statistics: fisher, stouffer or KS-test(default).
  --rankUse {st,pv}     Which criterion is used for ranking: 'st':statistics; 'pv':p-value(default)
  --SaveTest {0,1}      Whether significant test would be save. Default: 0 (not save)
  --RegionRankbyST {0,1}
                        Rank region of window according to the stastitcs of the region.
  --percentile PERCENTILE
                        The smallest percentile of p-values in the region is used for ranking. Only used when RegionRankbyST=1.
  --WindOvlp {0,1}      Whether two windows are overlapped. Default: 0 (not save). The overlapping size will be half of window. Only used when RegionRankbyST=1.
  --NA {,A,C,G,T}       The nucleotide type of interest. Default: ''(all 4 nucleotides). Only used  when RegionRankbyST=1.
```

### Example
```
For example,
    python NanoMod.py detect --wrkBase1 ctrl_oligo_SpeI_cut --wrkBase2 Nicotine_SpeI_cut
```

### Where to find the result:
The final results could be found in `--outFolder`/*`--FileID`*.txt/pdf. 




