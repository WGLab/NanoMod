

mresfolder_def = 'mRes/'

OUTPUT_DEBUG=0;
OUTPUT_INFO=1;
OUTPUT_WARNING=2;
OUTPUT_ERROR=3;

na_bp = {"A":"T", \
         "C":"G", \
         "G":"C", \
         "T":"A", \
         "a":"t", \
         "c":"g", \
         "g":"c", \
         "t":"a", \
         "N":"N", \
         "n":"n" \
         }

acgt = ('A', 'C', 'G', 'T', 'a', 'c', 'g', 't')
acgt = ('A', 'C', 'G', 'T', 'a', 'c', 'g', 't', 'N', 'n')

def getComplementary(na):
   com_na = []
   for i in range(len(na)):
      com_na.append(na_bp[na[i]])
   ''.join(com_na[::-1])

def format_last_letter_of_folder(cursub):
   if not cursub==None:
      if cursub[-1]=='/': return cursub;
      elif cursub[-1]=='\\': return cursub[:-1]+'/';
      else: return cursub+'/';

analyses_base = "Analyses"
basecall_1D_base = "Basecall_1D_000"
basecall_template_base = "BaseCalled_template"
basecall_events_base = "Events"
basecall_fastq_base = "Fastq"
basecall_move_base = "Move"

raw_base = 'Raw'
reads_base = "Reads"
signal_base = "Signal"

fast5_mo_correction = "NanomoCorrected_000"

rawGenomeCorrected_base = "RawGenomeCorrected_000"
################# for nanomod
rawGenomeCorrected_base = fast5_mo_correction
rawBaseCalled_template_base = "BaseCalled_template"
rawAlignment_base = "Alignment"
rawRead_segments_base = "read_segments"

raw_event_base = "Events"
rawGenome_alignment_base = "genome_alignment"
rawRead_alignment_base = "read_alignment"

map_chr_str = "mapped_chrom"
map_start_str = "mapped_start"
map_strand_str = "mapped_strand"
