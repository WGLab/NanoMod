# NanoMod: a computational tool to detect DNA modifications using Nanopore long-read sequencing data

NanoMod is a novel computational tool for the detection of DNA modifications using Nanopore long-read sequencing data. The evaluation on simulation data with different types of modifications and on a methylation data of E. coli suggested that NanoMod achieved better performance than other existing tools in detecting modifications without training data. 

## Methodology of NanoMod

NanoMod was designed for the detection of de novo DNA modifications (for example, synthetically introduced modifications). The inputs of NanoMod were a group of reads from a DNA sample with modification at specific bases and a group of reads from the matched non-modified sample. The nucleotide sequences for tested samples are assumed to be known, that is, the reference genome must be already known a priori. Currently, within NanoMod, we used albacore for basecalling, and then perform an indel error correction by aligning electric signals to a reference genome, similar to the procedure implemented in nanoraw. After that, two groups of electric signals for each genomic position were tested using Kolmogorov-Smirnov test in a per-base level to identify bases with significantly different distributions of signals between the two groups. Finally, weighted Stouffer’s method was used to combine the effects of neighboring bases since some modifications (especially bulky ones) may have strong neighbor effects that affect electric signals in neighboring non-modified bases.

## Inputs of NanoMod

The input of NanoMod is a dataset with two groups of reads: one from a sample with DNA modifications at specific positions and the other is the matched non-modified sample. A known sequence (or de novo assembly results) would be needed for indel correction.

## Usage

Please refer to [Usage](https://github.com/WGLab/NanoMod/blob/master/docs/Usage.md) for how to use NanoMod.

## Revision History

For release history, please visit [here](https://github.com/WGLab/NanoMod/releases). For details, please go [here](https://github.com/WGLab/NanoMod/blob/master/README.md).

## Contact

If you have any questions/issues/bugs, please post them on [GitHub](https://github.com/WGLab/NanoMod/issues). They would also be helpful to other users. 

## Reference

Liu Qian, Daniela C. Georgieva, Dietrich M. Egli, Kai Wang. NanoMod: a computational tool to detect DNA modifications using Nanopore long-read sequencing data.
