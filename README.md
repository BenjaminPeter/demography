# Simulator of modern human and Neanderthal demography

*Dependencies*: `numpy`, `pandas`, `msprime`, `jupyter`.

```
$ ./demography.py  --help
usage: demography.py [-h] [--eurA start,duration,rate1,rate2]
                     [--eurB start,duration,rate1,rate2]
                     [--afr start,duration,rate1,rate2]
                     [--neand start,duration,rate1,rate2]
                     [--Ne-nonafr NE_NONAFR] [--Ne-afr NE_AFR] --seq-len
                     SEQ_LEN
                     (--eur-ages EUR_AGES [EUR_AGES ...] | --neur NEUR)
                     [--nasn NASN] [--nafrA NAFRA] [--nafrB NAFRB]
                     [--introgression [from-to [from-to ...]]] [--snps]
                     [--writesnps]
                     [--stats {true_neand,asc_neand,indirect,direct,afr_f4} [{true_neand,asc_neand,indirect,direct,afr_f4} ...]]
                     [--output-prefix FILE] [--debug] [--chrom CHROM]
                     [--rec REC]

optional arguments:
  -h, --help            show this help message and exit
  --eurA start,duration,rate1,rate2
                        EUR <-> AFR_A geneflow
  --eurB start,duration,rate1,rate2
                        EUR <-> AFR_B geneflow
  --afr start,duration,rate1,rate2
                        AFR_A <-> AFR_B geneflow
  --neand start,duration,rate1,rate2
                        nonAfr <-> Neanderthal geneflow
  --Ne-nonafr NE_NONAFR
  --Ne-afr NE_AFR
  --seq-len SEQ_LEN     Sequence length
  --eur-ages EUR_AGES [EUR_AGES ...]
                        Ages of European samples [years BP]
  --neur NEUR           Number of European chromosomes to simulate
  --nasn NASN           Number of Asian chromosomes to simulate
  --nafrA NAFRA         Number of AfrA chromosomes to simulate
  --nafrB NAFRB         Number of AfrB chromosomes to simulate
  --introgression [from-to [from-to ...]]
                        Pairs of populations for introgression detection
  --snps                Save all SNPs to a file
  --writesnps           write all SNPs to a file in a memory efficient manner.
                        this will disable generation of the all_snps table
                        used for e.g. stats calculations
  --stats {true_neand,asc_neand,indirect,direct,afr_f4} [{true_neand,asc_neand,indirect,direct,afr_f4} ...]
                        Which statistics to calculate?
  --output-prefix FILE  Prefix of output files
  --debug               Debug info
  --chrom CHROM         chromosome id for admixfrog output
  --rec REC             file with recombination map
