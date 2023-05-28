# PooledPSMC
Auto-PSMC with pooled individuals

## usage

### input files
1. A bam_list file that contain all the bam file path to the individuals. For example:

```
/home/user1/project1/indiv1.bam
/home/user1/project1/indiv2.bam
/home/user1/project1/indiv3.bam
/home/user1/project1/indiv4.bam
```

2. A project name list with all project/individual names corresponding to the bam file. For example:

```
indiv1
indiv2
indiv3
indiv4
```

3. Reference genome

### to run:

simply call python PooledPSMC.py (or python PooledPSMC.py -h to see the help)

Make sure you have 
```
samtools
bcftools
fq2psmcfa
psmc
splitfa
vcfutils.pl
psmc_plot.pl
bamtools
pandas
matplotlib
numpy
```

in your environment.
