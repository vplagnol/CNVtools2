CONTENT of annotation_files/
############################
- candidate_outliers_sample_qc_vanguard_explained.txt
    Samples that failed qc explained
- CNV_annotation.txt
    CNV to probe annotation mapping: each CNV is associated with
    targetting probes via ProbeIndex.  This file will soon change to
    reflect the probe type on the column headers. The file will also
    contain the T1D regions as well as Novel insertions and INS VNTR
    mappings.
- probe_annotation_all_mappings_v2.txt
    Probe annotation file, containing all possible mappings of each
    probe to a genomic location.

- vanguard_sample_map.txt
    File mapping the sample name taken from PEDFILE to the 
For these files: 
 - immunochip-exclusions-2011-05-26.tab
 - ogt-vanguard-sample-subject-lookup-2011-05-26.tab
 - ogt-vanguard-slide-sample-lookup-2011-05-26.tab
 - Vanguard.987Samples.Gender.CR.DNAtype.txt
 - Vanguard.987Samples.PEDFILE.txt
see file: readme_sample_info.txt


CONTENT of log_ratio/
#####################


RData files split by chromosome

RData files to change in the future to keep track of probe type and 3
different types of probe summary(mean, median,pca)
log ratio is the log_10 ratio, will change to log_2 in future releases of the pipeline

Each .Rdata contains two variables: 

1) "probe_dat_sub"    
   Size:         NUM_PROBESx(NUMSAMPLES+3)
   Column IDs:   Chr StartCoord ProbeIndex   SAMPLE1 SAMPLE2 ... SAMPLE_M (M = NUMSAMPLES)
   Row:          chr, start,ProbeIndex, normalised intensities for samples 1 to M


2) "summary_data_std"
   Size:         1xNUMSAMPLES
   Column IDs:   SAMPLE1 SAMPLE2 ... SAMPLE_M  
   Row:          PCA summarised intensities for samples 1 to M    



Each chrom directory has a table summarising PCA summaries for each CNV in the chromosome: 
"all_per_CNV_summarised_data_chrZ.txt"
