Vanguard.987Samples.Gender.CR.DNAtype.txt
-----------------------------------------
Ped-file for the 987 Vanguard sample

File size = 987
Columns:  
            Sample.Name  = Sample ID
            DNA          = DNA type, Cell-line or Cell-Pack
            sex          = Gender, 1 for MALE, 2 for FEMALE
            CallRate     = Call rate on Immunochip

Example  
Sample.Name       DNA         sex  CallRate
10103601          Cell-Line   1    0.9851621
10103602          Cell-Line   2    0.9780638
10103603          Cell-Line   2    0.9788219
10103604          Cell-Line   1    0.9851367
10143901          Cell-Line   1    0.9853300
10143902          Cell-Line   2    0.9798142

Vanguard.987Samples.PEDFILE.txt
-------------------------------

File size = 987
Columns:   
           FamID        = Family ID     
           Sample.Name  = Sample ID
           fid          = Father ID
           mid          = Mother ID
           sex          = Gender, 1 for MALE, 2 for FEMALE
           t1d          = Affected, 0 for missing, 1 for 
                          unaffected, 2 for affected

FamID     Sample.Name      fid      mid       sex t1d
137260    13726002         0        0          2   0
143425    14342503         14342501 14342502   1   2
163303    16330303         16330301 16330302   2   2
163303    16330304         16330301 16330302   2   2
167880    16788007         16788001 16788002   1   2
285806    28580603         28580601 28580602   2   1


immunochip-exclusions-2011-05-26.tab
-------------------------------------
This file contains a lit of exclusions that UVA called in IMmunochip,
and include some of relevance in Vanguard data

File size = 7
family X.subjectid immunochip_exclusion_reason relevant_to_ogt
143425    14342503          MISS_RATE_5PERCENT              no
167880    16788007         SEXQC_MISCLASS_MALE             yes
411921    41192101                   DUPFAMILY             yes
411921    41192102                   DUPFAMILY             yes
411921    41192103                   DUPFAMILY             yes
411921    41192104                   DUPFAMILY             yes
167880    16788001          MISMATCH_PREV_DATA             yes


ogt-vanguard-sample-subject-lookup-2011-05-26.tab
-------------------------------------------------
File that links these samples to subjects, and their phenotypes.  
This includes one edit proposed in the Immunochip by UVA - the swapping 
of 40717103 and 40717104, which has no great distinction until one 
analyses age-at-onset, as both are female cases with the same parents;

File size = 988

#sampleid  #subjectid familyid member father mother    sex t1d onset  collection
10103601    10103601   101036      1      0      0   male  no    -1  T1DGC-AP
10103602    10103602   101036      2      0      0 female  no    -1  T1DGC-AP
10103603    10103603   101036      3      1      2 female yes    11  T1DGC-AP
10103604    10103604   101036      4      1      2   male yes     1  T1DGC-AP
10143901    10143901   101439      1      0      0   male  no    -1  T1DGC-AP
10143902    10143902   101439      2      0      0 female  no    -1  T1DGC-AP
  

ogt-vanguard-slide-sample-lookup-2011-05-26.tab
-----------------------------------------------
File that links the samples to the slides (and from which file 
names could be constructed automatically, give or take the issue of 
whether these are first or second scans);

File size = 1,008

#sampleid Slide_Barcode Location_on_Slide
23553803  253302210104               1_1
23553804  253302210104               1_2
23553801  253302210104               1_3
23553802  253302210104               1_4
22740202  253302210104               2_1
22740203  253302210104               2_2


