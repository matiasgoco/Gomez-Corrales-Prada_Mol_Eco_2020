Reanalysis of Dziedzic et al. (2019) O. faveolata data (<a href="https://www.ncbi.nlm.nih.gov//bioproject/PRJNA413258" class="uri">https://www.ncbi.nlm.nih.gov//bioproject/PRJNA413258</a>)
--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

Create and activate a conda environment
---------------------------------------

    conda create -n Orbicella
    conda activate Orbicella

Downloading files from SRA
--------------------------

First we install the sratoolkit using conda
-------------------------------------------

    conda install -c bioconda sra-tools 
    conda install -c bioconda bioconductor-sradb 
    conda install -c bioconda entrez-direct 
    conda install -c conda-forge parallel 

NCBI link to this BioProject <a href="https://www.ncbi.nlm.nih.gov/biosample?LinkName=bioproject_biosample_all&amp;from_uid=397196" class="uri">https://www.ncbi.nlm.nih.gov/biosample?LinkName=bioproject_biosample_all&amp;from_uid=397196</a>
------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

Downloading SRA sequences individually fron NCBI as I was not able to do it automatically with the below commands. This means I had to go to the BioProject and get every link individually, then convert it to a fastq file and then gzip it and rename it for dDocent pipeline
--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

List of samples and SRA available here: <a href="https://www.ncbi.nlm.nih.gov/sra" class="uri">https://www.ncbi.nlm.nih.gov/sra</a>
-----------------------------------------------------------------------------------------------------------------------------------

\#Sample 43 does not have a SRA link \#Sample 42 does not have a SRA
link \#Sample 41 does not have a SRA link

1) Sample Name: Ofav41
======================

    wget https://sra-download.ncbi.nlm.nih.gov/traces/sra29/SRR/008637/SRR8844782
    fastq-dump --gzip SRR8844782 
    mv  SRR8844782.fastq.gz PAN_Ofav41.F.fq.gz
    rm -i  SRR8844782
    y

2) Sample Name: Ofav40
======================

            
    wget https://sra-download.ncbi.nlm.nih.gov/traces/sra21/SRR/008637/SRR8844779
    fastq-dump --gzip SRR8844779 
    mv SRR8844779.fastq.gz  PAN_Ofav40.F.fq.gz
    rm -i  SRR8844779
    y

3)Sample Name:Ofav39
====================

    wget https://sra-download.ncbi.nlm.nih.gov/traces/sra30/SRR/008637/SRR8844780
    fastq-dump --gzip SRR8844780
    mv  SRR8844780.fastq.gz PAN_Ofav39.F.fq.gz
    rm -i  SRR8844780
    y 

4) Sample Name:Ofav38
=====================

    wget https://sra-download.ncbi.nlm.nih.gov/traces/sra46/SRR/008637/SRR8844777
    fastq-dump --gzip SRR8844777
    mv SRR8844777.fastq.gz PAN_Ofav38.F.fq.gz
    rm -i SRR8844777

5) Sample Name:Ofav37
=====================

    wget https://sra-download.ncbi.nlm.nih.gov/traces/sra7/SRR/008637/SRR8844778
    fastq-dump --gzip SRR8844778
    mv   SRR8844778.fastq.gz PAN_Ofav37.F.fq.gz
    rm -i SRR8844778

6) Sample Name:Ofav36
=====================

    wget https://sra-download.ncbi.nlm.nih.gov/traces/sra6/SRR/008637/SRR8844775
    fastq-dump --gzip SRR8844775
    mv   SRR8844775.fastq.gz PAN_Ofav36.F.fq.gz
    rm -i SRR8844775

7) Sample Name:PAN\_Ofav35
==========================

    wget https://sra-download.ncbi.nlm.nih.gov/traces/sra43/SRR/008637/SRR8844776
    fastq-dump --gzip SRR8844776
    mv   SRR8844776.fastq.gz PAN_Ofav35.F.fq.gz
    rm -i  SRR8844776

8) Sample Name:Ofav34
=====================

    wget https://sra-download.ncbi.nlm.nih.gov/sos/sra-pub-run-2/SRR8844804/SRR8844804.1
    fastq-dump --gzip SRR8844804.1
    mv   SRR8844804.1.fastq.gz PAN_Ofav34.F.fq.gz
    rm -i  SRR8844804.1
    y

9) Sample Name:Ofav33
=====================

    wget https://sra-download.ncbi.nlm.nih.gov/traces/sra19/SRR/008637/SRR8844803
    fastq-dump --gzip SRR8844803
    mv   SRR8844803.fastq.gz PAN_Ofav33.F.fq.gz
    rm -i  SRR8844803
    y

10) Sample Name: Ofav32
=======================

    wget https://sra-download.ncbi.nlm.nih.gov/traces/sra24/SRR/008637/SRR8844810
    fastq-dump --gzip SRR8844810
    mv   SRR8844810.fastq.gz PAN_Ofav32.F.fq.gz
    rm -i SRR8844810 
    y

11) Sample Name: Ofav31
=======================

    wget https://sra-download.ncbi.nlm.nih.gov/traces/sra17/SRR/008637/SRR8844809
    fastq-dump --gzip SRR8844809
    mv   SRR8844809.fastq.gz PAN_Ofav31.F.fq.gz
    rm -i SRR8844809 
    y

12) Sample Name: Ofav30
=======================

    wget https://sra-download.ncbi.nlm.nih.gov/traces/sra8/SRR/008637/SRR8844812
    fastq-dump --gzip SRR8844812
    mv   SRR8844812.fastq.gz PAN_Ofav30.F.fq.gz
    rm -i  SRR8844812
    y

13) Sample Name:Ofav29
======================

    wget https://sra-download.ncbi.nlm.nih.gov/traces/sra9/SRR/008637/SRR8844811
    fastq-dump --gzip SRR8844811
    mv   SRR8844811.fastq.gz PAN_Ofav29.F.fq.gz
    rm -i  SRR8844811
    y

14) Sample Name: Ofav28
=======================

    wget https://sra-download.ncbi.nlm.nih.gov/traces/sra33/SRR/008637/SRR8844806
    fastq-dump --gzip SRR8844806
    mv   SRR8844806.fastq.gz PAN_Ofav28.F.fq.gz
    rm -i  SRR8844806
    y

Sample 29 does not contain file, only for PCR
---------------------------------------------

15) Sample Name: Ofav26
=======================

    wget https://sra-download.ncbi.nlm.nih.gov/traces/sra40/SRR/008637/SRR8844805
    fastq-dump --gzip SRR8844805
    mv   SRR8844805.fastq.gz PAN_Ofav26.F.fq.gz
    rm -i  SRR8844805
    y

16 ) Sample Name: Ofav25
========================

    wget https://sra-download.ncbi.nlm.nih.gov/traces/sra13/SRR/008637/SRR8844808
    fastq-dump --gzip SRR8844808
    mv   SRR8844808.fastq.gz PAN_Ofav25.F.fq.gz
    rm -i  SRR8844808
    y

17) Sample Name: Ofav24
=======================

    wget https://sra-download.ncbi.nlm.nih.gov/traces/sra53/SRR/008637/SRR8844807
    fastq-dump --gzip SRR8844807
    mv SRR8844807.fastq.gz PAN_Ofav24.F.fq.gz
    rm -i  SRR8844807
    y

18) Sample Name:Ofav23
======================

    wget https://sra-download.ncbi.nlm.nih.gov/traces/sra13/SRR/008637/SRR8844791
    fastq-dump --gzip SRR8844791
    mv   SRR8844791.fastq.gz PAN_Ofav23.F.fq.gz
    rm -i  SRR8844791
    y

19) Sample Name:Ofav22
======================

    wget https://sra-download.ncbi.nlm.nih.gov/traces/sra30/SRR/008637/SRR8844792
    fastq-dump --gzip SRR8844792
    mv   SRR8844792.fastq.gz PAN_Ofav22.F.fq.gz
    rm -i  SRR8844792
    y

20) Sample Name: Ofav20
=======================

    wget https://sra-download.ncbi.nlm.nih.gov/sos/sra-pub-run-2/SRR8844787/SRR8844787.1
    fastq-dump --gzip SRR8844787.1
    mv   SRR8844787.1.fastq.gz PAN_Ofav20.F.fq.gz
    rm -i  SRR8844787.1
    y

21) Sample Name: Ofav19
=======================

    wget https://sra-download.ncbi.nlm.nih.gov/sos/sra-pub-run-2/SRR8844788/SRR8844788.1
    fastq-dump --gzip SRR8844788.1
    mv   SRR8844788.1.fastq.gz PAN_Ofav19.F.fq.gz
    rm -i  SRR8844788.1
    y

22) Sample Name: Ofav18
=======================

    wget https://sra-download.ncbi.nlm.nih.gov/traces/sra32/SRR/008637/SRR8844789
    fastq-dump --gzip SRR8844789
    mv   SRR8844789.fastq.gz PAN_Ofav18.F.fq.gz
    rm -i  SRR8844789
    y

23) Sample Name: Ofav17
=======================

    wget https://sra-download.ncbi.nlm.nih.gov/traces/sra53/SRR/008637/SRR8844790
    fastq-dump --gzip SRR8844790
    mv   SRR8844790.fastq.gz PAN_Ofav17.F.fq.gz
    rm -i SRR8844790 
    y

24) Sample Name: Ofav16
=======================

    wget https://sra-download.ncbi.nlm.nih.gov/traces/sra10/SRR/008637/SRR8844783
    fastq-dump --gzip SRR8844783
    mv  SRR8844783.fastq.gz PAN_Ofav16.F.fq.gz
    rm -i SRR8844783 
    y

25) Sample Name:Ofav15
======================

    wget https://sra-download.ncbi.nlm.nih.gov/traces/sra25/SRR/008637/SRR8844784
    fastq-dump --gzip SRR8844784
    mv   SRR8844784.fastq.gz PAN_Ofav15.F.fq.gz
    rm -i SRR8844784  
    y

26) Sample Name: Ofav14
=======================

    wget https://sra-download.ncbi.nlm.nih.gov/traces/sra48/SRR/008637/SRR8844785
    fastq-dump --gzip SRR8844785
    mv SRR8844785.fastq.gz PAN_Ofav14.F.fq.gz
    rm -i  SRR8844785
    y

27) Sample Name:PAN\_Ofav13
===========================

    wget https://sra-download.ncbi.nlm.nih.gov/traces/sra40/SRR/008637/SRR8844786
    fastq-dump --gzip SRR8844786
    mv   SRR8844786.fastq.gz PAN_Ofav13.F.fq.gz
    rm -i  SRR8844786
    y

28) Sample Name:Ofav12
======================

    wget https://sra-download.ncbi.nlm.nih.gov/traces/sra32/SRR/008637/SRR8844793
    fastq-dump --gzip SRR8844793
    mv   SRR8844793.fastq.gz PAN_Ofav12.F.fq.gz
    rm -i  SRR8844793
    y

29) Sample Name:Ofav11
======================

    wget https://sra-download.ncbi.nlm.nih.gov/sos/sra-pub-run-2/SRR8844794/SRR8844794.1
    fastq-dump --gzip SRR8844794.1
    mv   SRR8844794.1.fastq.gz PAN_Ofav11.F.fq.gz
    rm -i  SRR8844794.1
    y

30) Sample Name: Ofav10
=======================

    wget https://sra-download.ncbi.nlm.nih.gov/sos/sra-pub-run-2/SRR8844801/SRR8844801.1
    fastq-dump --gzip SRR8844801.1
    mv   SRR8844801.1.fastq.gz PAN_Ofav10.F.fq.gz
    rm -i  SRR8844801.1
    y

31) Sample Name: Ofav9
======================

    wget https://sra-download.ncbi.nlm.nih.gov/sos/sra-pub-run-2/SRR8844802/SRR8844802.1
    fastq-dump --gzip SRR8844802.1
    mv   SRR8844802.1.fastq.gz PAN_Ofav9.F.fq.gz
    rm -i  SRR8844802.1
    y

32) Sample Name:Ofav8
=====================

    wget https://sra-download.ncbi.nlm.nih.gov/traces/sra11/SRR/008637/SRR8844799
    fastq-dump --gzip SRR8844799
    mv   SRR8844799.fastq.gz PAN_Ofav8.F.fq.gz
    rm -i  SRR8844799
    y

33) Sample Name:Ofav7
=====================

    wget https://sra-download.ncbi.nlm.nih.gov/traces/sra39/SRR/008637/SRR8844800
    fastq-dump --gzip SRR8844800
    mv   SRR8844800.fastq.gz PAN_Ofav7.F.fq.gz
    rm -i  SRR8844800
    y

For some reason this sample (SRR8844797.1) did not download so I downloaded to my mac an then uploaded it manually using R
--------------------------------------------------------------------------------------------------------------------------

34) Sample Name:Ofav6
=====================

    wget https://sra-download.ncbi.nlm.nih.gov/sos/sra-pub-run-2/SRR8844797/SRR8844797.1
    fastq-dump --gzip SRR8844797.1
    mv   SRR8844797.1.fastq.gz PAN_Ofav6.F.fq.gz
    rm -i  SRR8844797.1
    y

35) Sample Name:Ofav5
=====================

    wget https://sra-download.ncbi.nlm.nih.gov/traces/sra23/SRR/008637/SRR8844798
    fastq-dump --gzip SRR8844798
    mv   SRR8844798.fastq.gz PAN_Ofav5.F.fq.gz
    rm -i  SRR8844798
    y

There is no sample 4
--------------------

36) Sample Name: Ofav3
======================

    wget https://sra-download.ncbi.nlm.nih.gov/traces/sra42/SRR/008637/SRR8844795
    fastq-dump --gzip SRR8844795
    mv   SRR8844795.fastq.gz PAN_Ofav3.F.fq.gz
    rm -i  SRR8844795
    y

There is no sample 2
--------------------

37) Sample Name: \_Ofav1
========================

    wget https://sra-download.ncbi.nlm.nih.gov/traces/sra0/SRR/008637/SRR8844796
    fastq-dump --gzip SRR8844796
    mv   SRR8844796.fastq.gz PAN_Ofav1.F.fq.gz
    rm -i  SRR8844796
    y

38) Sample Name: Ofav45
=======================

    wget https://sra-download.ncbi.nlm.nih.gov/traces/sra0/SRR/008637/SRR8844774
    fastq-dump --gzip SRR8844774 
    mv  SRR8844774.fastq.gz PAN_Ofav45.F.fq.gz
    rm -i  SRR8844774
    y

39) Sample Name: Ofav44
=======================

    wget https://sra-download.ncbi.nlm.nih.gov/sos/sra-pub-run-2/SRR8844781/SRR8844781.1
    fastq-dump --gzip SRR8844781.1 
    mv  SRR8844781.1.fastq.gz PAN_Ofav44.F.fq.gz
    rm -i  SRR8844781.1
    y

Now we need to copy our geneome reference to run dDocent pipeline or you can download it from <a href="ftp://ftp.ncbi.nlm.nih.gov/sra/wgs_aux/MZ/GG/MZGG01/MZGG01.1.fsa_nt.gz" class="uri">ftp://ftp.ncbi.nlm.nih.gov/sra/wgs_aux/MZ/GG/MZGG01/MZGG01.1.fsa_nt.gz</a>
---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

Copy O.faveolata genome to working directory
--------------------------------------------

    cp /RAID_STORAGE2/mgomez/Orbicella_Raw/Libraries_compiled/reference.fasta .

Running dDocent on demultiplex files Dziedzic\_2019\_Data
---------------------------------------------------------

    dDocent
    dDocent 2.7.7 

    Contact jpuritz@uri.edu with any problems 

     
    Checking for required software

    All required software is installed!

    dDocent run started Thu Aug 15 17:45:26 EDT 2019 

    39 individuals are detected. Is this correct? Enter yes or no and press [ENTER]
    yes

    Proceeding with 39 individuals
    dDocent detects 80 processors available on this system.
    Please enter the maximum number of processors to use for this analysis.
    20 
    dDocent detects 503 gigabytes of maximum memory available on this system.
    Please enter the maximum memory to use for this analysis in gigabytes
    For example, to limit dDocent to ten gigabytes, enter 10
    This option does not work with all distributions of Linux.  If runs are hanging at variant calling, enter 0
    Then press [ENTER]
    0
    Do you want to quality trim your reads?
    Type yes or no and press [ENTER]?
    yes 

    Do you want to perform an assembly?
    Type yes or no and press [ENTER].
    no

    Reference contigs need to be in a file named reference.fasta

    Do you want to map reads?  Type yes or no and press [ENTER]
    yes
    BWA will be used to map reads.  You may need to adjust -A -B and -O parameters for your taxa.
    Would you like to enter a new parameters now? Type yes or no and press [ENTER]
    no

    Proceeding with default values for BWA read mapping.
    Do you want to use FreeBayes to call SNPs?  Please type yes or no and press [ENTER]
    yes

    Please enter your email address.  dDocent will email you when it is finished running.
    Don`t worry; dDocent has no financial need to sell your email address to spammers.
    matias_gomez@uri.edu

    At this point, all configuration information has been entered and dDocent may take several hours to run.
    It is recommended that you move this script to a background operation and disable terminal input and output.
    All data and logfiles will still be recorded.
    To do this:
    Press control and Z simultaneously
    Type 'bg' without the quotes and press enter
    Type 'disown -h' again without the quotes and press enter

    Now sit back, relax, and wait for your analysis to finish

    Trimming reads
    Removing the _1 character and replacing with /1 in the name of every sequence
    ^Z
    [1]+  Stopped                 dDocent
    (Orbicella) [mgomez@KITT Dziedzic_2019_Data]$ bg
    [1]+ dDocent &
    (Orbicella) [mgomez@KITT Dziedzic_2019_Data]$ disown -h

dDocent has finished with an analysis in
/RAID\_STORAGE2/mgomez/Orbicella\_Raw/Dziedzic\_2019\_Data

dDocent started Thu Aug 15 17:45:26 EDT 2019

dDocent finished Fri Aug 16 08:13:08 EDT 2019

After filtering, kept 11611 out of a possible 383160 Sites

dDocent 2.7.7 The ‘d’ is silent, hillbilly.

SNP filtering (<a href="http://www.ddocent.com/filtering/" class="uri">http://www.ddocent.com/filtering/</a>)
-------------------------------------------------------------------------------------------------------------

First, let´s create a working dir and change to it, then copy the raw VCF file generetated by dDocent to start with the SNP filtering
-------------------------------------------------------------------------------------------------------------------------------------

    mkdir SNP_Filtering_Dziedzic_2019
    cd SNP_Filtering_Dziedzic_2019
    cp /RAID_STORAGE2/mgomez/Orbicella_Raw/Dziedzic_2019_Data/TotalRawSNPs.vcf .
