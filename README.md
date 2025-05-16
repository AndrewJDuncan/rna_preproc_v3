~~ v4 notes ~~  
For my own reference. Some of the below pasted from salmon docs.

**Pipeline v4 steps:**  
Initial stats  
PhiX removal  
UMI extraction  
Read cleaning  
Alignment  
Sorting/indexing  
UMI deduplication  
Convert BAM to FASTQ for Salmon  
Final stats  


# Salmon output --> counts table  
Salmon output is saved in a newly-created directory, salmon_quant. Each sample will have its own subdirectory. In each is the main quantification file, quant.sf, This has a single header row and 5 columns as below.  
Name / Length / EffectiveLength / TPM / NumReads  
  - Name — This is the name of the target transcript provided in the input transcript database (FASTA file)    
  - Length — This is the length of the target transcript in nucleotides  
  - EffectiveLength — This is the computed effective length of the target transcript. It takes into account all factors being modeled that will effect the probability of sampling fragments from this transcript, including the fragment length distribution and sequence-specific and gc-fragment bias (if they are being modeled).  
  - TPM — This is salmon’s estimate of the relative abundance of this transcript in units of Transcripts Per Million (TPM). *TPM is the recommended relative abundance measure to use for downstream analysis.*  
  - NumReads — This is salmon’s estimate of the number of reads mapping to each transcript that was quantified. It is an “estimate” insofar as it is the expected number of reads that have originated from each transcript given the structure of the uniquely mapping and multi-mapping reads and the relative abundance estimates for each transcript.


**Step 1**: Checking output 
```cd``` into main directory for sequencing run  
```ls -lah ./salmon_quant/*/quant.sf``` #see how many output files built  
```head ./salmon_quant/$SAMPLENAME/quant.sf``` #obvs put the name, right?  
Should show the correct output structure  
  
```tail -n +2 ./salmon_quant/$SAMPLENAME/quant.sf | head```  
- Tail -n +2 means start from the second row (omit header row)  
- Pipes to ```head``` to make sure it looks right  
  
We actually only want the TPM column (column 4):  
```tail -n +2 ./salmon_quant/$SAMPLENAME/quant.sf | head | cut -f4 > $SAMPLENAME.out.tab.count```  
- Pipes to ```cut``` then makes an out.tab with a single column of data: gene counts  

**Make counts files**
*ensure you're in the right directory, per part 1*  
This command:  
(1) makes directories for counts output
(2) loops through all the read count files, printing the names of each and assigning a file name to the ‘${sample}’ variable.  
(3)  runs the above command to extract col 4 and exclude row 1, creating a counts file for each sample
```cd /share/biocore/workshop/mrnaseq_workshop/$USER/rnaseq_example
mkdir Counts
mkdir Counts/tmp
for sample in `cat samples.txt`; do \
    echo ${sample}
    cat 02-STAR_alignment/${sample}/${sample}_ReadsPerGene.out.tab | tail -n +5 | cut -f4 > 03-Counts/tmp/${sample}.count
done```
