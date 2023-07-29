#########################################################################################################################
#	This shell script takes input as folder name of Normal and Tumor.						#
#	e.g WGS_pipeline_TumorVsNormal_hg38_v0.1.sh Blood Tumor 							#
#	Last Modified on 27/07/2023											#
#	Last modified: Annovar correction										#
#															#
#															#
#															#
#	Tools Used in this pipeline											#
#	1. Fastqc													#
#	2. Fastp													#
#	3. STAR														#
#	4. Markduplicate												#
#	5. cigar													#
#	6. AddOrReplaceReadGroups											#
#	7. bqsr														#
#	8. Star_Fusion													#
#	9. Fusion_Catcher												#
#	10. RSEM													#
#	11. Arriba													#
#	12. Mutect2													#
#	13. Varscan													#
#	14. HaplotypeCaller												#
#########################################################################################################################

# Taking user inputs and setting resource paths
input1="$1"	# Taking input of normal folder name as argument
workdir=$(pwd)	# storing path of Working director
raw_tumor="$workdir/$input1"	# Assigning path of tumor sample 
out="$workdir/Output_hg38_${input1}"	# Assigning path of output folder
db="/mnt/database"		# Assigning database folder path
hg38="$db/GRCH38.P14/hg38.fa"	# Assigning path of human genome reference file


# Removing pre-exiting files if any

rm $out/RNASeq_analysis.log
rm $out/parallel_fastqc
rm $out/parallel_fastp
rm $out/parallel_fastqc_after_fastp
rm $out/parallel_STAR

# Merging fastq files and Creating list file
mkdir -p $out	# Creating Output folder
echo "Merging fastq files for $input1 Started" >> $out/RNASeq_analysis.log	# Writing in log file
date >> $out/RNASeq_analysis.log					# Writing in log file

cd $input1

fn=$(ls | head -1| cut -f 1-4 -d "_")
echo "Merged file base name $fn"
zcat *R1* > ${fn}_All_R1_001.fastq
zcat *R2* > ${fn}_All_R2_001.fastq

ls *_All_R1_001* > ../list_${input1}

echo "Generated merged file size" >> $out/RNASeq_analysis.log	# Writing in log file
ls -trlh *All*  >> $out/RNASeq_analysis.log	# Writing in log file
echo "Merging fastq files for $input1 Completed" >> $out/RNASeq_analysis.log	# Writing in log file
date >> $out/RNASeq_analysis.log					# Writing in log file
echo "##############################" >> $out/RNASeq_analysis.log	# Writing in log file

cd $workdir

input2="list_${input1}"

#Dynamic allocation of cpus

t="$(nproc --all)"      # Fetching total number of cpus
tnp=$((`expr $t - 1`)) # Maximum number of cpus will be utilized
val=$((`expr $tnp / 2`)) # Arithmatic operation on Max cpu count
np=$((`printf "%.*f\n" 0 $val`)) # Number of shared cpu



echo "First input $input1 Second input $input2" # Print inputs as a test

# Creating folders for the intermediate results
mkdir -p $out/$input1/alignments_stats
mkdir -p $out/$input1/FastQC_output
mkdir -p $out/$input1/fastp_output
mkdir -p $out/$input1/FastQC_output_after_fastp
mkdir -p $out/$input1/STAR
mkdir -p $out/$input1/Markduplicate
mkdir -p $out/$input1/cigar
mkdir -p $out/$input1/AddOrReplaceReadGroups
mkdir -p $out/$input1/bqsr
mkdir -p $out/$input1/Star_Fusion
mkdir -p $out/$input1/Fusion_Catcher
mkdir -p $out/$input1/RSEM
mkdir -p $out/$input1/Arriba
mkdir -p $out/variantCaller/Mutect2
mkdir -p $out/variantCaller/Varscan
mkdir -p $out/variantCaller/HaplotypeCaller


eval "$(conda shell.bash hook)"		# Setting bash for conda environment
conda activate wgs_gatk4		# Activating conda environment

# Fastqc and Trimmomatics run for Normal
while IFS= read -r line			# While loop to add scripts for each lane in "parallel_fastqc_Normal" 
	do 				# "parallel_trimmomatic_Normal" and "parallel_after_trimmomatic_Normal"
		echo $line
		sm=$(echo $raw_tumor/$line | xargs -n 1 basename | cut -f 1-4 -d"_"); 	# Extracting file name from file 
		nam=$(echo $raw_tumor/$line | xargs -n 1 basename | cut -f 1-3 -d"_");	# list
		echo "sm base name $sm"
		echo "nam base name $nam"

# QC checking

		echo "fastqc $raw_tumor/${sm}_All_R1_001.fastq $raw_tumor/${sm}_All_R2_001.fastq -t $np -o $out/$input1/FastQC_output/" >> $out/parallel_fastqc

# Adaptor Trimming and cleaning

		echo "fastp --detect_adapter_for_pe --overrepresentation_analysis --correction --cut_right --thread $tnp --html Sample_Name.fastp.html --json Sample_Name.fastp.json -i $raw_tumor/${sm}_All_R1_001.fastq -o $out/$input1/fastp_output/${sm}_R1-trimmed.fastq -I $raw_tumor/${sm}_All_R2_001.fastq -O $out/$input1/fastp_output/${sm}_R2-trimmed.fastq" >> $out/parallel_fastp

# QC-rechecking after trimming

		echo "fastqc -t $tnp $out/$input1/fastp_output/${sm}_R1-trimmed.fastq $out/$input1/fastp_output/${sm}_R2-trimmed.fastq -o $out/$input1/FastQC_output_after_fastp/" >> $out/parallel_fastqc_after_fastp

	done < "$input2"


echo "Fastqc for $input1 Started" >> $out/RNASeq_analysis.log	# Writing in log file
date >> $out/RNASeq_analysis.log					# Writing in log

parallel -j 1 < $out/parallel_fastqc		# Running script in parallel

echo "Fastqc for $input1 and $input2 Completed" >> $out/RNASeq_analysis.log	# Writing in log file
date >> $out/RNASeq_analysis.log					# Writing in log file
echo "##############################" >> $out/RNASeq_analysis.log	# Writing in log file

echo "Trimmomatic for $input1 Started" >> $out/RNASeq_analysis.log	# Writing in log file
date >> $out/RNASeq_analysis.log					# Writing in log file

parallel -j 1 < $out/parallel_fastp		# Running script in parallel

echo "Trimmomatic for $input1 and $input2 Completed" >> $out/RNASeq_analysis.log	# Writing in log file
date >> $out/RNASeq_analysis.log						# Writing in log file
echo "##############################" >> $out/RNASeq_analysis.log		# Writing in log file

echo "Fastqc after trimmomatic for $input1 and $input2 Started" >> $out/RNASeq_analysis.log	# Writing in log file
date >> $out/RNASeq_analysis.log							# Writing in log file

parallel -j 1 < $out/parallel_fastqc_after_fastp		# Running script in parallel

echo "Fastqc after trimmomatic for $input1 Completed" >> $out/RNASeq_analysis.log	# Writing in log file
date >> $out/RNASeq_analysis.log							# Writing in log file
echo "##############################" >> $out/RNASeq_analysis.log			# Writing in log file



############# STAR for Tumor sample

while IFS= read -r line
do
		
	echo $line
	header=$(cat $raw_tumor/$line | head -n 1);	# Extracting header from fastq file. Will not work if fastq file
							# is not gunzip file. Change command to cat
	id=$(echo $header1 | cut -f 3-4 -d":" | sed 's/@//');	# Extracting run ID from fastq file
	sm=$(echo $raw_tumor/$line | xargs -n 1 basename | cut -f 1-4 -d"_");	# Extracting normal file name upto 3 place
										# delimited by '_'
	nam=$(echo $raw_tumor/$line | xargs -n 1 basename | cut -f 1-3 -d"_");  # Extracting normal file name upto
										# 2 places delimited by '_'
	echo "header contant" $header # printing header
	echo $id
	echo $sm
	echo $nam
	if [ -z "$(ls -A  $out/$input1/fastp_output/)" ]   # Checking folder if empty
	then

		echo "STAR --genomeDir $db/GRCH38.P14/Transcriptomic_db/Genecode_GRCH38.P13/Genomic_db/indexes_v2 --readFilesIn $raw_tumor/${sm}_All_R1_001.fastq $raw_tumor/${sm}_All_R2_001.fastq --outReadsUnmapped None --twopassMode Basic --readFilesCommand cat --outSAMstrandField intronMotif --outSAMunmapped Within --chimSegmentMin 12 --chimJunctionOverhangMin 8 --chimOutJunctionFormat 1 --alignSJDBoverhangMin 10 --alignMatesGapMax 100000 --alignIntronMax 100000 --alignSJstitchMismatchNmax 5 -1 5 5  --chimMultimapScoreRange 3 --chimScoreJunctionNonGTAG -4 --chimMultimapNmax 20 --chimNonchimScoreDropMin 10 --peOverlapNbasesMin 12 --peOverlapMMp 0.1 --alignInsertionFlush Right --alignSplicedMateMapLminOverLmate 0 --alignSplicedMateMapLmin 30 --outSAMtype BAM SortedByCoordinate --quantMode GeneCounts --sjdbGTFfile $db/GRCH38.P14/Transcriptomic_db/Genecode_GRCH38.P13/GTF_file/gencode.v38.annotation.gtf --sjdbOverhang 99 --runThreadN $tnp --outTmpDir $out/$input1/STAR/align/tmp --outFileNamePrefix $out/$input1/STAR/align/${sm}_" > $out/parallel_STAR   # Writing file for parallel run when trimmomatic_output is empty

	else
		
		echo "STAR --genomeDir $db/GRCH38.P14/Transcriptomic_db/Genecode_GRCH38.P13/Genomic_db/indexes_v2 --readFilesIn $out/$input1/fastp_output/${sm}_R1-trimmed.fastq $out/$input1/fastp_output/${sm}_R2-trimmed.fastq --outReadsUnmapped None --twopassMode Basic --readFilesCommand cat --outSAMstrandField intronMotif --outSAMunmapped Within --chimSegmentMin 12 --chimJunctionOverhangMin 8 --chimOutJunctionFormat 1 --alignSJDBoverhangMin 10 --alignMatesGapMax 100000 --alignIntronMax 100000 --alignSJstitchMismatchNmax 5 -1 5 5  --chimMultimapScoreRange 3 --chimScoreJunctionNonGTAG -4 --chimMultimapNmax 20 --chimNonchimScoreDropMin 10 --peOverlapNbasesMin 12 --peOverlapMMp 0.1 --alignInsertionFlush Right --alignSplicedMateMapLminOverLmate 0 --alignSplicedMateMapLmin 30 --outSAMtype BAM SortedByCoordinate --quantMode GeneCounts --sjdbGTFfile $db/GRCH38.P14/Transcriptomic_db/Genecode_GRCH38.P13/GTF_file/gencode.v38.annotation.gtf --sjdbOverhang 99 --runThreadN $tnp --outTmpDir $out/$input1/STAR/align/tmp --outFileNamePrefix $out/$input1/STAR/align/${sm}_" > $out/parallel_STAR   # Writing file for parallel run when trimmomatic_output is empty
		
	fi
	done < "$input2"	# Passing argument for list normal

eval "$(conda shell.bash hook)"         # Setting bash for conda environment
conda activate base 
echo "STAR for $input1 Started" >> $out/RNASeq_analysis.log	# Writing in log file
date >> $out/RNASeq_analysis.log						# Writing in log file

parallel -j 1 < $out/parallel_STAR	# Running script in parallel

echo "STAR for $input1 Completed" >> $out/RNASeq_analysis.log	# Writing in log file
date >> $out/RNASeq_analysis.log						# Writing in log file
echo "##############################" >> $out/RNASeq_analysis.log		# Writing in log file

eval "$(conda shell.bash hook)"		# Setting bash for conda environment
conda activate wgs_gatk4	# Activating conda environme

#markduplicates

echo "Markduplicate for $input1 Started" >> $out/RNASeq_analysis.log     # Writing in log file
date >> $out/RNASeq_analysis.log                                                # Writing in log file

gatk --java-options "-Xmx64G" MarkDuplicates -I $out/$input1/STAR/align/${sm}_Aligned.sortedByCoord.out.bam -O $out/$input1/Markduplicate/${sm}_Aligned.sortedByCoord.out.dedup.bam -M Metrics.txt --MAX_RECORDS_IN_RAM 1000000000

echo "Markduplicate for $input1 Completed" >> $out/RNASeq_analysis.log   # Writing in log file
date >> $out/RNASeq_analysis.log                                                # Writing in log file
echo "##############################" >> $out/RNASeq_analysis.log               # Writing in log file


#SplitNCigarReads

echo "SplitNCigarReads for $input1 Started" >> $out/RNASeq_analysis.log     # Writing in log file
date >> $out/RNASeq_analysis.log                                                # Writing in log file

gatk SplitNCigarReads -R /mnt/database/GRCH38.P14/Transcriptomic_db/Genecode_GRCH38.P13/Genomic_db/GRCh38.p13.genome.fa -I $out/$input1/Markduplicate/${sm}_Aligned.sortedByCoord.out.dedup.bam -O $out/$input1/cigar/${sm}_Aligned.sortedByCoord.out.dedup.cigar.bam

echo "SplitNCigarReads for $input1 Completed" >> $out/RNASeq_analysis.log   # Writing in log file
date >> $out/RNASeq_analysis.log                                                # Writing in log file
echo "##############################" >> $out/RNASeq_analysis.log               # Writing in log file

#AddOrReplaceReadGroup

echo "AddOrReplaceReadGroups for $input1 Started" >> $out/RNASeq_analysis.log     # Writing in log file
date >> $out/RNASeq_analysis.log                                                # Writing in log file

gatk --java-options "-Xmx64G" AddOrReplaceReadGroups -I $out/$input1/cigar/${sm}_Aligned.sortedByCoord.out.dedup.cigar.bam -O $out/$input1/cigar/${sm}_Aligned.sortedByCoord.out.dedup.cigar.RG.bam -SORT_ORDER coordinate -RGID 4 -RGLB lib1 -RGPL illumina -RGSM 20 -RGPU unit1 -CREATE_INDEX True --MAX_RECORDS_IN_RAM 100000000

echo "AddOrReplaceReadGroups for $input1 Completed" >> $out/RNASeq_analysis.log   # Writing in log file
date >> $out/RNASeq_analysis.log                                                # Writing in log file
echo "##############################" >> $out/RNASeq_analysis.log               # Writing in log file

#bqsr

echo "BaseRecalibrator for $input1 Started" >> $out/RNASeq_analysis.log	# Writing in log file
date >> $out/RNASeq_analysis.log                                                # Writing in log file

gatk --java-options "-Xmx64G" BaseRecalibrator -R /mnt/database/GRCH38.P14/Transcriptomic_db/Genecode_GRCH38.P13/Genomic_db/GRCh38.p13.genome.fa -I $out/$input1/cigar/${sm}_Aligned.sortedByCoord.out.dedup.cigar.RG.bam --known-sites /mnt/database/dbSNP_hg38/Homo_sapiens_assembly38.dbsnp138.vcf.gz --known-sites /mnt/database/dbSNP_hg38/Homo_sapiens_assembly38.known_indels.vcf -O $out/$input1/bqsr/${sm}_output_baserecal.table

#applybqsr

gatk --java-options "-Xmx64G" ApplyBQSR -R /mnt/database/GRCH38.P14/Transcriptomic_db/Genecode_GRCH38.P13/Genomic_db/GRCh38.p13.genome.fa -I $out/$input1/cigar/${sm}_Aligned.sortedByCoord.out.dedup.cigar.RG.bam -bqsr $out/$input1/bqsr/${sm}_output_baserecal.table -O $out/$input1/bqsr/${sm}_Aligned.sortedByCoord.out.dedup.cigar.RG.bqsr.bam


echo "BaseRecalibrator for $input1 Completed" >> $out/RNASeq_analysis.log   # Writing in log file
date >> $out/RNASeq_analysis.log                                                # Writing in log file
echo "##############################" >> $out/RNASeq_analysis.log               # Writing in log file

##VariantCaller

echo "HaplotypeCaller for $input1 Started" >> $out/RNASeq_analysis.log	# Writing in log file
date >> $out/RNASeq_analysis.log                                        # Writing in log file

pids=
for chr in chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX chrY

do	
               gatk --java-options "-Xmx5g" HaplotypeCaller -R /mnt/database/GRCH38.P14/Transcriptomic_db/Genecode_GRCH38.P13/Genomic_db/GRCh38.p13.genome.fa -I $out/$input1/bqsr/${sm}_Aligned.sortedByCoord.out.dedup.cigar.RG.bqsr.bam --dbsnp /mnt/database/dbSNP_hg38/Homo_sapiens_assembly38.dbsnp138.vcf.gz -O $out/variantCaller/HaplotypeCaller/$chr.vcf --standard-min-confidence-threshold-for-calling 20  --native-pair-hmm-threads 10 -L $chr &
                      
                        pids+=" $!"
                done;

gatk --java-options "-Xmx16g" MergeVcfs \
      I=$out/variantCaller/HaplotypeCaller/chr1.vcf \
      I=$out/variantCaller/HaplotypeCaller/chr2.vcf \
      I=$out/variantCaller/HaplotypeCaller/chr3.vcf \
      I=$out/variantCaller/HaplotypeCaller/chr4.vcf \
      I=$out/variantCaller/HaplotypeCaller/chr5.vcf \
      I=$out/variantCaller/HaplotypeCaller/chr6.vcf \
      I=$out/variantCaller/HaplotypeCaller/chr7.vcf \
      I=$out/variantCaller/HaplotypeCaller/chr8.vcf \
      I=$out/variantCaller/HaplotypeCaller/chr9.vcf \
      I=$out/variantCaller/HaplotypeCaller/chr10.vcf \
      I=$out/variantCaller/HaplotypeCaller/chr11.vcf \
      I=$out/variantCaller/HaplotypeCaller/chr12.vcf \
      I=$out/variantCaller/HaplotypeCaller/chr13.vcf \
      I=$out/variantCaller/HaplotypeCaller/chr14.vcf \
      I=$out/variantCaller/HaplotypeCaller/chr15.vcf \
      I=$out/variantCaller/HaplotypeCaller/chr16.vcf \
      I=$out/variantCaller/HaplotypeCaller/chr17.vcf \
      I=$out/variantCaller/HaplotypeCaller/chr18.vcf \
      I=$out/variantCaller/HaplotypeCaller/chr19.vcf \
      I=$out/variantCaller/HaplotypeCaller/chr20.vcf \
      I=$out/variantCaller/HaplotypeCaller/chr21.vcf \
      I=$out/variantCaller/HaplotypeCaller/chr22.vcf \
      I=$out/variantCaller/HaplotypeCaller/chrX.vcf \
      I=$out/variantCaller/HaplotypeCaller/chrY.vcf \
      O=$out/variantCaller/HaplotypeCaller/${sm}_HaplotypeCaller.vcf

echo "Rscript mutect2Reformat.R $out/variantCaller/HaplotypeCaller/${sm}_HaplotypeCaller.vcf" >> $out/parallel_Rscript

echo "HaplotypeCaller for $input1 Completed" >> $out/RNASeq_analysis.log   # Writing in log file
date >> $out/RNASeq_analysis.log                                                # Writing in log file
echo "##############################" >> $out/RNASeq_analysis.log               # Writing in log file


#Mutect2

echo "Mutect2 for $input1 Started" >> $out/RNASeq_analysis.log   # Writing in log file
date >> $out/RNASeq_analysis.log                                                # Writing in log file

eval "$(conda shell.bash hook)"         # Setting bash for conda environment
conda activate wgs_gatk4                # Activating conda environment


pids=
for chr in chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX chrY

do
                gatk --java-options "-Xmx5g" Mutect2 -R /mnt/database/GRCH38.P14/Transcriptomic_db/Genecode_GRCH38.P13/Genomic_db/GRCh38.p13.genome.fa -I /mnt/data/Abhishek/Abhishek_RNASeq/bqsr/${sm}_Aligned.sortedByCoord.out.dedup.cigar.RG.bqsr.bam --germline-resource /mnt/database/GRCH38.P14/gnomAD/af-only-gnomad.hg38.vcf.gz --panel-of-normals /mnt/database/GRCH38.P14/panel-of-normals/1000g_pon.hg38.vcf.gz -O $out/variantCaller/Mutect2/$chr.vcf --native-pair-hmm-threads 10 -L $chr &

                        pids+=" $!"
                done;

gatk --java-options "-Xmx16g" MergeVcfs \
	I=$out/variantCaller/Mutect2/chr1.vcf \
	I=$out/variantCaller/Mutect2/chr2.vcf \
	I=$out/variantCaller/Mutect2/chr3.vcf \
	I=$out/variantCaller/Mutect2/chr4.vcf \
	I=$out/variantCaller/Mutect2/chr5.vcf \
	I=$out/variantCaller/Mutect2/chr6.vcf \
	I=$out/variantCaller/Mutect2/chr7.vcf \
	I=$out/variantCaller/Mutect2/chr8.vcf \
	I=$out/variantCaller/Mutect2/chr9.vcf \
	I=$out/variantCaller/Mutect2/chr10.vcf \
	I=$out/variantCaller/Mutect2/chr11.vcf \
	I=$out/variantCaller/Mutect2/chr12.vcf \
	I=$out/variantCaller/Mutect2/chr13.vcf \
	I=$out/variantCaller/Mutect2/chr14.vcf \
	I=$out/variantCaller/Mutect2/chr15.vcf \
	I=$out/variantCaller/Mutect2/chr16.vcf \
	I=$out/variantCaller/Mutect2/chr17.vcf \
	I=$out/variantCaller/Mutect2/chr18.vcf \
	I=$out/variantCaller/Mutect2/chr19.vcf \
	I=$out/variantCaller/Mutect2/chr20.vcf \
	I=$out/variantCaller/Mutect2/chr21.vcf \
	I=$out/variantCaller/Mutect2/chr22.vcf \
	I=$out/variantCaller/Mutect2/chrX.vcf \
	I=$out/variantCaller/Mutect2/chrY.vcf \
	O=$out/variantCaller/Mutect2/${sm}_Mutect2.vcf

echo "Rscript mutect2Reformat.R $out/variantCaller/Mutect2/${sm}_Mutect2.vcf" >> $out/parallel_Rscript

parallel -j 1 < $out/parallel_Rscript

echo "Mutect2 for $input1 Completed" >> $out/RNASeq_analysis.log   # Writing in log file
date >> $out/RNASeq_analysis.log                                                # Writing in log file
echo "##############################" >> $out/RNASeq_analysis.log               # Writing in log file


#Varscan

echo "Varscan for $input1 Started" >> $out/RNASeq_analysis.log   # Writing in log file
date >> $out/RNASeq_analysis.log                                 # Writing in log file

samtools mpileup -B -f /mnt/database/GRCH38.P14/Transcriptomic_db/Genecode_GRCH38.P13/Genomic_db/GRCh38.p13.genome.fa $out/$input1/bqsr/${sm}_Aligned.sortedByCoord.out.dedup.cigar.RG.bqsr.bam > $out/variantCaller/Varscan/${sm}_Aligned.sortedByCoord.out.dedup.cigar.RG.bqsr.mpileup.bam

java -XX:+UseParallelGC -XX:ParallelGCThreads=64 -Xms64g -Xmx64g -jar /home/genomics/anaconda3/envs/wgs_gatk4/share/varscan-2.3.7-4/VarScan.jar mpileup2cns $out/variantCaller/Varscan/${sm}_Aligned.sortedByCoord.out.dedup.cigar.RG.bqsr.mpileup.bam --min-coverage 10 --min-reads2 5 --output-vcf 1 --variants > $out/variantCaller/Varscan/${sm}_Varscan.vcf

echo "Varscan for $input1 Completed" >> $out/RNASeq_analysis.log   # Writing in log file
date >> $out/RNASeq_analysis.log                                                # Writing in log file
echo "##############################" >> $out/RNASeq_analysis.log               # Writing in log file


#Star-Fusion

echo "Star-Fusion for $input1 Started" >> $out/RNASeq_analysis.log   # Writing in log file
date >> $out/RNASeq_analysis.log                                 # Writing in log file

/home/genomics/software/STAR-Fusion-v1.12.0/./STAR-Fusion --left_fq $out/$input1/fastp_output/${sm}_R1-trimmed.fastq --right_fq $out/$input1/fastp_output/${sm}_R2-trimmed.fastq --genome_lib_dir /mnt/database/GRCH38.P14/Transcriptomic_db/Genecode_GRCH38.P13/GRCh38_gencode_v37_CTAT_lib_Mar012021_plug-n-play/GRCh38_gencode_v37_CTAT_lib_Mar012021.plug-n-play/ctat_genome_lib_build_dir --CPU $tnp --output_dir $out/$input1/Star_Fusion/

echo "Star-Fusion for $input1 Completed" >> $out/RNASeq_analysis.log   # Writing in log file
date >> $out/RNASeq_analysis.log                                                # Writing in log file
echo "##############################" >> $out/RNASeq_analysis.log               # Writing in log file


#RSEM

echo "RSEM for $input1 Started" >> $out/RNASeq_analysis.log   # Writing in log file
date >> $out/RNASeq_analysis.log                                 # Writing in log file

rsem-calculate-expression -p $tnp --paired-end  --estimate-rspd --append-names --output-genome-bam  $out/$input1/fastp_output/${sm}_R1-trimmed.fastq $out/$input1/fastp_output/${sm}_R2-trimmed.fastq /mnt/database/GRCH38.P14/Transcriptomic_db/Genecode_GRCH38.P13/RSEM/ref/human_gencode $out/$input1/RSEM/$sm

echo "RSEM for $input1 Completed" >> $out/RNASeq_analysis.log   # Writing in log file
date >> $out/RNASeq_analysis.log                                                # Writing in log file
echo "##############################" >> $out/RNASeq_analysis.log               # Writing in log file

#Fusion_catcher

echo "Fusion_catcher for $input1 Started" >> $out/RNASeq_analysis.log   # Writing in log file
date >> $out/RNASeq_analysis.log                                                # Writing in log file


eval "$(conda shell.bash hook)"
conda activate fusioncatcher

fusioncatcher -d /mnt/database/fusioncatcher_v2/ -i $out/$input1/fastp_output/${sm}_R1-trimmed.fastq,$out/$input1/fastp_output/${sm}_R2-trimmed.fastq  --Xmx=32G -p $tnp -o $out/$input1/Fusioncatcher_Output --aligners blat,star,bowtie2

echo "Fusion_catcher for $input1 Completed" >> $out/RNASeq_analysis.log   # Writing in log file
date >> $out/RNASeq_analysis.log                                                # Writing in log file
echo "##############################" >> $out/RNASeq_analysis.log               # Writing in log file


#Arriba

eval "$(conda shell.bash hook)"
conda activate base

pigz -k -p $tnp $out/$input1/fastp_output/*fastq

cd  $out/$input1/Arriba

echo "Arriba for $input1 Started" >> $out/RNASeq_analysis.log   # Writing in log file
date >> $out/RNASeq_analysis.log                                                # Writing in log file

/home/genomics/software/arriba_v2.4.0/./run_arriba.sh /mnt/database/GRCH38.P14/Transcriptomic_db/Genecode_GRCH38.P13/Genomic_db/indexes_v2 /mnt/database/GRCH38.P14/Transcriptomic_db/Genecode_GRCH38.P13/GTF_file/gencode.v38.annotation.gtf /mnt/database/GRCH38.P14/Transcriptomic_db/Genecode_GRCH38.P13/Genomic_db/GRCh38.p13.genome.fa /home/genomics/software/arriba_v2.4.0/database/blacklist_hg38_GRCh38_v2.4.0.tsv /home/genomics/software/arriba_v2.4.0/database/known_fusions_hg38_GRCh38_v2.4.0.tsv /home/genomics/software/arriba_v2.4.0/database/protein_domains_hg38_GRCh38_v2.4.0.gff3 $tnp   $out/$input1/fastp_output/${sm}_R1-trimmed.fastq.gz  $out/$input1/fastp_output/${sm}_R2-trimmed.fastq.gz

cd $workdir

echo "Arriba for $input1 Started" >> $out/RNASeq_analysis.log   # Writing in log file
date >> $out/RNASeq_analysis.log                                                # Writing in log file
echo "##############################" >> $out/RNASeq_analysis.log               # Writing in log file
