ANNOVAR_DATABASE='/lustre/scratch116/casm/cgp/users/ms44/src/hg19/annovar/'
CONVERT_2_ANNOVAR='/lustre/scratch116/casm/cgp/users/ms44/tools/annovar/convert2annovar.pl'
TABLE_ANNOVAR='/lustre/scratch116/casm/cgp/users/ms44/tools/annovar/table_annovar.pl'
VCF_REFORMAT='/lustre/scratch116/casm/cgp/users/ms44/tools/javaApps/SamtoolsVCFReformat/'
ANNOTATE_BAM_STATISTICS='/lustre/scratch116/casm/cgp/users/ms44/tools/annotateBamStatistics8/annotateBamStatistics8'

#---------------------------------------------
# Input arguments
#---------------------------------------------

PREFIX=`echo $1 | sed "s:\.vcf::"`
TUMOR_BAM_FILE=$2
NORMAL_BAM_FILE=$3
THREADS=$4

#---------------------------------------------
# Fix INDEL annotations
#---------------------------------------------

echo -e 'Info: Fix INDEL annotations'
java -Xmx1G -cp $VCF_REFORMAT samtoolsVCFReformat ${PREFIX}.vcf ${PREFIX}.txt

#---------------------------------------------
# Convert VCF to Annovar file
#---------------------------------------------

echo -e 'Info: Convert VCF to Annovar file'

$CONVERT_2_ANNOVAR --format vcf4old ${PREFIX}.txt > ${PREFIX}.annovar

if [[ $? -ne 0 ]]
then
	echo 'Error: could not convert to Annovar file'
	rm ${PREFIX}.annovar
	return -1
fi

#---------------------------------------------
# Annotate the variant file with Annovar
#---------------------------------------------

echo -e 'Info: Annotating the variant file with Annovar'

$TABLE_ANNOVAR ${PREFIX}.annovar $ANNOVAR_DATABASE --remove --buildver hg19 --thread $THREADS --protocol refGene,ljb26_all,popfreq_all_20150413,avsnp147,cosmic70 --operation g,f,f,f,f --nastring ""
if [[ $? -ne 0 ]]
then
	echo -e 'Error: Could not annotate the variant file'
	rm ${PREFIX}.annovar
	return -1
fi

rm ${PREFIX}.annovar
mv ${PREFIX}.annovar.hg19_multianno.txt ${PREFIX}.txt

#---------------------------------------------
# Annotate variant file with BAM statistics
#---------------------------------------------

echo -e 'Info: Annotating the variant file with BAM statistics'

$ANNOTATE_BAM_STATISTICS --threads $THREADS --pileup-regions --min-alignment-score 40 --annovar-file ${PREFIX}.txt --bam-files ${TUMOR_BAM_FILE},${NORMAL_BAM_FILE} > ${PREFIX}.temp

if [[ $? -ne 0 ]]
then
	echo -e 'Error: An error has been encountered with AnnotateBamStatistics'
	#rm ${PREFIX}.temp
	return -1
fi

mv ${PREFIX}.temp ${PREFIX}.txt

exit 0
