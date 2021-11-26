#!/bin/bash

#--------------------------------------------------------------------------------------------------------------------------
# Name		: runLCMFiltering_25112021.sh
# Author	: Mathijs A. Sanders
# Version	: 1.00
# Copyright	: None
# Description	: runLCMFiltering_25112021.sh -pt project_IDs -st samples_tumour_csv -pn project_IDs -sn samples_normal_csv -o output_directory [-t 10] [-f 4] [-r] [-s snp_database] [-g reference_genome] [-h]
#--------------------------------------------------------------------------------------------------------------------------

#--------------------------------------------------------------------------------------------------------------------------
# Annotate variants
#--------------------------------------------------------------------------------------------------------------------------

annotateVariants() {
	
	local PT_DIR=$1
	local ST_ID=$2
	local PN_DIR=$3
	local SN_ID=$4
	local REFERENCE=$5
	local SNP_DB=$6
	local OUTPUT_DIR=$7	
	local THREADS=$8
	
	VCF=$(find /nfs/cancer_ref01/nst_links/live/${PT_DIR}/${ST_ID}/ -mindepth 1 -maxdepth 1 -name "*caveman_c.annot.vcf.gz")

	if [[ -z $VCF ]]
	then
		echo "Error: Could not find the VCF file belonging to sample $ST_ID in project directory $PT_DIR. Skipping this sample."
		return 0
	fi
		
	NM=$(basename $VCF .caveman_c.annot.vcf.gz)_complete.vcf

	(zgrep '^#' $VCF; zgrep -v '^#' $VCF | grep 'PASS' | grep -P 'ASMD=1[4-9][0-9]' | grep -P 'CLPM=0') > ${OUTPUT_DIR}/${NM}

	BAM=${VCF/.caveman_c.annot.vcf.gz/.sample.dupmarked.bam}
	
	if [[ -z $BAM ]]
	then
		echo "Error: Could not find the BAM file belonging to sample $ST_ID in project directory $PT_DIR. Skipping this sample."
		return 0
	fi

	CONTROL_BAM=$(find /nfs/cancer_ref01/nst_links/live/${PN_DIR}/${SN_ID}/ -mindepth 1 -maxdepth 1 -type l -name "*.sample.dupmarked.bam")

	if [[ -z $CONTROL_BAM ]]
	then
		echo "Error: Could not fint the control BAM file belonging to sample $SN_ID in project directory $PN_DIR. Skipping this sample."
		return 0
	fi

	OUTPUT_NM=${NM/.vcf/.txt}
	TNM=${OUTPUT_NM%_complete.txt}
	FNM=${OUTPUT_NM%.txt}

	bsub -q normal -J "${TNM}1" -o ${OUTPUT_DIR}/log/log1.%J -e ${OUTPUT_DIR}/log/err1.%J -n${THREADS} -R"select[mem>10000] span[hosts=1] order[-slots] rusage[mem=10000]" -M10000 bash -e /lustre/scratch116/casm/cgp/users/ms44/scripts/FixedScriptFarm5/annotateVariantsTumorNormal_25112021.sh ${OUTPUT_DIR}/${NM} ${BAM} ${CONTROL_BAM} $THREADS
	
	sleep 5
        cd ${OUTPUT_DIR}
	
	bsub -w "done(${TNM}1)" -q normal -J "${TNM}2" -o ${OUTPUT_DIR}/log_sec/log2.%J -e ${OUTPUT_DIR}/log_sec/err2.%J -n${THREADS} -R"select[mem>20000] span[hosts=1] order[io] rusage[mem=20000]" -M20000 java -Xmx20G -jar /lustre/scratch116/casm/cgp/users/ms44/tools/javaApps/additionalBamStatistics_speedup/additionalBamStatistics_speedup_second.jar --input-annovar-file $OUTPUT_NM --input-bam-file $BAM --output-file ${FNM}_second_filter.txt --reference $REFERENCE --snp-database $SNP_DB --threads $THREADS &

	sleep 5

	bsub -w "done(${TNM}1)&&done(${TNM}2)" -K -q normal -J "${TNM}3" -o ${OUTPUT_DIR}/log_sec/log3.%J -e ${OUTPUT_DIR}/log_sec/err3.%J -n1 bash -e 'sleep 5'

	wait

	return 0
}
export -f annotateVariants

display_usage() {
cat << HEREDOC
Usage: runLCMFiltering_25112021.sh -pt project_tumour_IDS_csv -st sample_tumour_ids_csv -pn project_normal_ids -sn sample_normal_ids -o output_directory [-f 4][-retain] [-g reference_genome]

Required arguments:

	-pt/--projects_tumour <comma_separated_project_IDs>		Comma-separated project IDs for each individual sample (as used in CanApps).
	-st/--samples_tumour <comma_separated_sample_IDs>		Comma-separated sample IDs (required to be present in the matched project directory).
	-pn/--projects_normal <comma_separated_project_IDs>		Comma-separated project IDs for each individual matched control sample (as used in CanApps).
	-sn/--samples_normal <comma_separated_sample_IDs>		Comma-separated sample IDs (required to be present int he matched project directory).
	-o/--output <output_directory>					Output directory for storing all intermediate and final results.
	
Optional arguments:

	-t/--threads <threads>			Number of threads. Default: 10
	-f/--fragment <fragment_threshold>	Fragment count threshold used in the fragment filter. Default: 4
	-r/--retain      			Retain all files up until to a certain point (depending on other provided arguments). Default: Overwrite
	-s/--snp <snp_database>			dbSNP database containing common SNPs. Default: dbSNP147
	-g/--genome <reference_genome>		Reference genome used for aligning and variant calling. Default: hg19
	-h/--help				Get help information.

HEREDOC
}

#------------------------------------------------------------------------------------
# Get input parameters
#------------------------------------------------------------------------------------

THREADS=10
FRAG_THRESHOLD=4
REFERENCE='/lustre/scratch117/core/sciops_repository/references/Homo_sapiens/1000Genomes_hs37d5/all/fasta/hs37d5.fa'
SNP_DB='/lustre/scratch117/casm/team154/ms44/src/dbsnp/common_all_20180423.vcf.gz'

[ ! $(module is-loaded samtools) ] && module load samtools

POSITIONAL=()
while [[ $# -gt 0 ]]
do
	key="$1"

	case $key in
		-pt|--projects_tumour)
		PROJECT_TUMOUR_IDS="$2"
    		shift
    		shift
    		;;
		-st|--samples_tumour)
		SAMPLE_TUMOUR_IDS="$2"
		shift
		shift
		;;
		-pn|--projects_normal)
		PROJECT_NORMAL_IDS=$2
		shift
		shift
		;;
		-sn|--samples_normal)
		SAMPLE_NORMAL_IDS=$2
		shift
		shift
		;;
    		-o|--output)
    		OUTPUT_DIR="$2"
    		shift
    		shift
    		;;
		-t|--threads)
		THREADS="$2"
		shift
		shift
		;;
		-f|--fragment)
		FRAG_THRESHOLD="$2"
		shift
		shift
		;;
    		-r|--retain)
    		RETAIN=YES
    		shift
    		;;
		-s|--snp)
		SNP_DB="$2"
		shift
		shift
		;;
		-g|--genome)
		GENOME="$2"
		shift
		shift
		;;
		-h|--help)
		display_usage
		exit 0
		shift
		;;
    		*)
    		display_usage
		exit -1
    		shift
    		;;
	esac
done
set -- "${POSITIONAL[@]}"

if [[ -z $PROJECT_TUMOUR_IDS ]]
then
	echo "Error: Please provide project directories for the samples of interest"
	exit -1
fi

if [[ -z $SAMPLE_TUMOUR_IDS ]]
then
        echo "Error: Please provide sample ids for the samples of interest matching the project ids listed under -pt"
        exit -1
fi

if [[ -z $PROJECT_NORMAL_IDS ]]
then
        echo "Error: Please provide project directories for the controls"
        exit -1
fi

if [[ -z $SAMPLE_NORMAL_IDS ]]
then
        echo "Error: Please provide sample ids for the controls matching the projects id listed under -pn"
        exit -1
fi


OLD=$IFS
IFS=','
PT_ARRAY=(${PROJECT_TUMOUR_IDS})
ST_ARRAY=(${SAMPLE_TUMOUR_IDS})
PN_ARRAY=(${PROJECT_NORMAL_IDS})
SN_ARRAY=(${SAMPLE_NORMAL_IDS})
IFS=$OLD

for t in "${PT_ARRAY[@]}"
do
	IDIR=/nfs/cancer_ref01/nst_links/live/${t}/
	if [[ ! -d $IDIR ]]
	then
		echo "Error: The provided project - $IDIR - does not exist on the iRods system"
		exit -1
	fi
done

if [[ ${#PT_ARRAY[@]} -ne ${#ST_ARRAY[@]} ]]
then
	echo "Error: The number of project directories and samples of interest are not equal. Please review provided project IDs and samples IDs for cases for which variant filtering is wished"
	exit -1
fi

for ((i=0; i<${#PT_ARRAY[@]}; i++))
do
	IDIR=/nfs/cancer_ref01/nst_links/live/${PT_ARRAY[$i]}/${ST_ARRAY[$i]}
	if [[ ! -d $IDIR ]]
	then
		echo "Error: The provided sample - ${ST_ARRAY[$i]} - does not exist in ${PT_ARRAY[$i]}"
		exit -1
	fi
done

for t in "${PN_ARRAY[@]}"
do
	IDIR=/nfs/cancer_ref01/nst_links/live/${t}/
	if [[ ! -d $IDIR ]]
	then
		echo "Error: The provided directory - $IDIR - does not exist on the iRods system"
		exit -1
	fi
done

if [[ ${#PN_ARRAY[@]} -ne ${#SN_ARRAY[@]} ]]
then
	echo "Error: The number of project directories and control samples are not equal. Please review the provided project IDs and samples IDs for use as controls"
	exit -1
fi

for ((i=0; i<${#PN_ARRAY[@]}; i++))
do
	IDIR=/nfs/cancer_ref01/nst_links/live/${PN_ARRAY[$i]}/${SN_ARRAY[$i]}
	if [[ ! -d $IDIR ]]
	then
		echo "Error: The control sample - ${SN_ARRAY[$i]} - does not exist in the project directory ${PN_ARRAY[$i]}"
		exit -1
	fi
done

if [[ ${#PT_ARRAY[@]} -ne ${#PN_ARRAY[@]} ]]
then
	echo "Error: The number of samples of interest does not match the number of control samples. Please review these variables"
	exit -1
fi

if [[ -z $OUTPUT_DIR ]]
then
	echo "Error: Please provide a output directory (-h for help)"
	exit -1
fi

OUTPUT_DIR="${OUTPUT_DIR/#\~/$HOME}"

regex='^[0-9]+$'

if ! [[ $THREADS =~ $regex ]]
then
	echo "Error: Please provide an integer value for the number of threads"
	exit -1
fi

if ! [[ $FRAG_THRESHOLD =~ $regex ]]
then
	echo "Error: Please provide an integer value for the fragment threshold"
	exit -1
fi

if [[ ! -f $SNP_DB ]]
then
	echo "Error: The provided dbSNP database file does not exist. Please provide the correct path to the dbSNP database."
	exit -1
fi
if [[ ! -f $REFERENCE ]]
then
	echo "Error: The provided reference sequence file does not exist. Please provde the correct path to the reference FASTA file."
	exit -1
fi

if [[ -z $RETAIN ]]
then
	mkdir -p ${OUTPUT_DIR}
	mkdir -p ${OUTPUT_DIR}/log
	mkdir -p ${OUTPUT_DIR}/log_sec

	if [[ ! -d $OUTPUT_DIR ]]
	then
        	echo "Error: Unable to make the directory. Please check permissions."
        	exit -1
	fi
	echo 'Starting jobs...'
	awk -v STR_PT=$PROJECT_TUMOUR_IDS -v STR_ST=$SAMPLE_TUMOUR_IDS -v STR_PN=$PROJECT_NORMAL_IDS -v STR_SN=$SAMPLE_NORMAL_IDS -v OUTPUT_DIR=$OUTPUT_DIR -v REF=$REFERENCE -v SNP=$SNP_DB -v THR=$THREADS 'BEGIN{n=split(STR_PT, apt, ","); split(STR_ST,ast , ","); split(STR_PN, apn, ","); split(STR_SN , asn, ","); for(i=1; i<=n ; i++){print apt[i] " " ast[i] " " apn[i] " " asn[i] " " REF " " SNP " " OUTPUT_DIR " " THR}}' | xargs -n1 -P0 -I {} bash -c 'annotateVariants {}'

	cd ${OUTPUT_DIR}

else
	cd ${OUTPUT_DIR}
fi

for f in *_second_filter.txt; do NM=${f%_second_filter.txt}; (grep '^#' ${NM}.vcf; awk -v ft=$FRAG_THRESHOLD -F $'\t' '(NR==FNR && $63 >= ft && $120 >= ft && (($112 == "NA" && $118 > 2) || ($117 == "NA" && $113 > 2) || (($110 > 1 && $113 > 2) || ($115 > 1 && $118 > 2))) && (($110 <= 1 && $115 > 1 && (($116/$115) <= 0.9 || ($117 > 0 && $118 >= 4))) || ($115 <= 1 && $110 > 1 && (($111/$110)<=0.9 || ($112 > 0 && $113 >= 4))) || ($110 > 1 && $115 > 1 && (($111/$110)<=0.9 || ($110 > 2 && $112 > 2) || ($115 > 1 && $118 > 10)) && (($116/$115)<=0.9 || ($115 > 2 && $117>2)|| ($110 > 1 && $113 > 10))))){a[$1"_"$2"_"$4"_"$5]++} NR != FNR&&a[$1"_"$2"_"$4"_"$5]{print}' <(tail -n+2 $f) <(grep -v '^#' ${NM}.vcf)) > ${NM}_final_retained_${FRAG_THRESHOLD}.vcf; done

for f in *_second_filter.txt; do NM=${f%_second_filter.txt}; (grep '^#' ${NM}.vcf; awk -v ft=$FRAG_THRESHOLD -F $'\t' '(NR==FNR && $63 >= ft && $120 >= ft && (($112 == "NA" && $118 > 2) || ($117 == "NA" && $113 > 2) || (($110 > 1 && $113 > 2) || ($115 > 1 && $118 > 2))) && (($110 <= 1 && $115 > 1 && (($116/$115) <= 0.9 || ($117 > 0 && $118 >= 4))) || ($115 <= 1 && $110 > 1 && (($111/$110)<=0.9 || ($112 > 0 && $113 >= 4))) || ($110 > 1 && $115 > 1 && (($111/$110)<=0.9 || ($110 > 2 && $112 > 2) || ($115 > 1 && $118 > 10)) && (($116/$115)<=0.9 || ($115 > 2 && $117>2)|| ($110 > 1 && $113 > 10))))){a[$1"_"$2"_"$4"_"$5]++} NR != FNR&&!a[$1"_"$2"_"$4"_"$5]{print}' <(tail -n+2 $f) <(grep -v '^#' ${NM}.vcf)) > ${NM}_final_removed_${FRAG_THRESHOLD}.vcf; done

for f in *_second_filter.txt; do NM=${f%_second_filter.txt}; (head -1 $f; awk -v ft=$FRAG_THRESHOLD -v OFS=$'\t' -F$'\t' '(NR > 1 && $63 >= ft && $120 >= ft && (($112 == "NA" && $118 > 2) || ($117 == "NA" && $113 > 2) || (($110 > 1 && $113 > 2) || ($115 > 1 && $118 > 2))) && (($110 <= 1 && $115 > 1 && (($116/$115) <= 0.9 || ($117 > 0 && $118 >= 4))) || ($115 <= 1 && $110 > 1 && (($111/$110)<=0.9 || ($112 > 0 && $113 >= 4))) || ($110 > 1 && $115 > 1 && (($111/$110)<=0.9 || ($110 > 2 && $112 > 2) || ($115 > 1 && $118 > 10)) && (($116/$115)<=0.9 || ($115 > 2 && $117>2)|| ($110 > 1 && $113 > 10)))))' $f) > ${NM}_final_allinfo_${FRAG_THRESHOLD}.txt; done

awk 'BEGIN{print "Sample ID\tChr\tPos\t‘PM-Tum’"}'> ./VAF_longlist_${FRAG_THRESHOLD}.txt

(echo -e 'Sample ID\tChr\tPos\t‘PM-Tum’'; for f in *_final_allinfo_${FRAG_THRESHOLD}.txt; do NM=${f%_complete_final_allinfo_${FRAG_THRESHOLD}.txt}; awk -F$'\t' -v OFS=$'\t' -v x=${NM} 'NR>1{print x,$1,$2,$64}' $f; done) > ./VAF_longlist_${FRAG_THRESHOLD}.txt

exit 0
