#!/bin/bash
# This script automatically maps sequencing reads (RNAseq, WGS e.t.c.) to a reference genome
# There is no need to download FastQs, omicstreamer tool streams sequence datasets in the cloud
# Dependencies: magicblast

display_help () {

cat << EOF
Usage: ./mapreads.sh [OPTIONS] -f [ACCESSION_FILE]...
./mapreads.sh [-h] [-r REFERENCE_GENOME_FILE] -f [ACCESSION_FILE] [-p PHENOTYPE_NAME] -g [GENOME_REFNAME] -n [CPU_CORES]
e.g. ./mapreads.sh -r data/sarscov2/sarscov2.fasta -f data/covid19/covid19_Acc_List3.txt -p covid19 -g sarscov2 -n 4
Takes Accession IDs listed in a text file and automatically maps to the specified reference genome
Mandatory arguments include -r REFERENCE_GENOME_FILE -f PATH_TO_ACCESSION_FILE
Where:
  -r    typically use -r <path to reference genome file>
  -f    file containing the list of ACCESSIONS
  -p	  phenotype or disease name
  -g	  genome reference name
  -n	  number of threads (integer)
  -h    give this help
Report bugs to <Olaitan I. Awe - laitanawe@gmail.com>
EOF

}

while getopts r:hf:p:g:n: option
do
case "${option}"
in
   h) display_help
      exit
      ;;
   r) REFERENCE=${OPTARG}
      ;;
   f) FILE_DLOAD=$OPTARG
      ;;
   p) PHENO=${OPTARG}
      ;;
   g) REFNAME=$OPTARG
      ;;
   n) NUM_THREADS=$OPTARG
      ;;
   :) printf "missing argument for -%s\n" "$OPTARG" >&2
      echo "$usage" >&2
      exit 1
      ;;
  \?) printf "illegal option: -%s\n" "$OPTARG" >&2
      echo "$usage" >&2
      exit 1
      ;;
  esac
done
shift $((OPTIND - 1))


#Specify the data directory, this will be the parent directory for the read-mapping results
DATADIR="data"
OUTNAME="results"
OUTPUT_DB="${DATADIR}/${REFNAME}/${PHENO}_${REFNAME}_db"
RESULTS_DIR="${DATADIR}/${OUTNAME}"
SAM_DIR="sam"
BAM_DIR="bam"
SORTEDBAM_DIR="sorted_bam"

counturlsam=0
counturlbam=0
counturlsortedbam=0
counturlindexbam=0
refdb_exists=0

# If we do have a reference genome
if [ $REFERENCE == "${DATADIR}/${REFNAME}/${REFNAME}.fasta" ]
  then
     echo "Reference Format Used: FASTA"
elif [ $REFERENCE -eq "${DATADIR}/${REFNAME}/hg38.fasta" ]
  then
     echo "Using Default Reference Format: Hg38"
else echo "Reference Path: " $REFERENCE
fi

if [[ -f "${OUTPUT_DB}.nsq" ]];
  then
     ((refdb_exists = 1))
     echo "We don't need to create a Reference Genome DB. RefDb Exists already!"
else
     echo "REFERENCE GENOME CREATION IS IN PROGRESS"
     makeblastdb -in $REFERENCE -out "${OUTPUT_DB}" -parse_seqids -dbtype nucl
fi

# Align FASTQs to our reference genome
samfiles=(`find "${RESULTS_DIR}/${SAM_DIR}/" -maxdepth 1 -name "*.sam"`)
if [[ "${#samfiles[@]}" == 0 ]];
  then
     echo "SAM files do not exist, magicblast is needed!"
     # Create all SAM results subdirectory
     mkdir -pv "${RESULTS_DIR}/${SAM_DIR}"

     # Pick up all the resulting Run ID's and put it into an array
     run_ids_array=(`awk '{print $1}' ${FILE_DLOAD}`)

     echo "Starting alignment of FASTQ files to reference genome."

     # We have one line for each Accession and we are runing them sequentially
     for run_id in "${run_ids_array[@]}"

       do
       ((counturlsam++))
       echo "Currently processing Accession"$counturlsam ": " $run_id

       magicblast -sra "${run_id}"  -db "${OUTPUT_DB}" -no_discordant \
           -paired -num_threads "${NUM_THREADS}" -no_unaligned -out "${RESULTS_DIR}/${SAM_DIR}/${run_id}_${REFNAME}.sam"
       done

       echo "Completed alignment map creation successfully."

else
       echo ${#samfiles[@]}" SAM files exist, we do not need to generate SAM files!"
       for index in ${!samfiles[@]}; do
         echo "SAM file "$((index+1))/${#samfiles[@]} = "${samfiles[index]}"
       done
fi

bamfiles=(`find "${RESULTS_DIR}/${BAM_DIR}/" -maxdepth 1 -name "*.bam"`)
if [[ "${#bamfiles[@]}" == 0 ]];
  then
     echo "BAM files do not exist, conversion of SAM to BAM is needed!"
     # Create all BAM results sub-directory
     mkdir -pv "${RESULTS_DIR}/${BAM_DIR}"

     # We have one line for each Accession and we are doing the SAM to BAM conversion, sorting and indexing sequentially
     # Pick up all the resulting text-based alignment maps from the SAM directory and put it into an array

     run_ids_array_aligned=(`ls ${RESULTS_DIR}/${SAM_DIR}`)
     for run_id_aligned in "${run_ids_array_aligned[@]}"

       do
       ((counturlbam++))
       echo "Starting Conversion of SAM to BAM."
       echo "Currently generating BAM"$counturlbam "from :" $run_id_aligned

       samtools view -S -b "${RESULTS_DIR}/${SAM_DIR}/"$run_id_aligned > "${RESULTS_DIR}/${BAM_DIR}/${run_id_aligned}.bam"

       echo "Completed generation of BAM file"$counturlbam " successfully."

     done

else
       echo ${#bamfiles[@]}" BAM files exist, conversion of SAM to BAM is not needed!"
       for index in ${!bamfiles[@]}; do
         echo "BAM file "$((index+1))/${#bamfiles[@]} = "${bamfiles[index]}"
       done
fi


## To do: Fix automated sorting and indexing from this script. I need some sleep so as to be productive next day.
sortedbamfiles=(`find "${RESULTS_DIR}/${SORTEDBAM_DIR}/" -maxdepth 1 -name "*.bam"`)
if [[ "${#sortedbamfiles[@]}" == 0 ]];
  then

     echo "sorted BAM files do not exist, sorting of BAM file is needed!"

     # Create all sorted BAM results sub-directory
     mkdir -pv "${RESULTS_DIR}/${SORTEDBAM_DIR}"

     # Pick up all the resulting Run binary alignment maps from the BAM directory, put it into an array and sort them
     run_ids_array_aligned=(`ls ${RESULTS_DIR}/${BAM_DIR}`)
     for run_id_aligned in "${run_ids_array_aligned[@]}"

       do
       ((counturlsortedbam++))

       echo "Starting the sorting of Alignment Map in genome order."
       echo "Currently sorting BAM"$counturlsortedbam ": "$run_id_aligned

       echo "${RESULTS_DIR}/${BAM_DIR}/"$run_id_aligned

       samtools sort "${RESULTS_DIR}/${BAM_DIR}/"$run_id_aligned -o "${RESULTS_DIR}/${SORTEDBAM_DIR}/${run_id_aligned}.sorted.bam"

       echo "Completed the sorting of Alignment Map in genome order successfully!"

     done

else
       echo ${#sortedbamfiles[@]}" sorted BAM files exist, sorting of BAM files not needed!"
       for index in ${!sortedbamfiles[@]}; do
         echo "Sorted BAM file "$((index+1))/${#sortedbamfiles[@]} = "${sortedbamfiles[index]}"
       done
fi

indexbamfiles=(`find "${RESULTS_DIR}/${SORTEDBAM_DIR}/" -maxdepth 1 -name "*.bam.bai"`)
if [[ "${#indexbamfiles[@]}" == 0 ]];
  then

     echo "indexed BAM files do not exist, BAM file indexing is needed!"

     # Pick up all the sorted binary alignment maps from the BAM directory, put it into an array and sort them
     run_ids_array_aligned=(`ls ${RESULTS_DIR}/${SORTEDBAM_DIR}`)
     for run_id_aligned in "${run_ids_array_aligned[@]}"

       do
       ((counturlindexbam++))

       echo "Starting the Indexing of genome sorted BAM file, so we can quickly extract alignments overlapping genomic regions"
       echo "Currently indexing SORTED-BAM"$counturlindexbam ": "$run_id_aligned
       samtools index "${RESULTS_DIR}/${SORTEDBAM_DIR}/"${run_id_aligned}

       echo "Completed the Indexing of genome sorted BAM file"

     done

else
       echo ${#indexbamfiles[@]}" indexed BAM files exist. Indexing of sorted BAM files not needed!"
       for index in ${!indexbamfiles[@]}; do
         echo "index for BAM file "$((index+1))/${#indexbamfiles[@]} = "${indexbamfiles[index]}"
       done

       echo "SAM, BAM, sorted BAM and indexing files are in their respective directories!"

fi

echo "Processing of datasets fully completed successfully!"
