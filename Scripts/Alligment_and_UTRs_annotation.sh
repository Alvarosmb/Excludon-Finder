cd ..

### Ask for file names

echo Introduce fasta file, please

read var1

echo Introduce fastq file, please

read var2

echo Introduce GFF file, please

read var3
## Replace ID for gene_id if neccessary. gff file from NCBI is annotated as ID but the hole script uses the gene_id annotation
sed -i 's/ID/gene_id/g' $var3

### Concatenate FASTQs
#mkdir -p results/concatenated_fastq

#cat /home/vant/Navarrabiomed/Identify_Operons/Data/fastq_pass/*.fastq.gz >/home/vant/Navarrabiomed/Identify_Operons/results/concatenated_fastq/${Sample_ID}.fastq.gz

####THIS PREVIOUS STEP IS NOT NECCESSARY SINCE THE FASTQ IS PROVIDED. I ONLY USED THE SAMPLE ERR2300649.fastq

string1=$var1
string2=$var2
IFS='/' read -ra split <<< "$string1"
input_file1=${split[-1]}

IFS='/' read -ra split <<< "$string2"
input_file2=${split[-1]}

string3=$input_file1
IFS='.' read -ra split <<< "$string3"
genome=${split[0]}

string4=$input_file2
IFS='.' read -ra split <<< "$string4"
sample=${split[0]}

echo  input files: $input_file2, $input_file2



echo "################ MAKNG ALIGMENT #####################"


mkdir -p results/minimap

minimap2 -p 0.99 -ax splice -k14 --MD -uf \
  $var1 \
  $var2| \
  samtools sort -O bam - > results/minimap/${sample}_sorted.bam
 
samtools index results/minimap/${sample}_sorted.bam

echo " ### DONE ### "


echo "################ GET DEPTH AT EACH POSITION IN THE GENOME #####################"

mkdir -p results/coverage_data


samtools view -F 0x10 -h -b results/minimap/${sample}_sorted.bam | \
  samtools depth -a -d 0 -o results/coverage_data/${sample}_plus_depth.txt -
echo ". . . "

samtools view -f 0x10 -h -b results/minimap/${sample}_sorted.bam | \
  samtools depth -a -d 0 -o results/coverage_data/${sample}_minus_depth.txt -
 
 echo " ### DONE ### "
var4="results/minimap/${sample}_sorted.bam"


echo " ################## GET COUNTS OF MAPPED READS ######################## " 
mkdir -p results/operons

Rscript Scripts/identify_Overlapping_reads_paralelized.R $var1  $var3 $var4  $sample
  
echo "### The analysis is finished ###"

 
