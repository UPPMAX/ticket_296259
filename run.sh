#!/bin/bash -l
#SBATCH -A staff 
#SBATCH -t 1:00:00 
#SBATCH -p devcore

module load bioinfo-tools
module load samtools
module load vartrix
module load java/OpenJDK_17+35

# Added by Richel
module load picard

# Added by Richel
echo "Skip the first part"

# Added by Richel
if [[ 1 == 2 ]]; then

    # Commented out by Richel
    # <<'_SKIP_'

    # Original bam
    cp ../8_EXPRESSION.D/output.d/01_out/05.0_rsem/SS2_19_037/H13/SS2_19_037-H13.genome.sorted.bam ./H13/SS2_19_037-H13.genome.sorted-original.bam
    cp ../8_EXPRESSION.D/output.d/01_out/05.0_rsem/SS2_19_037/H13/SS2_19_037-H13.genome.sorted.bam.bai .H13/SS2_19_037-H13.genome.sorted-original.bam
    bamO="SS2_19_037-H13.genome.sorted-original.bam"
    bai="SS2_19_037-H13.genome.sorted-original.bam.bai"

    # Split to chr2
    ## aligment
    samtools view -b ${bamO} chr2 > SS2_19_037-H13_chr2.bam



    ## fasta reference
    fastafile="/castor/project/proj/maria.d/8_EXPRESSION.D/09.0_vartrix/data-d/ref_fasta-d/00.0_chrom_seq_removed_GRCh38.primary_assembly.genome-nochrY_ERCC92.fa"
    cp ${fastafile} .
    samtools faidx 00.0_chrom_seq_removed_GRCh38.primary_assembly.genome-nochrY_ERCC92.fa chr2 >  ./00.0_chrom_seq_removed_GRCh38.primary_assembly.genome-nochrY_ERCC92-chr2.fa
    samtools faidx ./00.0_chrom_seq_removed_GRCh38.primary_assembly.genome-nochrY_ERCC92-chr2.fa
    rm 00.0_chrom_seq_removed_GRCh38.primary_assembly.genome-nochrY_ERCC92.fa.fai

    # Commented out by Richel
    # _SKIP_
fi

# Added by Richel
echo "Add missing tag"

# Add missing tag
tooldir="/home/mararc/bin/jvarkit/dist/"
LINE="./SS2_19_037-H13_chr2.bam"
outdir1="./"
cell="SS2_19_037-H13"

# Added by Richel
echo "Running jvarkit"
if [ -f ${tooldir}jvarkit.jar ]; then

    # Added by Richel
    echo "Running jvarkit in user's way"

    java -jar ${tooldir}jvarkit.jar samjdk -e 'String c=record.getReadName(); int h=0; int s=21; record.setAttribute("CB",c.substring(h,s));return record;' ${LINE}  > ${outdir1}${cell}-CB_chr2.sam
else

    # Added by Richel
    echo "Check if jvarkit can be run in Richel's way"

    if [ ! -f jvarkit.sif ]; then
        echo "ERROR: 'jvarkit.sif' not found. See https://docs.uppmax.uu.se/software/jvarkit how to create it"
        exit 42
    fi

    # Added by Richel
    echo "Running jvarkit in Richel's way"

    ./jvarkit.sif java -jar /opt/jvarkit/dist/jvarkit.jar samjdk -e 'String c=record.getReadName(); int h=0; int s=21; record.setAttribute(\"CB\",c.substring(h,s));return record;' ${LINE}  > ${outdir1}${cell}-CB_chr2.sam

    # Added by Richel
    echo "Done running jvarkit in Richel's way"
  
fi

# Added by Richel
echo "Convert sam to bam"


# convert back from sam to bam -- I met with Richel @UPPMAX-SUPPORT and he pointed out the file jvarkit spits out is a text, he's right. I actually get a sam format and not bam
# Therefore, I will now convert back 
samtools view -bS ./SS2_19_037-H13-CB_chr2.sam > ./SS2_19_037-H13_chr2_bam.bam

# Added by Richel
echo "Indexing"

## and index
samtools index SS2_19_037-H13_chr2_bam.bam

# Added by Richel
echo "Done indexing"

################# Added by Richel, start
echo "----------------------------------------------"
echo "Fix problem here"
# The file, SS2_19_037-H13-CB_chr2.sam, is already ready:
head SS2_19_037-H13-CB_chr2.sam
echo "----------------------------------------------"
################# Added by Richel, end


# Added by Richel
echo "Setting vcffile"

# Run - it works!! The only thing left to do is add the tag to the bam header sth like @CB:ZXXXXX
# reference vcf file
vcffile="./ALK_chr.vcf"

# Added by Richel
echo "Setting fastafile"

# fasta file for the reference genome used
fastafile="./00.0_chrom_seq_removed_GRCh38.primary_assembly.genome-nochrY_ERCC92-chr2.fa"

# Added by Richel
echo "Setting LINE"

LINE="./"

# Added by Richel
echo "Setting cell"

cell="SS2_19_037-H13"

# Added by Richel
echo "Running VarTrix"

vartrix --vcf ${vcffile} --bam ${LINE}${cell}_chr2_bam.bam --fasta ${fastafile} --cell-barcodes ${LINE}barcodes --out-matrix ./out_mat --out-variants ./out_var &> log


samtools view -S -b sample.sam > sample.bam
