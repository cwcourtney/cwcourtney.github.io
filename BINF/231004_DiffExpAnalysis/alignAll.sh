#!/usr/bin/env bash
# alignAll.sh
outDir='quant/'
# sample='Aip02' # TODO: update to loop over all Aip## samples
fastq_path='/work/courses/BINF6309/AiptasiaMiSeq/fastq/'
leftSuffix='.R1.fastq'
rightSuffix='.R2.fastq'
function align {
    for leftInFile in $fastq_path"Aip"*$leftSuffix
    do
        #Remove the path from the filename and assign to pathRemoved
        pathRemoved="${leftInFile/$fastq_path/}"
        #Remove the left-read suffix from $pathRemoved and assign to suffixRemoved
        suffixRemoved="${pathRemoved/$leftSuffix/}"
        #echo $suffixRemoved
        salmon quant -l IU \
        -1 /work/courses/BINF6309/AiptasiaMiSeq/fastq/${suffixRemoved}.R1.fastq \
        -2 /work/courses/BINF6309/AiptasiaMiSeq/fastq/${suffixRemoved}.R2.fastq \
        -i AipIndex \
        --validateMappings \
        -o ${outDir}${suffixRemoved}
    done
}

align 
