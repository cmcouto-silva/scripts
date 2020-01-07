#!/bin/bash
set -e

echo ""
echo " ---------------------"
echo " | Liftover Plink |"
echo " | Patrick Deelen |"
echo " | UMCG |"
echo " ---------------------"
echo "
 ---------------------
 Modified by:
 Cainã Max Couto-Silva
 ---------------------
"

if [ "$#" -ne 3 ]; then
echo "Usage:"
echo " - Arg 1: input plink data (without extentions)"
echo " - Arg 2: output (without extentions)"
echo " - Arg 3: chain file"
echo ""
exit 1
fi

input=$1
output=$2
chain=$3
echo "Input: ${input}"
echo "Output: ${output}"
echo "Chain: ${chain}"
echo ""

bed_IN="${output}.tmpIN.bed"
bed_OUT="${output}.tmpOUT.bed"
unmapped_OUT="${output}.tmp.unmapped.txt"
unmapped_OUT_2="${output}.tmp.unmapped2.txt"
mappingUpdate="${output}.tmp.mappingUpdate.txt"
output_unSorted="${output}_OUT_unsorted.tmp"

#Create bed
awk '{print "chr"$1,$4,$4+1,$2}' OFS="\t" "${input}.bim" > ${bed_IN}

#Get updated mappings
liftOver -bedPlus=4 ${bed_IN} ${chain} ${bed_OUT} ${unmapped_OUT}

#Parse unmapped SNPs
awk '/^[^#]/ {print $4}' ${unmapped_OUT} > ${unmapped_OUT_2}

#Create mapping update list used by Plink
awk '{print $4, $2}' OFS="\t" ${bed_OUT} > ${mappingUpdate}

#Update plink mappings
plink \
--bfile ${input} \
--make-bed \
--exclude ${unmapped_OUT_2} \
--update-map ${mappingUpdate} \
--keep-allele-order \
--allow-no-sex \
--out ${output_unSorted}

#No we have to again create a plink file to make sure the implied order is correct after liftover.
plink --bfile ${output_unSorted} --keep-allele-order --allow-no-sex --make-bed --out ${output}

rm ${bed_IN}
rm ${bed_OUT}
rm ${unmapped_OUT}
rm ${unmapped_OUT_2}
rm ${mappingUpdate}
rm "${output_unSorted}.bed"
rm "${output_unSorted}.bim"
rm "${output_unSorted}.fam"
rm "${output_unSorted}.nosex"
rm "${output_unSorted}.log"
