set -e

echo ""
echo " --------------------"
echo " | Liftover Plink |"
echo " | Patrick Deelen |"
echo " | UMCG |"
echo " --------------------"
echo ""

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

bed_b38="${output}.tmpB38.bed"
bed_b37="${output}.tmpB37.bed"
unmapped_b37="${output}.tmp.unmapped.txt"
unmapped_b37_2="${output}.tmp.unmapped2.txt"
mappingUpdate="${output}.tmp.mappingUpdate.txt"
output_unSorted="${output}_b37_unsorted.tmp"

#Create bed
awk '{print "chr"$1,$4,$4+1,$2}' OFS="\t" "${input}.bim" > ${bed_b38}

#Get updated mappings
liftOver -bedPlus=4 ${bed_b38} ${chain} ${bed_b37} ${unmapped_b37}

#Parse unmapped SNPs
awk '/^[^#]/ {print $4}' ${unmapped_b37} > ${unmapped_b37_2}

#Create mapping update list used by Plink
awk '{print $4, $2}' OFS="\t" ${bed_b37} > ${mappingUpdate}

#Update plink mappings
plink \
--bfile ${input} \
--make-bed \
--exclude ${unmapped_b37_2} \
--update-map ${mappingUpdate} \
--keep-allele-order \
--allow-no-sex \
--out ${output_unSorted}

#No we have to again create a plink file to make sure the implied order is correct after liftover.
plink --bfile ${output_unSorted} --keep-allele-order --allow-no-sex --make-bed --out ${output}

rm ${bed_b38}
rm ${bed_b37}
rm ${unmapped_b37}
rm ${unmapped_b37_2}
rm ${mappingUpdate}
rm "${output_unSorted}.bed"
rm "${output_unSorted}.bim"
rm "${output_unSorted}.fam"
rm "${output_unSorted}.nosex"
rm "${output_unSorted}.log"
