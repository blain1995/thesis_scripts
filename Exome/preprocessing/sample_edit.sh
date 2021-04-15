# A script to identify vcf files that require column rearrangement before entry to the bash script for extracting VAFS
# Written by Alex Blain, January 2021
# Run this is the vcf folder of your analysis

mkdir vcf_R

for i in *.vcf;do
	grep  -v '^##' $i > vcf_R/$i;
done
cd vcf_R
mkdir headers

for i in *.vcf;do
	head -n 1 $i >  headers/$i.header;
done

cat headers/*.header | \
awk '$11 ~ /Cons/' | cut -f 10  > samples_to_edit.txt

mkdir original_samples

for i in $(cat samples_to_edit.txt);do
	mv *$i* original_samples;
done
cd original_samples

for i in *.vcf;do
        awk 'BEGIN {FS="\t"; OFS="\t"} {print $1,$2,$3,$4,$5,$6,$7,$8,$9,$11,$10}' $i > temp.vcf && mv temp.vcf ../$i;
done
cd ../


echo "=================="
echo "Format: chr       pos ref_allele  alt_allele      ref_reads       var_reads total_depth allel_freq"

for FILE in *.vcf;do
        B_NAME=`basename $FILE .vcf`
        echo "----------------"
        echo "Extract variants from $B_NAME"
        grep -v "#" $B_NAME.vcf | \
        awk 'BEGIN {OFS="\t"} {print $1,$2,$4,$5,$11}'| \
        awk -F "[\t '':,]" 'BEGIN {OFS="\t"}{print $1,$2,$3,$4,$6,$7,$9,$8*100}' | \
        sed 's/%//g' | \
        sed 's/chr//g' > $B_NAME.tsv
        echo "Pasting header for $B_NAME"
        gsed -i "1s/^/chr\tpos\tref_allele\talt_allele\tref_reads\tvar_reads\ttotal_depth\tallel_freq\n/" $B_NAME.tsv
        echo "Done"
done
