#!/bin/bash
#SBATCH --job-name=getENCODE
#SBATCH --partition=gc
#SBATCH -t 1-00:00
#SBATCH --ntasks=2
#SBATCH --mem-per-cpu=4000

##### get the exon boundaries #####
wget ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_human/release_19/gencode.v19.annotation.gtf.gz
gunzip gencode.v19.annotation.gtf.gz
# get the start, stop, name of transcript, exon number
awk '{if ($3 == "exon"){print $12, $4, $5, $26}}' gencode.v19.annotation.gtf | tr -d '";' > exon_boundaries_per_transcript

#### Get missense ####
# get the transcript and the position of the point mutation
sed 's/:.*//' ../exacdb.nonpsych_missense | awk '{print $6,$2}' > missense_mutations_in_transcript

# get transcript, point mutation position, and exon number
join <(sort -k1,1 missense_mutations_in_transcript) <(sort -k1,1 exon_boundaries_per_transcript) | awk '{if ($2 >= $3 && $2 <= $4) {print $1,$2,$5}}' > transcript_w_missense_mutation_exon

# get the transcript, exon number, and number of missense per exon
python sum_mutations_per_exon.py -i transcript_w_missense_mutation_exon -o summed_missense_per_exon

# write in lines for exons with 0 missense mutations
python write_zeros.py -i summed_missense_per_exon -b exon_boundaries_per_transcript -o summed_missense_per_exon_zeros

# some transcripts only exist in synonymous file. Add these transcripts to the missense file with '0' as missense counts.
comm -23 <(cut -f1 summed_synonymous_per_exon_zeros | sort | uniq) <(cut -f1 summed_missense_per_exon_zeros | sort | uniq) > only_in_synonymous
grep -f only_in_synonymous exon_boundaries_per_transcript | cut -f1,4 -d ' ' | awk '{print $1,"\t",$2,"\t",0}' > ADD_missing_in_missense
cat summed_missense_per_exon_zeros ADD_missing_in_missense | sort -k1,1 -k2,2n > counts_missense_exon

#### Get synonymous ####
# get the transcript and the position of the point mutation
sed 's/:.*//' ../exacdb.nonpsych_synonymous | awk '{print $6,$2}' > synonymous_mutations_in_transcript

# get transcript, point mutation position, and exon number
join <(sort -k1,1 synonymous_mutations_in_transcript) <(sort -k1,1 exon_boundaries_per_transcript) | awk '{if ($2 >= $3 && $2 <= $4) {print $1,$2,$5}}' > transcript_w_synonymous_mutation_exon

# get the transcript, exon number, and number of synonymous per exon
python sum_mutations_per_exon.py -i transcript_w_synonymous_mutation_exon -o summed_synonymous_per_exon

# write in lines for exons with 0 synonymous mutations
python write_zeros.py -i summed_synonymous_per_exon -b exon_boundaries_per_transcript -o summed_synonymous_per_exon_zeros

### Fill in empty transcripts with 0 for both missense and synonymous ###
# some transcripts only exist in synonymous file. Add these transcripts to the missense file with '0' as missense counts.
comm -23 <(cut -f1 summed_synonymous_per_exon_zeros | sort | uniq) <(cut -f1 summed_missense_per_exon_zeros | sort | uniq) > only_in_synonymous
grep -f only_in_synonymous exon_boundaries_per_transcript | cut -f1,4 -d ' ' | awk '{print $1,"\t",$2,"\t",0}' > ADD_missing_in_missense
cat summed_missense_per_exon_zeros ADD_missing_in_missense | sort -k1,1 -k2,2n > counts_missense_exon

# some transcripts only exist in missense file. Add these transcripts to the synonymous file with '0' as synonymous counts.
comm -23 <(cut -f1 summed_missense_per_exon_zeros | sort | uniq) <(cut -f1 summed_synonymous_per_exon_zeros | sort | uniq) > only_in_missense
grep -f only_in_missense exon_boundaries_per_transcript | cut -f1,4 -d ' ' | awk '{print $1,"\t",$2,"\t",0}' > ADD_missing_in_synonymous
cat summed_synonymous_per_exon_zeros ADD_missing_in_synonymous | sort -k1,1 -k2,2n > counts_synonymous_exon

#### paste together total mutation counts, calculate score of (miss/syn) / ((total miss - miss)/(total syn - syn)) ####
paste counts_missense_exon <(cut -f3 counts_synonymous_exon) > counts_miss_and_syn
paste <(cut -f1,3 ../exacdb.nonpsych_missense_total) <(cut -f3 ../exacdb.nonpsych_synonymous_total) > counts_total_mutations
# shows total missense, synonymous per exon and per transcript
join -t $'\t' <(sort -k1,1 counts_miss_and_syn) <(sort -k1,1 counts_total_mutations) > all_mutations
awk '{if($4 == 0 || ($6-$4) == 0 || ($5-$3) == 0){print $0"\tNaN"}else{print $0"\t"($3/$4)/(($5-$3)/($6-$4))}}' all_mutations | sort -k1,1 -k2,2n > scored_exons

#### calculate p-values and adjust ####
Rscript pvalue_exons.R scored_exons pvalue_exons
Rscript fdr_exons.R pvalue_exons fdr_exons
join -t $'\t' -a 1 <(sort -k1,1 fdr_exons) <(sort -k1,1 ../Transcript2GeneName) | sort -k1,1 -k2,2n > fdr_exons_annot
awk '{if($9<=0.05){print $0}}' fdr_exons_annot > fdr_exons_annot_top 

###############################################

#### filter MAF < 0.1 and repeat #####
sed 's/CSQ=/\t/g' ../exacdb.nonpsych | cut -f1,2,4,7,8,9- | awk -F';' '{print $1,$NF}' | sed 's/ //g' | sed 's/,/\t/g' | sed 's/AC=//g' |  perl -ane '@l=split(/\t/); chomp(@l); $i=5; while($i<scalar(@l)){print "$l[0]\t$l[1]\t$l[2]\t$l[3]\t$l[4]\t$l[$i]\n"; $i++;}' | awk '{if ($5 <=123){print $0}}' > exacdb.nonpsych_MAF_filter
grep "missense_variant" exacdb.nonpsych_MAF_filter | awk '{if($3=="A" || $3=="G" || $3=="T" || $3=="C"){print $0}}' | grep "PASS" | sed 's/|/\t/g' | cut -f1,2,3,4,6,35 > exacdb.nonpsych_MAF_missense
grep "synonymous_variant" exacdb.nonpsych_MAF_filter | awk '{if($3=="A" || $3=="G" || $3=="T" || $3=="C"){print $0}}' | grep "PASS" | sed 's/|/\t/g' | cut -f1,2,3,4,6,35 > exacdb.nonpsych_MAF_synonymous

#### Get missense ####
# get the transcript and the position of the point mutation
sed 's/:.*//' exacdb.nonpsych_MAF_missense | awk '{print $6,$2}' > MAF_missense_mutations_in_transcript

# get transcript, point mutation position, and exon number
join <(sort -k1,1 MAF_missense_mutations_in_transcript) <(sort -k1,1 exon_boundaries_per_transcript) | awk '{if ($2 >= $3 && $2 <= $4) {print $1,$2,$5}}' > MAF_transcript_w_missense_mutation_exon

# get the transcript, exon number, and number of missense per exon
python sum_mutations_per_exon.py -i MAF_transcript_w_missense_mutation_exon -o MAF_summed_missense_per_exon

# write in lines for exons with 0 missense mutations
python write_zeros.py -i MAF_summed_missense_per_exon -b exon_boundaries_per_transcript -o MAF_summed_missense_per_exon_zeros

#### Get synonymous ####
# get the transcript and the position of the point mutation
sed 's/:.*//' exacdb.nonpsych_MAF_synonymous | awk '{print $6,$2}' > MAF_synonymous_mutations_in_transcript

# get transcript, point mutation position, and exon number
join <(sort -k1,1 MAF_synonymous_mutations_in_transcript) <(sort -k1,1 exon_boundaries_per_transcript) | awk '{if ($2 >= $3 && $2 <= $4) {print $1,$2,$5}}' > MAF_transcript_w_synonymous_mutation_exon

# get the transcript, exon number, and number of synonymous per exon
python sum_mutations_per_exon.py -i MAF_transcript_w_synonymous_mutation_exon -o MAF_summed_synonymous_per_exon

# write in lines for exons with 0 synonymous mutations
python write_zeros.py -i MAF_summed_synonymous_per_exon -b exon_boundaries_per_transcript -o MAF_summed_synonymous_per_exon_zeros

### Fill in empty transcripts with 0 for both missense and synonymous ###
# some transcripts only exist in synonymous file. Add these transcripts to the missense file with '0' as missense counts.
comm -23 <(cut -f1 MAF_summed_synonymous_per_exon_zeros | sort | uniq) <(cut -f1 MAF_summed_missense_per_exon_zeros | sort | uniq) > MAF_only_in_synonymous
grep -f MAF_only_in_synonymous exon_boundaries_per_transcript | cut -f1,4 -d ' ' | awk '{print $1,"\t",$2,"\t",0}' > MAF_ADD_missing_in_missense
cat MAF_summed_missense_per_exon_zeros MAF_ADD_missing_in_missense | sort -k1,1 -k2,2n > MAF_counts_missense_exon

# some transcripts only exist in missense file. Add these transcripts to the synonymous file with '0' as synonymous counts.
comm -23 <(cut -f1 MAF_summed_missense_per_exon_zeros | sort | uniq) <(cut -f1 MAF_summed_synonymous_per_exon_zeros | sort | uniq) > MAF_only_in_missense
grep -f MAF_only_in_missense exon_boundaries_per_transcript | cut -f1,4 -d ' ' | awk '{print $1,"\t",$2,"\t",0}' > MAF_ADD_missing_in_synonymous
cat MAF_summed_synonymous_per_exon_zeros MAF_ADD_missing_in_synonymous | sort -k1,1 -k2,2n > MAF_counts_synonymous_exon

#### paste together total mutation counts, calculate score of (miss/syn) / ((total miss - miss)/(total syn - syn)) ####
paste MAF_counts_missense_exon <(cut -f3 MAF_counts_synonymous_exon) > MAF_counts_miss_and_syn
join -t $'\t' <(sort -k1,1 MAF_counts_miss_and_syn) <(sort -k1,1 counts_total_mutations) > MAF_all_mutations
awk '{if($4 == 0 || ($6-$4) == 0 || ($5-$3) == 0){print $0"\tNaN"}else{print $0"\t"($3/$4)/(($5-$3)/($6-$4))}}' MAF_all_mutations | sort -k1,1 -k2,2n > MAF_scored_exons

#### calculate p-values and adjust ####
Rscript pvalue_exons.R MAF_scored_exons MAF_pvalue_exons
Rscript fdr_exons.R MAF_pvalue_exons MAF_fdr_exons
join -t $'\t' -a 1 <(sort -k1,1 MAF_fdr_exons) <(sort -k1,1 ../Transcript2GeneName) | sort -k1,1 -k2,2n > MAF_fdr_exons_annot
awk '{if($9<=0.05){print $0}}' MAF_fdr_exons_annot > MAF_fdr_exons_annot_top 

#### get regions found by Samocha but NOT by MAF ####
# chi-sq > 10.8 being significant according to Samocha
# 148353-3.csv was retrieved from https://www.biorxiv.org/content/early/2017/06/12/148353.figures-only
# header is transcript, gene, start, stop, observed expected missense, chi-sq
cut -d, -f1,2,5,6,7,8,10 148353-3.csv | awk -F, '{if ($7>=10.8) {print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7}}' > chisq_gt_10.8
# getting exon boundaries only for those transcripts relevant to Samocha
grep -f <(cut -f1 chisq_gt_10.8) exon_boundaries_per_transcript > exons_chisq_gt_10.8
# output those specific exons of transcripts in samocha with chi-sq > 10.8
python get_exon_ranges_in_samocha.py -i chisq_gt_10.8 -b exons_chisq_gt_10.8 -o exon_ranges_in_samocha
sort -k1,1 exon_ranges_in_samocha > exon_ranges_in_samocha_1
rm exon_ranges_in_samocha
mv exon_ranges_in_samocha_1 exon_ranges_in_samocha
# compress the transcript and exon into a single column
awk '{print $1"_"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10}' MAF_fdr_exons_annot > MAF_fdr_exons_annot_combined_transcript_exon
# show the sig regions of Samocha in context of all transcripts available after MAF filtering
join -t $'\t' -a1 -a2 <(sort -k1,1 MAF_fdr_exons_annot_combined_transcript_exon) <(sort -k1,1 exon_ranges_in_samocha)| sed 's/_/\'$'\t/' | sort -k1,1 -k2,2n  > show_sig_regions_of_samocha

#### get pLI of genes of what we picked vs Samocha ####
# Note that some genes do not have pLI values, which is why genes_in_top_MAF_not_Samocha is larger than genes_in_top_MAF_not_Samocha_w_pli
cut -f10 MAF_fdr_exons_annot_top | sort | uniq > MAF_fdr_exons_annot_top_genes
cut -d, -f2 148353-3.csv | sort | uniq > samocha_gene_2700_sorted.txt 
comm -12 <(sort MAF_fdr_exons_annot_top_genes | uniq) <(sort samocha_gene_2700_sorted.txt | uniq) > genes_common_Samocha_MAF_top
comm -23 <(sort MAF_fdr_exons_annot_top_genes | uniq) <(sort samocha_gene_2700_sorted.txt | uniq) > genes_in_top_MAF_not_Samocha
comm -23 <(sort samocha_gene_2700_sorted.txt | uniq) <(sort MAF_fdr_exons_annot_top_genes | uniq) > genes_in_Samocha_not_top_MAF
join -t $'\t' <(sort -k1,1 genes_common_Samocha_MAF_top) <(sort -k1,1 ../Gene2pLI) | awk '{print $0"\tS,Us"}' > genes_common_Samocha_top_MAF_w_pli
join -t $'\t' <(sort -k1,1 genes_in_top_MAF_not_Samocha) <(sort -k1,1 ../Gene2pLI) | awk '{print $0"\tUs"}' > genes_in_top_MAF_not_Samocha_w_pli
join -t $'\t' <(sort -k1,1 genes_in_Samocha_not_top_MAF) <(sort -k1,1 ../Gene2pLI) | awk '{print $0"\tS"}' > genes_in_Samocha_not_top_MAF_w_pli
cat genes_common_Samocha_MAF_w_pli genes_in_top_MAF_not_Samocha_w_pli genes_in_Samocha_not_top_MAF_w_pli | sort -k1,1 > genes_in_Samocha_or_top_MAF_w_pli
awk '{if ($2 >= 0.90) {print $0}}' genes_in_Samocha_or_top_MAF_w_pli > genes_in_Samocha_or_top_MAF_w_0.90_pli

# Just get pLI of MAF_fdr_exons_annot_top
join -t $'\t' -a 1 -1 10 -2 1 <(sort -k10,10 MAF_fdr_exons_annot_top) <(sort -k1,1 ../Gene2pLI) | sort -k2,2 -k3,3n | sed 's/ //g' > MAF_fdr_exons_annot_top_pli
awk '{if ($11 > 0.90) {print $0}}' MAF_fdr_exons_annot_top_pli > MAF_fdr_exons_annot_top_pli_90
cut -f1 MAF_fdr_exons_annot_top_pli_90 | sort | uniq > MAF_fdr_exons_annot_top_pli_90_genes

# Just get pLI for samocha's tops (chi-squared >=10.8)
join -t $'\t' -a 1 -1 2 -2 1 <(sort -k2,2 chisq_gt_10.8) <(sort -k1,1 ../Gene2pLI) | sort -k2,2 -k3,3 > samocha_gt_10.8_pli
awk '{if ($8 > 0.90) {print $0}}' samocha_gt_10.8_pli > samocha_gt_10.8_pli_90
cut -f1 samocha_gt_10.8_pli_90 | sort | uniq > samocha_gt_10.8_pli_90_genes

# Just get pLI for sliding window
awk '{if ($14 > 0.90) {print $0}}' ../exacdb.nonpsych_merged_pli_top | cut -f1 | sort | uniq >  exacdb.nonpsych_merged_top_pli_gt90_genes

#### Kruskal-wallis + Wilcoxon rank sum. average score of Samocha-identified regions comparison ####
cut -d, -f1 148353-3.csv | sort | uniq | tail -n +2 > samocha_transcripts
grep -f samocha_transcripts show_sig_regions_of_samocha | cut -f1,2,7,10,13 | awk '{if (NF==5) {print $1,$2,$3,$4,$5} else {print $1,$2,$3,$4}}' | sort -k1,1 -k2,2n -k5,5 > samocha_transcripts_in_context
python get_score_vectors.py -i samocha_transcripts_in_context -o scores_per_samocha_region_context -n scores_per_samocha_region_context_no_nan
python kruskal_wallis_format.py -i scores_per_samocha_region_context -o kruskal_format
sed -i '/ENST00000326558.5/d' kruskal_format_1 # remove that transcript that had two regions with exactly same scores
Rscript kruskal_wallis.R kruskal_format kruskal_pvalue
awk '{if ($7 < 0.05) {print $0}}' kruskal_pvalue > kruskal_top
python wilcoxon_format.py -i kruskal_top -o wilcoxon_format
Rscript wilcoxon.R wilcoxon_format wilcoxon_pvalue
Rscript wilcoxon_fdr.R wilcoxon_pvalue wilcoxon_fdr
awk '{if ($6 < 0.05) {print $0}}' wilcoxon_fdr > wilcoxon_top
cut -f2 wilcoxon_top | sort | uniq > wilcoxon_top_genes

#### MISSENSE: denovodb de novo mutation count per exon, for top MAF only ####
python exon_boundary_from_1.py -i MAF_fdr_exons_annot_top_genes -o exon_boundary_from_1
join -t $'\t' <(sort -k1,1 <(awk '{print $1"_"$4"\t"$2"\t"$3}' exon_boundary_from_1)) <(sort -k1,1 <(awk '{print $1"_"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10}' MAF_fdr_exons_annot_top))| sed 's/_/\'$'\t/g' | awk '{print $12"\t"$1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10"\t"$11}' > MAF_fdr_exons_annot_top_from_1
python count_mutations_denovo.py --regions MAF_fdr_exons_annot_top_from_1 --denovo_counts ../denovodb.missense_diseases_counts --output MAF_fdr_exons_annot_top_denovo_counts
awk '{if ($13 > 0) {print $0}}' MAF_fdr_exons_annot_top_denovo_counts > MAF_fdr_exons_annot_top_denovo_counts_gt0
# attach the pLI
join -t $'\t' -a 1 -1 1 -2 1 <(sort -k1,1 MAF_fdr_exons_annot_top_denovo_counts_gt0) <(sort -k1,1 ../Gene2pLI) | sort -k2,2 -k3,3n > MAF_fdr_exons_annot_top_denovo_counts_gt0_pli
awk '{if ( $14 > 0.90) {print $0}}' MAF_fdr_exons_annot_top_denovo_counts_gt0_pli > MAF_fdr_exons_annot_top_denovo_counts_gt0_pli_90
cut -f1 MAF_fdr_exons_annot_top_denovo_counts_gt0_pli_90 | sort | uniq > MAF_fdr_exons_annot_top_denovo_counts_gt0_pli_90_genes

#### MISSENSE: denovodb de novo mutation count per region, for S ####
python exon_boundary_from_1_for_S.py -i exons_chisq_gt_10 -c chisq_gt_10 -o samocha_regions_from_1
python count_mutations_denovo.py --regions samocha_regions_from_1 --denovo_counts ../denovodb.missense_diseases_counts --output samocha_regions_denovo_counts
awk '{if ($8 > 0) {print $0}}' samocha_regions_denovo_counts > samocha_regions_denovo_counts_gt0
# attach the pLI
join -t $'\t' -a 1 -1 1 -2 1 <(sort -k1,1 samocha_regions_denovo_counts_gt0) <(sort -k1,1 ../Gene2pLI) | sort -k2,2 -k3,3n > samocha_regions_denovo_counts_gt0_pli
awk '{if ( $9 > 0.90) {print $0}}' samocha_regions_denovo_counts_gt0_pli > samocha_regions_denovo_counts_gt0_pli_90
cut -f1 samocha_regions_denovo_counts_gt0_pli_90 | sort | uniq > samocha_regions_denovo_counts_gt0_pli_90_genes

#### MISSENSE: denovodb de novo mutation in missense depleted region, from sliding window
awk '{if ($14 > 0) {print $0}}' ../exacdb.nonpsych_merged_pli_top_denovo > sliding_window_merged_top_denovo_counts_gt0
# attach the pLI
join -t $'\t' -a 1 -1 1 -2 1 <(sort -k1,1 sliding_window_merged_top_denovo_counts_gt0) <(sort -k1,1 ../Gene2pLI) | sort -k2,2 -k3,3n > sliding_window_merged_top_denovo_counts_gt0_pli
awk '{if ( $15 > 0.90) {print $0}}' sliding_window_merged_top_denovo_counts_gt0_pli > sliding_window_merged_top_denovo_counts_gt0_pli_90
cut -f1 sliding_window_merged_top_denovo_counts_gt0_pli_90 | sort | uniq > sliding_window_merged_top_denovo_counts_gt0_pli_90_genes

#### SYNONYMOUS: denovodb de novo mutation count per exon, for top MAF only ####
grep -w "synonymous" ../denovo-db.variants.v.1.5.tsv | grep -v -w "control" | cut -f19,23 | sed 's/c\.//g' | sed 's/.>/\t&/g' | sort -u > denovodb.synonymous_diseases
join -t $'\t' -1 1 -2 2 <(sort -k1,1 denovodb.synonymous_diseases) <(sort -k2,2 ../Ensembl2RefSeq)| cut -f2- | sort -u | awk '{print $3"\t"$1}' | grep -v -w "NA" | sort | uniq -c | sed 's/^ *//' | sed 's/ /\t/' | awk '{print $2"\t"$3"\t"$1}' > denovodb.synonymous_diseases_counts # Warning: RefSeq conversion file it is not including all de novo mutations in de novo db 
python count_mutations_denovo.py --regions MAF_fdr_exons_annot_top_from_1 --denovo_counts denovodb.synonymous_diseases_counts --output MAF_fdr_exons_annot_top_denovo_counts_synonymous
awk '{if ($13 > 0) {print $0}}' MAF_fdr_exons_annot_top_denovo_counts_synonymous > MAF_fdr_exons_annot_top_denovo_counts_synonymous_gt0

#### SYNONYMOUS: denovodb de novo mutation count per window, for top sliding window ####
python ../count_mutations_denovo.py --regions ../exacdb.nonpsych_merged_pli_top --denovo_counts denovodb.synonymous_diseases_counts --output exacdb.nonpsych_merged_pli_top_denovo_synonymous
awk '{if ($14 > 0) {print $0}}' exacdb.nonpsych_merged_pli_top_denovo_synonymous > exacdb.nonpsych_merged_pli_top_denovo_synonymous_gt0

#### significance of protein protein interactions ####
python get_ppi_degrees.py -p /share/hormozdiarilab/Data/Networks/PPI/StringNew_HPRD -a MAF_fdr_exons_annot_top_genes -b samocha_gene_2700_sorted -x maf_ppis.csv -y sam_ppis.csv -z maf_split_sam_ppis.csv
python get_ppi_degrees.py -p /share/hormozdiarilab/Data/Networks/PPI/StringNew_HPRD -a MAF_fdr_exons_annot_top_genes -b sliding_window_top_genes -x maf_ppis.csv -y sw_ppis.csv -z maf_split_sw_ppis.csv
python get_ppi_degrees.py -p /share/hormozdiarilab/Data/Networks/PPI/StringNew_HPRD -a samocha_gene_2700_sorted -b sliding_window_top_genes -x sam_ppis.csv -y sw_ppis.csv -z sam_split_sw_ppis.csv

#### uncorrected or corrected p-value MAF exon boundary approach vs. Samocha: 4 comparison groups ####

# attach the gene name to the transcript, pLI
join -t $'\t' -a 1 -1 1 -2 1  <(join -t $'\t' -a 1 <(sort -k1,1 counts_total_mutations) <(sort -k1,1 ../Transcript2GeneName) | awk '{if ($3 == 0) {print $4"\t"$1"\t"$2"\t"$3"\tNaN"} else {print $4"\t"$1"\t"$2"\t"$3"\t"$2/$3}}' | sort -k1,1) <(sort -k1,1 ../Gene2pLI) > mis_syn_counts_pli_per_all_20k_genes
# incorporating score into the 4 comparison groups #
python get_avg_score_per_transcript.py -i scored_exons -o avg_scores_per_transcript
join -t $'\t' -1 1 -2 2  <(sort -k1,1 avg_scores_per_transcript) <(sort -k2,2 mis_syn_counts_pli_per_all_20k_genes) | awk '{print $3"\t"$1"\t"$4"\t"$5"\t"$6"\t"$2"\t"$7}' | sort -k1,1 > mis_syn_counts_pli_score_per_all_20k_genes
# get the average of missense, synonymous, m/s per gene
python get_average_mis_syn_per_gene.py -i mis_syn_counts_pli_score_per_all_20k_genes -o avg_mis_syn_counts_pli_score_per_all_20k_genes 
sort -k1,1 avg_mis_syn_counts_pli_score_per_all_20k_genes > avg_mis_syn_counts_pli_score_per_all_20k_genes_sorted

# uncorrected
awk '{if ($8 < 0.05) {print $0}}' MAF_fdr_exons_annot > MAF_UNCORRECTED_annot_top
cut -f10 MAF_UNCORRECTED_annot_top > MAF_UNCORRECTED_annot_top_genes
comm -23 <(sort MAF_UNCORRECTED_annot_top_genes) <(sort samocha_gene_2700_sorted) > only_in_uncorr_not_sam_3833_genes
comm -23 <(sort samocha_gene_2700_sorted) <(sort MAF_UNCORRECTED_annot_top_genes) > only_in_sam_not_uncorr_1107_genes
comm -12 <(sort MAF_UNCORRECTED_annot_top_genes) <(sort samocha_gene_2700_sorted) > both_in_sam_and_uncorr_1593_genes
# get the genes not found by any approach
comm -23 <(cut -f1 avg_mis_syn_counts_pli_score_per_all_20k_genes_sorted | sort | uniq) <(cat only_in_uncorr_not_sam_3833_genes only_in_sam_not_uncorr_1107_genes both_in_sam_and_uncorr_1593_genes | sort | uniq) > NOT_FOUND_by_any_13687_genes
#uncorrected MAF < 0.05. retrieve lines from avg_mis_syn_counts_pli_score_per_all_20k_genes_sorted contaning only the 3833, 1107, 1593 genes
join -t $'\t' <(sort -k1,1 only_in_uncorr_not_sam_3833_genes) <(sort -k1,1 avg_mis_syn_counts_pli_score_per_all_20k_genes_sorted) > avg_mis_syn_counts_score_pli_only_in_uncorr_not_sam_3833
join -t $'\t' <(sort -k1,1 only_in_sam_not_uncorr_1107_genes) <(sort -k1,1 avg_mis_syn_counts_pli_score_per_all_20k_genes_sorted) > avg_mis_syn_counts_score_pli_only_in_sam_not_uncorr_1107  # PHF15 and PHF16 do not exist according to our exacdb filtering. they are dropped
join -t $'\t' <(sort -k1,1 both_in_sam_and_uncorr_1593_genes) <(sort -k1,1 avg_mis_syn_counts_pli_score_per_all_20k_genes_sorted) > avg_mis_syn_counts_score_pli_both_in_sam_and_uncorr_1593  
join -t $'\t' <(sort -k1,1 NOT_FOUND_by_any_13687_genes) <(sort -k1,1 avg_mis_syn_counts_pli_score_per_all_20k_genes_sorted) > avg_NOT_FOUND_score_by_any_13687_genes
# by gene, get the total average #miss, #syn, m/s, pLI for above files
python get_total_avg_mis_syn_pli_group_genes.py -i avg_mis_syn_counts_score_pli_only_in_uncorr_not_sam_3833
python get_total_avg_mis_syn_pli_group_genes.py -i avg_mis_syn_counts_score_pli_only_in_sam_not_uncorr_1107
python get_total_avg_mis_syn_pli_group_genes.py -i avg_mis_syn_counts_score_pli_both_in_sam_and_uncorr_1593
python get_total_avg_mis_syn_pli_group_genes.py -i avg_NOT_FOUND_score_by_any_13687_genes

#### corrected p-value MAF exon boundary approach vs. Samocha: 4 comparison groups ####
comm -23 <(sort MAF_fdr_exons_annot_top_genes) <(sort samocha_gene_2700_sorted) > only_in_corr_not_sam_494_genes
comm -23 <(sort samocha_gene_2700_sorted) <(sort MAF_fdr_exons_annot_top_genes) > only_in_sam_not_corr_2322_genes
comm -12 <(sort MAF_fdr_exons_annot_top_genes) <(sort samocha_gene_2700_sorted) > both_in_sam_and_corr_378_genes
comm -23 <(cut -f1 avg_mis_syn_counts_pli_score_per_all_20k_genes_sorted | sort | uniq) <(cat only_in_corr_not_sam_494_genes only_in_sam_not_corr_2322_genes both_in_sam_and_corr_378_genes | sort | uniq) > NOT_FOUND_by_any_corr_17026_genes
#corrected MAF < 0.05 retrieve lines from avg_mis_syn_counts_pli_score_per_all_20k_genes_sorted containing only the 494, 2322, 378 genes
join -t $'\t' <(sort -k1,1 only_in_corr_not_sam_494_genes) <(sort -k1,1 avg_mis_syn_counts_pli_score_per_all_20k_genes_sorted) > avg_mis_syn_counts_score_pli_only_in_corr_not_sam_494
join -t $'\t' <(sort -k1,1 only_in_sam_not_corr_2322_genes) <(sort -k1,1 avg_mis_syn_counts_pli_score_per_all_20k_genes_sorted) > avg_mis_syn_counts_score_pli_only_in_sam_not_corr_2322  # JHDM1D, PHF15, PHF16, and CXorf61 do not exist according to our exacdb filtering. they are dropped
join -t $'\t' <(sort -k1,1 both_in_sam_and_corr_378_genes) <(sort -k1,1 avg_mis_syn_counts_pli_score_per_all_20k_genes_sorted) > avg_mis_syn_counts_score_pli_both_in_sam_and_corr_378  
join -t $'\t' <(sort -k1,1 NOT_FOUND_by_any_corr_17026_genes) <(sort -k1,1 avg_mis_syn_counts_pli_score_per_all_20k_genes_sorted) > avg_NOT_FOUND_by_score_any_corr_17026_genes
# by gene, get the total average #miss, #syn, m/s, pLI for above files
python get_total_avg_mis_syn_pli_group_genes.py -i avg_mis_syn_counts_score_pli_only_in_corr_not_sam_494
python get_total_avg_mis_syn_pli_group_genes.py -i avg_mis_syn_counts_score_pli_only_in_sam_not_corr_2322
python get_total_avg_mis_syn_pli_group_genes.py -i avg_mis_syn_counts_score_pli_both_in_sam_and_corr_378
python get_total_avg_mis_syn_pli_group_genes.py -i avg_NOT_FOUND_by_score_any_corr_17026_genes

#### looking closely at regions where Samocha maybe wasn't right ####
# find regions in context that we found, but S did not 
awk '{if ($9 < 0.05 && $11 == "") {print $0}}' show_sig_regions_of_samocha | cut -f1 | sort | uniq > show_sig_regions_only_us_1536_transcripts
join -t $'\t' -1 1 -2 1 <(sort -k1,1 show_sig_regions_only_us_1536_transcripts) <(sort -k1,1 show_sig_regions_of_samocha) | sort -k10,10 -k1,1 -k2,2n | awk '{print $10"\t"$1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$11"\t"$12"\t"$13}' > show_sig_regions_only_us_in_context

# find regions in context that S found, but we did not
awk '{if($9 > 0.05 && $11 != "") {print $0}}' show_sig_regions_of_samocha | cut -f1 | sort | uniq > show_sig_regions_only_S_2627_transcripts
join -t $'\t' -1 1 -2 1 <(sort -k1,1 show_sig_regions_only_S_2627_transcripts) <(sort -k1,1 show_sig_regions_of_samocha) | sort -k10,10 -k1,1 -k2,2n | awk '{print $10"\t"$1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$11"\t"$12"\t"$13}' > show_sig_regions_only_S_in_context

# randomly choose 5 genes from show_sig_regions_only_us_in_context   
cut -f1 show_sig_regions_only_us_in_context | sort | uniq | sort -R | head -n 5
cut -f1 show_sig_regions_only_S_in_context | sort | uniq | sort -R | head -n 5

python get_total_avg_mis_syn_pli_group_genes.py -i <(join -t $'\t' <(sort -k1,1 <(cut -f1 show_sig_regions_only_us_in_context | sort | uniq | sort -R | head -n 5 )) <(sort -k1,1 avg_mis_syn_counts_pli_score_per_all_20k_genes_sorted))
python get_total_avg_mis_syn_pli_group_genes.py -i <(join -t $'\t' <(sort -k1,1 <(cut -f1 show_sig_regions_only_S_in_context | sort | uniq | sort -R | head -n 5 )) <(sort -k1,1 avg_mis_syn_counts_pli_score_per_all_20k_genes_sorted))

# use the handpicked 5 genes for averages table
python get_total_avg_mis_syn_pli_group_genes.py -i <(join -t $'\t' <(sort -k1,1 handpicked_only_us) <(sort -k1,1 avg_mis_syn_counts_pli_score_per_all_20k_genes_sorted))
python get_total_avg_mis_syn_pli_group_genes.py -i <(join -t $'\t' <(sort -k1,1 handpicked_only_S) <(sort -k1,1 avg_mis_syn_counts_pli_score_per_all_20k_genes_sorted))

#### Get only the signficant exons that we picked ####
# Get the table of averages for MAF_fdr_exons_annot_top_pli
join -t $'\t' <(sort genes_in_top_MAF_not_Samocha) <(sort -k1,1 MAF_fdr_exons_annot_top_pli) > only_us_top_corrected_MAF_exons
python get_table_avgs_for_only_sig_reg.py -i only_us_top_corrected_MAF_exons

join -t $'\t' -1 10 -2 1 <(sort -k10,10 <(join -t $'\t' -1 1 -2 1 <(sort -k1,1 <(join -t $'\t' -a1 -1 1 -2 1 <(awk '{print $1"_"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7}' scored_exons | sort -k1,1 ) <(sort -k1,1 exon_ranges_in_samocha) | awk '{if ($7 != "") {print $0"\tp\tq"}}' | sed 's/_/\t/g'| cut -f1,2,3,4,5,6,7,11,12)) <(sort -k1,1 ../Transcript2GeneName) | join -t $'\t' -1 1 -2 1 <(sort -k1,1 <(join -t $'\t' -a1 -1 1 -2 1 <(awk '{print $1"_"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7}' scored_exons | sort -k1,1 ) <(sort -k1,1 exon_ranges_in_samocha) | awk '{if ($7 != "") {print $0"\tp\tq"}}' | sed 's/_/\t/g'| cut -f1,2,3,4,5,6,7,11,12)) <(sort -k1,1 ../Transcript2GeneName))) <(sort -k1,1 ../Gene2pLI) > samocha_top_exons
comm -23 <(cut -f2,3 samocha_top_exons | awk '{print $1"_"$2}' | sort) <(cut -f1,2 MAF_fdr_exons_annot_top | awk '{print $1"_"$2}' | sort) > trans_ex_only_in_sam
join -t $'\t' -1 1 -2 2 <(sort -k1,1 trans_ex_only_in_sam) <(sort -k2,2 <(awk '{print $1"\t"$2"_"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10"\t"$11}' samocha_top_exons)) | sort -u | awk '{print $2"\t"$1"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10}' | sed 's/_/\t/g' > only_in_S_not_corr_sig_exons
python get_table_avgs_for_only_sig_reg.py -i only_in_S_not_corr_sig_exons

# Get the table of averages for what we both picked- only sig exons, overlapping both top corrected MAF and top Samocha
join -t $'\t' <(sort genes_common_Samocha_MAF_top) <(sort -k1,1 MAF_fdr_exons_annot_top_pli) > both_us_S_top_corrected_exons
python get_table_avgs_for_only_sig_reg.py -i both_us_S_top_corrected_exons


