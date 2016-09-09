#!/bin/bash
#single end 300 454 reads, converted from sff to fastq - put reads in /sff folder:
python ~/s/convert_folder.py reads/sff sff fasta
#put all reads onto single file
cat reads/*.fastq > reads.fastq
#remove 4bp leader from 454 reads
python clipleft.py reads.fastq 4 fastq && mv reads.fastq.clip reads.fastq
#convert reads to upper case so that barcodes can be matched
upper.py reads.fastq fastq
mv upper.fastq reads.fastq
#relabel reads by barcode label in barcodes.fa

#Uparse pipeline http://drive5.com/usearch/manual/upp_454.html
#Usearch version used for following: usearch v8.1.1861_i86linux64
python fastq_strip_barcode_relabel2.py reads.fastq CATGCTGCCTCCCGTAGGAGT barcodes.fa fixed > fixed.fastq
#examine read lengths to find reasonable trim length
python read-histogram.py reads.fastq fastq 100
#chose 300bp
#trim reads
usearch -fastq_filter fixed.fastq -fastq_trunclen 300 -fastaout trimmed.fa
#quality filter and trim to separate file
usearch -fastq_filter fixed.fastq -fastq_maxee 1.0 -fastq_trunclen 300 -fastaout filtered.fa
#cluster identical reads
usearch -derep_fulllength filtered.fa -sizeout -fastaout uniques.fa
#cluster otus to 3% radius, discarding singletons
usearch -cluster_otus uniques.fa -minsize 2 -otus otus.fa -relabel Otu
#assign taxonomy using RDP database, http://drive5.com/utax/data/utax_rdp_16s_tainset15.tar.gz
#make udb if required, uncomment below:
#usearch -makeudb_utax ~/db/utax/16s_ref.fa -output ~/db/utax/16s.udb
usearch -utax otus.fa -db ~/db/utax/16s.udb -strand both -strand both -fastaout otus_tax.fa -utax_cutoff 0.9 -log utax.log
#optional, make alignment file to show taxonomy assignment
usearch -usearch_global otus.fa -db ~/db/utax/16s.udb -strand both -id 0.97 -alnout otus_ref.aln -userout otus_ref.user -userfields query+target+id
#map trimmed reads to otus to generate abundance profile
usearch -usearch_global trimmed.fa -db otus_tax.fa -strand both -id 0.97 -log make_otutab.log -otutabout otutab.txt -biomout otutab.biom
#optional, depreplicate taxa in the otu table, taking sums when taxonomy is the same, n.b. at higher taxonomy this can mean very different otus are agglomerated.
python ~/s/sum-tax-reps.py otutab.txt otu-noreps.txt

#QIIME analysis of uparse results
#QIIME version 1.8
#get otu sequences for tree and calculate unifrac distance matrix
cut -f1 otu-noreps.txt | sed -n '1!p' > otu.list
biom convert -i otu-noreps.txt -o otu.biom --to-json --table-type "OTU table" -m mapping.txt
#n.b.
beta_diversity.py -i otu.biom -m weighted_unifrac -o ./ -t reps.tree
#calculate principal coordiantes on unifrac
principal_coordinates.py -i weighted_unifrac_otu.txt -o pco.out
#normalise the otu table with css for pearson and anova:
normalize_table.py -i otutab.biom -o otu-css.biom
biom convert -i otu-css.biom -o otu-css.txt --to-tsv --header-key taxonomy
group_significance.py -i otu-css.biom -m qiime_mapping.txt -c Description -o otu-css.anova -s ANOVA

#PICRUSt - can ony use QIIME output with greengenes taxon codes. Therefore otu table also made using QIIME:
#convert reads to qual - put all .sff files in a /sff folder
python convert_folder.py reads/sff sff qual
python convert_folder.py reads/sff sff fasta

cat reads/*.fasta > reads/reads.fasta
cat reads/*.qual > reads/reads.qual
#convert sequence to uppercase
tr '[:lower:]' '[:upper:]' < reads.fasta >reads2.fasta
#qiime label reads
split_libraries.py -b 16 -m qiime_mapping.txt -f reads/reads2.fasta -q reads/reads.qual -o ./
#clip reads at reverse primer
python clip_at_primer.py reads/reads2.fasta reads/seqs_clip.fna TGACTGAGCGGGC
#pick closed reference otus
pick_closed_reference_otus.py -i seqs_clip.fna -o output -r ~/db/gg_13_8_otus/rep_set/97_otus.fasta -t ~/db/gg_13_8_otus/taxonomy/97_otu_taxonomy.txt -f  -p params.txt
biom convert -i output/otu_table.biom -o output/otu_table.txt --to-tsv

#PICRUSt v1.0.0
normalize_by_copy_number.py -i output/otu_table.biom -o output/otu-norm.biom
predict_metagenomes.py -i output/otu-norm.biom -o output/picrust.biom
categorize_by_function.py -i output/picrust.biom -c KEGG_Pathways -l 3 -o l3-kegg.biom
biom convert -i l3-kegg.biom -o l3-kegg.txt --to-tsv

#0.5% abundance filter (QIIME):
filter_otus_from_otu_table.py -i l3-kegg.biom -o kegg.filt.biom --min_count_fraction 0.005
biom convert -i l3-kegg.biom -o kegg.txt --to-tsv
#normalise:
normalize_table.py -i l3-kegg.biom -o kegg-norm.biom
#ANOVA
group_significance.py -i kegg-norm.biom -m qiime_mapping.txt -c Description -o kegg-norm.anova.txt -s ANOVA


