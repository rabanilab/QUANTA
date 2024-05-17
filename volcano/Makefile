# ---------------------------------------------------------------------------------
# 3'UTR annotations
#
# Required files (Downloads):
#  1. Downloads/biomart_ENS_to_Refseq.txt
#  2. Downloads/ucsc.refseq_curated.fa
#  3. Downloads/biomart_GRCz11_3UTR_ensembl103.fa
# ---------------------------------------------------------------------------------
MIN_LEN = 10

# select longest 3'UTR annotation per gene
filter_length:
	cat $(f).txt \
	| perl -e 'while(<STDIN>){chomp $$_;@l=split("\t",$$_);if(length($$h{$$l[0]})<length($$l[3])){$$h{$$l[0]}=$$l[3];$$u{$$l[0]}=$$l[1]."\t".$$l[2]}}foreach $$k(keys %h){print "$$k\t$$u{$$k}\t$$h{$$k}\n"}' \
	> $(f).longest.txt;


# DB specific annotations
gencode_3utr:
	cat Downloads/ucsc.gencode.fa \
	| fasta2stab.pl \
	| sed 's/ .*\t/\t/g' \
	| sed 's/.*_knownGene_//g' \
	| sort \
	| uniq \
	| sed 's/\.[0-9]*\t/\t/g' \
	| join.pl -o NA - Downloads/biomart_ENS_to_geneid.txt -2 2 \
	| cut.pl -f 3,1,5,2 \
	> 3utr_seqs.gencode.txt; \
	make filter_length f=3utr_seqs.gencode;

refseq_3utr:
	cat Downloads/biomart_ENS_to_Refseq.txt \
	| body.pl 2 -1 \
	| cut.pl -f 3,1,2 \
	| grep -E -v "^\t" \
	| sort \
	| uniq \
	> c1; \
	cat Downloads/biomart_ENS_to_Refseq.txt \
	| body.pl 2 -1 \
	| cut.pl -f 4,1,2 \
	| grep -E -v "^\t" \
	| sort \
	| uniq \
	> c2; \
	cat c1 > conv.txt; \
	join.pl -neg c2 c1 -1 2,3 -2 2,3 | cut.pl -f 3,1,2 >> conv.txt; \
	rm -rf c1 c2; \
	\
	cat Downloads/ucsc.refseq_curated.fa \
	| fasta2stab.pl \
	| sed 's/ .*\t/\t/g' \
	| sed 's/.*_ncbiRefSeqCurated_//g' \
	| sort \
	| uniq \
	| sed 's/\.[0-9]*\t/\t/g' \
	| join.pl conv.txt - \
	| cut.pl -f 2,3,1,4 \
	> 3utr_seqs.refseq.txt; \
	make filter_length f=3utr_seqs.refseq;

ensembl_3utr:
	cat Downloads/biomart_*.fa \
	| fasta2stab.pl \
	| grep -v SEQUENCEUNAVAILABLE \
	| tr "|" "\t" \
	| cut -f 1,3,5,6 \
	| sort \
	| uniq \
	> 3utr_seqs.ensembl.txt; \
	cat Downloads/ucsc.ensembl.fa \
	| fasta2stab.pl \
	| sed 's/ .*\t/\t/g' \
	| sed 's/.*_ensGene_//g' \
	| sort \
	| uniq \
	| sed 's/\.[0-9]*\t/\t/g' \
	| join.pl -o NA - Downloads/biomart_ENS_to_geneid.txt -2 2 \
	| cut.pl -f 3,1,5,2 \
	>> 3utr_seqs.ensembl.txt; \
	cat 3utr_seqs.ensembl.txt \
	| join.pl - Downloads/biomart_ENS_to_geneid.txt -1 1,2 -2 1,2 -o NA \
	| cut.pl -f 1,2,5,3 \
	> w; mv w 3utr_seqs.ensembl.txt; \
	make filter_length f=3utr_seqs.ensembl;

# merge annotations
3utr_merge_refseq_ensembl:
	cat 3utr_seqs.ensembl.longest.txt \
	| cut -f 1,3 \
	| join.pl 3utr_seqs.refseq.longest.txt - -o NA \
	| cut.pl -f 1,3,5,4 \
	> 3utr_seqs.all.txt; \
	join.pl -neg 3utr_seqs.ensembl.longest.txt 3utr_seqs.refseq.longest.txt \
	| cut -f 1,3,4 \
	| add_column.pl -c 1 -s NA \
	>> 3utr_seqs.all.txt; \
	\
	cat 3utr_seqs.all.txt \
	| cut -f 4 \
	| tr -d "N" \
	| paste 3utr_seqs.all.txt - \
	| cut -f 1,5 \
	| stab2length.pl \
	| filter.pl -c 1 -min $(MIN_LEN) \
	| cut -f 1 \
	| join.pl - 3utr_seqs.all.txt \
	| cut -f 1,4 \
	> 3utr_seqs.genes.txt;

3utr_merge_refseq_gencode_ensembl:
	cat 3utr_seqs.refseq.longest.txt \
	| cut -f 1,4 \
	> tmp.all; \
	cat 3utr_seqs.gencode.longest.txt \
	| cut -f 1,4 \
	| join.pl -neg - 3utr_seqs.all.txt \
	>> tmp.all; \
	cat 3utr_seqs.ensembl.longest.txt \
	| cut -f 1,4 \
	| join.pl -neg - 3utr_seqs.all.txt \
	>> tmp.all; \
	\
	cat tmp.all \
	| cut -f 1 \
	| join.pl - Downloads/biomart_ENS_to_Refseq.txt -o NA \
	| cut -f 1,3 \
	| add_column.pl -s NA \
	| sed 's/\t\t/\t/g' \
	| cut -f 1,2 \
	| join.pl - Downloads/biomart_ENS_to_geneid.txt -o NA \
	| cut -f 1,2,5 \
	| add_column.pl -s NA \
	| sed 's/\t\t/\t/g' \
	| cut -f 1,2,3 \
	> tmp.ids; \
	join.pl tmp.ids tmp.all \
	> 3utr_seqs.all.txt; \
	\
	cat 3utr_seqs.all.txt \
	| cut -f 4 \
	| tr -d "N" \
	| paste 3utr_seqs.all.txt - \
	| cut -f 1,5 \
	| stab2length.pl \
	| filter.pl -c 1 -min $(MIN_LEN) \
	| cut -f 1 \
	| join.pl - 3utr_seqs.all.txt \
	| cut -f 1,4 \
	> 3utr_seqs.genes.txt; \
	rm -rf tmp.all tmp.genes tmp.ids;

3utr_stat:
	wc -l 3utr_seqs.* \
	| grep -v total \
	| sed 's/^  *//g;s/  */\t/g;s/.txt//g' \
	> stat.txt;
