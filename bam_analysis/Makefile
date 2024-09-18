# ----------------------------------------------------------
# RNA-seq analysis
# Requirements:
# 1. SAMTOOLS
# 2. CUFFLINKS
# ----------------------------------------------------------
work_dir = cufflinks
bam_dir = bam_files
cufflinks_type = fr-firststrand
mem = 64
time = 24:0:0

# zebrafish
ZTX_REF = gtf_files/zebrafish_GRCz11.gtf
zpre_flag = -t -o DAR

# frog
FTX_REF = gtf_files/frog_XENTR_10.0_xenbase.new.gtf
fpre_flag = -x

# mouse
MTX_REF = gtf_files/Mus_musculus.GRCm39.108.pre.gtf
mpre_flag = -t -o MUS

# human
HTX_REF = gtf_files/Homo_sapiens.GRCh38.110.pre.gtf
hpre_flag = -t -o ""



# run cufflinks on BAM files
process_cufflinks:
	mkdir $(work_dir)/$(id).cuff; \
	cufflinks -v -q --no-update-check -o $(work_dir)/$(id).cuff -G $(TX_REF) --library-type $(cufflinks_type) $(bam_dir)/$(id).bam; \
	cat $(work_dir)/$(id).cuff/genes.fpkm_tracking \
	| cut -f 1,4,5,10 \
	> $(work_dir)/$(id).g.fpkm; \
	cat $(work_dir)/$(id).cuff/isoforms.fpkm_tracking \
	| cut -f 1,4,5,10 \
	> $(work_dir)/$(id).i.fpkm; \
	echo "gene_id#gene_short_name#Mature_FPKM#Pre_FPKM" \
	| tr "#" "\t" \
	> $(work_dir)/$(id).p.fpkm; \
	cat $(work_dir)/$(id).i.fpkm \
	| ./body.pl 2 -1 \
	| ./cufflinks_fpkm_collect.pl $(pre_flag) \
	>> $(work_dir)/$(id).p.fpkm; \
	rm -rf $(work_dir)/$(id).cuff;

rnaseq_cuff:
	mkdir $(work_dir); \
	rm -rf cuff_array.txt; \
	$(foreach id, $(shell /bin/ls -1 $(bam_dir)/*$(o).bam | sed 's/.*\///g;s/.bam.*//g'), \
		samtools view $(bam_dir)/$(id).bam | head; \
		echo "make process_cufflinks id=$(id) TX_REF=$(TX_REF) pre_flag=$(pre_flag)" >> cuff_array.txt; \
	) \
	./sbatch.pl -t $(time) -R $(mem) -a cuff_array.txt;


# specific command for each organism
rnaseq_cuff_zebrafish:
	make rnaseq_cuff TX_REF=$(ZTX_REF) pre_flag=$(zpre_flag);

rnaseq_cuff_frog:
	make rnaseq_cuff TX_REF=$(FTX_REF) pre_flag=$(fpre_flag);

rnaseq_cuff_human:
	make make rnaseq_cuff TX_REF=$(HTX_REF) pre_flag=$(hpre_flag);

rnaseq_cuff_mouse:
	make rnaseq_cuff TX_REF=$(MTX_REF) pre_flag=$(mpre_flag);


# test running the analysis
test_run:
	make rnaseq_cuff TX_REF=$(ZTX_REF) pre_flag=$(zpre_flag) o=zebrafish;
