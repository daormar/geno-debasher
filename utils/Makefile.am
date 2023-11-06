bin_SCRIPTS= genop_bam_analysis genop_proc_dataset filter_vcf		\
genop_gen_bed_for_genome genop_check_ref_and_bam_contigs		\
genop_get_bam_contigs genop_convert_snppos_to_snpgcc			\
genop_create_snv_pos_ascat genop_query_ega_metadata			\
genop_query_icgc_metadata genop_download_gdc_metadata			\
genop_query_gdc_metadata genop_create_genref_for_bam			\
filter_contig_from_genref genop_map_vcf_contnames			\
genop_get_two_elem_combs genop_filter_nondiscordant_pheno_entries

libexec_SCRIPTS= genop_genref_processes genop_bam_common_lib		\
genop_bam_download_processes genop_bam_manip_processes			\
genop_cleaning_processes genop_bam_analysis_processes			\
genop_bam_summariz_processes genop_bam_filter_processes			\
genop_bam_facets_processes genop_bam_ascat_processes genop_run_facets	\
genop_run_sequenza genop_run_ascat genop_file_shuffle			\
genop_get_entrez_fasta genop_get_filtered_sam_align			\
genop_get_filtered_sam_header genop_reorder_fa_seqs			\
genop_reorder_fa_seqs_lowmem genop_convert_allele_counts		\
genop_map_contnames

CLEANFILES = $(bin_SCRIPTS) $(libexec_SCRIPTS)

EXTRA_DIST = genop_bam_common_lib.sh genop_genref_processes.sh		\
genop_bam_download_processes.sh genop_bam_manip_processes.sh		\
genop_cleaning_processes.sh genop_bam_analysis_processes.sh		\
genop_bam_summariz_processes.sh genop_bam_analysis.sh			\
genop_bam_filter_processes.sh genop_bam_facets_processes.sh		\
genop_bam_ascat_processes.sh genop_proc_dataset.sh			\
genop_file_shuffle.sh filter_vcf.py genop_gen_bed_for_genome.sh		\
genop_check_ref_and_bam_contigs.sh genop_get_bam_contigs.sh		\
genop_get_entrez_fasta.py genop_convert_snppos_to_snpgcc.sh		\
genop_create_snv_pos_ascat.py genop_query_ega_metadata.py		\
genop_download_gdc_metadata.sh genop_query_gdc_metadata.py		\
genop_query_icgc_metadata.py genop_create_genref_for_bam.sh		\
filter_contig_from_genref.py genop_get_filtered_sam_header.py		\
genop_get_filtered_sam_align.py genop_map_contnames.py			\
genop_map_vcf_contnames.py genop_run_facets.R genop_run_sequenza.R	\
genop_convert_allele_counts.R genop_run_ascat.R				\
genop_get_two_elem_combs.py genop_filter_nondiscordant_pheno_entries.py	\
genop_reorder_fa_seqs.py genop_reorder_fa_seqs_lowmem.sh

# Build rules for PanPipe modules
bam_analysis: bam_analysis.sh
	rm -f $@
	echo "genopanpipe_bindir=\"$(bindir)\"" >> $@
	echo "genopanpipe_libexecdir=\"$(libexecdir)\"" >> $@
	echo "HEAD=\"$(HEAD)\"" >> $@
	echo "TAIL=\"$(TAIL)\"" >> $@
	echo "PASTE=\"$(PASTE)\"" >> $@
	echo "CUT=\"$(CUT)\"" >> $@
	echo "AWK=\"$(AWK)\"" >> $@
	echo "GREP=\"$(GREP)\"" >> $@
	echo "SED=\"$(SED)\"" >> $@
	echo "FIND=\"$(FIND)\"" >> $@
	echo "CAT=\"$(CAT)\"" >> $@
	echo "ZCAT=\"$(ZCAT)\"" >> $@
	echo "GZIP=\"$(GZIP)\"" >> $@
	echo "BASENAME=\"$(BASENAME)\"" >> $@
	echo "DIRNAME=\"$(DIRNAME)\"" >> $@
	echo "PLATYPUS_HOME_DIR=\"$(PLATYPUS_HOME_DIR)\"" >> $@
	echo "LUMPY_HOME_DIR=\"$(LUMPY_HOME_DIR)\"" >> $@
	echo "ICGCSCOR_HOME_DIR=\"$(ICGCSCOR_HOME_DIR)\"" >> $@
	echo "ASPERA_HOME_DIR=\"$(ASPERA_HOME_DIR)\"" >> $@
	echo "EGADECRYPT_HOME_DIR=\"$(EGADECRYPT_HOME_DIR)\"" >> $@
	echo "DECSINGLE_HOME_DIR=\"$(DECSINGLE_HOME_DIR)\"" >> $@
	echo "DECSINGLE_JAVA_HOME_DIR=\"$(DECSINGLE_JAVA_HOME_DIR)\"" >> $@
	cat bam_analysis.sh >> $@

bam_facets: bam_facets.sh
	rm -f $@
	echo "genopanpipe_bindir=\"$(bindir)\"" >> $@
	echo "genopanpipe_libexecdir=\"$(libexecdir)\"" >> $@
	echo "HEAD=\"$(HEAD)\"" >> $@
	echo "TAIL=\"$(TAIL)\"" >> $@
	echo "PASTE=\"$(PASTE)\"" >> $@
	echo "CUT=\"$(CUT)\"" >> $@
	echo "AWK=\"$(AWK)\"" >> $@
	echo "GREP=\"$(GREP)\"" >> $@
	echo "SED=\"$(SED)\"" >> $@
	echo "ZCAT=\"$(ZCAT)\"" >> $@
	echo "GZIP=\"$(GZIP)\"" >> $@
	echo "BASENAME=\"$(BASENAME)\"" >> $@
	echo "DIRNAME=\"$(DIRNAME)\"" >> $@
	echo "PLATYPUS_HOME_DIR=\"$(PLATYPUS_HOME_DIR)\"" >> $@
	echo "LUMPY_HOME_DIR=\"$(LUMPY_HOME_DIR)\"" >> $@
	echo "ICGCSCOR_HOME_DIR=\"$(ICGCSCOR_HOME_DIR)\"" >> $@
	echo "ASPERA_HOME_DIR=\"$(ASPERA_HOME_DIR)\"" >> $@
	echo "EGADECRYPT_HOME_DIR=\"$(EGADECRYPT_HOME_DIR)\"" >> $@
	echo "DECSINGLE_HOME_DIR=\"$(DECSINGLE_HOME_DIR)\"" >> $@
	echo "DECSINGLE_JAVA_HOME_DIR=\"$(DECSINGLE_JAVA_HOME_DIR)\"" >> $@
	cat bam_facets.sh >> $@

bam_ascat: bam_ascat.sh
	rm -f $@
	echo "genopanpipe_bindir=\"$(bindir)\"" >> $@
	echo "genopanpipe_libexecdir=\"$(libexecdir)\"" >> $@
	echo "HEAD=\"$(HEAD)\"" >> $@
	echo "TAIL=\"$(TAIL)\"" >> $@
	echo "PASTE=\"$(PASTE)\"" >> $@
	echo "CUT=\"$(CUT)\"" >> $@
	echo "AWK=\"$(AWK)\"" >> $@
	echo "GREP=\"$(GREP)\"" >> $@
	echo "SED=\"$(SED)\"" >> $@
	echo "CAT=\"$(CAT)\"" >> $@
	echo "ZCAT=\"$(ZCAT)\"" >> $@
	echo "GZIP=\"$(GZIP)\"" >> $@
	echo "BASENAME=\"$(BASENAME)\"" >> $@
	echo "DIRNAME=\"$(DIRNAME)\"" >> $@
	echo "PLATYPUS_HOME_DIR=\"$(PLATYPUS_HOME_DIR)\"" >> $@
	echo "LUMPY_HOME_DIR=\"$(LUMPY_HOME_DIR)\"" >> $@
	echo "ICGCSCOR_HOME_DIR=\"$(ICGCSCOR_HOME_DIR)\"" >> $@
	echo "ASPERA_HOME_DIR=\"$(ASPERA_HOME_DIR)\"" >> $@
	echo "EGADECRYPT_HOME_DIR=\"$(EGADECRYPT_HOME_DIR)\"" >> $@
	echo "DECSINGLE_HOME_DIR=\"$(DECSINGLE_HOME_DIR)\"" >> $@
	echo "DECSINGLE_JAVA_HOME_DIR=\"$(DECSINGLE_JAVA_HOME_DIR)\"" >> $@
	cat bam_ascat.sh >> $@

# Build using suffix rules (portable)
SUFFIXES= .sh .py .r .R

.sh :
	rm -f $@
	echo "#! " $(BASH) > $@
	echo "genopanpipe_pkgname=\"$(PACKAGE_NAME)\"" >> $@
	echo "genopanpipe_version=\"$(PACKAGE_VERSION)\"" >> $@
	echo "genopanpipe_bindir=\"$(bindir)\"" >> $@
	echo "genopanpipe_libexecdir=\"$(libexecdir)\"" >> $@
	echo "genopanpipe_datadir=\"$(datadir)/$(PACKAGE_NAME)\"" >> $@
	echo "BASH=\"$(BASH)\"" >> $@
	echo "AWK=\"$(AWK)\"" >> $@
	echo "SORT=\"$(SORT)\"" >> $@
	echo "sortpars=\"$(sortpars)\"" >> $@
	echo "sortT=\"$(sortT)\"" >> $@
	echo "UNIQ=\"$(UNIQ)\"" >> $@
	echo "GREP=\"$(GREP)\"" >> $@
	echo "SED=\"$(SED)\"" >> $@
	echo "BASENAME=\"$(BASENAME)\"" >> $@
	echo "DIRNAME=\"$(DIRNAME)\"" >> $@
	echo "SPLIT=\"$(SPLIT)\"" >> $@
	echo "CUT=\"$(CUT)\"" >> $@
	echo "HEAD=\"$(HEAD)\"" >> $@
	echo "TAIL=\"$(TAIL)\"" >> $@
	echo "PASTE=\"$(PASTE)\"" >> $@
	echo "DIFF=\"$(DIFF)\"" >> $@
	echo "FIND=\"$(FIND)\"" >> $@
	echo "WC=\"$(WC)\"" >> $@
	echo "MKTEMP=\"$(MKTEMP)\"" >> $@
	echo "FLOCK=\"$(FLOCK)\"" >> $@
	echo "CAT=\"$(CAT)\"" >> $@
	echo "ZCAT=\"$(ZCAT)\"" >> $@
	echo "GZIP=\"$(GZIP)\"" >> $@
	echo "WGET=\"$(WGET)\"" >> $@
	echo "PANDOC=\"$(PANDOC)\"" >> $@
	echo "PYTHON=\"$(PYTHON)\"" >> $@
	echo "RSCRIPT=\"$(RSCRIPT)\"" >> $@
	echo "JAVA=\"$(JAVA)\"" >> $@
	echo "PANPIPE_HOME_DIR=\"$(PANPIPE_HOME_DIR)\"" >> $@
	echo "PLATYPUS_HOME_DIR=\"$(PLATYPUS_HOME_DIR)\"" >> $@
	echo "LUMPY_HOME_DIR=\"$(LUMPY_HOME_DIR)\"" >> $@
	echo "ICGCSCOR_HOME_DIR=\"$(ICGCSCOR_HOME_DIR)\"" >> $@
	echo "ASPERA_HOME_DIR=\"$(ASPERA_HOME_DIR)\"" >> $@
	echo "EGADECRYPT_HOME_DIR=\"$(EGADECRYPT_HOME_DIR)\"" >> $@
	echo "DECSINGLE_HOME_DIR=\"$(DECSINGLE_HOME_DIR)\"" >> $@
	echo "DECSINGLE_JAVA_HOME_DIR=\"$(DECSINGLE_JAVA_HOME_DIR)\"" >> $@
	cat $< >> $@
	chmod ugo+x $@

.py:
	rm -f $@
	echo "#! " $(PYTHON) > $@
	cat $< >> $@
	chmod ugo+x $@

.r:
	rm -f $@
	echo "#! " $(RSCRIPT) > $@
	cat $< >> $@
	chmod ugo+x $@

.R:
	rm -f $@
	echo "#! " $(RSCRIPT) > $@
	cat $< >> $@
	chmod ugo+x $@
