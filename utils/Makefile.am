bin_SCRIPTS= genodb_bam_analysis genodb_proc_dataset genodb_filter_vcf	\
genodb_gen_bed_for_genome genodb_check_ref_and_bam_contigs		\
genodb_get_bam_contigs genodb_convert_snppos_to_snpgcc			\
genodb_create_snv_pos_ascat genodb_query_ega_metadata			\
genodb_query_icgc_metadata genodb_download_gdc_metadata			\
genodb_query_gdc_metadata genodb_create_genref_for_bam			\
genodb_map_vcf_contnames genodb_get_two_elem_combs			\
genodb_filter_nondiscordant_pheno_entries

libexec_SCRIPTS= genodb_genref_processes genodb_bam_common_lib	\
genodb_bam_download_processes genodb_bam_manip_processes	\
genodb_cleaning_processes genodb_bam_analysis_processes		\
genodb_bam_summariz_processes genodb_bam_filter_processes	\
genodb_bam_facets_processes genodb_bam_ascat_processes		\
genodb_run_facets genodb_run_sequenza genodb_run_ascat		\
genodb_file_shuffle genodb_get_entrez_fasta			\
genodb_get_filtered_sam_align genodb_get_filtered_sam_header	\
genodb_reorder_fa_seqs genodb_reorder_fa_seqs_lowmem		\
genodb_convert_allele_counts genodb_map_contnames		\
genodb_filter_contig_from_genref

CLEANFILES = $(bin_SCRIPTS) $(libexec_SCRIPTS)

EXTRA_DIST = genodb_bam_common_lib.sh genodb_genref_processes.sh	\
genodb_bam_download_processes.sh genodb_bam_manip_processes.sh		\
genodb_cleaning_processes.sh genodb_bam_analysis_processes.sh		\
genodb_bam_summariz_processes.sh genodb_bam_analysis.sh			\
genodb_bam_filter_processes.sh genodb_bam_facets_processes.sh		\
genodb_bam_ascat_processes.sh genodb_proc_dataset.sh			\
genodb_file_shuffle.sh genodb_filter_vcf.py				\
genodb_gen_bed_for_genome.sh genodb_check_ref_and_bam_contigs.sh	\
genodb_get_bam_contigs.sh genodb_get_entrez_fasta.py			\
genodb_convert_snppos_to_snpgcc.sh genodb_create_snv_pos_ascat.py	\
genodb_query_ega_metadata.py genodb_download_gdc_metadata.sh		\
genodb_query_gdc_metadata.py genodb_query_icgc_metadata.py		\
genodb_create_genref_for_bam.sh genodb_filter_contig_from_genref.py	\
genodb_get_filtered_sam_header.py genodb_get_filtered_sam_align.py	\
genodb_map_contnames.py genodb_map_vcf_contnames.py genodb_run_facets.R	\
genodb_run_sequenza.R genodb_convert_allele_counts.R genodb_run_ascat.R	\
genodb_get_two_elem_combs.py						\
genodb_filter_nondiscordant_pheno_entries.py genodb_reorder_fa_seqs.py	\
genodb_reorder_fa_seqs_lowmem.sh

# Build rules for Debasher modules
bam_analysis: bam_analysis.sh
	rm -f $@
	echo "genodebasher_bindir=\"$(bindir)\"" >> $@
	echo "genodebasher_libexecdir=\"$(libexecdir)\"" >> $@
	echo "ECHO=\"$(ECHO)\"" >> $@
	echo "HEAD=\"$(HEAD)\"" >> $@
	echo "TAIL=\"$(TAIL)\"" >> $@
	echo "PASTE=\"$(PASTE)\"" >> $@
	echo "CUT=\"$(CUT)\"" >> $@
	echo "AWK=\"$(AWK)\"" >> $@
	echo "GREP=\"$(GREP)\"" >> $@
	echo "SED=\"$(SED)\"" >> $@
	echo "FIND=\"$(FIND)\"" >> $@
	echo "CAT=\"$(CAT)\"" >> $@
	echo "RM=\"$(RM_CMD)\"" >> $@
	echo "MV=\"$(MV)\"" >> $@
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
	echo "genodebasher_bindir=\"$(bindir)\"" >> $@
	echo "genodebasher_libexecdir=\"$(libexecdir)\"" >> $@
	echo "ECHO=\"$(ECHO)\"" >> $@
	echo "HEAD=\"$(HEAD)\"" >> $@
	echo "TAIL=\"$(TAIL)\"" >> $@
	echo "PASTE=\"$(PASTE)\"" >> $@
	echo "CUT=\"$(CUT)\"" >> $@
	echo "AWK=\"$(AWK)\"" >> $@
	echo "GREP=\"$(GREP)\"" >> $@
	echo "SED=\"$(SED)\"" >> $@
	echo "CAT=\"$(CAT)\"" >> $@
	echo "RM=\"$(RM_CMD)\"" >> $@
	echo "MV=\"$(MV)\"" >> $@
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
	echo "genodebasher_bindir=\"$(bindir)\"" >> $@
	echo "genodebasher_libexecdir=\"$(libexecdir)\"" >> $@
	echo "ECHO=\"$(ECHO)\"" >> $@
	echo "HEAD=\"$(HEAD)\"" >> $@
	echo "TAIL=\"$(TAIL)\"" >> $@
	echo "PASTE=\"$(PASTE)\"" >> $@
	echo "CUT=\"$(CUT)\"" >> $@
	echo "AWK=\"$(AWK)\"" >> $@
	echo "GREP=\"$(GREP)\"" >> $@
	echo "SED=\"$(SED)\"" >> $@
	echo "CAT=\"$(CAT)\"" >> $@
	echo "RM=\"$(RM_CMD)\"" >> $@
	echo "MV=\"$(MV)\"" >> $@
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
	echo "genodebasher_pkgname=\"$(PACKAGE_NAME)\"" >> $@
	echo "genodebasher_version=\"$(PACKAGE_VERSION)\"" >> $@
	echo "genodebasher_bindir=\"$(bindir)\"" >> $@
	echo "genodebasher_libexecdir=\"$(libexecdir)\"" >> $@
	echo "genodebasher_datadir=\"$(datadir)/$(PACKAGE_NAME)\"" >> $@
	echo "BASH=\"$(BASH)\"" >> $@
	echo "ECHO=\"$(ECHO)\"" >> $@
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
	echo "CAT=\"$(CAT)\"" >> $@
	echo "RM=\"$(RM_CMD)\"" >> $@
	echo "MV=\"$(MV)\"" >> $@
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
	echo "DEBASHER_HOME_DIR=\"$(DEBASHER_HOME_DIR)\"" >> $@
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
