SCRIPTS_TO_BE_BUILT= bam_common_lib genref_steps bam_download_steps	\
bam_manip_steps cleaning_steps bam_analysis_steps bam_analysis		\
bam_filter bam_facets bam_ascat bam_ascatngs analyze_dataset		\
file_shuffle filter_vcf gen_bed_for_genome check_ref_and_bam_contigs	\
get_bam_contigs get_entrez_fasta gen_wisecondorx_ref			\
convert_snppos_to_snpgcc create_snv_pos_ascat query_ega_metadata	\
query_icgc_metadata create_genref_for_bam filter_contig_from_genref	\
get_filtered_sam_header get_filtered_sam_align map_contnames		\
map_vcf_contnames run_facets run_sequenza convert_allele_counts		\
run_ascat get_two_elem_combs

bin_SCRIPTS = $(SCRIPTS_TO_BE_BUILT)

CLEANFILES = $(SCRIPTS_TO_BE_BUILT)

EXTRA_DIST = bam_common_lib.sh genref_steps.sh bam_download_steps.sh	\
bam_manip_steps.sh cleaning_steps.sh bam_analysis_steps.sh		\
bam_analysis.sh bam_filter.sh bam_facets.sh bam_ascat.sh		\
bam_ascatngs.sh analyze_dataset.sh file_shuffle.sh filter_vcf.py	\
gen_bed_for_genome.sh check_ref_and_bam_contigs.sh get_bam_contigs.sh	\
get_entrez_fasta.py gen_wisecondorx_ref.sh convert_snppos_to_snpgcc.sh	\
create_snv_pos_ascat.py query_ega_metadata.py query_icgc_metadata.py	\
create_genref_for_bam.sh filter_contig_from_genref.py			\
get_filtered_sam_header.py get_filtered_sam_align.py map_contnames.py	\
map_vcf_contnames.py run_facets.R run_sequenza.R			\
convert_allele_counts.R run_ascat.R get_two_elem_combs.py

# Build rules for libraries and bash files implementing steps
bam_common_lib: bam_common_lib.sh
	rm -f $@
	cat bam_common_lib.sh >> $@

genref_steps: genref_steps.sh
	rm -f $@
	cat genref_steps.sh >> $@

bam_download_steps: bam_download_steps.sh
	rm -f $@
	cat bam_download_steps.sh >> $@

bam_manip_steps: bam_manip_steps.sh
	rm -f $@
	cat bam_manip_steps.sh >> $@

bam_analysis_steps: bam_analysis_steps.sh
	rm -f $@
	cat bam_analysis_steps.sh >> $@

cleaning_steps: cleaning_steps.sh
	rm -f $@
	cat cleaning_steps.sh >> $@

# Build rules for PanPipe modules
bam_analysis: bam_analysis.sh
	rm -f $@
	echo "biopanpipe_bindir="$(prefix)/bin/ >> $@
	echo "HEAD="$(HEAD) >> $@
	echo "TAIL="$(TAIL) >> $@
	echo "PASTE="$(PASTE) >> $@
	echo "CUT="$(CUT) >> $@
	echo "AWK="$(AWK) >> $@
	echo "GREP="$(GREP) >> $@
	echo "SED="$(SED) >> $@
	echo "FIND="$(FIND) >> $@
	echo "ZCAT="$(ZCAT) >> $@
	echo "GZIP="$(GZIP) >> $@
	echo "BASENAME="$(BASENAME) >> $@
	echo "DIRNAME="$(DIRNAME) >> $@
	echo "PLATYPUS_HOME_DIR="$(PLATYPUS_HOME_DIR) >> $@
	echo "LUMPY_HOME_DIR="$(LUMPY_HOME_DIR) >> $@
	echo "ICGCSCOR_HOME_DIR="$(ICGCSCOR_HOME_DIR) >> $@
	echo "ASPERA_HOME_DIR="$(ASPERA_HOME_DIR) >> $@
	echo "EGADECRYPT_HOME_DIR="$(EGADECRYPT_HOME_DIR) >> $@
	echo "DECSINGLE_HOME_DIR="$(DECSINGLE_HOME_DIR) >> $@
	echo "DECSINGLE_JAVA_HOME_DIR="$(DECSINGLE_JAVA_HOME_DIR) >> $@
	cat bam_analysis.sh >> $@

bam_facets: bam_facets.sh
	rm -f $@
	echo "biopanpipe_bindir="$(prefix)/bin/ >> $@
	echo "HEAD="$(HEAD) >> $@
	echo "TAIL="$(TAIL) >> $@
	echo "PASTE="$(PASTE) >> $@
	echo "CUT="$(CUT) >> $@
	echo "AWK="$(AWK) >> $@
	echo "GREP="$(GREP) >> $@
	echo "SED="$(SED) >> $@
	echo "ZCAT="$(ZCAT) >> $@
	echo "GZIP="$(GZIP) >> $@
	echo "BASENAME="$(BASENAME) >> $@
	echo "DIRNAME="$(DIRNAME) >> $@
	echo "PLATYPUS_HOME_DIR="$(PLATYPUS_HOME_DIR) >> $@
	echo "LUMPY_HOME_DIR="$(LUMPY_HOME_DIR) >> $@
	echo "ICGCSCOR_HOME_DIR="$(ICGCSCOR_HOME_DIR) >> $@
	echo "ASPERA_HOME_DIR="$(ASPERA_HOME_DIR) >> $@
	echo "EGADECRYPT_HOME_DIR="$(EGADECRYPT_HOME_DIR) >> $@
	echo "DECSINGLE_HOME_DIR="$(DECSINGLE_HOME_DIR) >> $@
	echo "DECSINGLE_JAVA_HOME_DIR="$(DECSINGLE_JAVA_HOME_DIR) >> $@
	cat bam_facets.sh >> $@

bam_ascat: bam_ascat.sh
	rm -f $@
	echo "biopanpipe_bindir="$(prefix)/bin/ >> $@
	echo "HEAD="$(HEAD) >> $@
	echo "TAIL="$(TAIL) >> $@
	echo "PASTE="$(PASTE) >> $@
	echo "CUT="$(CUT) >> $@
	echo "AWK="$(AWK) >> $@
	echo "GREP="$(GREP) >> $@
	echo "SED="$(SED) >> $@
	echo "ZCAT="$(ZCAT) >> $@
	echo "GZIP="$(GZIP) >> $@
	echo "BASENAME="$(BASENAME) >> $@
	echo "DIRNAME="$(DIRNAME) >> $@
	echo "PLATYPUS_HOME_DIR="$(PLATYPUS_HOME_DIR) >> $@
	echo "LUMPY_HOME_DIR="$(LUMPY_HOME_DIR) >> $@
	echo "ICGCSCOR_HOME_DIR="$(ICGCSCOR_HOME_DIR) >> $@
	echo "ASPERA_HOME_DIR="$(ASPERA_HOME_DIR) >> $@
	echo "EGADECRYPT_HOME_DIR="$(EGADECRYPT_HOME_DIR) >> $@
	echo "DECSINGLE_HOME_DIR="$(DECSINGLE_HOME_DIR) >> $@
	echo "DECSINGLE_JAVA_HOME_DIR="$(DECSINGLE_JAVA_HOME_DIR) >> $@
	cat bam_ascat.sh >> $@

bam_ascatngs: bam_ascatngs.sh
	rm -f $@
	echo "biopanpipe_bindir="$(prefix)/bin/ >> $@
	echo "CUT="$(CUT) >> $@
	echo "AWK="$(AWK) >> $@
	echo "GREP="$(GREP) >> $@
	echo "SED="$(SED) >> $@
	echo "ZCAT="$(ZCAT) >> $@
	echo "GZIP="$(GZIP) >> $@
	echo "BASENAME="$(BASENAME) >> $@
	echo "DIRNAME="$(DIRNAME) >> $@
	echo "PLATYPUS_HOME_DIR="$(PLATYPUS_HOME_DIR) >> $@
	echo "LUMPY_HOME_DIR="$(LUMPY_HOME_DIR) >> $@
	echo "ICGCSCOR_HOME_DIR="$(ICGCSCOR_HOME_DIR) >> $@
	echo "ASPERA_HOME_DIR="$(ASPERA_HOME_DIR) >> $@
	echo "EGADECRYPT_HOME_DIR="$(EGADECRYPT_HOME_DIR) >> $@
	echo "DECSINGLE_HOME_DIR="$(DECSINGLE_HOME_DIR) >> $@
	echo "DECSINGLE_JAVA_HOME_DIR="$(DECSINGLE_JAVA_HOME_DIR) >> $@
	cat bam_ascatngs.sh >> $@

# Build using suffix rules (portable)
SUFFIXES= .sh .py .r .R

.sh :
	rm -f $@
	echo "#! " $(BASH) > $@
	echo "biopanpipe_pkgname=\""$(PACKAGE_NAME)"\"" >> $@
	echo "biopanpipe_version=\""$(PACKAGE_VERSION)"\"" >> $@
	echo "biopanpipe_bindir="$(prefix)/bin/ >> $@
	echo "biopanpipe_datadir="$(datadir)/$(PACKAGE_NAME) >> $@
	echo "BASH="$(BASH) >> $@
	echo "AWK="$(AWK) >> $@
	echo "SORT="$(SORT) >> $@
	echo "sortpars=\""$(sortpars)"\"" >> $@
	echo "sortT=\""$(sortT)"\"" >> $@
	echo "UNIQ="$(UNIQ) >> $@
	echo "GREP="$(GREP) >> $@
	echo "SED="$(SED) >> $@
	echo "BASENAME="$(BASENAME) >> $@
	echo "DIRNAME="$(DIRNAME) >> $@
	echo "SPLIT="$(SPLIT) >> $@
	echo "CUT="$(CUT) >> $@
	echo "HEAD="$(HEAD) >> $@
	echo "TAIL="$(TAIL) >> $@
	echo "PASTE="$(PASTE) >> $@
	echo "DIFF="$(DIFF) >> $@
	echo "FIND="$(FIND) >> $@
	echo "WC="$(WC) >> $@
	echo "MKTEMP="$(MKTEMP) >> $@
	echo "FLOCK="$(FLOCK) >> $@
	echo "ZCAT="$(ZCAT) >> $@
	echo "GZIP="$(GZIP) >> $@
	echo "WGET="$(WGET) >> $@
	echo "PANDOC="$(PANDOC) >> $@
	echo "PYTHON="$(PYTHON) >> $@
	echo "RSCRIPT="$(RSCRIPT) >> $@
	echo "JAVA="$(JAVA) >> $@
	echo "PANPIPE_HOME_DIR="$(PANPIPE_HOME_DIR) >> $@
	echo "PLATYPUS_HOME_DIR="$(PLATYPUS_HOME_DIR) >> $@
	echo "LUMPY_HOME_DIR="$(LUMPY_HOME_DIR) >> $@
	echo "ICGCSCOR_HOME_DIR="$(ICGCSCOR_HOME_DIR) >> $@
	echo "ASPERA_HOME_DIR="$(ASPERA_HOME_DIR) >> $@
	echo "EGADECRYPT_HOME_DIR="$(EGADECRYPT_HOME_DIR) >> $@
	echo "DECSINGLE_HOME_DIR="$(DECSINGLE_HOME_DIR) >> $@
	echo "DECSINGLE_JAVA_HOME_DIR="$(DECSINGLE_JAVA_HOME_DIR) >> $@
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
