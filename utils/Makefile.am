SCRIPTS_TO_BE_BUILT= bam_analysis bam_facets bam_ascat_r bam_ascatngs	\
analyze_dataset gen_bed_for_genome check_ref_and_bam_contigs		\
gen_wisecondorx_ref convert_snppos_to_snpgcc create_snv_pos_ascat	\
query_ega_metadata query_icgc_metadata filter_contig_from_genref	\
run_facets run_sequenza convert_allele_counts run_ascat

bin_SCRIPTS = $(SCRIPTS_TO_BE_BUILT)

CLEANFILES = $(SCRIPTS_TO_BE_BUILT)

EXTRA_DIST = bam_analysis.sh bam_facets.sh bam_ascat_r.sh	\
bam_ascatngs.sh analyze_dataset.sh gen_bed_for_genome.sh	\
check_ref_and_bam_contigs.sh gen_wisecondorx_ref.sh		\
convert_snppos_to_snpgcc.sh create_snv_pos_ascat.py		\
query_ega_metadata.py query_icgc_metadata.py			\
filter_contig_from_genref.py run_facets.R run_sequenza.R	\
convert_allele_counts.R run_ascat.R

# Build rules for libraries
bam_analysis: bam_analysis.sh
	rm -f $@
	echo "biopanpipe_bindir="$(prefix)/bin/ >> $@
	echo "AWK="$(AWK) >> $@
	echo "GREP="$(GREP) >> $@
	echo "SED="$(SED) >> $@
	echo "GZIP="$(GZIP) >> $@
	echo "BASENAME="$(BASENAME) >> $@
	echo "DIRNAME="$(DIRNAME) >> $@
	echo "PLATYPUS_HOME_DIR="$(PLATYPUS_HOME_DIR) >> $@
	echo "FACETS_HOME_DIR="$(FACETS_HOME_DIR) >> $@
	echo "ICGCSTOR_HOME_DIR="$(ICGCSTOR_HOME_DIR) >> $@
	echo "ASPERA_HOME_DIR="$(ASPERA_HOME_DIR) >> $@
	echo "EGADECRYPT_HOME_DIR="$(EGADECRYPT_HOME_DIR) >> $@
	cat bam_analysis.sh >> $@

bam_facets: bam_facets.sh
	rm -f $@
	echo "biopanpipe_bindir="$(prefix)/bin/ >> $@
	echo "AWK="$(AWK) >> $@
	echo "GREP="$(GREP) >> $@
	echo "SED="$(SED) >> $@
	echo "GZIP="$(GZIP) >> $@
	echo "BASENAME="$(BASENAME) >> $@
	echo "DIRNAME="$(DIRNAME) >> $@
	echo "PLATYPUS_HOME_DIR="$(PLATYPUS_HOME_DIR) >> $@
	echo "FACETS_HOME_DIR="$(FACETS_HOME_DIR) >> $@
	echo "ICGCSTOR_HOME_DIR="$(ICGCSTOR_HOME_DIR) >> $@
	echo "ASPERA_HOME_DIR="$(ASPERA_HOME_DIR) >> $@
	echo "EGADECRYPT_HOME_DIR="$(EGADECRYPT_HOME_DIR) >> $@
	cat bam_facets.sh >> $@

bam_ascat_r: bam_ascat_r.sh
	rm -f $@
	echo "biopanpipe_bindir="$(prefix)/bin/ >> $@
	echo "AWK="$(AWK) >> $@
	echo "GREP="$(GREP) >> $@
	echo "SED="$(SED) >> $@
	echo "GZIP="$(GZIP) >> $@
	echo "BASENAME="$(BASENAME) >> $@
	echo "DIRNAME="$(DIRNAME) >> $@
	echo "PLATYPUS_HOME_DIR="$(PLATYPUS_HOME_DIR) >> $@
	echo "FACETS_HOME_DIR="$(FACETS_HOME_DIR) >> $@
	echo "ICGCSTOR_HOME_DIR="$(ICGCSTOR_HOME_DIR) >> $@
	echo "ASPERA_HOME_DIR="$(ASPERA_HOME_DIR) >> $@
	echo "EGADECRYPT_HOME_DIR="$(EGADECRYPT_HOME_DIR) >> $@
	cat bam_ascat_r.sh >> $@

bam_ascat_ngs: bam_ascat_r.sh
	rm -f $@
	echo "biopanpipe_bindir="$(prefix)/bin/ >> $@
	echo "AWK="$(AWK) >> $@
	echo "GREP="$(GREP) >> $@
	echo "SED="$(SED) >> $@
	echo "GZIP="$(GZIP) >> $@
	echo "BASENAME="$(BASENAME) >> $@
	echo "DIRNAME="$(DIRNAME) >> $@
	echo "PLATYPUS_HOME_DIR="$(PLATYPUS_HOME_DIR) >> $@
	echo "FACETS_HOME_DIR="$(FACETS_HOME_DIR) >> $@
	echo "ICGCSTOR_HOME_DIR="$(ICGCSTOR_HOME_DIR) >> $@
	echo "ASPERA_HOME_DIR="$(ASPERA_HOME_DIR) >> $@
	echo "EGADECRYPT_HOME_DIR="$(EGADECRYPT_HOME_DIR) >> $@
	cat bam_ascat_r.sh >> $@

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
	echo "HEAD="$(HEAD) >> $@
	echo "TAIL="$(TAIL) >> $@
	echo "DIFF="$(DIFF) >> $@
	echo "FIND="$(FIND) >> $@
	echo "WC="$(WC) >> $@
	echo "MKTEMP="$(MKTEMP) >> $@
	echo "FLOCK="$(FLOCK) >> $@
	echo "GZIP="$(GZIP) >> $@
	echo "WGET="$(WGET) >> $@
	echo "PANDOC="$(PANDOC) >> $@
	echo "PYTHON="$(PYTHON) >> $@
	echo "RSCRIPT="$(RSCRIPT) >> $@
	echo "JAVA="$(JAVA) >> $@
	echo "PANPIPE_HOME_DIR="$(PANPIPE_HOME_DIR) >> $@
	echo "PLATYPUS_HOME_DIR="$(PLATYPUS_HOME_DIR) >> $@
	echo "FACETS_HOME_DIR="$(FACETS_HOME_DIR) >> $@
	echo "ICGCSTOR_HOME_DIR="$(ICGCSTOR_HOME_DIR) >> $@
	echo "ASPERA_HOME_DIR="$(ASPERA_HOME_DIR) >> $@
	echo "EGADECRYPT_HOME_DIR="$(EGADECRYPT_HOME_DIR) >> $@
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
