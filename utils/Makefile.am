SCRIPTS_TO_BE_BUILT= bam_utils_lib submit_bam_analysis analyze_dataset	\
get_analysis_status gen_bed_for_genome gen_wisecondorx_ref		\
convert_snppos_to_snpgcc create_snv_pos_ascat query_ega_metadata	\
query_icgc_metadata filter_contig_from_genref run_facets

bin_SCRIPTS = $(SCRIPTS_TO_BE_BUILT)

CLEANFILES = $(SCRIPTS_TO_BE_BUILT)

EXTRA_DIST = bam_utils_lib.sh submit_bam_analysis.sh analyze_dataset.sh	\
get_analysis_status.sh gen_bed_for_genome.sh gen_wisecondorx_ref.sh	\
convert_snppos_to_snpgcc.sh create_snv_pos_ascat.py			\
query_ega_metadata.py query_icgc_metadata.py				\
filter_contig_from_genref.py run_facets.R

# Build rule for shell library
bam_utils_lib: bam_utils_lib.sh
	rm -f $@
	$(PYTHON) -c "import os; t=os.path.getmtime('bam_utils_lib.sh'); print '# Lib time stamp:',int(t)" > $@
	cat bam_utils_lib.sh >> $@

# Build using suffix rules (portable)
SUFFIXES= .sh .py .r .R

.sh :
	rm -f $@
	echo "#! " $(BASH) > $@
	echo "pkgname=\""$(PACKAGE_NAME)"\"" >> $@
	echo "version=\""$(PACKAGE_VERSION)"\"" >> $@
	echo "bindir="$(prefix)/bin/ >> $@
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
	echo "MKTEMP="$(MKTEMP) >> $@
	echo "QSUB=\""$(QSUB)"\"" >> $@
	echo "QSUB_WORKS=\""$(QSUB_WORKS)"\"" >> $@
	echo "QSUB_TERSE_OPT=\""$(QSUB_TERSE_OPT)"\"" >> $@
	echo "QSTAT="$(QSTAT) >> $@
	echo "QSTAT_J_OPT=\""$(QSTAT_J_OPT)"\"" >> $@
	echo "QRLS="$(QRLS) >> $@
	echo "SRUN="$(SRUN) >> $@
	echo "SBATCH="$(SBATCH) >> $@
	echo "WGET="$(WGET) >> $@
	echo "PANDOC="$(PANDOC) >> $@
	echo "PYTHON="$(PYTHON) >> $@
	echo "JAVA="$(JAVA) >> $@
	echo "SAMTOOLS="$(SAMTOOLS) >> $@
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
