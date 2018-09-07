SCRIPTS_TO_BE_BUILT= submit_bam_analysis get_analysis_status		\
convert_snppos_to_snpgcc create_snv_pos_ascat query_ega_metadata

bin_SCRIPTS = bam_utils_lib.sh $(SCRIPTS_TO_BE_BUILT)

CLEANFILES = $(SCRIPTS_TO_BE_BUILT)

EXTRA_DIST = bam_utils_lib.sh submit_bam_analysis.sh	\
get_analysis_status.sh convert_snppos_to_snpgcc.sh	\
create_snv_pos_ascat.py query_ega_metadata.py

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
	echo "SPLIT="$(SPLIT) >> $@
	echo "HEAD="$(HEAD) >> $@
	echo "TAIL="$(TAIL) >> $@
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
	echo "SAMTOOLS="$(SAMTOOLS) >> $@
	echo "PLATYPUS_BUILD_DIR="$(PLATYPUS_BUILD_DIR) >> $@
	echo "PYEGA_BUILD_DIR="$(PYEGA_BUILD_DIR) >> $@
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
