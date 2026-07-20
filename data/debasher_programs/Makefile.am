genodbfilesdir = $(datadir)/$(PACKAGE)/debasher_programs
genodbfiles_DATA= genodb_bam_analysis.sh genodb_genref_processes.sh	\
genodb_bam_download_processes.sh genodb_bam_manip_processes.sh		\
genodb_cleaning_processes.sh genodb_bam_analysis_processes.sh		\
genodb_bam_summariz_processes.sh genodb_bam_filter_processes.sh		\
genodb_bam_facets_processes.sh genodb_bam_ascat_processes.sh		\
genodb_nodownload_basic.sh genodb_nodownload_extended.sh		\
genodb_nodownload_xgs.sh genodb_ega_basic.sh genodb_ega_asp_basic.sh	\
genodb_aws_basic.sh genodb_collab_basic.sh genodb_ega_asp_extended.sh	\
genodb_ega_decsingle_extended.sh genodb_copy_extended.sh		\
genodb_copy_harmoniz.sh genodb_copy_debug.sh genodb_gdc_extended.sh	\
genodb_tcga_xgs.sh

EXTRA_DIST= $(genodbfiles_DATA)

genodb_bam_analysis.sh: genodb_bam_analysis_pred.sh
	echo "genodebasher_pkglibdir=\"$(libdir)/$(PACKAGE_NAME)\"" >> $@
	echo "genodebasher_debasher_programs_dir=\"$(datadir)/$(PACKAGE_NAME)/debasher_programs\"" >> $@
	cat $< >> $@
