# *- bash -*

#################
# CFG FUNCTIONS #
#################

########
bam_ascat_shared_dirs()
{
    define_shared_dir ${DATADIR_BASENAME}
}

########
bam_ascat_fifos()
{
    :
}

###############
# ASCAT STEPS #
###############

########
extract_snp_ids_from_locis()
{
    local locis=$1
    ${AWK} '{print $1}' ${locis}
}

########
remove_snp_ids_from_locis()
{
    local locis=$1
    ${AWK} '{
              for(i=2;i<=NF;++i)
              {
               printf"%s",$i 
               if(i!=NF) printf"\t"
              } 
              printf"\n"
            }' $locis
}

########
preproc_allelecounter_locis()
{
    local locis=$1
    local contig_mapping=$2

    if [ ${contig_mapping} = ${NOFILE} ]; then
        remove_snp_ids_from_locis ${locis}
    else
        remove_snp_ids_from_locis ${locis} | ${biopanpipe_bindir}/map_contnames -m ${contig_mapping} -c 0
    fi
}

########
postproc_allelecounter_output()
{
    local allelecounterfile=$1
    local contig_mapping=$2

    if [ ${contig_mapping} = ${NOFILE} ]; then
        cat ${allelecounterfile}        
    else
        ${biopanpipe_bindir}/map_contnames -m ${contig_mapping} -f ${allelecounterfile} -c 0 --invert
    fi
}

########
allele_counter_norm_explain_cmdline_opts()
{
    # -l option
    description="Loci (SNP position) file. IMPORTANT: Chromosome ids should not contain the 'chr' string prefix, first field represents the SNP id"
    explain_cmdline_req_opt "-l" "<string>" "$description"

    # -r option
    description="Reference genome file"
    explain_cmdline_req_opt "-r" "<string>" "$description"

    # -n option
    description="Normal bam file (required if no downloading steps have been defined)"
    explain_cmdline_opt "-n" "<string>" "$description"

    # -ma option
    description="File containing a mapping between standard contig names expected by ASCAT (without the 'chr' prefix) and bam contig names"
    explain_cmdline_opt "-ma" "<string>" "$description"
}

########
allele_counter_norm_define_opts()
{
    # Initialize variables
    local cmdline=$1
    local jobspec=$2
    local optlist=""

    # Define the -step-outd option, the output directory for the step
    local step_outd=`get_step_outdir_given_stepspec "$stepspec"`
    define_opt "-step-outd" ${step_outd} optlist || exit 1

    # -l option
    define_cmdline_infile_opt "$cmdline" "-l" optlist || exit 1

    # -r option
    define_cmdline_infile_opt "$cmdline" "-r" optlist || exit 1

    # -normalbam option
    local normalbam
    normalbam=`get_normal_bam_filename "$cmdline"` || exit 1
    define_opt "-normalbam" $normalbam optlist || exit 1

    # -ma option
    define_cmdline_infile_nonmand_opt "$cmdline" "-ma" ${NOFILE} optlist || exit 1

    # Save option list
    save_opt_list optlist
}

########
allele_counter_norm()
{
    # Initialize variables
    local step_outd=`read_opt_value_from_line "$*" "-step-outd"`
    local locis=`read_opt_value_from_line "$*" "-l"`
    local ref=`read_opt_value_from_line "$*" "-r"`
    local normalbam=`read_opt_value_from_line "$*" "-normalbam"`
    local contig_mapping=`read_opt_value_from_line "$*" "-ma"`

    # Extract SNP ids to a separate file
    logmsg "* Extracting SNP ids..."
    extract_snp_ids_from_locis ${locis} > ${step_outd}/snpids
    
    # Create file for alleleCounter, mapping contigs if required
    logmsg "* Preprocessing locis..."
    preproc_allelecounter_locis ${locis} ${contig_mapping} > ${step_outd}/allele_counter_norm.preproc_locis
    
    # Activate conda environment    
    logmsg "* Activating conda environment..."
    conda activate allelecount 2>&1 || exit 1

    # Execute alleleCounter
    logmsg "* Executing alleleCounter..."
    alleleCounter -l ${step_outd}/allele_counter_norm.preproc_locis -r ${ref} -b ${normalbam} -o ${step_outd}/allele_counter_norm.csv 2>&1 || exit 1
    
    # Deactivate conda environment
    logmsg "* Dectivating conda environment..."
    conda deactivate 2>&1

    # Postprocess alleleCounter output
    logmsg "* Postprocessing alleleCounter output..."
    postproc_allelecounter_output ${step_outd}/allele_counter_norm.csv ${contig_mapping} > ${step_outd}/allele_counter_norm_postproc.csv
}

########
allele_counter_norm_conda_envs()
{
    define_conda_env allelecount allelecount.yml
}

########
allele_counter_tumor_explain_cmdline_opts()
{
    # -l option
    description="Loci (SNP position) file. IMPORTANT: Chromosome ids should not contain the 'chr' string prefix, first field represents the SNP id"
    explain_cmdline_req_opt "-l" "<string>" "$description"

    # -r option
    description="Reference genome file"
    explain_cmdline_req_opt "-r" "<string>" "$description"

    # -t option
    description="Tumor bam file (required if no downloading steps have been defined)"
    explain_cmdline_opt "-t" "<string>" "$description"

    # -ma option
    description="File containing a mapping between standard contig names expected by ASCAT (without the 'chr' prefix) and bam contig names"
    explain_cmdline_opt "-ma" "<string>" "$description"
}

########
allele_counter_tumor_define_opts()
{
    # Initialize variables
    local cmdline=$1
    local jobspec=$2
    local optlist=""

    # Define the -step-outd option, the output directory for the step
    local step_outd=`get_step_outdir_given_stepspec "$stepspec"`
    define_opt "-step-outd" ${step_outd} optlist || exit 1

    # -l option
    define_cmdline_infile_opt "$cmdline" "-l" optlist || exit 1

    # -r option
    define_cmdline_infile_opt "$cmdline" "-r" optlist || exit 1

    # -tumorbam option
    local normalbam
    tumorbam=`get_tumor_bam_filename "$cmdline"` || exit 1
    define_opt "-tumorbam" $tumorbam optlist || exit 1

    # -ma option
    define_cmdline_infile_nonmand_opt "$cmdline" "-ma" ${NOFILE} optlist || exit 1

    # Save option list
    save_opt_list optlist
}

########
allele_counter_tumor()
{
    # Initialize variables
    local step_outd=`read_opt_value_from_line "$*" "-step-outd"`
    local locis=`read_opt_value_from_line "$*" "-l"`
    local ref=`read_opt_value_from_line "$*" "-r"`
    local tumorbam=`read_opt_value_from_line "$*" "-tumorbam"`
    local contig_mapping=`read_opt_value_from_line "$*" "-ma"`

    # Extract SNP ids to a separate file
    logmsg "* Extracting SNP ids..."
    extract_snp_ids_from_locis ${locis} > ${step_outd}/snpids
    
    # Create file for alleleCounter, mapping contigs if required
    logmsg "* Preprocessing locis..."
    preproc_allelecounter_locis ${locis} ${contig_mapping} > ${step_outd}/allele_counter_tumor.preproc_locis

    # Activate conda environment
    logmsg "* Activating conda environment..."
    conda activate allelecount 2>&1 || exit 1
    
    # Execute alleleCounter
    logmsg "* Executing alleleCounter..."
    alleleCounter -l ${step_outd}/allele_counter_tumor.preproc_locis -r ${ref} -b ${tumorbam} -o ${step_outd}/allele_counter_tumor.csv 2>&1 || exit 1

    # Deactivate conda environment
    logmsg "* Dectivating conda environment..."
    conda deactivate 2>&1

    # Postprocess alleleCounter output
    logmsg "* Postprocessing alleleCounter output..."
    postproc_allelecounter_output ${step_outd}/allele_counter_tumor.csv ${contig_mapping} > ${step_outd}/allele_counter_tumor_postproc.csv
}

########
allele_counter_tumor_conda_envs()
{
    define_conda_env allelecount allelecount.yml
}

########
ascat_explain_cmdline_opts()
{
    # -acn option
    description="alleleCounter normal file"
    explain_cmdline_opt "-acn" "<string>" "$description"

    # -act option
    description="alleleCounter tumor file"
    explain_cmdline_opt "-act" "<string>" "$description"

    # -snpids option
    description="File with snp ids for alleleCounter output"
    explain_cmdline_opt "-snpids" "<string>" "$description"

    # -g option
    description="Sample gender: XX|XY"
    explain_cmdline_req_opt "-g" "<string>" "$description"

    # -sg option
    description="SNP GC correction file. IMPORTANT: Chromosome ids should not contain the 'chr' string prefix, first field represents the SNP id"
    explain_cmdline_req_opt "-sg" "<string>" "$description"
}

########
ascat_define_opts()
{
    # Initialize variables
    local cmdline=$1
    local jobspec=$2
    local optlist=""
    
    # Define the -step-outd option, the output directory for the step
    local step_outd=`get_step_outdir_given_stepspec "$stepspec"`
    define_opt "-step-outd" ${step_outd} optlist || exit 1

    # Define -acn option or retrieve dependency
    local allelecountnorm_dep=`find_dependency_for_step "${jobspec}" allele_counter_norm`
    if [ ${allelecountnorm_dep} = ${DEP_NOT_FOUND} ]; then
        define_cmdline_infile_opt "${cmdline}" "-acn"  optlist || exit 1
    else
        local acnorm_outd=`get_outd_for_dep "${allelecountnorm_dep}"`
        local allelecount_norm_file=${acnorm_outd}/allele_counter_norm_postproc.csv
        define_opt "-acn" ${allelecount_norm_file} optlist || exit 1
    fi
    
    # Define -act option or retrieve dependency
    local allelecounttumor_dep=`find_dependency_for_step "${jobspec}" allele_counter_tumor`
    if [ ${allelecounttumor_dep} = ${DEP_NOT_FOUND} ]; then
        define_cmdline_infile_opt "${cmdline}" "-act"  optlist || exit 1
    else
        local actumor_outd=`get_outd_for_dep "${allelecounttumor_dep}"`
        local allelecount_tumor_file=${acnorm_outd}/allele_counter_tumor_postproc.csv
        define_opt "-act" ${allelecount_norm_file} optlist || exit 1
    fi

    # Define -snpids option or retrieve dependency
    local allelecountnorm_dep=`find_dependency_for_step "${jobspec}" allele_counter_norm`
    if [ ${allelecountnorm_dep} = ${DEP_NOT_FOUND} ]; then
        define_cmdline_infile_opt "${cmdline}" "-snpids"  optlist || exit 1
    else
        local acnorm_outd=`get_outd_for_dep "${allelecountnorm_dep}"`
        local snpids_file=${acnorm_outd}/snpids
        define_opt "-snpids" ${snpids_file} optlist || exit 1
    fi

    # -g option
    define_cmdline_opt "$cmdline" "-g" optlist || exit 1

    # -sg option
    define_cmdline_infile_opt "$cmdline" "-sg" optlist || exit 1
 
    # Save option list
    save_opt_list optlist
}

########
add_snpids_to_convert_allele_counts_outfile()
{
    local snpids=$1
    local allc_outfile=$2

    ${HEAD} -1 ${allc_outfile} | ${AWK} '{printf"\t%s",$0}'
    ${PASTE} ${snpids} <(${TAIL} -n +2 ${allc_outfile})
}

########
ascat()
{
    # Initialize variables
    local step_outd=`read_opt_value_from_line "$*" "-step-outd"`
    local allelecounternormal=`read_opt_value_from_line "$*" "-acn"`
    local allelecountertumor=`read_opt_value_from_line "$*" "-act"`
    local snpids=`read_opt_value_from_line "$*" "-snpids"`
    local gender=`read_opt_value_from_line "$*" "-g"`
    local snpgccorr=`read_opt_value_from_line "$*" "-sg"`

    # Activate conda environment
    logmsg "* Activating conda environment..."
    conda activate ascat 2>&1 || exit 1
    
    # Convert allele counts
    logmsg "* Executing convert_allele_counts..."
    Rscript ${biopanpipe_bindir}/convert_allele_counts "tumor" ${allelecountertumor} "normal" ${allelecounternormal} ${gender} ${step_outd} || exit 1

    # Add SNP ids information to convert_allele_counts output
    logmsg "* Add SNP ids information to convert_allele_counts output..."
    add_snpids_to_convert_allele_counts_outfile ${snpids} ${step_outd}/tumor.BAF > ${step_outd}/tumor_snpids.BAF
    add_snpids_to_convert_allele_counts_outfile ${snpids} ${step_outd}/tumor.LogR > ${step_outd}/tumor_snpids.LogR
    add_snpids_to_convert_allele_counts_outfile ${snpids} ${step_outd}/normal.BAF > ${step_outd}/normal_snpids.BAF
    add_snpids_to_convert_allele_counts_outfile ${snpids} ${step_outd}/normal.LogR > ${step_outd}/normal_snpids.LogR
    
    # Run ascat
    logmsg "* Executing run_ascat..."
    Rscript ${biopanpipe_bindir}/run_ascat --tumor_baf="${step_outd}/tumor_snpids.BAF" --tumor_logr="${step_outd}/tumor_snpids.LogR" --normal_baf="${step_outd}/normal_snpids.BAF" --normal_logr="${step_outd}/tumor_snpids.LogR" --tumor_name="tumor" --gc_correction=${snpgccorr} --out_dir="${step_outd}/" || exit 1

    # Deactivate conda environment
    logmsg "* Deactivating conda environment..."
    conda deactivate 2>&1
}

########
ascat_conda_envs()
{
    define_conda_env ascat ascat.yml
}
