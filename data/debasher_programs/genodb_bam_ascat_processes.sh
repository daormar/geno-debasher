# Geno-DeBasher package
# Copyright (C) 2019-2024 Daniel Ortiz-Mart\'inez and Jose
# Espinosa-Carrasco (initial version)
#
# This library is free software; you can redistribute it and/or
# modify it under the terms of the GNU Lesser General Public License
# as published by the Free Software Foundation; either version 3
# of the License, or (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with this program; If not, see <http://www.gnu.org/licenses/>.

# *- bash -*

###################
# ASCAT PROCESSES #
###################

########
extract_snp_ids_from_locis()
{
    local locis=$1
    "${AWK}" '{print $1}' "${locis}"
}

########
remove_snp_ids_from_locis()
{
    local locis=$1
    "${AWK}" '{
              for(i=2;i<=NF;++i)
              {
               printf"%s",$i
               if(i!=NF) printf"\t"
              }
              printf"\n"
            }' "$locis"
}

########
preproc_allelecounter_locis()
{
    local locis=$1
    local contig_mapping=$2

    if [ "${contig_mapping}" = ${NOFILE} ]; then
        remove_snp_ids_from_locis "${locis}"
    else
        remove_snp_ids_from_locis "${locis}" | "${genodebasher_libexecdir}"/genodb_map_contnames -m "${contig_mapping}" -c 0
    fi
}

########
postproc_allelecounter_output()
{
    local allelecounterfile=$1
    local contig_mapping=$2

    if [ "${contig_mapping}" = ${NOFILE} ]; then
        cat "${allelecounterfile}"
    else
        "${genodebasher_libexecdir}"/genodb_map_contnames -m "${contig_mapping}" -f "${allelecounterfile}" -c 0 --invert
    fi
}

########
allele_counter_norm_explain_opts()
{
    # -l option
    description="Loci (SNP position) file. IMPORTANT: Chromosome ids should not contain the 'chr' string prefix, first field represents the SNP id"
    explain_opt "-l" "<file>" "$description"

    # -r option
    description="Reference genome file"
    explain_opt "-r" "<file>" "$description"

    # -normalbam option
    description="Normal bam file (required if no downloading processes have been defined)"
    explain_opt "-normalbam" "<file>" "$description"

    # -ma option
    description="File containing a mapping from standard contig names expected by ASCAT (without the 'chr' prefix) to bam contig names"
    explain_opt "-ma" "<file>" "$description"

    # -process-outd option
    local description="output directory"
    explain_opt "-process-outd" "<file>" "$description"

    # -outcsv option
    description="Postprocessed alleleCounter output for the normal sample"
    explain_opt "-outcsv" "<file>" "$description"

    # -outsnpids option
    description="File with snp ids for alleleCounter output"
    explain_opt "-outsnpids" "<file>" "$description"
}

########
allele_counter_norm_identify_cmdline_opts()
{
    opt_is_cmdline "-l"
    opt_is_cmdline "-r"
    opt_is_cmdline "-normalbam"
    opt_is_cmdline "-ma"
}

########
allele_counter_norm_define_opts()
{
    # Initialize variables
    local cmdline=$1
    local process_spec=$2
    local process_name=$3
    local process_outdir=$4
    local optlist=""

    # Define the -process-outd option, the output directory for the process
    define_opt "-process-outd" "${process_outdir}" optlist || exit 1

    # -l option
    define_cmdline_infile_opt "$cmdline" "-l" optlist || exit 1

    # -r option
    define_cmdline_infile_opt "$cmdline" "-r" optlist || exit 1

    # -normalbam option
    local normalbam
    normalbam=`genodb_bam_common::get_normal_bam_filename "$cmdline"` || exit 1
    define_opt "-normalbam" "$normalbam" optlist || exit 1

    # -ma option
    define_cmdline_infile_opt_if_given "$cmdline" "-ma" optlist || exit 1

    # -outcsv option
    local outcsv="${process_outdir}"/allele_counter_norm_postproc.csv
    define_opt "-outcsv" "$outcsv" optlist || exit 1

    # -outsnpids option
    local outsnpids="${process_outdir}"/snpids
    define_opt "-outsnpids" "$outsnpids" optlist || exit 1

    # Save option list
    save_opt_list optlist
}

########
allele_counter_norm()
{
    # Initialize variables
    local process_outd=`read_opt_value_from_func_args "-process-outd" "$@"`
    local locis=`read_opt_value_from_func_args "-l" "$@"`
    local ref=`read_opt_value_from_func_args "-r" "$@"`
    local normalbam=`read_opt_value_from_func_args "-normalbam" "$@"`
    local contig_mapping=`read_opt_value_from_func_args "-ma" "$@"`
    if [ "$contig_mapping" = "${DEBASHER_OPT_NOT_FOUND}" ]; then
        contig_mapping=${NOFILE}
    fi
    local outcsv=`read_opt_value_from_func_args "-outcsv" "$@"`
    local outsnpids=`read_opt_value_from_func_args "-outsnpids" "$@"`

    # Extract SNP ids to a separate file
    logmsg "* Extracting SNP ids..."
    extract_snp_ids_from_locis "${locis}" > "${outsnpids}"

    # Create file for alleleCounter, mapping contigs if required
    logmsg "* Preprocessing locis..."
    preproc_allelecounter_locis "${locis}" "${contig_mapping}" > "${process_outd}"/allele_counter_norm.preproc_locis

    # Activate conda environment
    logmsg "* Activating conda environment..."
    conda activate allelecount 2>&1 || exit 1

    # Execute alleleCounter
    logmsg "* Executing alleleCounter..."
    alleleCounter -l "${process_outd}"/allele_counter_norm.preproc_locis -r "${ref}" -b "${normalbam}" -o "${process_outd}"/allele_counter_norm.csv 2>&1 || exit 1

    # Deactivate conda environment
    logmsg "* Dectivating conda environment..."
    conda deactivate 2>&1

    # Postprocess alleleCounter output
    logmsg "* Postprocessing alleleCounter output..."
    postproc_allelecounter_output "${process_outd}"/allele_counter_norm.csv "${contig_mapping}" > "${outcsv}"
}

########
allele_counter_norm_conda_envs()
{
    define_conda_env allelecount allelecount.yml
}

########
allele_counter_tumor_explain_opts()
{
    # -l option
    description="Loci (SNP position) file. IMPORTANT: Chromosome ids should not contain the 'chr' string prefix, first field represents the SNP id"
    explain_opt "-l" "<file>" "$description"

    # -r option
    description="Reference genome file"
    explain_opt "-r" "<file>" "$description"

    # -tumorbam option
    description="Tumor bam file (required if no downloading processes have been defined)"
    explain_opt "-tumorbam" "<file>" "$description"

    # -ma option
    description="File containing a mapping from standard contig names expected by ASCAT (without the 'chr' prefix) to bam contig names"
    explain_opt "-ma" "<file>" "$description"

    # -process-outd option
    local description="output directory"
    explain_opt "-process-outd" "<file>" "$description"

    # -outcsv option
    description="Postprocessed alleleCounter output for the tumor sample"
    explain_opt "-outcsv" "<file>" "$description"
}

########
allele_counter_tumor_identify_cmdline_opts()
{
    opt_is_cmdline "-l"
    opt_is_cmdline "-r"
    opt_is_cmdline "-tumorbam"
    opt_is_cmdline "-ma"
}

########
allele_counter_tumor_define_opts()
{
    # Initialize variables
    local cmdline=$1
    local process_spec=$2
    local process_name=$3
    local process_outdir=$4
    local optlist=""

    # Define the -process-outd option, the output directory for the process
    define_opt "-process-outd" "${process_outdir}" optlist || exit 1

    # -l option
    define_cmdline_infile_opt "$cmdline" "-l" optlist || exit 1

    # -r option
    define_cmdline_infile_opt "$cmdline" "-r" optlist || exit 1

    # -tumorbam option
    local tumorbam
    tumorbam=`genodb_bam_common::get_tumor_bam_filename "$cmdline"` || exit 1
    define_opt "-tumorbam" "$tumorbam" optlist || exit 1

    # -ma option
    define_cmdline_infile_opt_if_given "$cmdline" "-ma" optlist || exit 1

    # -outcsv option
    local outcsv="${process_outdir}"/allele_counter_tumor_postproc.csv
    define_opt "-outcsv" "$outcsv" optlist || exit 1

    # Save option list
    save_opt_list optlist
}

########
allele_counter_tumor()
{
    # Initialize variables
    local process_outd=`read_opt_value_from_func_args "-process-outd" "$@"`
    local locis=`read_opt_value_from_func_args "-l" "$@"`
    local ref=`read_opt_value_from_func_args "-r" "$@"`
    local tumorbam=`read_opt_value_from_func_args "-tumorbam" "$@"`
    local contig_mapping=`read_opt_value_from_func_args "-ma" "$@"`
    if [ "$contig_mapping" = "${DEBASHER_OPT_NOT_FOUND}" ]; then
        contig_mapping=${NOFILE}
    fi
    local outcsv=`read_opt_value_from_func_args "-outcsv" "$@"`

    # Extract SNP ids to a separate file
    logmsg "* Extracting SNP ids..."
    extract_snp_ids_from_locis "${locis}" > "${process_outd}"/snpids

    # Create file for alleleCounter, mapping contigs if required
    logmsg "* Preprocessing locis..."
    preproc_allelecounter_locis "${locis}" "${contig_mapping}" > "${process_outd}"/allele_counter_tumor.preproc_locis

    # Activate conda environment
    logmsg "* Activating conda environment..."
    conda activate allelecount 2>&1 || exit 1

    # Execute alleleCounter
    logmsg "* Executing alleleCounter..."
    alleleCounter -l "${process_outd}"/allele_counter_tumor.preproc_locis -r "${ref}" -b "${tumorbam}" -o "${process_outd}"/allele_counter_tumor.csv 2>&1 || exit 1

    # Deactivate conda environment
    logmsg "* Dectivating conda environment..."
    conda deactivate 2>&1

    # Postprocess alleleCounter output
    logmsg "* Postprocessing alleleCounter output..."
    postproc_allelecounter_output "${process_outd}"/allele_counter_tumor.csv "${contig_mapping}" > "${outcsv}"
}

########
allele_counter_tumor_conda_envs()
{
    define_conda_env allelecount allelecount.yml
}

########
ascat_explain_opts()
{
    # -acn option
    description="alleleCounter normal file"
    explain_opt "-acn" "<file>" "$description"

    # -act option
    description="alleleCounter tumor file"
    explain_opt "-act" "<file>" "$description"

    # -snpids option
    description="File with snp ids for alleleCounter output"
    explain_opt "-snpids" "<file>" "$description"

    # -g option
    description="Sample gender: XX|XY"
    explain_opt "-g" "<string>" "$description"

    # -sg option
    description="SNP GC correction file. IMPORTANT: Chromosome ids should not contain the 'chr' string prefix, first field represents the SNP id"
    explain_opt "-sg" "<file>" "$description"

    # -process-outd option
    local description="output directory"
    explain_opt "-process-outd" "<file>" "$description"
}

########
ascat_identify_cmdline_opts()
{
    opt_is_cmdline "-acn"
    opt_is_cmdline "-act"
    opt_is_cmdline "-snpids"
    opt_is_cmdline "-g"
    opt_is_cmdline "-sg"
}

########
ascat_define_opts()
{
    # Initialize variables
    local cmdline=$1
    local process_spec=$2
    local process_name=$3
    local process_outdir=$4
    local optlist=""

    # Define the -process-outd option, the output directory for the process
    define_opt "-process-outd" "${process_outdir}" optlist || exit 1

    # -acn option
    define_opt_from_proc_out "-acn" "allele_counter_norm" "-outcsv" optlist || exit 1

    # -act option
    define_opt_from_proc_out "-act" "allele_counter_tumor" "-outcsv" optlist || exit 1

    # -snpids option
    define_opt_from_proc_out "-snpids" "allele_counter_norm" "-outsnpids" optlist || exit 1

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

    "${HEAD}" -1 ${allc_outfile} | "${AWK}" '{printf"\t%s\n",$0}'
    "${PASTE}" ${snpids} <("${TAIL}" -n +2 "${allc_outfile}")
}

########
ascat()
{
    # Initialize variables
    local process_outd=`read_opt_value_from_func_args "-process-outd" "$@"`
    local allelecounternormal=`read_opt_value_from_func_args "-acn" "$@"`
    local allelecountertumor=`read_opt_value_from_func_args "-act" "$@"`
    local snpids=`read_opt_value_from_func_args "-snpids" "$@"`
    local gender=`read_opt_value_from_func_args "-g" "$@"`
    local snpgccorr=`read_opt_value_from_func_args "-sg" "$@"`

    # Activate conda environment
    logmsg "* Activating conda environment..."
    conda activate ascat 2>&1 || exit 1

    # Convert allele counts
    logmsg "* Executing convert_allele_counts..."
    Rscript "${genodebasher_libexecdir}"/genodb_convert_allele_counts "tumor" "${allelecountertumor}" "normal" "${allelecounternormal}" ${gender} "${process_outd}" || exit 1

    # Add SNP ids information to convert_allele_counts output
    logmsg "* Add SNP ids information to convert_allele_counts output..."
    add_snpids_to_convert_allele_counts_outfile "${snpids}" "${process_outd}"/tumor.BAF > "${process_outd}"/tumor_snpids.BAF
    add_snpids_to_convert_allele_counts_outfile "${snpids}" "${process_outd}"/tumor.LogR > "${process_outd}"/tumor_snpids.LogR
    add_snpids_to_convert_allele_counts_outfile "${snpids}" "${process_outd}"/normal.BAF > "${process_outd}"/normal_snpids.BAF
    add_snpids_to_convert_allele_counts_outfile "${snpids}" "${process_outd}"/normal.LogR > "${process_outd}"/normal_snpids.LogR

    # Run ascat
    logmsg "* Executing ascat..."
    Rscript "${genodebasher_libexecdir}"/genodb_run_ascat --tumor_baf="${process_outd}/tumor_snpids.BAF" --tumor_logr="${process_outd}/tumor_snpids.LogR" --normal_baf="${process_outd}/normal_snpids.BAF" --normal_logr="${process_outd}/tumor_snpids.LogR" --tumor_name="sample" --gc_correction=${snpgccorr} --out_dir="${process_outd}/" || exit 1

    # Deactivate conda environment
    logmsg "* Deactivating conda environment..."
    conda deactivate 2>&1
}

########
ascat_conda_envs()
{
    define_conda_env ascat ascat.yml
}
