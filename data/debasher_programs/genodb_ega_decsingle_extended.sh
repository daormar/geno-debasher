# *- bash -*

load_debasher_module "genodb_bam_analysis"

########
genodb.index_norm_bam_define_opts()
{
    # Initialize variables
    local cmdline=$1
    local process_spec=$2
    local optlist=""

    # -normalbam option
    define_opt_from_proc_out "-normalbam" "genodb.decsingle_ega_norm_bam" "-out-nb" optlist || return 1

    # -out-nbidx option
    local abs_datadir=$(get_absolute_shdirname "data")
    define_opt "-out-nbidx" "${abs_datadir}/normal.bam.bai" optlist || return 1

    # -out-nb option (republished once indexed, so downstream processes
    # can connect to it and depend on indexing having completed)
    define_opt "-out-nb" "${abs_datadir}/normal.bam" optlist || return 1

    # Save option list
    save_opt_list optlist
}

########
genodb.index_tum_bam_define_opts()
{
    # Initialize variables
    local cmdline=$1
    local process_spec=$2
    local optlist=""

    # -tumorbam option
    define_opt_from_proc_out "-tumorbam" "genodb.decsingle_ega_tum_bam" "-out-tb" optlist || return 1

    # -out-tbidx option
    local abs_datadir=$(get_absolute_shdirname "data")
    define_opt "-out-tbidx" "${abs_datadir}/tumor.bam.bai" optlist || return 1

    # -out-tb option (republished once indexed, so downstream processes
    # can connect to it and depend on indexing having completed)
    define_opt "-out-tb" "${abs_datadir}/tumor.bam" optlist || return 1

    # Save option list
    save_opt_list optlist
}

########
genodb.bedtools_genomecov_norm_bam_identify_cmdline_opts()
{
    :
}

########
genodb.bedtools_genomecov_norm_bam_define_opts()
{
    # Initialize variables
    local cmdline=$1
    local process_spec=$2
    local process_name=$3
    local process_outdir=$4
    local optlist=""

    # Define the -out-processdir option, the output directory for the process
    define_opt "-out-processdir" "${process_outdir}" optlist || return 1

    # -normalbam option
    define_opt_from_proc_out "-normalbam" "genodb.index_norm_bam" "-out-nb" optlist || return 1

    # Save option list
    save_opt_list optlist
}

########
genodb.bedtools_genomecov_tum_bam_identify_cmdline_opts()
{
    :
}

########
genodb.bedtools_genomecov_tum_bam_define_opts()
{
    # Initialize variables
    local cmdline=$1
    local process_spec=$2
    local process_name=$3
    local process_outdir=$4
    local optlist=""

    # Define the -out-processdir option, the output directory for the process
    define_opt "-out-processdir" "${process_outdir}" optlist || return 1

    # -tumorbam option
    define_opt_from_proc_out "-tumorbam" "genodb.index_tum_bam" "-out-tb" optlist || return 1

    # Save option list
    save_opt_list optlist
}

########
genodb.create_genref_for_bam_identify_cmdline_opts()
{
    opt_is_cmdline "-br"
    opt_is_cmdline "-cm"
    opt_is_cmdline "-fbr"
}

########
genodb.create_genref_for_bam_define_opts()
{
    # Initialize variables
    local cmdline=$1
    local process_spec=$2
    local process_name=$3
    local process_outdir=$4
    local optlist=""

    # Define the -out-processdir option, the output directory for the process
    define_opt "-out-processdir" "${process_outdir}" optlist || return 1

    # -br option
    define_cmdline_infile_opt "$cmdline" "-br" optlist || return 1

    # -bam option
    define_opt_from_proc_out "-bam" "genodb.index_norm_bam" "-out-nb" optlist || return 1

    # -cm option
    define_cmdline_infile_opt_if_given "$cmdline" "-cm" optlist || return 1

    # -fbr option
    define_cmdline_infile_opt_if_given "$cmdline" "-fbr" optlist || return 1

    # Get data directory
    local abs_datadir=$(get_absolute_shdirname "data")

    # -outfile option
    local outfile="${abs_datadir}"/genref.fa
    define_opt "-outfile" "$outfile" optlist || return 1

    # Save option list
    save_opt_list optlist
}

########
genodb.parallel_split_norm_bam_identify_cmdline_opts()
{
    opt_is_cmdline "-lc"
}

########
genodb.parallel_split_norm_bam_define_opts()
{
    # Initialize variables
    local cmdline=$1
    local process_spec=$2
    local process_name=$3
    local process_outdir=$4

    # Obtain splitdir directory
    local abs_splitdir=$(get_absolute_shdirname "split")

    # Get name of contig list file
    local clist
    clist=$(read_opt_value_from_line "$cmdline" "-lc") || { errmsg "Error: -lc option not found"; return 1; }

    # Array of contigs to process, one task per contig
    array=( $(get_contig_list_from_file $clist) ) || return 1

    for idx in "${!array[@]}"; do
        local optlist=""

        # -normalbam option
        define_opt_from_proc_out "-normalbam" "genodb.index_norm_bam" "-out-nb" optlist || return 1

        # -contig option
        local contig=${array[$idx]}
        define_opt "-contig" "$contig" optlist || return 1

        # -outfile option
        local outfile="${abs_splitdir}"/normal_${contig}.bam
        define_opt "-outfile" "$outfile" optlist || return 1

        # Save option list
        save_opt_list optlist
    done
}

########
genodb.parallel_split_tum_bam_identify_cmdline_opts()
{
    opt_is_cmdline "-lc"
}

########
genodb.parallel_split_tum_bam_define_opts()
{
    # Initialize variables
    local cmdline=$1
    local process_spec=$2
    local process_name=$3
    local process_outdir=$4

    # Obtain splitdir directory
    local abs_splitdir=$(get_absolute_shdirname "split")

    # Get name of contig list file
    local clist
    clist=$(read_opt_value_from_line "$cmdline" "-lc") || { errmsg "Error: -lc option not found"; return 1; }

    # Array of contigs to process, one task per contig
    array=( $(get_contig_list_from_file $clist) ) || return 1

    for idx in "${!array[@]}"; do
        local optlist=""

        # -tumorbam option
        define_opt_from_proc_out "-tumorbam" "genodb.index_tum_bam" "-out-tb" optlist || return 1

        # -contig option
        local contig=${array[$idx]}
        define_opt "-contig" "$contig" optlist || return 1

        # -outfile option
        local outfile="${abs_splitdir}"/tumor_${contig}.bam
        define_opt "-outfile" "$outfile" optlist || return 1

        # Save option list
        save_opt_list optlist
    done
}

########
genodb.parallel_samtools_mpileup_norm_bam_identify_cmdline_opts()
{
    opt_is_cmdline "-mpb"
    opt_is_cmdline "-lc"
}

########
genodb.parallel_samtools_mpileup_norm_bam_define_opts()
{
    # Initialize variables
    local cmdline=$1
    local process_spec=$2
    local process_name=$3
    local process_outdir=$4

    # Get name of contig list file
    local clist
    clist=$(read_opt_value_from_line "$cmdline" "-lc") || { errmsg "Error: -lc option not found"; return 1; }

    # Array of contigs to process, one task per contig
    array=( $(get_contig_list_from_file $clist) ) || return 1

    for idx in "${!array[@]}"; do
        local optlist=""

        # -out-processdir option, the output directory for the process
        define_opt "-out-processdir" "${process_outdir}" optlist || return 1

        # -r option
        define_opt_from_proc_out "-r" "genodb.create_genref_for_bam" "-outfile" optlist || return 1

        # -mpb option
        define_cmdline_opt_if_given "$cmdline" "-mpb" optlist

        # -normalbam option
        local contig=${array[$idx]}
        define_opt_from_proc_task_out "-normalbam" "genodb.parallel_split_norm_bam" "${idx}" "-outfile" optlist || return 1

        # -contig option
        define_opt "-contig" "$contig" optlist || return 1

        # -outfile option
        local outfile="${process_outdir}"/normal_${contig}.pileup.gz
        define_opt "-outfile" "$outfile" optlist || return 1

        # Save option list
        save_opt_list optlist
    done
}

########
genodb.parallel_samtools_mpileup_tum_bam_identify_cmdline_opts()
{
    opt_is_cmdline "-mpb"
    opt_is_cmdline "-lc"
}

########
genodb.parallel_samtools_mpileup_tum_bam_define_opts()
{
    # Initialize variables
    local cmdline=$1
    local process_spec=$2
    local process_name=$3
    local process_outdir=$4

    # Get name of contig list file
    local clist
    clist=$(read_opt_value_from_line "$cmdline" "-lc") || { errmsg "Error: -lc option not found"; return 1; }

    # Array of contigs to process, one task per contig
    array=( $(get_contig_list_from_file $clist) ) || return 1

    for idx in "${!array[@]}"; do
        local optlist=""

        # -out-processdir option, the output directory for the process
        define_opt "-out-processdir" "${process_outdir}" optlist || return 1

        # -r option
        define_opt_from_proc_out "-r" "genodb.create_genref_for_bam" "-outfile" optlist || return 1

        # -mpb option
        define_cmdline_opt_if_given "$cmdline" "-mpb" optlist

        # -tumorbam option
        local contig=${array[$idx]}
        define_opt_from_proc_task_out "-tumorbam" "genodb.parallel_split_tum_bam" "${idx}" "-outfile" optlist || return 1

        # -contig option
        define_opt "-contig" "$contig" optlist || return 1

        # -outfile option
        local outfile="${process_outdir}"/tumor_${contig}.pileup.gz
        define_opt "-outfile" "$outfile" optlist || return 1

        # Save option list
        save_opt_list optlist
    done
}

########
genodb.gen_sequenza_gcc_identify_cmdline_opts()
{
    :
}

########
genodb.gen_sequenza_gcc_define_opts()
{
    # Initialize variables
    local cmdline=$1
    local process_spec=$2
    local process_name=$3
    local process_outdir=$4
    local optlist=""

    # -r option
    define_opt_from_proc_out "-r" "genodb.create_genref_for_bam" "-outfile" optlist || return 1

    # Get data directory
    local abs_datadir=$(get_absolute_shdirname "data")

    # -outfile option
    define_opt "-outfile" "${abs_datadir}"/sequenza_gccfile.txt.gz optlist || return 1

    # Save option list
    save_opt_list optlist
}

########
genodb.parallel_bam2seqz_identify_cmdline_opts()
{
    opt_is_cmdline "-lc"
}

########
genodb.parallel_bam2seqz_define_opts()
{
    # Initialize variables
    local cmdline=$1
    local process_spec=$2
    local process_name=$3
    local process_outdir=$4

    # Get name of contig list file
    local clist
    clist=$(read_opt_value_from_line "$cmdline" "-lc") || { errmsg "Error: -lc option not found"; return 1; }

    # Array of contigs to process, one task per contig
    array=( $(get_contig_list_from_file $clist) ) || return 1

    for idx in "${!array[@]}"; do
        local optlist=""

        # -out-processdir option, the output directory for the process
        define_opt "-out-processdir" "${process_outdir}" optlist || return 1

        # -gcc option
        define_opt_from_proc_out "-gcc" "genodb.gen_sequenza_gcc" "-outfile" optlist || return 1

        # -npileup option
        define_opt_from_proc_task_out "-npileup" "genodb.parallel_samtools_mpileup_norm_bam" "${idx}" "-outfile" optlist || return 1

        # -tpileup option
        define_opt_from_proc_task_out "-tpileup" "genodb.parallel_samtools_mpileup_tum_bam" "${idx}" "-outfile" optlist || return 1

        # -contig option
        local contig=${array[$idx]}
        define_opt "-contig" "$contig" optlist || return 1

        # -outfile option
        local outfile="${process_outdir}"/${contig}_seqz.gz
        define_opt "-outfile" "$outfile" optlist || return 1

        # Save option list
        save_opt_list optlist
    done
}

########
genodb.parallel_delly_identify_cmdline_opts()
{
    opt_is_cmdline "-dx"
    opt_is_cmdline "-lc"
}

########
genodb.parallel_delly_define_opts()
{
    # Initialize variables
    local cmdline=$1
    local process_spec=$2
    local process_name=$3
    local process_outdir=$4

    # Get name of contig list file
    local clist
    clist=$(read_opt_value_from_line "$cmdline" "-lc") || { errmsg "Error: -lc option not found"; return 1; }

    # Array of contigs to process, one task per contig
    array=( $(get_contig_list_from_file $clist) ) || return 1

    for idx in "${!array[@]}"; do
        local optlist=""

        # -out-processdir option, the output directory for the process
        define_opt "-out-processdir" "${process_outdir}" optlist || return 1

        # -r option
        define_opt_from_proc_out "-r" "genodb.create_genref_for_bam" "-outfile" optlist || return 1

        # -dx option
        define_cmdline_infile_opt_if_given "$cmdline" "-dx" optlist || return 1

        # -normalbam option
        define_opt_from_proc_task_out "-normalbam" "genodb.parallel_split_norm_bam" "${idx}" "-outfile" optlist || return 1

        # -tumorbam option
        define_opt_from_proc_task_out "-tumorbam" "genodb.parallel_split_tum_bam" "${idx}" "-outfile" optlist || return 1

        # -contig option
        local contig=${array[$idx]}
        define_opt "-contig" "$contig" optlist || return 1

        # Save option list
        save_opt_list optlist
    done
}

########
genodb.parallel_lumpy_identify_cmdline_opts()
{
    opt_is_cmdline "-lc"
    opt_is_cmdline "-lx"
}

########
genodb.parallel_lumpy_define_opts()
{
    # Initialize variables
    local cmdline=$1
    local process_spec=$2
    local process_name=$3
    local process_outdir=$4

    # Get name of contig list file
    local clist
    clist=$(read_opt_value_from_line "$cmdline" "-lc") || { errmsg "Error: -lc option not found"; return 1; }

    # Array of contigs to process, one task per contig
    array=( $(get_contig_list_from_file $clist) ) || return 1

    for idx in "${!array[@]}"; do
        local optlist=""

        # -out-processdir option, the output directory for the process
        define_opt "-out-processdir" "${process_outdir}" optlist || return 1

        # -lx option
        define_cmdline_infile_opt_if_given "$cmdline" "-lx" optlist || return 1

        # -normalbam option
        local contig=${array[$idx]}
        define_opt_from_proc_task_out "-normalbam" "genodb.parallel_split_norm_bam" "${idx}" "-outfile" optlist || return 1

        # -tumorbam option
        define_opt_from_proc_task_out "-tumorbam" "genodb.parallel_split_tum_bam" "${idx}" "-outfile" optlist || return 1

        # -contig option
        define_opt "-contig" "$contig" optlist || return 1

        # -outfile option
        local outfile="${process_outdir}"/out${contig}.vcf
        define_opt "-outfile" "$outfile" optlist || return 1

        # Save option list
        save_opt_list optlist
    done
}

########
genodb.parallel_svtyper_identify_cmdline_opts()
{
    opt_is_cmdline "-lc"
}

########
genodb.parallel_svtyper_define_opts()
{
    # Initialize variables
    local cmdline=$1
    local process_spec=$2
    local process_name=$3
    local process_outdir=$4

    # Get name of contig list file
    local clist
    clist=$(read_opt_value_from_line "$cmdline" "-lc") || { errmsg "Error: -lc option not found"; return 1; }

    # Array of contigs to process, one task per contig
    array=( $(get_contig_list_from_file $clist) ) || return 1

    for idx in "${!array[@]}"; do
        local optlist=""

        # -normalbam option (the full bam, not the per-contig split: svtyper
        # needs whole-genome read evidence for breakpoint genotyping)
        define_opt_from_proc_out "-normalbam" "genodb.index_norm_bam" "-out-nb" optlist || return 1

        # -tumorbam option (full bam, same reason)
        define_opt_from_proc_out "-tumorbam" "genodb.index_tum_bam" "-out-tb" optlist || return 1

        # -contig option
        local contig=${array[$idx]}
        define_opt "-contig" "$contig" optlist || return 1

        # -vcf option
        define_opt_from_proc_task_out "-vcf" "genodb.parallel_lumpy" "${idx}" "-outfile" optlist || return 1

        # -outfile option
        local outfile="${process_outdir}"/out${contig}.vcf
        define_opt "-outfile" "$outfile" optlist || return 1

        # Save option list
        save_opt_list optlist
    done
}

########
genodb.strelka_germline_identify_cmdline_opts()
{
    opt_is_cmdline "-cr"
}

########
genodb.strelka_germline_define_opts()
{
    # Initialize variables
    local cmdline=$1
    local process_spec=$2
    local process_name=$3
    local process_outdir=$4
    local optlist=""

    # Define the -out-processdir option, the output directory for the process
    define_opt "-out-processdir" "${process_outdir}" optlist || return 1

    # -r option
    define_opt_from_proc_out "-r" "genodb.create_genref_for_bam" "-outfile" optlist || return 1

    # -normalbam option
    define_opt_from_proc_out "-normalbam" "genodb.index_norm_bam" "-out-nb" optlist || return 1

    # -cr option
    define_cmdline_infile_opt_if_given "$cmdline" "-cr" optlist || return 1

    # -out-summarydir option
    define_opt_from_shared_dir "-out-summarydir" "summary/germline_snvs" optlist || return 1

    # -cpus option
    define_procspec_opt "${process_spec}" "-cpus" "cpus" optlist || return 1

    # Save option list
    save_opt_list optlist
}

########
genodb.manta_somatic_identify_cmdline_opts()
{
    opt_is_cmdline "-cr"
}

########
genodb.manta_somatic_define_opts()
{
    # Initialize variables
    local cmdline=$1
    local process_spec=$2
    local process_name=$3
    local process_outdir=$4
    local optlist=""

    # Define the -out-processdir option, the output directory for the process
    define_opt "-out-processdir" "${process_outdir}" optlist || return 1

    # -r option
    define_opt_from_proc_out "-r" "genodb.create_genref_for_bam" "-outfile" optlist || return 1

    # -normalbam option
    define_opt_from_proc_out "-normalbam" "genodb.index_norm_bam" "-out-nb" optlist || return 1

    # -tumorbam option
    define_opt_from_proc_out "-tumorbam" "genodb.index_tum_bam" "-out-tb" optlist || return 1

    # -cr option
    define_cmdline_infile_opt_if_given "$cmdline" "-cr" optlist || return 1

    # -cpus option
    define_procspec_opt "${process_spec}" "-cpus" "cpus" optlist || return 1

    # Save option list
    save_opt_list optlist
}

########
genodb.strelka_somatic_identify_cmdline_opts()
{
    opt_is_cmdline "-cr"
}

########
genodb.strelka_somatic_define_opts()
{
    # Initialize variables
    local cmdline=$1
    local process_spec=$2
    local process_name=$3
    local process_outdir=$4
    local optlist=""

    # Define the -out-processdir option, the output directory for the process
    define_opt "-out-processdir" "${process_outdir}" optlist || return 1

    # -r option
    define_opt_from_proc_out "-r" "genodb.create_genref_for_bam" "-outfile" optlist || return 1

    # -normalbam option
    define_opt_from_proc_out "-normalbam" "genodb.index_norm_bam" "-out-nb" optlist || return 1

    # -tumorbam option
    define_opt_from_proc_out "-tumorbam" "genodb.index_tum_bam" "-out-tb" optlist || return 1

    # -manta-outd option
    define_opt_from_proc_out "-manta-outd" "genodb.manta_somatic" "-out-processdir" optlist || return 1

    # -cr option
    define_cmdline_infile_opt_if_given "$cmdline" "-cr" optlist || return 1

    # -cpus option
    define_procspec_opt "${process_spec}" "-cpus" "cpus" optlist || return 1

    # Save option list
    save_opt_list optlist
}

########
genodb.msisensor_pro_identify_cmdline_opts()
{
    :
}

########
genodb.msisensor_pro_define_opts()
{
    # Initialize variables
    local cmdline=$1
    local process_spec=$2
    local process_name=$3
    local process_outdir=$4
    local optlist=""

    # Define the -out-processdir option, the output directory for the process
    define_opt "-out-processdir" "${process_outdir}" optlist || return 1

    # -r option
    define_opt_from_proc_out "-r" "genodb.create_genref_for_bam" "-outfile" optlist || return 1

    # -normalbam option
    define_opt_from_proc_out "-normalbam" "genodb.index_norm_bam" "-out-nb" optlist || return 1

    # -tumorbam option
    define_opt_from_proc_out "-tumorbam" "genodb.index_tum_bam" "-out-tb" optlist || return 1

    # -cpus option
    define_procspec_opt "${process_spec}" "-cpus" "cpus" optlist || return 1

    # Save option list
    save_opt_list optlist
}

########
genodb.snp_pileup_plus_facets_identify_cmdline_opts()
{
    opt_is_cmdline "-sv"
    opt_is_cmdline "-md"
}

########
genodb.snp_pileup_plus_facets_define_opts()
{
    # Initialize variables
    local cmdline=$1
    local process_spec=$2
    local process_name=$3
    local process_outdir=$4
    local optlist=""

    # Define the -out-processdir option, the output directory for the process
    define_opt "-out-processdir" "${process_outdir}" optlist || return 1

    # -sv option
    define_cmdline_infile_opt "$cmdline" "-sv" optlist || return 1

    # -md option
    define_cmdline_opt_if_given "$cmdline" "-md" optlist || return 1

    # -normalbam option
    define_opt_from_proc_out "-normalbam" "genodb.index_norm_bam" "-out-nb" optlist || return 1

    # -tumorbam option
    define_opt_from_proc_out "-tumorbam" "genodb.index_tum_bam" "-out-tb" optlist || return 1

    # Save option list
    save_opt_list optlist
}

########
genodb.cnvkit_identify_cmdline_opts()
{
    :
}

########
genodb.cnvkit_define_opts()
{
    # Initialize variables
    local cmdline=$1
    local process_spec=$2
    local process_name=$3
    local process_outdir=$4
    local optlist=""

    # Define the -out-processdir option, the output directory for the process
    define_opt "-out-processdir" "${process_outdir}" optlist || return 1

    # -r option
    define_opt_from_proc_out "-r" "genodb.create_genref_for_bam" "-outfile" optlist || return 1

    # -normalbam option
    define_opt_from_proc_out "-normalbam" "genodb.index_norm_bam" "-out-nb" optlist || return 1

    # -tumorbam option
    define_opt_from_proc_out "-tumorbam" "genodb.index_tum_bam" "-out-tb" optlist || return 1

    # -cpus option
    define_procspec_opt "${process_spec}" "-cpus" "cpus" optlist || return 1

    # Save option list
    save_opt_list optlist
}

########
genodb.platypus_germline_identify_cmdline_opts()
{
    :
}

########
genodb.platypus_germline_define_opts()
{
    # Initialize variables
    local cmdline=$1
    local process_spec=$2
    local process_name=$3
    local process_outdir=$4
    local optlist=""

    # Define the -out-processdir option, the output directory for the process
    define_opt "-out-processdir" "${process_outdir}" optlist || return 1

    # -r option
    define_opt_from_proc_out "-r" "genodb.create_genref_for_bam" "-outfile" optlist || return 1

    # -normalbam option
    define_opt_from_proc_out "-normalbam" "genodb.index_norm_bam" "-out-nb" optlist || return 1

    # -out-summarydir option
    define_opt_from_shared_dir "-out-summarydir" "summary/germline_snvs" optlist || return 1

    # -cpus option
    define_procspec_opt "${process_spec}" "-cpus" "cpus" optlist || return 1

    # Save option list
    save_opt_list optlist
}

genodb_ega_decsingle_extended_program()
{
    add_debasher_process "genodb.decsingle_ega_norm_bam"      "cpus=1  mem=2048    time=24:00:00" "processdeps=none"
    add_debasher_process "genodb.decsingle_ega_tum_bam"       "cpus=1  mem=2048    time=24:00:00" "processdeps=none"
    add_debasher_process "genodb.index_norm_bam"              "cpus=1  mem=1024    time=4:00:00"
    add_debasher_process "genodb.index_tum_bam"               "cpus=1  mem=1024    time=4:00:00"
    add_debasher_process "genodb.bedtools_genomecov_norm_bam" "cpus=1  mem=2G      time=5:00:00"
    add_debasher_process "genodb.bedtools_genomecov_tum_bam"  "cpus=1  mem=2G      time=5:00:00"
    add_debasher_process "genodb.create_genref_for_bam"       "cpus=1  mem=8G      time=4:00:00"
    add_debasher_process "genodb.parallel_split_norm_bam"     "cpus=1  mem=2G      time=5:00:00,10:00:00  throttle=16"
    add_debasher_process "genodb.parallel_split_tum_bam"      "cpus=1  mem=2G      time=5:00:00,10:00:00  throttle=16"
    add_debasher_process "genodb.parallel_samtools_mpileup_norm_bam" "cpus=1  mem=4G  time=8:00:00,16:00:00 throttle=24"
    add_debasher_process "genodb.parallel_samtools_mpileup_tum_bam"  "cpus=1  mem=4G  time=8:00:00,16:00:00 throttle=24"
    add_debasher_process "genodb.gen_sequenza_gcc"            "cpus=1  mem=1G      time=01:00:00"
    add_debasher_process "genodb.parallel_bam2seqz"           "cpus=1  mem=2G      time=5:00:00  throttle=16"
    add_debasher_process "genodb.seqzmerge_plus_sequenza"     "cpus=1  mem=10G     time=5:00:00,10:00:00"
    add_debasher_process "genodb.parallel_delly"              "cpus=1  mem=10G,30G time=5:00:00  throttle=16"
    add_debasher_process "genodb.parallel_lumpy"              "cpus=2  mem=10G,30G time=5:00:00,10:00:00  throttle=16"
    add_debasher_process "genodb.parallel_svtyper"            "cpus=1  mem=8G,16G  time=6:00:00,24:00:00  throttle=16"
    add_debasher_process "genodb.strelka_germline"            "cpus=8  mem=6G      time=6:00:00,12:00:00"
    add_debasher_process "genodb.manta_somatic"               "cpus=8  mem=8G      time=8:00:00,16:00:00"
    add_debasher_process "genodb.strelka_somatic"             "cpus=8  mem=6G      time=8:00:00,16:00:00"
    add_debasher_process "genodb.msisensor_pro"               "cpus=2  mem=10G,20G time=8:00:00,24:00:00"
    add_debasher_process "genodb.snp_pileup_plus_facets"      "cpus=1  mem=25G     time=8:00:00,24:00:00"
    add_debasher_process "genodb.cnvkit"                      "cpus=8  mem=30G     time=6:00:00,24:00:00"
    add_debasher_process "genodb.platypus_germline"           "cpus=1  mem=4096    time=8:00:00,16:00:00"
    add_debasher_process "genodb.clear_datadir"               "cpus=1  mem=1024    time=0:10:00"  "processdeps=afterok:genodb.parallel_delly,afterok:genodb.parallel_lumpy,afterok:genodb.parallel_svtyper,afterok:genodb.manta_somatic,afterok:genodb.strelka_germline,afterok:genodb.strelka_somatic,afterok:genodb.msisensor_pro,afterok:genodb.cnvkit,afterok:genodb.snp_pileup_plus_facets,afterok:genodb.platypus_germline,afterok:genodb.seqzmerge_plus_sequenza"
}
