# *- bash -*

load_debasher_module "genodb_bam_analysis"

########
parallel_samtools_mpileup_norm_bam_identify_cmdline_opts()
{
    opt_is_cmdline "-mpb"
    opt_is_cmdline "-lc"
}

########
parallel_samtools_mpileup_norm_bam_define_opts()
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

        # -out-processdir option, the output directory for the process
        define_opt "-out-processdir" "${process_outdir}" optlist || return 1

        # -r option
        define_opt_from_proc_out "-r" "create_genref_for_bam" "-outfile" optlist || return 1

        # -mpb option
        define_cmdline_opt_if_given "$cmdline" "-mpb" optlist

        # -normalbam option
        local contig=${array[$idx]}
        local normalbam="${abs_splitdir}"/normal_${contig}.bam
        define_opt "-normalbam" "$normalbam" optlist || return 1

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
parallel_samtools_mpileup_tum_bam_identify_cmdline_opts()
{
    opt_is_cmdline "-mpb"
    opt_is_cmdline "-lc"
}

########
parallel_samtools_mpileup_tum_bam_define_opts()
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

        # -out-processdir option, the output directory for the process
        define_opt "-out-processdir" "${process_outdir}" optlist || return 1

        # -r option
        define_opt_from_proc_out "-r" "create_genref_for_bam" "-outfile" optlist || return 1

        # -mpb option
        define_cmdline_opt_if_given "$cmdline" "-mpb" optlist

        # -tumorbam option
        local contig=${array[$idx]}
        local tumorbam="${abs_splitdir}"/tumor_${contig}.bam
        define_opt "-tumorbam" "$tumorbam" optlist || return 1

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
gen_sequenza_gcc_identify_cmdline_opts()
{
    :
}

########
gen_sequenza_gcc_define_opts()
{
    # Initialize variables
    local cmdline=$1
    local process_spec=$2
    local process_name=$3
    local process_outdir=$4
    local optlist=""

    # -r option
    define_opt_from_proc_out "-r" "create_genref_for_bam" "-outfile" optlist || return 1

    # Get data directory
    local abs_datadir=$(get_absolute_shdirname "data")

    # -outfile option
    define_opt "-outfile" "${abs_datadir}"/sequenza_gccfile.txt.gz optlist || return 1

    # Save option list
    save_opt_list optlist
}

########
parallel_delly_identify_cmdline_opts()
{
    opt_is_cmdline "-dx"
    opt_is_cmdline "-lc"
}

########
parallel_delly_define_opts()
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

        # -out-processdir option, the output directory for the process
        define_opt "-out-processdir" "${process_outdir}" optlist || return 1

        # -r option
        define_opt_from_proc_out "-r" "create_genref_for_bam" "-outfile" optlist || return 1

        # -dx option
        define_cmdline_infile_opt_if_given "$cmdline" "-dx" optlist || return 1

        # -normalbam option
        local contig=${array[$idx]}
        local normalbam="${abs_splitdir}"/normal_${contig}.bam
        define_opt "-normalbam" "$normalbam" optlist || return 1

        # -tumorbam option
        local tumorbam="${abs_splitdir}"/tumor_${contig}.bam
        define_opt "-tumorbam" "$tumorbam" optlist || return 1

        # -contig option
        define_opt "-contig" "$contig" optlist || return 1

        # Save option list
        save_opt_list optlist
    done
}

########
strelka_germline_identify_cmdline_opts()
{
    opt_is_cmdline "-normalbam"
    opt_is_cmdline "-cr"
}

########
strelka_germline_define_opts()
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
    define_opt_from_proc_out "-r" "create_genref_for_bam" "-outfile" optlist || return 1

    # -normalbam option
    define_cmdline_infile_opt "$cmdline" "-normalbam" optlist || return 1

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
manta_somatic_identify_cmdline_opts()
{
    opt_is_cmdline "-normalbam"
    opt_is_cmdline "-tumorbam"
    opt_is_cmdline "-cr"
}

########
manta_somatic_define_opts()
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
    define_opt_from_proc_out "-r" "create_genref_for_bam" "-outfile" optlist || return 1

    # -normalbam option
    define_cmdline_infile_opt "$cmdline" "-normalbam" optlist || return 1

    # -tumorbam option
    define_cmdline_infile_opt "$cmdline" "-tumorbam" optlist || return 1

    # -cr option
    define_cmdline_infile_opt_if_given "$cmdline" "-cr" optlist || return 1

    # -cpus option
    define_procspec_opt "${process_spec}" "-cpus" "cpus" optlist || return 1

    # Save option list
    save_opt_list optlist
}

########
strelka_somatic_identify_cmdline_opts()
{
    opt_is_cmdline "-normalbam"
    opt_is_cmdline "-tumorbam"
    opt_is_cmdline "-cr"
}

########
strelka_somatic_define_opts()
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
    define_opt_from_proc_out "-r" "create_genref_for_bam" "-outfile" optlist || return 1

    # -normalbam option
    define_cmdline_infile_opt "$cmdline" "-normalbam" optlist || return 1

    # -tumorbam option
    define_cmdline_infile_opt "$cmdline" "-tumorbam" optlist || return 1

    # -manta-outd option
    define_opt_from_proc_out "-manta-outd" "manta_somatic" "-out-processdir" optlist || return 1

    # -cr option
    define_cmdline_infile_opt_if_given "$cmdline" "-cr" optlist || return 1

    # -cpus option
    define_procspec_opt "${process_spec}" "-cpus" "cpus" optlist || return 1

    # Save option list
    save_opt_list optlist
}

########
msisensor_pro_identify_cmdline_opts()
{
    opt_is_cmdline "-normalbam"
    opt_is_cmdline "-tumorbam"
}

########
msisensor_pro_define_opts()
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
    define_opt_from_proc_out "-r" "create_genref_for_bam" "-outfile" optlist || return 1

    # -normalbam option
    define_cmdline_infile_opt "$cmdline" "-normalbam" optlist || return 1

    # -tumorbam option
    define_cmdline_infile_opt "$cmdline" "-tumorbam" optlist || return 1

    # -cpus option
    define_procspec_opt "${process_spec}" "-cpus" "cpus" optlist || return 1

    # Save option list
    save_opt_list optlist
}

########
cnvkit_identify_cmdline_opts()
{
    opt_is_cmdline "-normalbam"
    opt_is_cmdline "-tumorbam"
}

########
cnvkit_define_opts()
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
    define_opt_from_proc_out "-r" "create_genref_for_bam" "-outfile" optlist || return 1

    # -normalbam option
    define_cmdline_infile_opt "$cmdline" "-normalbam" optlist || return 1

    # -tumorbam option
    define_cmdline_infile_opt "$cmdline" "-tumorbam" optlist || return 1

    # -cpus option
    define_procspec_opt "${process_spec}" "-cpus" "cpus" optlist || return 1

    # Save option list
    save_opt_list optlist
}

########
platypus_germline_identify_cmdline_opts()
{
    opt_is_cmdline "-normalbam"
}

########
platypus_germline_define_opts()
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
    define_opt_from_proc_out "-r" "create_genref_for_bam" "-outfile" optlist || return 1

    # -normalbam option
    define_cmdline_infile_opt "$cmdline" "-normalbam" optlist || return 1

    # -out-summarydir option
    define_opt_from_shared_dir "-out-summarydir" "summary/germline_snvs" optlist || return 1

    # -cpus option
    define_procspec_opt "${process_spec}" "-cpus" "cpus" optlist || return 1

    # Save option list
    save_opt_list optlist
}

genodb_nodownload_extended_program()
{
    add_debasher_process "create_genref_for_bam"     "cpus=1  mem=8G      time=4:00:00"  "processdeps=none"
    add_debasher_process "parallel_split_norm_bam"   "cpus=1  mem=2G      time=5:00:00,10:00:00  throttle=16" "processdeps=none"
    add_debasher_process "parallel_split_tum_bam"    "cpus=1  mem=2G      time=5:00:00,10:00:00  throttle=16" "processdeps=none"
    add_debasher_process "parallel_samtools_mpileup_norm_bam" "cpus=1  mem=4G  time=8:00:00,16:00:00 throttle=24" "processdeps=aftercorr:parallel_split_norm_bam"
    add_debasher_process "parallel_samtools_mpileup_tum_bam"  "cpus=1  mem=4G  time=8:00:00,16:00:00 throttle=24" "processdeps=aftercorr:parallel_split_tum_bam"
    add_debasher_process "gen_sequenza_gcc"          "cpus=1  mem=1G      time=1:00:00"
    add_debasher_process "parallel_bam2seqz"         "cpus=1  mem=2G      time=5:00:00  throttle=16" "processdeps=afterok:gen_sequenza_gcc,aftercorr:parallel_samtools_mpileup_norm_bam,aftercorr:parallel_samtools_mpileup_tum_bam"
    add_debasher_process "seqzmerge_plus_sequenza"   "cpus=1  mem=10G     time=5:00:00,10:00:00"  "processdeps=afterok:parallel_bam2seqz"
    add_debasher_process "parallel_delly"            "cpus=1  mem=10G,30G time=5:00:00  throttle=16" "processdeps=aftercorr:parallel_split_norm_bam,aftercorr:parallel_split_tum_bam"
    add_debasher_process "parallel_lumpy"            "cpus=2  mem=10G,30G time=5:00:00,10:00:00  throttle=16" "processdeps=aftercorr:parallel_split_norm_bam,aftercorr:parallel_split_tum_bam"
    add_debasher_process "parallel_svtyper"          "cpus=1  mem=8G,16G  time=6:00:00,24:00:00  throttle=16" "processdeps=aftercorr:parallel_lumpy"
    add_debasher_process "strelka_germline"          "cpus=8  mem=6G      time=6:00:00,12:00:00"
    add_debasher_process "manta_somatic"             "cpus=8  mem=8G      time=8:00:00,16:00:00"
    add_debasher_process "strelka_somatic"           "cpus=8  mem=6G      time=8:00:00,16:00:00"
    add_debasher_process "msisensor_pro"             "cpus=2  mem=10G,20G time=8:00:00,24:00:00"
    add_debasher_process "snp_pileup_plus_facets"    "cpus=1  mem=25G     time=8:00:00,24:00:00"
    add_debasher_process "cnvkit"                    "cpus=8  mem=30G     time=6:00:00,24:00:00"
    add_debasher_process "platypus_germline"         "cpus=1  mem=4096    time=8:00:00,16:00:00"
}
