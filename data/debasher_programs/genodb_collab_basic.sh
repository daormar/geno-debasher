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
    define_opt_from_proc_out "-normalbam" "genodb.download_collab_norm_bam" "-out-nb" optlist || return 1

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
    define_opt_from_proc_out "-tumorbam" "genodb.download_collab_tum_bam" "-out-tb" optlist || return 1

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
genodb.manta_somatic_identify_cmdline_opts()
{
    opt_is_cmdline "-r"
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
    define_cmdline_infile_opt "$cmdline" "-r" optlist || return 1

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
    opt_is_cmdline "-r"
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
    define_cmdline_infile_opt "$cmdline" "-r" optlist || return 1

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
    opt_is_cmdline "-r"
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
    define_cmdline_infile_opt "$cmdline" "-r" optlist || return 1

    # -normalbam option
    define_opt_from_proc_out "-normalbam" "genodb.index_norm_bam" "-out-nb" optlist || return 1

    # -tumorbam option
    define_opt_from_proc_out "-tumorbam" "genodb.index_tum_bam" "-out-tb" optlist || return 1

    # -cpus option
    define_procspec_opt "${process_spec}" "-cpus" "cpus" optlist || return 1

    # Save option list
    save_opt_list optlist
}

genodb_collab_basic_program()
{
    add_debasher_process "genodb.download_collab_norm_bam" "cpus=1  mem=2048 time=48:00:00" "processdeps=none"
    add_debasher_process "genodb.download_collab_tum_bam"  "cpus=1  mem=2048 time=48:00:00" "processdeps=none"
    add_debasher_process "genodb.index_norm_bam"           "cpus=1  mem=1024 time=48:00:00"
    add_debasher_process "genodb.index_tum_bam"            "cpus=1  mem=1024 time=48:00:00"
    add_debasher_process "genodb.manta_somatic"            "cpus=8  mem=8G   time=48:00:00"
    add_debasher_process "genodb.strelka_somatic"          "cpus=8  mem=4096 time=48:00:00"
    add_debasher_process "genodb.msisensor_pro"            "cpus=2  mem=10G  time=48:00:00"
    add_debasher_process "genodb.clear_datadir"            "cpus=1  mem=1024 time=48:00:00" "processdeps=afterok:genodb.manta_somatic,afterok:genodb.strelka_somatic"
}
