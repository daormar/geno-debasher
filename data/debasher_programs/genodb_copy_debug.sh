# *- bash -*

load_debasher_module "genodb_bam_analysis"

########
index_norm_bam_define_opts()
{
    # Initialize variables
    local cmdline=$1
    local process_spec=$2
    local optlist=""

    # -normalbam option
    define_opt_from_proc_out "-normalbam" "copy_norm_bam" "-out-nb" optlist || return 1

    # -out-nbidx option
    local abs_datadir=`get_absolute_shdirname "data"`
    define_opt "-out-nbidx" "${abs_datadir}/normal.bam.bai" optlist || return 1

    # -out-nb option (republished once indexed, so downstream processes
    # can connect to it and depend on indexing having completed)
    define_opt "-out-nb" "${abs_datadir}/normal.bam" optlist || return 1

    # Save option list
    save_opt_list optlist
}

########
index_tum_bam_define_opts()
{
    # Initialize variables
    local cmdline=$1
    local process_spec=$2
    local optlist=""

    # -tumorbam option
    define_opt_from_proc_out "-tumorbam" "copy_tum_bam" "-out-tb" optlist || return 1

    # -out-tbidx option
    local abs_datadir=`get_absolute_shdirname "data"`
    define_opt "-out-tbidx" "${abs_datadir}/tumor.bam.bai" optlist || return 1

    # -out-tb option (republished once indexed, so downstream processes
    # can connect to it and depend on indexing having completed)
    define_opt "-out-tb" "${abs_datadir}/tumor.bam" optlist || return 1

    # Save option list
    save_opt_list optlist
}

########
create_genref_for_bam_define_opts()
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
    define_opt_from_proc_out "-bam" "index_norm_bam" "-out-nb" optlist || return 1

    # -cm option
    define_cmdline_infile_opt_if_given "$cmdline" "-cm" optlist || return 1

    # -fbr option
    define_cmdline_infile_opt_if_given "$cmdline" "-fbr" optlist || return 1

    # Get data directory
    local abs_datadir=`get_absolute_shdirname "data"`

    # -outfile option
    local outfile="${abs_datadir}"/genref.fa
    define_opt "-outfile" "$outfile" optlist || return 1

    # Save option list
    save_opt_list optlist
}

genodb_copy_debug_program()
{
    add_debasher_process "copy_norm_bam"         "cpus=1  mem=2048    time=24:00:00"
    add_debasher_process "copy_tum_bam"          "cpus=1  mem=2048    time=24:00:00"
    add_debasher_process "index_norm_bam"        "cpus=1  mem=1024    time=4:00:00"
    add_debasher_process "index_tum_bam"         "cpus=1  mem=1024    time=4:00:00"
    add_debasher_process "create_genref_for_bam" "cpus=1  mem=8G      time=4:00:00"
}
