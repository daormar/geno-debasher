# *- bash -*

load_debasher_module "genodb_bam_analysis"

copy_norm_bam_define_opts()
{
    # Initialize variables
    local cmdline=$1
    local process_spec=$2
    local optlist=""

    # -extn option
    define_cmdline_opt "$cmdline" "-extn" optlist || return 1

    # -out-nb option
    local abs_datadir=`get_absolute_shdirname "${DATADIR_BASENAME}"`
    local normalbam="${abs_datadir}"/copy_normal.bam
    define_opt "-out-nb" "$normalbam" optlist || return 1

    # Save option list
    save_opt_list optlist
}

copy_tum_bam_define_opts()
{
    # Initialize variables
    local cmdline=$1
    local process_spec=$2
    local optlist=""


    # -extt option
    define_cmdline_opt "$cmdline" "-extt" optlist || return 1

    # -tumorbam option
    local abs_datadir=`get_absolute_shdirname "${DATADIR_BASENAME}"`
    local tumorbam="${abs_datadir}"/copy_tumor.bam
    define_opt "-out-tb" "$tumorbam" optlist || return 1

    # Save option list
    save_opt_list optlist
}

norm_bam_to_ubam_define_opts()
{
    # Initialize variables
    local cmdline=$1
    local process_spec=$2
    local process_name=$3
    local process_outdir=$4
    local optlist=""

    # Define the -out-processdir option, the output directory for the process
    define_opt "-out-processdir" "${process_outdir}" optlist || return 1

    # Get data directory
    local abs_datadir=`get_absolute_shdirname "${DATADIR_BASENAME}"`

    # -normalbam option
    define_opt_from_proc_out "-normalbam" "copy_norm_bam" "-out-nb" optlist || return 1

    # -mrec option
    define_cmdline_opt_if_given "$cmdline" "-mrec" optlist || return 1

    # -outfile option
    local outfile="${abs_datadir}"/normal_unmapped.bam
    define_opt "-outfile" "$outfile" optlist || return 1

    # Save option list
    save_opt_list optlist
}

tum_bam_to_ubam_define_opts()
{
    # Initialize variables
    local cmdline=$1
    local process_spec=$2
    local process_name=$3
    local process_outdir=$4
    local optlist=""

    # Define the -out-processdir option, the output directory for the process
    define_opt "-out-processdir" "${process_outdir}" optlist || return 1

    # Get data directory
    local abs_datadir=`get_absolute_shdirname "${DATADIR_BASENAME}"`

    # -tumorbam option
    define_opt_from_proc_out "-tumorbam" "copy_tum_bam" "-out-tb" optlist || return 1

    # -mrec option
    define_cmdline_opt_if_given "$cmdline" "-mrec" optlist || return 1

    # -outfile option
    local outfile="${abs_datadir}"/tumor_unmapped.bam
    define_opt "-outfile" "$outfile" optlist || return 1

    # Save option list
    save_opt_list optlist
}

align_norm_ubam_define_opts()
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
    local genref
    genref=`get_ref_filename "$cmdline"` || return 1
    define_opt "-r" "$genref" optlist || return 1

    # Get data directory
    local abs_datadir=`get_absolute_shdirname "${DATADIR_BASENAME}"`

    # -normalbam option
    define_opt_from_proc_out "-normalbam" "norm_bam_to_ubam" "-outfile" optlist || return 1

    # -mrec option
    define_cmdline_opt_if_given "$cmdline" "-mrec" optlist || return 1

    # -cpus option
    local cpus
    cpus=`extract_cpus_from_process_spec "$process_spec"` || return 1
    define_opt "-cpus" $cpus optlist

    # -outfile option
    local outfile="${abs_datadir}"/normal.bam
    define_opt "-outfile" "$outfile" optlist || return 1

    # Save option list
    save_opt_list optlist
}

align_tum_ubam_define_opts()
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
    local genref
    genref=`get_ref_filename "$cmdline"` || return 1
    define_opt "-r" "$genref" optlist || return 1

    # Get data directory
    local abs_datadir=`get_absolute_shdirname "${DATADIR_BASENAME}"`

    # -tumorbam option
    define_opt_from_proc_out "-tumorbam" "tum_bam_to_ubam" "-outfile" optlist || return 1

    # -mrec option
    define_cmdline_opt_if_given "$cmdline" "-mrec" optlist || return 1

    # -cpus option
    local cpus
    cpus=`extract_cpus_from_process_spec "$process_spec"` || return 1
    define_opt "-cpus" $cpus optlist

    # -outfile option
    local outfile="${abs_datadir}"/tumor.bam
    define_opt "-outfile" "$outfile" optlist || return 1

    # Save option list
    save_opt_list optlist
}

genodb_copy_harmoniz_program()
{
    add_debasher_process "copy_norm_bam"         "cpus=1  mem=2048 time=48:00:00"
    add_debasher_process "copy_tum_bam"          "cpus=1  mem=2048 time=48:00:00"
    add_debasher_process "norm_bam_to_ubam"      "cpus=1  mem=4096 time=48:00:00"
    add_debasher_process "align_norm_ubam"       "cpus=6  mem=8192 time=48:00:00"
    add_debasher_process "tum_bam_to_ubam"       "cpus=1  mem=4096 time=48:00:00"
    add_debasher_process "align_tum_ubam"        "cpus=6  mem=8192 time=48:00:00"
    add_debasher_process "index_norm_bam"        "cpus=1  mem=1024 time=48:00:00"
    add_debasher_process "index_tum_bam"         "cpus=1  mem=1024 time=48:00:00"
    add_debasher_process "strelka_germline"      "cpus=4  mem=6G   time=48:00:00"
    add_debasher_process "platypus_germline"     "cpus=1  mem=4096 time=48:00:00"
    add_debasher_process "gatk_haplotypecaller"  "cpus=4  mem=4096 time=48:00:00"
    add_debasher_process "manta_somatic"         "cpus=2  mem=4096 time=48:00:00"
    add_debasher_process "strelka_somatic"       "cpus=4  mem=6G   time=48:00:00"
    add_debasher_process "mutect2_somatic"       "cpus=4  mem=8G   time=48:00:00"
    add_debasher_process "lofreq_somatic"        "cpus=4  mem=4096 time=48:00:00"
    add_debasher_process "msisensor_pro"         "cpus=2  mem=4096 time=48:00:00"
    add_debasher_process "concat_germline_snvs"  "cpus=1  mem=2048 time=48:00:00"
}
