# *- bash -*

load_panpipe_module "genop_bam_analysis"

genop_copy_harmoniz_pipeline()
{
    add_panpipe_process "copy_norm_bam"         "cpus=1  mem=2048 time=48:00:00"
    add_panpipe_process "copy_tum_bam"          "cpus=1  mem=2048 time=48:00:00"
    add_panpipe_process "norm_bam_to_ubam"      "cpus=1  mem=4096 time=48:00:00"
    add_panpipe_process "align_norm_ubam"       "cpus=6  mem=8192 time=48:00:00"
    add_panpipe_process "tum_bam_to_ubam"       "cpus=1  mem=4096 time=48:00:00"
    add_panpipe_process "align_tum_ubam"        "cpus=6  mem=8192 time=48:00:00"
    add_panpipe_process "index_norm_bam"        "cpus=1  mem=1024 time=48:00:00"
    add_panpipe_process "index_tum_bam"         "cpus=1  mem=1024 time=48:00:00"
    add_panpipe_process "strelka_germline"      "cpus=4  mem=6G   time=48:00:00"
    add_panpipe_process "platypus_germline"     "cpus=1  mem=4096 time=48:00:00"
    add_panpipe_process "gatk_haplotypecaller"  "cpus=4  mem=4096 time=48:00:00"
    add_panpipe_process "manta_somatic"         "cpus=2  mem=4096 time=48:00:00"
    add_panpipe_process "strelka_somatic"       "cpus=4  mem=6G   time=48:00:00"
    add_panpipe_process "mutect2_somatic"       "cpus=4  mem=8G   time=48:00:00"
    add_panpipe_process "lofreq_somatic"        "cpus=4  mem=4096 time=48:00:00"
    add_panpipe_process "msisensor_pro"         "cpus=2  mem=4096 time=48:00:00"
    add_panpipe_process "concat_germline_snvs"  "cpus=1  mem=2048 time=48:00:00"
}

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
    define_opt "-normalbam" "${abs_datadir}"/copy_normal.bam optlist || return 1

    # -mrec option
    define_cmdline_nonmandatory_opt "$cmdline" "-mrec" ${DEFAULT_MAX_RECORDS_IN_RAM_GATK} optlist || return 1

    # -outfile option
    define_opt "-outfile" "${abs_datadir}"/normal_unmapped.bam optlist || return 1

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
    define_opt "-tumorbam" "${abs_datadir}"/copy_tumor.bam optlist || return 1

    # -mrec option
    define_cmdline_nonmandatory_opt "$cmdline" "-mrec" ${DEFAULT_MAX_RECORDS_IN_RAM_GATK} optlist || return 1

    # -outfile option
    define_opt "-outfile" "${abs_datadir}"/tumor_unmapped.bam optlist || return 1

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
    define_opt "-normalbam" "${abs_datadir}"/normal_unmapped.bam optlist || return 1

    # -mrec option
    define_cmdline_nonmandatory_opt "$cmdline" "-mrec" ${DEFAULT_MAX_RECORDS_IN_RAM_GATK} optlist || return 1

    # -cpus option
    local cpus
    cpus=`extract_cpus_from_process_spec "$process_spec"` || return 1
    define_opt "-cpus" $cpus optlist

    # -outfile option
    define_opt "-outfile" "${abs_datadir}"/normal.bam optlist || return 1

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
    define_opt "-tumorbam" "${abs_datadir}"/tumor_unmapped.bam optlist || return 1

    # -mrec option
    define_cmdline_nonmandatory_opt "$cmdline" "-mrec" ${DEFAULT_MAX_RECORDS_IN_RAM_GATK} optlist || return 1

    # -cpus option
    local cpus
    cpus=`extract_cpus_from_process_spec "$process_spec"` || return 1
    define_opt "-cpus" $cpus optlist

    # -outfile option
    define_opt "-outfile" "${abs_datadir}"/tumor.bam optlist || return 1

    # Save option list
    save_opt_list optlist
}
