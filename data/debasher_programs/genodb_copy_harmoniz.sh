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

########
norm_bam_to_ubam_identify_cmdline_opts()
{
    opt_is_cmdline "-mrec"
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

########
tum_bam_to_ubam_identify_cmdline_opts()
{
    opt_is_cmdline "-mrec"
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

########
align_norm_ubam_identify_cmdline_opts()
{
    opt_is_cmdline "-r"
    opt_is_cmdline "-mrec"
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
    define_cmdline_infile_opt "$cmdline" "-r" optlist || return 1

    # Get data directory
    local abs_datadir=`get_absolute_shdirname "${DATADIR_BASENAME}"`

    # -normalbam option
    define_opt_from_proc_out "-normalbam" "norm_bam_to_ubam" "-outfile" optlist || return 1

    # -mrec option
    define_cmdline_opt_if_given "$cmdline" "-mrec" optlist || return 1

    # -cpus option
    define_procspec_opt "${process_spec}" "-cpus" "cpus" optlist || return 1

    # -outfile option
    local outfile="${abs_datadir}"/normal.bam
    define_opt "-outfile" "$outfile" optlist || return 1

    # Save option list
    save_opt_list optlist
}

########
align_tum_ubam_identify_cmdline_opts()
{
    opt_is_cmdline "-r"
    opt_is_cmdline "-mrec"
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
    define_cmdline_infile_opt "$cmdline" "-r" optlist || return 1

    # Get data directory
    local abs_datadir=`get_absolute_shdirname "${DATADIR_BASENAME}"`

    # -tumorbam option
    define_opt_from_proc_out "-tumorbam" "tum_bam_to_ubam" "-outfile" optlist || return 1

    # -mrec option
    define_cmdline_opt_if_given "$cmdline" "-mrec" optlist || return 1

    # -cpus option
    define_procspec_opt "${process_spec}" "-cpus" "cpus" optlist || return 1

    # -outfile option
    local outfile="${abs_datadir}"/tumor.bam
    define_opt "-outfile" "$outfile" optlist || return 1

    # Save option list
    save_opt_list optlist
}

########
index_norm_bam_define_opts()
{
    # Initialize variables
    local cmdline=$1
    local process_spec=$2
    local optlist=""

    # -normalbam option
    define_opt_from_proc_out "-normalbam" "align_norm_ubam" "-outfile" optlist || return 1

    # -out-nbidx option
    local abs_datadir=`get_absolute_shdirname "${DATADIR_BASENAME}"`
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
    define_opt_from_proc_out "-tumorbam" "align_tum_ubam" "-outfile" optlist || return 1

    # -out-tbidx option
    local abs_datadir=`get_absolute_shdirname "${DATADIR_BASENAME}"`
    define_opt "-out-tbidx" "${abs_datadir}/tumor.bam.bai" optlist || return 1

    # -out-tb option (republished once indexed, so downstream processes
    # can connect to it and depend on indexing having completed)
    define_opt "-out-tb" "${abs_datadir}/tumor.bam" optlist || return 1

    # Save option list
    save_opt_list optlist
}

########
strelka_germline_identify_cmdline_opts()
{
    opt_is_cmdline "-r"
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
    define_cmdline_infile_opt "$cmdline" "-r" optlist || return 1

    # -normalbam option
    define_opt_from_proc_out "-normalbam" "index_norm_bam" "-out-nb" optlist || return 1

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
platypus_germline_identify_cmdline_opts()
{
    opt_is_cmdline "-r"
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
    define_cmdline_infile_opt "$cmdline" "-r" optlist || return 1

    # -normalbam option
    define_opt_from_proc_out "-normalbam" "index_norm_bam" "-out-nb" optlist || return 1

    # -out-summarydir option
    define_opt_from_shared_dir "-out-summarydir" "summary/germline_snvs" optlist || return 1

    # -cpus option
    define_procspec_opt "${process_spec}" "-cpus" "cpus" optlist || return 1

    # Save option list
    save_opt_list optlist
}

########
gatk_haplotypecaller_identify_cmdline_opts()
{
    opt_is_cmdline "-r"
    opt_is_cmdline "-sample-name"
}

########
gatk_haplotypecaller_define_opts()
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
    define_opt_from_proc_out "-normalbam" "index_norm_bam" "-out-nb" optlist || return 1

    # -sample-name option
    define_cmdline_opt "$cmdline" "-sample-name" optlist || return 1

    # -cpus option
    define_procspec_opt "${process_spec}" "-cpus" "cpus" optlist || return 1

    # -mem option
    define_procspec_opt "${process_spec}" "-mem" "mem" optlist || return 1

    # Save option list
    save_opt_list optlist
}

########
manta_somatic_identify_cmdline_opts()
{
    opt_is_cmdline "-r"
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
    define_cmdline_infile_opt "$cmdline" "-r" optlist || return 1

    # -normalbam option
    define_opt_from_proc_out "-normalbam" "index_norm_bam" "-out-nb" optlist || return 1

    # -tumorbam option
    define_opt_from_proc_out "-tumorbam" "index_tum_bam" "-out-tb" optlist || return 1

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
    opt_is_cmdline "-r"
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
    define_cmdline_infile_opt "$cmdline" "-r" optlist || return 1

    # -normalbam option
    define_opt_from_proc_out "-normalbam" "index_norm_bam" "-out-nb" optlist || return 1

    # -tumorbam option
    define_opt_from_proc_out "-tumorbam" "index_tum_bam" "-out-tb" optlist || return 1

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
mutect2_somatic_identify_cmdline_opts()
{
    opt_is_cmdline "-r"
    opt_is_cmdline "-norm-sample-name"
    opt_is_cmdline "-panel-of-normals"
}

########
mutect2_somatic_define_opts()
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
    define_opt_from_proc_out "-normalbam" "index_norm_bam" "-out-nb" optlist || return 1

    # -tumorbam option
    define_opt_from_proc_out "-tumorbam" "index_tum_bam" "-out-tb" optlist || return 1

    # -norm-sample-name option
    define_cmdline_opt "$cmdline" "-norm-sample-name" optlist || return 1

    # -panel-of-normals option
    define_cmdline_opt "$cmdline" "-panel-of-normals" optlist || return 1

    # -cpus option
    define_procspec_opt "${process_spec}" "-cpus" "cpus" optlist || return 1

    # -mem option
    define_procspec_opt "${process_spec}" "-mem" "mem" optlist || return 1

    # Save option list
    save_opt_list optlist
}

########
lofreq_somatic_identify_cmdline_opts()
{
    opt_is_cmdline "-r"
}

########
lofreq_somatic_define_opts()
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
    define_opt_from_proc_out "-normalbam" "index_norm_bam" "-out-nb" optlist || return 1

    # -tumorbam option
    define_opt_from_proc_out "-tumorbam" "index_tum_bam" "-out-tb" optlist || return 1

    # -cpus option
    define_procspec_opt "${process_spec}" "-cpus" "cpus" optlist || return 1

    # Save option list
    save_opt_list optlist
}

########
msisensor_pro_identify_cmdline_opts()
{
    opt_is_cmdline "-r"
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
    define_cmdline_infile_opt "$cmdline" "-r" optlist || return 1

    # -normalbam option
    define_opt_from_proc_out "-normalbam" "index_norm_bam" "-out-nb" optlist || return 1

    # -tumorbam option
    define_opt_from_proc_out "-tumorbam" "index_tum_bam" "-out-tb" optlist || return 1

    # -cpus option
    define_procspec_opt "${process_spec}" "-cpus" "cpus" optlist || return 1

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
