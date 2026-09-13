# *- bash -*

load_debasher_module "genodb_bam_analysis"

########
genodb.strelka_germline_identify_cmdline_opts()
{
    opt_is_cmdline "-normalbam"
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
genodb.manta_somatic_identify_cmdline_opts()
{
    opt_is_cmdline "-normalbam"
    opt_is_cmdline "-tumorbam"
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
genodb.strelka_somatic_identify_cmdline_opts()
{
    opt_is_cmdline "-normalbam"
    opt_is_cmdline "-tumorbam"
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
    define_cmdline_infile_opt "$cmdline" "-normalbam" optlist || return 1

    # -tumorbam option
    define_cmdline_infile_opt "$cmdline" "-tumorbam" optlist || return 1

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
    opt_is_cmdline "-normalbam"
    opt_is_cmdline "-tumorbam"
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
    define_cmdline_infile_opt "$cmdline" "-normalbam" optlist || return 1

    # -tumorbam option
    define_cmdline_infile_opt "$cmdline" "-tumorbam" optlist || return 1

    # -cpus option
    define_procspec_opt "${process_spec}" "-cpus" "cpus" optlist || return 1

    # Save option list
    save_opt_list optlist
}

########
genodb.platypus_germline_identify_cmdline_opts()
{
    opt_is_cmdline "-normalbam"
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
    define_cmdline_infile_opt "$cmdline" "-normalbam" optlist || return 1

    # -out-summarydir option
    define_opt_from_shared_dir "-out-summarydir" "summary/germline_snvs" optlist || return 1

    # -cpus option
    define_procspec_opt "${process_spec}" "-cpus" "cpus" optlist || return 1

    # Save option list
    save_opt_list optlist
}

genodb_nodownload_xgs_program()
{
    add_debasher_process "genodb.create_genref_for_bam"       "cpus=1  mem=8G      time=4:00:00"  "processdeps=none"
    add_debasher_process "genodb.bedtools_genomecov_norm_bam" "cpus=1  mem=1024    time=4:00:00"  "processdeps=none"
    add_debasher_process "genodb.bedtools_genomecov_tum_bam"  "cpus=1  mem=1024    time=4:00:00"  "processdeps=none"
    add_debasher_process "genodb.strelka_germline"            "cpus=8  mem=6G      time=6:00:00,12:00:00"
    add_debasher_process "genodb.manta_somatic"               "cpus=8  mem=8G      time=8:00:00,16:00:00"
    add_debasher_process "genodb.strelka_somatic"             "cpus=8  mem=6G      time=8:00:00,16:00:00"
    add_debasher_process "genodb.msisensor_pro"               "cpus=2  mem=8G      time=8:00:00,24:00:00"
    add_debasher_process "genodb.snp_pileup_plus_facets"      "cpus=1  mem=8G      time=8:00:00,24:00:00"  "processdeps=afterok:genodb.create_genref_for_bam"
    add_debasher_process "genodb.platypus_germline"           "cpus=1  mem=4096    time=8:00:00,16:00:00"
}
