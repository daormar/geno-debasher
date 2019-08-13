# *- bash -*

##################
# CLEANING STEPS #
##################

########
delete_bam_files_explain_cmdline_opts()
{
    :
}

########
delete_bam_files_define_opts()
{
    # Initialize variables
    local cmdline=$1
    local stepspec=$2
    local optlist=""

    # Define the -step-outd option, the output directory for the step
    local step_outd=`get_step_outdir_given_stepspec "$stepspec"`
    define_opt "-step-outd" ${step_outd} optlist || exit 1

    # -datadir option
    abs_datadir=`get_absolute_shdirname ${DATADIR_BASENAME}`
    define_opt "-datadir" ${abs_datadir} optlist || exit 1

    # Save option list
    save_opt_list optlist
}

########
delete_bam_files()
{
    # Initialize variables
    local step_outd=`read_opt_value_from_line "$*" "-step-outd"`
    local abs_datadir=`read_opt_value_from_line "$*" "-datadir"`

    # Delete bam files
    rm -f ${abs_datadir}/*.bam || exit 1
}

########
clear_datadir_explain_cmdline_opts()
{
    :
}

########
clear_datadir_define_opts()
{
    # Initialize variables
    local cmdline=$1
    local stepspec=$2
    local optlist=""

    # Define the -step-outd option, the output directory for the step
    local step_outd=`get_step_outdir_given_stepspec "$stepspec"`
    define_opt "-step-outd" ${step_outd} optlist || exit 1

    # -datadir option
    abs_datadir=`get_absolute_shdirname ${DATADIR_BASENAME}`
    define_opt "-datadir" ${abs_datadir} optlist || exit 1

    # Save option list
    save_opt_list optlist
}

########
clear_datadir()
{
    # Initialize variables
    local step_outd=`read_opt_value_from_line "$*" "-step-outd"`
    local abs_datadir=`read_opt_value_from_line "$*" "-datadir"`

    # Delete bam files
    rm -rf ${abs_datadir}/* || exit 1

    # Print README.txt file
    echo "NOTE: This directory was cleared by means of the 'clear_datadir' step" > ${abs_datadir}/README.txt || exit 1
}
