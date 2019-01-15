# *- bash -*

#############
# CONSTANTS #
#############

###################################
# PIPELINE SOFTWARE TESTING STEPS #
###################################

########

########
step_a_explain_cmdline_opts()
{
    :
}

########
step_a_define_opts()
{
    # Initialize variables
    local cmdline=$1
    local jobspec=$2
    optlist=""
    
    # Define the -step-outd option, the output directory for the step,
    # which will have the same name of the step
    define_default_step_outd_opt "$cmdline" "$jobspec" optlist || exit 1

    # Save option list
    save_opt_list optlist
}

########
step_a()
{
    display_begin_step_message

    # Initialize variables
    local step_outd=`read_opt_value_from_line "$*" "-step-outd"`

    # sleep some time
    sleep 10

    # Signal that step execution was completed
    signal_step_completion ${step_outd}

    display_end_step_message
}

########
step_b_explain_cmdline_opts()
{
    :
}

########
step_b_define_opts()
{
    # Initialize variables
    local cmdline=$1
    local jobspec=$2
    optlist=""
    
    # Define the -step-outd option, the output directory for the step,
    # which will have the same name of the step
    define_default_step_outd_opt "$cmdline" "$jobspec" optlist || exit 1

    # Save option list
    save_opt_list optlist
}

########
step_b()
{
    display_begin_step_message

    # Initialize variables
    local step_outd=`read_opt_value_from_line "$*" "-step-outd"`

    # sleep some time
    sleep 10

    # Signal that step execution was completed
    signal_step_completion ${step_outd}

    display_end_step_message
}

########
step_c_explain_cmdline_opts()
{
    :
}

########
step_c_define_opts()
{
    # Initialize variables
    local cmdline=$1
    local jobspec=$2
    optlist=""
    
    # Define the -step-outd option, the output directory for the step,
    # which will have the same name of the step
    define_default_step_outd_opt "$cmdline" "$jobspec" optlist || exit 1

    # Save option list
    save_opt_list optlist
}

########
step_c()
{
    display_begin_step_message

    # Initialize variables
    local step_outd=`read_opt_value_from_line "$*" "-step-outd"`

    # sleep some time
    sleep 10
    
    # Signal that step execution was completed
    signal_step_completion ${step_outd}

    display_end_step_message
}
