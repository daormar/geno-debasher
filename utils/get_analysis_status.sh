# *- bash -*

# INCLUDE BASH LIBRARY
. ${bindir}/bam_utils_lib

########
print_desc()
{
    echo "get_analysis_status get status of analysis steps"
    echo "type \"submit_bam_analysis --help\" to get usage information"
}

########
usage()
{
    echo "get_analysis_status       -d <string> [-s <string>]"
    echo "                          [--help]"
    echo ""
    echo "-d <string>               Directory where the analysis steps are stored"
    echo "-s <string>               Step name whose status should be determined"
    echo "--help                    Display this help and exit"
}

########
read_pars()
{
    d_given=0
    s_given=0
    debug=0
    while [ $# -ne 0 ]; do
        case $1 in
            "--help") usage
                      exit 1
                      ;;
            "--version") version
                         exit 1
                         ;;
            "-d") shift
                  if [ $# -ne 0 ]; then
                      adir=$1
                      d_given=1
                  fi
                  ;;
            "-s") shift
                  if [ $# -ne 0 ]; then
                      given_stepname=$1
                      s_given=1
                  fi
                  ;;
        esac
        shift
    done   
}

########
check_pars()
{
    if [ ${d_given} -eq 0 ]; then
        echo "Error! -d parameter not given!" >&2
        exit 1
    else
        if [ ! -d ${adir} ]; then
            echo "Error! analysis directory does not exist" >&2 
            exit 1
        fi

        if [ ! -f ${adir}/command_line.sh ]; then
            echo "Error! ${adir}/command_line.sh file is missing" >&2 
            exit 1
        fi
    fi
}

########
get_orig_workdir()
{
    local command_line_file=$1
    local workdir=`$HEAD -1 ${command_line_file} | $AWK '{print $2}'` || return 1
    echo $workdir
}

########
get_cmdline()
{
    local command_line_file=$1
    local cmdline=`$TAIL -1 ${command_line_file}`
    echo $cmdline
}

########
get_afile()
{
    local command_line_file=$1
    local cmdline=`$TAIL -1 ${command_line_file}`
    local afile=`read_opt_value_from_line "$cmdline" "-a"` || return 1
    echo $afile
}

########
process_status_for_afile()
{
    local dirname=$1
    command_line_file=$dirname/command_line.sh
    
    # Extract information from command_line.sh file
    local orig_workdir=`get_orig_workdir ${command_line_file}` || return 1
    local cmdline=`get_cmdline ${command_line_file}` || return 1
    local afile=`get_afile ${command_line_file}` || return 1

    # Change directory
    cd ${orig_workdir}

    # Load pipeline modules
    load_pipeline_modules $afile 2>/dev/null || return 1
        
    # Read information about the steps to be executed
    lineno=1
    analysis_finished=1
    while read jobspec; do
        local jobspec_comment=`analysis_jobspec_is_comment "$jobspec"`
        local jobspec_ok=`analysis_jobspec_is_ok "$jobspec"`
        if [ ${jobspec_comment} = "no" -a ${jobspec_ok} = "yes" ]; then
            # Extract step information
            local stepname=`extract_stepname_from_jobspec "$jobspec"`

            # If s option was given, continue to next iteration if step
            # name does not match with the given one
            if [ ${s_given} -eq 1 -a "${given_stepname}" != $stepname ]; then
                continue
            fi

            local script_define_opts_funcname=`get_script_define_opts_funcname ${stepname}`
            ${script_define_opts_funcname} ${cmdline} ${jobspec} || return 1

            # Check step status
            local status=`get_step_status ${dirname} ${stepname}`

            # Print status
            echo "STEP: $stepname ; STATUS: $status"

            # Revise value of analysis_finished variable
            if [ "${status}" != "${FINISHED_STEP_STATUS}" ]; then
                analysis_finished=0
            fi
        else
            if [ ${jobspec_ok} = "no" ]; then
                echo "Error: incorrect job specification at line $lineno" >&2
                return 1
            fi
        fi
        
        # Increase lineno
        lineno=`expr $lineno + 1`
        
    done < ${afile}

    # Return error if analysis is not finished
    if [ ${analysis_finished} -eq 0 ]; then
        return 1
    fi
}

########

if [ $# -eq 0 ]; then
    print_desc
    exit 1
fi

read_pars $@ || exit 1

check_pars || exit 1

process_status_for_afile ${adir} ${afile} || exit 1
