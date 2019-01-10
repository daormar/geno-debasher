# *- bash -*

# INCLUDE BASH LIBRARY
. ${bindir}/bam_utils_lib

########
print_desc()
{
    echo "exec_pipeline_batch executes a batch of pipelines"
    echo "type \"exec_pipeline_batch --help\" to get usage information"
}

########
usage()
{
    echo "exec_pipeline_batch       -f <string> -m <int> [-o <string>]"
    echo "                          [--help]"
    echo ""
    echo "-f <string>               File with a set of exec_pipeline commands (one"
    echo "                          per line)"
    echo "-m <string>               Maximum number of pipelines executed simultaneously"
    echo "-o <string>               Output directory where the pipeline output should be"
    echo "                          moved (if not given, the output directories are"
    echo "                          provided by the exec_pipeline commands)"
    echo "--help                    Display this help and exit"
}

########
read_pars()
{
    f_given=0
    m_given=0
    o_given=0
    while [ $# -ne 0 ]; do
        case $1 in
            "--help") usage
                      exit 1
                      ;;
            "--version") version
                         exit 1
                         ;;
            "-f") shift
                  if [ $# -ne 0 ]; then
                      file=$1
                      f_given=1
                  fi
                  ;;
            "-m") shift
                  if [ $# -ne 0 ]; then
                      maxp=$1
                      m_given=1
                  fi
                  ;;
            "-o") shift
                  if [ $# -ne 0 ]; then
                      outd=$1
                      o_given=1
                  fi
                  ;;
        esac
        shift
    done   
}

########
wait_simul_exec_reduction()
{
    local maxp=$1
    local -n assoc_array=$2
    local SLEEP_TIME=10
    end=0
    
    while [ ${end} -eq 0 ] ; do
        # Iterate over active pipelines
        num_finished_pipelines=0
        num_unfinished_pipelines=0
        for pipeline_outd in "${!assoc_array[@]}"; do
            # Check if pipeline has finished execution
            ${bindir}/get_analysis_status -d ${pipeline_outd} > /dev/null 2>&1
            local exit_code=$?

            case ${exit_code} in
                ${ANALYSIS_FINISHED_EXIT_CODE})
                    num_finished_pipelines=`expr ${num_finished_pipelines} + 1`
                    ;;
                ${ANALYSIS_UNFINISHED_EXIT_CODE})
                    num_unfinished_pipelines=`expr ${num_unfinished_pipelines} + 1`                    
                    ;;
            esac
        done

        num_active_pipelines=${#assoc_array[@]} 

        # Sanity check: if all pipelines are unfinished, then it is not
        # possible to continue execution
        if [ ${num_unfinished_pipelines} -eq ${num_active_pipelines} ]; then
            echo "Error: all active pipelines are unfinished" >&2
            return 1
        fi
        
        # Obtain number of pending pipelines
        pending_pipelines=`expr ${num_active_pipelines} - ${num_finished_pipelines}`

        # Wait if number of pending pipelines is equal or greater than
        # maximum
        if [ ${pending_pipelines} -ge ${maxp} ]; then
            sleep ${SLEEP_TIME}
        else
            end=1
        fi
    done
}

########
update_active_pipelines()
{
    local outd=$1
    local -n assoc_array=$2

    # Iterate over active pipelines
    for pipeline_outd in "${!assoc_array[@]}"; do
        # Check if pipeline has finished execution
        ${bindir}/get_analysis_status -d ${pipeline_outd} > /dev/null 2>&1
        local exit_code=$?
        
        if [ ${exit_code} -eq ${ANALYSIS_FINISHED_EXIT_CODE} ]; then
            # Remove pipeline from array of active pipelines
            unset assoc_array[${pipeline_outd}]
            
            # Move directory if requested
            if [ ! -z "${outd}" ]; then
                mv ${pipeline_outd} ${outd}
            fi            
        fi
    done
}

########
add_cmd_to_assoc_array()
{
    local cmd=$1
    local -n assoc_array=$2

    # Extract output directory from command
    local dir=`read_opt_value_from_line "-o" ${cmd}`

    # Add command to associative array if directory was sucessfully retrieved
    if [ ${dir} != ${OPT_NOT_FOUND} ]; then
        assoc_array[${dir}]=${cmd}
    fi
}

########
execute_batches()
{
    # Read file with exec_pipeline commands
    lineno=1
    declare -A PIPELINE_COMMANDS
    while exec_pipeline_cmd; do

        # Wait until number of simultaneous executions is below the given maximum
        wait_simul_exec_reduction ${maxp} "PIPELINE_COMMANDS" || return 1
            
        # Update array of active pipelines
        update_active_pipelines ${outd} "PIPELINE_COMMANDS"

        # Execute command
        ${exec_pipeline_cmd} || return 1

        # Add command to associative array
        add_cmd_to_assoc_array "${exec_pipeline_cmd}" "PIPELINE_COMMANDS"
        
        # Increase lineno
        lineno=`expr $lineno + 1`

    done < ${file}
}

########
check_pars()
{
    if [ ${f_given} -eq 0 ]; then
        echo "Error! -f parameter not given!" >&2
        exit 1
    else
        if [ ! -f ${file} ]; then
            echo "Error! file ${file} does not exist" >&2 
            exit 1
        fi
    fi

    if [ ${m_given} -eq 0 ]; then
        echo "Error! -m parameter not given!" >&2
        exit 1
    fi

    if [ ${o_given} -eq 1 ]; then
        if [ ! -d ${outd} ]; then
            echo "Error! output directory does not exist" >&2 
            exit 1
        fi
    fi
}

########

if [ $# -eq 0 ]; then
    print_desc
    exit 1
fi

read_pars $@ || exit 1

check_pars || exit 1

execute_batches || exit 1
