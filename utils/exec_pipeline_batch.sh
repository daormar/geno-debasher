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
    # TBD
}

########
move_outdirs()
{
    # TBD
}

########
add_cmd_to_array()
{
    # TBD
}

########
execute_batches()
{
    # Read file with exec_pipeline commands
    lineno=1
    declare -A PIPELINE_COMMANDS
    while exec_pipeline_cmd; do

        # Wait until number of simultaneous executions is below the given maximum
        if [ ${#PIPELINE_COMMANDS[@]} -ge ${maxp} ]; then
            wait_simul_exec_reduction "PIPELINE_COMMANDS"
        fi
        
        # Move output directories if requested
        if [ ${o_given} -eq 1 ]; then
            move_outdirs ${outd} "PIPELINE_COMMANDS"
        fi

        # Execute command
        ${exec_pipeline_cmd}

        # Add command to associative array
        add_cmd_to_array "${exec_pipeline_cmd}" "PIPELINE_COMMANDS"
        
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
