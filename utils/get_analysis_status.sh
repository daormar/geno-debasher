# *- bash -*

########
print_desc()
{
    echo "get_analysis_status get status of analysis steps"
    echo "type \"submit_bam_analysis --help\" to get usage information"
}

########
usage()
{
    echo "get_analysis_status       -d <string> [-s <string>| -a <string>]"
    echo "                          [--help]"
    echo ""
    echo "-d <string>               Directory where the analysis steps are stored."
    echo "-s <string>               Step name whose status should be get"
    echo "-a <string>               File with steps to be performed. Expected format:"
    echo "                          <stepname> <cpus> <mem> <time>"
    echo "--help                    Display this help and exit."
}

########
read_pars()
{
    d_given=0
    s_given=0
    a_given=0
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
                      stepname=$1
                      s_given=1
                  fi
                  ;;
            "-a") shift
                  if [ $# -ne 0 ]; then
                      afile=$1
                      a_given=1
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
            echo "Warning! analysis directory does not exist" >&2 
        fi
    fi

    if [ ${a_given} -eq 0 -a ${s_given} -eq 0 ]; then
        echo "Error! -a or -s parameter not given!" >&2
        exit 1        
    fi

    if [ ${a_given} -eq 1 -a ${s_given} -eq 1 ]; then
        echo "Error! -a and -s parameters cannot be given simultaneously!" >&2
        exit 1        
    fi

    if [ ${a_given} -eq 1 ]; then   
        if [ ! -f ${afile} ]; then
            echo "Error! file ${afile} does not exist" >&2
            exit 1
        fi
    fi
}

########
get_step_dirname()
{
    local_stepname=$1
    echo ${adir}/${local_stepname}
}

########
get_step_status()
{
    local_stepname=$1
    local_stepdirname=`get_step_dirname ${local_stepname}`
    
    if [ -d ${local_stepdirname} ]; then
        if [ -f ${local_stepdirname}/finished ]; then
            echo "FINISHED"
        else
            echo "UNFINISHED"
        fi
    else
        echo "TO-DO"
    fi
}

########
entry_is_ok()
{
    local_entry=$1
    echo "${local_entry}" | $AWK '{if(NF>=4) print"yes\n"; else print"no\n"}'
}

########
extract_stepname_from_entry()
{
    local_entry=$1
    echo "${local_entry}" | $AWK '{print $1}'
}

########
get_status_for_afile()
{
    local_afile=$1
    
    # Read information about the steps to be executed
    while read entry; do
        entry_ok=`entry_is_ok "$entry"`
        if [ ${entry_ok} = "yes" ]; then
            # Extract entry information
            stepname=`extract_stepname_from_entry "$entry"`

            # Check step status
            status=`get_step_status $stepname`

            # Print status
            echo "STEP: $stepname ; STATUS: $status"
        fi
    done < ${local_afile}
}

########
process_pars()
{
    if [ ${s_given} -eq 1 ]; then
        get_step_status $stepname
    fi

    if [ ${a_given} -eq 1 ]; then
        get_status_for_afile ${afile}
    fi
}

########

if [ $# -eq 0 ]; then
    print_desc
    exit 1
fi

read_pars $@ || exit 1

check_pars || exit 1

process_pars
