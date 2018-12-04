# *- bash -*

# INCLUDE BASH LIBRARY
. ${bindir}/bam_utils_lib

########
print_desc()
{
    echo "exec_pipeline executes general purpose pipelines"
    echo "type \"exec_pipeline --help\" to get usage information"
}

########
usage()
{
    echo "exec_pipeline        -a <string> -o <string> [--showopts]"
    echo "                     [-debug] [--version] [--help]"
    echo ""
    echo "-a <string>          File with analysis steps to be performed."
    echo "                     Expected format:"
    echo "                     <stepname> <account> <partition> <cpus> <mem> <time> <jobdeps=stepname1:...>"
    echo "-o <string>          Output directory"
    echo "--showopts           Show pipeline options (-a option should be provided)"
    echo "--debug              Do not execute pipeline, only print status"
    echo "                     and input parameters"
    echo "--version            Display version information and exit"
    echo "--help               Display this help and exit"
}

########
save_command_line()
{
    input_pars="$*"
    command_name=$0
}

########
read_pars()
{
    a_given=0
    showopts_given=0
    debug=0
    while [ $# -ne 0 ]; do
        case $1 in
            "--help") usage
                      exit 1
                      ;;
            "--version") version
                         exit 1
                         ;;
            "-a") shift
                  if [ $# -ne 0 ]; then
                      afile=$1
                      a_given=1
                  fi
                  ;;
            "-o") shift
                  if [ $# -ne 0 ]; then
                      outd=$1
                      o_given=1
                  fi
                  ;;
            "--showopts") showopts_given=1
                  ;;
            "--debug") debug=1
                      debug_opt="-debug"
                      ;;
        esac
        shift
    done   
}

########
check_pars()
{
    if [ ${a_given} -eq 0 ]; then   
        echo "Error! -a parameter not given!" >&2
        exit 1
    else
        if [ ! -f ${afile} ]; then
            echo "Error! file ${afile} does not exist" >&2
            exit 1
        fi
    fi

    if [ ${o_given} -eq 0 ]; then
        echo "Error! -o parameter not given!" >&2
        exit 1
    else
        if [ -d ${outd} ]; then
            echo "Warning! output directory does exist" >&2 
        fi
    fi
}

########
absolutize_file_paths()
{
    if [ ${a_given} -eq 1 ]; then   
        afile=`get_absolute_path ${afile}`
    fi
}

########
show_pipeline_opts()
{
    # Read input parameters
    local cmdline=$1
    local dirname=$2
    local afile=$3

    # Load pipeline modules
    load_pipeline_modules $afile || return 1
        
    # Read information about the steps to be executed
    while read jobspec; do
        local jobspec_comment=`analysis_jobspec_is_comment "$jobspec"`
        local jobspec_ok=`analysis_jobspec_is_ok "$jobspec"`
        if [ ${jobspec_comment} = "no" -a ${jobspec_ok} = "yes" ]; then
            # Extract step information
            local stepname=`extract_stepname_from_jobspec "$jobspec"`
            local script_explain_cmdline_opts_funcname=`get_script_explain_cmdline_opts_funcname ${stepname}`
            ${script_explain_cmdline_opts_funcname}
        fi
    done < ${afile}

    # Print parameters
    echo "* Pipeline parameters:"
    print_pipeline_opts
}

########
check_pipeline_pars()
{
    echo "* Checking pipeline parameters..." >&2
    
    # Read input parameters
    local cmdline=$1
    local dirname=$2
    local afile=$3

    # Load pipeline modules
    load_pipeline_modules $afile 2>/dev/null || return 1
        
    # Read information about the steps to be executed
    while read jobspec; do
        local jobspec_comment=`analysis_jobspec_is_comment "$jobspec"`
        local jobspec_ok=`analysis_jobspec_is_ok "$jobspec"`
        if [ ${jobspec_comment} = "no" -a ${jobspec_ok} = "yes" ]; then
            # Extract step information
            local stepname=`extract_stepname_from_jobspec "$jobspec"`
            local script_define_opts_funcname=`get_script_define_opts_funcname ${stepname}`
            local script_opts
            script_opts=`${script_define_opts_funcname} ${cmdline} ${jobspec}` || return 1
        fi
    done < ${afile}
}

########
create_dirs()
{
    mkdir -p ${outd} || { echo "Error! cannot create output directory" >&2; return 1; }

    mkdir -p ${outd}/scripts || { echo "Error! cannot create scripts directory" >&2; return 1; }

    # Create shared directories required by the pipeline steps
    # IMPORTANT NOTE: the following function can only be executed after
    # executing check_pipeline_pars
    create_pipeline_shdirs
}

########
print_command_line()
{
    echo ${command_line} > ${outd}/command_line.sh
}

########
get_jobdeps_from_detailed_spec()
{
    local jobdeps_spec=$1
    local jdeps=""

    # Iterate over the elements of the job specification: type1:stepname1,...,typen:stepnamen
    for dep_spec in `echo ${jobdeps_spec} | $SED 's/,/ /g'`; do
        local deptype=`echo ${dep_spec} | $AWK -F ":" '{print $1}'`
        local step=`echo ${dep_spec} | $AWK -F ":" '{print $2}'`
        
        # Check if there is a jid for the step
        local step_jid=${step}_jid
        if [ ! -z "${!step_jid}" ]; then
            if [ -z "${jdeps}" ]; then
                jdeps=${deptype}":"${!step_jid}
            else
                jdeps=${jdeps}","${deptype}":"${!step_jid}
            fi
        fi
    done
    echo ${jdeps}
}

########
get_jobdeps()
{
    local jobdeps_spec=$1
    case ${jobdeps_spec} in
            "afterok:all") apply_deptype_to_jobids ${step_jids} afterok
                    ;;
            "none") echo ""
                    ;;
            *) get_jobdeps_from_detailed_spec ${jobdeps_spec}
               ;;
    esac
}

########
archive_script()
{
    local script_filename=$1
        
    # Archive script with date info
    local curr_date=`date '+%Y_%m_%d'`
    cp ${script_filename} ${script_filename}.${curr_date}
}

# ########
# check_script_is_older_than_lib()
# {
#     local script_filename=$1
#     # Check if script exists
#     if [ -f ${script_filename} ]; then
#         # script exists
#         lib_timestamp=`get_lib_timestamp`
#         script_timestamp=`get_file_timestamp ${script_filename}`
#         if [ ${script_timestamp} -lt ${lib_timestamp} ]; then
#             echo 1
#         else
#             echo 0
#         fi
#     else
#         # script does not exist
#         echo 0
#     fi
# }

########
execute_step()
{
    # Initialize variables
    local cmdline=$1
    local dirname=$2
    local stepname=$3
    local jobspec=$4
    local step_outd=`get_step_dirname ${outd} ${stepname}`
    
    # Execute step

    ## Obtain step status
    local status=`${bindir}/get_analysis_status -d ${dirname} -s "${stepname}"`
    echo "STEP: ${stepname} ; STATUS: ${status} ; JOBSPEC: ${jobspec}" >&2

    ## Decide whether the step should be executed
    if [ "${status}" != "FINISHED" ]; then
        # Initialize script variables
        local script_filename=`get_script_filename ${stepname}`
        local step_function=`get_step_function ${stepname}`
        local script_define_opts_funcname=`get_script_define_opts_funcname ${stepname}`
        local script_opts
        script_opts=`${script_define_opts_funcname} ${cmdline} ${jobspec}` || return 1
        echo "-> ${stepname} options: ${script_opts}" >&2
        
        # Create script
        create_script ${script_filename} ${step_function} "${script_opts}"

        # Archive script
        archive_script ${script_filename}

        # Execute script
        reset_outdir_for_step ${dirname} ${stepname} || return 1
        local jobdeps_spec=`extract_jobdeps_spec_from_jobspec "$jobspec"`
        local jobdeps="`get_jobdeps ${jobdeps_spec}`"
        local stepname_jid=${stepname}_jid
        launch ${script_filename} ${jobspec} "${jobdeps}" ${stepname_jid} || return 1
        
        # Update variables storing jids
        step_jids="${step_jids}:${!stepname_jid}"
    else
        local script_filename=`get_script_filename ${stepname}`
        prev_exec_script_older_than_lib=`check_script_is_older_than_lib ${script_filename}`
        if [ ${prev_exec_script_older_than_lib} -eq 1 ]; then
            echo "Warning: last execution of this script used an outdated shell library">&2
        fi
    fi
}

########
debug_step()
{
    # Initialize variables
    local cmdline=$1
    local dirname=$2
    local stepname=$3
    local jobspec=$4
    local step_outd=`get_step_dirname ${outd} ${stepname}`
    
    # Debug step

    ## Obtain step status
    local status=`${bindir}/get_analysis_status -d ${dirname} -s "${stepname}"`
    local step_function=`get_step_function ${stepname}`
    echo "STEP: ${stepname} ; STATUS: ${status} ; JOBSPEC: ${jobspec}" >&2

    ## Obtain step options
    local script_define_opts_funcname=`get_script_define_opts_funcname ${stepname}`
    local script_opts
    script_opts=`${script_define_opts_funcname} ${cmdline} ${jobspec}` || return 1
    echo "-> ${stepname} options: ${script_opts}" >&2
}

########
execute_pipeline_steps()
{
    echo "* Executing pipeline steps..." >&2

    # Read input parameters
    local cmdline=$1
    local dirname=$2
    local afile=$3

    # Load pipeline modules
    load_pipeline_modules $afile || return 1
    
    # step_jids will store the job ids of the analysis steps
    step_jids=""
    
    # Read information about the steps to be executed
    while read jobspec; do
        local jobspec_comment=`analysis_jobspec_is_comment "$jobspec"`
        local jobspec_ok=`analysis_jobspec_is_ok "$jobspec"`
        if [ ${jobspec_comment} = "no" -a ${jobspec_ok} = "yes" ]; then
            # Extract step name
            local stepname=`extract_stepname_from_jobspec "$jobspec"`

            # Decide whether to execute or debug step
            if [ $debug -eq 0 ]; then
                execute_step ${cmdline} ${dirname} ${stepname} ${jobspec} || return 1
            else
                debug_step ${cmdline} ${dirname} ${stepname} ${jobspec} || return 1                
            fi
        fi
    done < ${afile}
}

########

if [ $# -eq 0 ]; then
    print_desc
    exit 1
fi

# Save command line
command_line="$0 $*"

read_pars $@ || exit 1

check_pars || exit 1

absolutize_file_paths || exit 1

if [ ${showopts_given} -eq 1 ]; then
    show_pipeline_opts "${command_line}" ${outd} ${afile} || exit 1
else
    check_pipeline_pars "${command_line}" ${outd} ${afile} || exit 1

    create_dirs || exit 1

    print_command_line || exit 1

    execute_pipeline_steps "${command_line}" ${outd} ${afile} || exit 1
fi
