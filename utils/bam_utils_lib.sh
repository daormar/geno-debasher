# *- bash -*

#############
# CONSTANTS #
#############

NOFILE="_NONE_"
OPT_NOT_FOUND="_OPT_NOT_FOUND_"
DEP_NOT_FOUND="_DEP_NOT_FOUND_"
INVALID_JID="_INVALID_JID_"
VOID_VALUE="_VOID_VALUE_"

####################
# GLOBAL VARIABLES #
####################

# Declare associative array to store help about pipeline options
declare -A PIPELINE_OPT_DESC
declare -A PIPELINE_OPT_TYPE

# Declare associative array to store names of directories
declare -A PIPELINE_STEPDIRS
declare -A PIPELINE_SHDIRS

# Declare associative array to memoize command line options
declare -A MEMOIZED_OPTS

# Declare string variable to store last processed command line
LAST_PROC_LINE_MEMOPTS=""

# Declare variable used to save option lists for scripts
declare SCRIPT_OPT_LIST

#####################
# GENERAL FUNCTIONS #
#####################

init_bash_shebang_var()
{
    echo "#!${BASH}"
}

########
is_absolute_path()
{
    case $1 in
        /*) echo 1 ;;
       *) echo 0 ;;
    esac
}

########
get_absolute_path()
{
    local file=$1
    # Check if an absolute path was given
    local absolute=`is_absolute_path $file`
    if [ $absolute -eq 1 ]; then
        echo $file
    else
        local oldpwd=$PWD
        local basetmp=`$BASENAME $PWD/$file`
        local dirtmp=`$DIRNAME $PWD/$file`
        cd $dirtmp
        local result=${PWD}/${basetmp}
        cd $oldpwd
        echo $result
    fi
}

########
exclude_readonly_vars()
{
    $AWK -F "=" 'BEGIN{
                         readonlyvars["BASHOPTS"]=1
                         readonlyvars["BASH_VERSINFO"]=1
                         readonlyvars["EUID"]=1
                         readonlyvars["PPID"]=1
                         readonlyvars["SHELLOPTS"]=1
                         readonlyvars["UID"]=1
                        }
                        {
                         if(!($1 in readonlyvars)) printf"%s\n",$0
                        }'
}

########
exclude_bashisms()
{
    $AWK '{if(index($1,"=(")==0) printf"%s\n",$0}'
}

########
create_script()
{
    # Init variables
    local name=$1
    local command=$2
    local script_opts=$3
    
    # Write bash shebang
    local BASH_SHEBANG=`init_bash_shebang_var`
    echo ${BASH_SHEBANG} > ${name} || return 1

    # Write SLURM commands
    echo "#SBATCH --job-name=${command}" >> ${name} || return 1
    echo "#SBATCH --output=${name}.slurm_out" >> ${name} || return 1

    # Write environment variables
    set | exclude_readonly_vars | exclude_bashisms >> ${name} || return 1
    
    # Write command to be executed
    echo "${command} ${script_opts}" >> ${name} || return 1

    # Give execution permission
    chmod u+x ${name} || return 1
}

########
get_file_timestamp()
{
    local file=$1
    ${PYTHON} -c "import os; t=os.path.getmtime('${file}'); print int(t)"
}

########
get_account_opt()
{
    local account=$1

    if [ -z "${account}" ]; then
        echo ""
    else
        echo "-A ${account}"
    fi
}

########
get_partition_opt()
{
    local partition=$1

    if [ -z "${partition}" ]; then
        echo ""
    else
        echo "--partition=${partition}"
    fi
}

########
get_script_filename() 
{
    local stepname=$1
    
    echo ${dirname}/scripts/${stepname}
}

########
remove_suffix_from_stepname()
{
    local stepname=$1
    
    echo ${stepname} | $AWK '{if(index($1,"__")==0){print $1} else{printf "%s\n",substr($1,1,index($1,"__")-1)}}'
}

########
get_step_function()
{
    local stepname=$1

    local stepname_wo_suffix=`remove_suffix_from_stepname ${stepname}`
    
    echo "${stepname_wo_suffix}"
}

########
get_script_explain_cmdline_opts_funcname()
{
    local stepname=$1

    local stepname_wo_suffix=`remove_suffix_from_stepname ${stepname}`

    echo ${stepname_wo_suffix}_explain_cmdline_opts
}

########
get_script_define_opts_funcname()
{
    local stepname=$1

    local stepname_wo_suffix=`remove_suffix_from_stepname ${stepname}`

    echo ${stepname_wo_suffix}_define_opts
}

########
find_dependency_for_step()
{
    local jobspec=$1
    local stepname_part=$2

    jobdeps=`extract_jobdeps_from_jobspec "$jobspec"`
    prevIFS=$IFS
    IFS=','
    for local_dep in ${jobdeps}; do
        local stepname_part_in_dep=`echo ${local_dep} | $AWK -F ":" '{print $2}'`
        if [ ${stepname_part_in_dep} = ${stepname_part} ]; then
            echo ${local_dep}
            IFS=${prevIFS}
            return 0
        fi
    done
    IFS=${prevIFS}
    echo ${DEP_NOT_FOUND}
    return 1
}

########
get_default_outd_for_dep()
{
    local outd=$1
    local dep=$2

    if [ -z "${dep}" ]; then
        echo ""
    else
        local stepname_part=`echo ${dep} | $AWK -F ":" '{print $2}'`
        get_step_dirname ${outd} ${stepname_part}
    fi
}

########
apply_deptype_to_jobids()
{
    # Initialize variables
    local jids=$1
    local deptype=$2

    # Apply deptype
    result=""
    prevIFS=$IFS
    IFS=','
    for local_jid in ${jids}; do
        if [ -z "" ]; then
            result=${deptype}:${local_jid}
        else
            result=${result}","${deptype}:${local_jid}
        fi
    done
    IFS=${prevIFS}

    echo $result
}

########
get_slurm_dependency_opt()
{
    local jobdeps=$1

    # Create dependency option
    if [ -z "${jobdeps}" ]; then
        echo ""
    else
        echo "--dependency=${jobdeps}"
    fi
}

########
launch()
{
    # Initialize variables
    local file=$1
    local jobspec=$2
    local jobdeps=$3
    local outvar=$4

    # Launch file
    if [ -z "${SBATCH}" ]; then
        ## Launch without using any scheduler
        ${file} > ${file}.log 2>&1 || return 1
        eval "${outvar}=\"\""
    else
        ## Launch using slurm
        # Retrieve specification
        local account=`extract_account_from_jobspec "$jobspec"`
        local partition=`extract_partition_from_jobspec "$jobspec"`
        local cpus=`extract_cpus_from_jobspec "$jobspec"`
        local mem=`extract_mem_from_jobspec "$jobspec"`
        local time=`extract_time_from_jobspec "$jobspec"`

        # Define options for sbatch
        local account_opt=`get_account_opt ${account}`
        local partition_opt=`get_partition_opt ${partition}`
        local dependency_opt=`get_slurm_dependency_opt "${jobdeps}"`

        # Submit job
        local jid=$($SBATCH --cpus-per-task=${cpus} --mem=${mem} --time ${time} --parsable ${account_opt} ${partition_opt} ${dependency_opt} ${file})

        # Check for errors
        if [ -z "$jid" ]; then
            jid=${INVALID_JID}
        fi

        # Assign output variable to job identifier
        eval "${outvar}='${jid}'"
    fi
}

########
launch_step()
{
    # Initialize variables
    local stepname=$1
    local jobspec=$2
    local jobdeps=$3
    local script_opts=$4
    local jid=$5

    # Create script
    create_script ${tmpdir}/scripts/${stepname} ${stepname} "${script_opts}" || return 1

    # Launch script
    launch ${tmpdir}/scripts/${stepname} "${jobspec}" ${jobdeps} ${jid} || return 1
}

########
get_step_info()
{
    local stepname=$1
    local infofile=$2

    $AWK -v stepname=${stepname} '{if($1==stepname) print $0}' ${infofile}
}

########
analysis_jobspec_is_comment()
{
    local jobspec=$1
    echo "${jobspec}" | $AWK '{if(index($1,"#")==1) print"yes\n"; else print"no\n"}'
}

########
analysis_jobspec_is_ok()
{
    local jobspec=$1
    echo "${jobspec}" | $AWK '{if(NF>=4) print"yes\n"; else print"no\n"}'
}

########
extract_stepname_from_jobspec()
{
    local jobspec=$1
    echo "${jobspec}" | $AWK '{print $1}'
}

########
extract_account_from_jobspec()
{
    local jobspec=$1
    echo "${jobspec}" | $AWK '{print $2}'
}

########
extract_partition_from_jobspec()
{
    local jobspec=$1
    echo "${jobspec}" | $AWK '{print $3}'
}

########
extract_cpus_from_jobspec()
{
    local jobspec=$1
    echo "${jobspec}" | $AWK '{print $4}'
}

########
extract_mem_from_jobspec()
{
    local jobspec=$1
    echo "${jobspec}" | $AWK '{print $5}'
}

########
extract_time_from_jobspec()
{
    local jobspec=$1
    echo "${jobspec}" | $AWK '{print $6}'
}

########
extract_jobdeps_from_jobspec()
{
    local jobspec=$1
    echo "${jobspec}" | $AWK '{print substr($7,9)}'
}

########
get_pipeline_modules()
{
    local afile=$1
    local modules=`$AWK '{if($1=="#import") {$1=""; printf "%s ",$0}}' $afile | $AWK '{for(i=1;i<=NF;++i) printf"%s",$i}'`
    echo ${modules}
}

########
determine_full_module_name()
{
    local module=$1
    local absolute=`is_absolute_path $module`
    if [ $absolute -eq 1 ]; then
        fullmodname=${module}
    else
        fullmodname=${bindir}${module}
    fi

    echo $fullmodname
}

########
load_pipeline_module()
{
    local module=$1

    # Determine full module name
    fullmodname=`determine_full_module_name $module`

    echo "Loading module $module (${fullmodname})..." >&2

    # Check that module file exists
    if [ -f ${fullmodname} ]; then
        . ${fullmodname}
    else
        return 1
    fi
}

########
load_pipeline_modules()
{
    local afile=$1

    file_exists $afile || { echo "Error: file $afile does not exist" >&2 ; return 1; }
    
    comma_sep_modules=`get_pipeline_modules $afile`
    
    if [ -z "${comma_sep_modules}" ]; then
        echo "Error: no pipeline modules were given" >&2
        return 1
    else
        # Load modules
        prevIFS=$IFS
        IFS=','
        for mod in ${comma_sep_modules}; do
            load_pipeline_module $mod || { "Error: module ${fullmodname} does not exist" >&2 ; return 1 }
        done
        IFS=${prevIFS}
    fi
}

########
get_pipeline_fullmodnames()
{
    local afile=$1

    file_exists $afile || { echo "Error: file $afile does not exist" >&2 ; return 1; }
    
    comma_sep_modules=`get_pipeline_modules $afile`
    
    if [ -z "${comma_sep_modules}" ]; then
        echo "Warning: no pipeline modules were given" >&2
    else
        # Get names
        local fullmodnames
        prevIFS=$IFS
        IFS=','
        for mod in ${comma_sep_modules}; do
            local fullmodname=`determine_full_module_name $mod`
            if [ -z "${fullmodnames}" ]; then
                fullmodnames=${fullmodname}
            else
                fullmodnames="${fullmodnames} ${fullmodname}"
            fi
        done
        IFS=${prevIFS}
        echo "${fullmodnames}"
    fi
}

########
get_default_step_dirname()
{
    local dirname=$1
    local stepname=$2
    echo ${dirname}/${stepname}
}

########
get_step_dirname()
{
    local dirname=$1
    local stepname=$2

    if [ -z "${PIPELINE_STEPDIRS[${stepname}]}" ]; then
        echo "Warning! directory for step ${stepname} not defined, assuming default name" >&2
        get_default_step_dirname $dirname $stepname
    else
        echo ${PIPELINE_STEPDIRS[${stepname}]}
    fi
}

########
reset_outdir_for_step() 
{
    local dirname=$1
    local stepname=$2
    local outd=`get_step_dirname ${dirname} ${stepname}`

    if [ -d ${outd} ]; then
        echo "Warning: ${stepname} output directory already exists but analysis was not finished, directory content will be removed">&2
        rm -rf ${outd}/* || { echo "Error! could not clear output directory" >&2; return 1; }
    else
        mkdir ${outd} || { echo "Error! cannot create output directory" >&2; return 1; }
    fi
}

########
get_step_status()
{
    local dirname=$1
    local stepname=$2
    local stepdirname=`get_step_dirname ${dirname} ${stepname}`
    
    if [ -d ${stepdirname} ]; then
        if [ -f ${stepdirname}/finished ]; then
            echo "FINISHED"
        else
            echo "UNFINISHED"
        fi
    else
        echo "TO-DO"
    fi
}

########
display_begin_step_message()
{
    if [ -z "${SLURM_JOB_ID}" ]; then
        echo "Step started at `date`" >&2
    else
        echo "Step started at `date` (SLURM_JOB_ID= ${SLURM_JOB_ID})" >&2
    fi
}

########
display_end_step_message()
{
    echo "Step finished at `date`" >&2
}

########
errmsg()
{
    local msg=$1
    echo "$msg" >&2
}

########
logmsg()
{
    local msg=$1
    echo "$msg" >&2
}

########
file_exists()
{
    local file=$1
    if [ -f $file ]; then
        return 0
    else
        return 1
    fi
}

########
dir_exists()
{
    local dir=$1
    if [ -d $dir ]; then
        return 0
    else
        return 1
    fi
}

########
signal_step_completion()
{
    local step_outd=$1
    touch ${step_outd}/finished
}

########
memoize_opts()
{
    shift
    while [ $# -ne 0 ]; do
        # Check if argument is an option
        if [ ${1:0:1} = "-" -o ${1:0:2} = "--"  ]; then
            # Argument is option
            local opt=$1
            shift
            if [ $# -ne 0 ]; then
                # Check if next argument is option
                if [ ${1:0:1} = "-" -o ${1:0:2} = "--"  ]; then
                    # Next argument is option
                    MEMOIZED_OPTS[$opt]=${VOID_VALUE}
                else
                    # Next argument is value
                    value=$1
                    MEMOIZED_OPTS[$opt]=$value
                    shift
                fi
            else
                # There are no more arguments
                MEMOIZED_OPTS[$opt]=${VOID_VALUE}
            fi
        else
            echo "Warning: unexpected argument ($1), skipping..." >&2
            shift
        fi
    done
}

########
check_opt_given()
{
    line=$1
    opt=$2
    # Convert string to array
    local array
    IFS=' ' read -r -a array <<< $line
    # Scan array
    i=0
    while [ $i -lt ${#array[@]} ]; do
        if [ ${array[$i]} = "${opt}" ]; then
            return 0
        fi
        i=$((i+1))
    done

    # Option not given
    return 1
}

########
check_memoized_opt()
{
    local opt=$1

    # Check if option was not given
    if [ -z "${MEMOIZED_OPTS[$opt]}" ]; then
        return 1
    else
        return 0
    fi
}

########
check_opt_given_memoiz()
{
    line=$1
    opt=$2

    if [ "${LAST_PROC_LINE_MEMOPTS}" = "$line" ]; then
        # Given line was previously processed, return memoized result
        check_memoized_opt $opt || return 1
    else
        # Process not memoized line
        memoize_opts $line
        
        # Store processed line
        LAST_PROC_LINE_MEMOPTS="$line"

        # Return result
        check_memoized_opt $opt || return 1
    fi    
}

########
read_opt_value_from_line()
{
    line=$1
    opt=$2
    
    # Convert string to array
    local array
    IFS=' ' read -r -a array <<< $line
    # Scan array
    i=0
    while [ $i -lt ${#array[@]} ]; do
        if [ ${array[$i]} = "${opt}" ]; then
            i=$((i+1))
            if [ $i -lt ${#array[@]} ]; then
                echo ${array[$i]}
                return 0
            fi
        fi
        i=$((i+1))
    done

    # Option not given
    echo ${OPT_NOT_FOUND}
    return 1
}

########
read_memoized_opt_value()
{
    local opt=$1

    # Check if option was not given or it had void value
    if [ -z "${MEMOIZED_OPTS[$opt]}" -o "${MEMOIZED_OPTS[$opt]}" = ${VOID_VALUE} ]; then
        echo ${OPT_NOT_FOUND}
        return 1
    else
        echo ${MEMOIZED_OPTS[$opt]}
        return 0
    fi
}

########
read_opt_value_from_line_memoiz()
{
    line=$1
    opt=$2

    if [ "${LAST_PROC_LINE_MEMOPTS}" = "$line" ]; then
        # Given line was previously processed, return memoized result
        _OPT_VALUE_=`read_memoized_opt_value $opt` || return 1
    else
        # Process not memoized line
        memoize_opts $line
        
        # Store processed line
        LAST_PROC_LINE_MEMOPTS="$line"

        # Return result
        _OPT_VALUE_=`read_memoized_opt_value $opt` || return 1
    fi
}

########
explain_cmdline_opt()
{
    local opt=$1
    local type=$2
    local desc=$3

    # Store option in associative array
    PIPELINE_OPT_TYPE[$opt]=$type
    PIPELINE_OPT_DESC[$opt]=$desc
}

########
explain_cmdline_opt_wo_value()
{
    local opt=$1
    local desc=$2

    # Store option in associative array
    PIPELINE_OPT_TYPE[$opt]=""
    PIPELINE_OPT_DESC[$opt]=$desc
}

########
print_pipeline_opts()
{
    for opt in ${!PIPELINE_OPT_TYPE[@]}; do
        if [ -z ${PIPELINE_OPT_TYPE[$opt]} ]; then
            echo "${opt} ${PIPELINE_OPT_DESC[$opt]}"
        else
            echo "${opt} ${PIPELINE_OPT_TYPE[$opt]} ${PIPELINE_OPT_DESC[$opt]}"
        fi
    done
}

########
define_cmdline_opt()
{
    local cmdline=$1
    local opt=$2
    local varname=$3

    # Get value for option
    # local value
    # value=`read_opt_value_from_line "$cmdline" $opt` || { errmsg "$opt option not found" ; return 1; }
    read_opt_value_from_line_memoiz "$cmdline" $opt || { errmsg "$opt option not found" ; return 1; }
    local value=${_OPT_VALUE_}
    
    # Add option
    define_opt $opt $value $varname
}

########
define_cmdline_opt_wo_value()
{
    local cmdline=$1
    local opt=$2
    local varname=$3

    # Get value for option
    check_opt_given "$cmdline" $opt || { errmsg "$opt option not found" ; return 1; }

    # Add option
    define_opt_wo_value $opt $varname
}

########
define_cmdline_nonmandatory_opt()
{
    local cmdline=$1
    local opt=$2
    local default_value=$3
    local varname=$4

    # Get value for option
    # local value
    # value=`read_opt_value_from_line "$cmdline" $opt`
    read_opt_value_from_line_memoiz "$cmdline" $opt
    local value=${_OPT_VALUE_}

    if [ $value = ${OPT_NOT_FOUND} ]; then
        value=${default_value}
    fi
    
    # Add option
    define_opt $opt $value $varname    
}

########
define_cmdline_infile_opt()
{
    local cmdline=$1
    local opt=$2
    local varname=$3

    # Get value for option
    # local value
    # value=`read_opt_value_from_line "$cmdline" $opt` || { errmsg "$opt option not found" ; return 1; }
    read_opt_value_from_line_memoiz "$cmdline" $opt || { errmsg "$opt option not found" ; return 1; }
    local value=${_OPT_VALUE_}

    if [ $value != ${NOFILE} ]; then
        # Check if file exists
        file_exists $value || { errmsg "file $value does not exist ($opt option)" ; return 1; }
    fi

    # Absolutize path
    value=`get_absolute_path ${value}`
    
    # Add option
    define_opt $opt $value $varname
}

########
define_cmdline_opt_shdir()
{
    local cmdline=$1
    local opt=$2
    local varname=$3

    # Get value for option
    # local value
    # value=`read_opt_value_from_line "$cmdline" $opt` || { errmsg "$opt option not found" ; return 1; }
    read_opt_value_from_line_memoiz "$cmdline" $opt || { errmsg "$opt option not found" ; return 1; }
    local value=${_OPT_VALUE_}

    # Add option
    define_opt $opt $value $varname

    # Store shared directory name in associative array
    PIPELINE_SHDIRS["-$opt"]=$value
}

########
define_cmdline_nonmandatory_opt_shdir()
{
    local cmdline=$1
    local opt=$2
    local default_value=$3
    local varname=$4

    # Get value for option
    # local value
    # value=`read_opt_value_from_line "$cmdline" $opt`
    read_opt_value_from_line_memoiz "$cmdline" $opt
    local value=${_OPT_VALUE_}

    if [ $value = ${OPT_NOT_FOUND} ]; then
        value=${default_value}
    fi

    # Store shared directory name in associative array
    PIPELINE_SHDIRS["-$opt"]=$value
}

########
define_step_outd_opt()
{
    local stepname=$1
    local step_outd=$2
    local varname=$3

    # Add option
    define_opt "-step-outd" ${step_outd} $varname

    # Store directory name for step in associative array
    PIPELINE_STEPDIRS[${stepname}]=${step_outd}
}

########
define_default_step_outd_opt()
{
    local cmdline=$1
    local jobspec=$2
    local varname=$3

    # Get full path of directory
    # local outd
    # outd=`read_opt_value_from_line "$cmdline" "-o"` || return 1
    read_opt_value_from_line_memoiz "$cmdline" "-o" || return 1
    local outd=${_OPT_VALUE_}

    outd=`get_absolute_path ${outd}`
    local stepname=`extract_stepname_from_jobspec ${jobspec}`
    local step_outd=`get_default_step_dirname ${outd} ${stepname}`

    # Add option
    define_step_outd_opt ${stepname} ${step_outd} ${varname}
}

########
define_opt()
{
    local opt=$1
    local value=$2
    local varname=$3

    eval "${varname}='${!varname} ${opt} ${value}'"
}

########
define_opt_wo_value()
{
    local opt=$1
    local varname=$2

    eval "${varname}='${!varname} ${opt}'"
}

########
define_infile_opt()
{
    local opt=$1
    local value=$2
    local varname=$3

    # Check if file exists
    file_exists $value || { errmsg "file $value does not exist ($opt option)" ; return 1; }

    # Absolutize path
    value=`get_absolute_path ${value}`

    eval "${varname}='${!varname} ${opt} ${value}'"
}

########
define_indir_opt()
{
    local opt=$1
    local value=$2
    local varname=$3

    # Check if file exists
    dir_exists $value || { errmsg "directory $value does not exist ($opt option)" ; return 1; }

    # Absolutize path
    value=`get_absolute_path ${value}`

    eval "${varname}='${!varname} ${opt} ${value}'"
}

########
create_pipeline_shdirs()
{
    local outd=$1
    
    for dirname in "${PIPELINE_SHDIRS[@]}"; do
        absdir=`get_absolute_shdirname $outd $dirname`
        if [ ! -d ${absdir} ]; then
           mkdir ${absdir} || exit 1
        fi
    done
}

########
get_absolute_shdirname()
{
    local outd=$1
    local shdirname=$2
    echo ${outd}/${shdirname}
}

########
get_default_shdirname()
{
    local cmdline=$1
    local shdiropt=$2

    # Get full path of directory
    # local outd
    # outd=`read_opt_value_from_line "$cmdline" "-o"` || return 1
    read_opt_value_from_line_memoiz "$cmdline" "-o" || return 1
    local outd=${_OPT_VALUE_}
    outd=`get_absolute_path ${outd}`

    # Get name of shared dir
    # local shdir
    # shdir=`read_opt_value_from_line "$cmdline" "${shdiropt}"` || return 1
    read_opt_value_from_line_memoiz "$cmdline" "${shdiropt}" || return 1
    local shdir=${_OPT_VALUE_}

    get_absolute_shdirname $outd $shdir
}

########
save_opt_list()
{
    local optlist_varname=$1
    SCRIPT_OPT_LIST=${!optlist_varname}
}
