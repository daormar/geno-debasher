# *- bash -*

# INCLUDE BASH LIBRARY
. ${bindir}/bam_utils_lib

########
print_desc()
{
    echo "submit_bam_analysis performs analyses given normal and tumor bam files"
    echo "type \"submit_bam_analysis --help\" to get usage information"
}

########
usage()
{
    echo "submit_bam_analysis  -r <string>"
    echo "                     -n <string>|-extn <string>"
    echo "                     -t <string>|-extt <string>"
    echo "                     -g <string> -a <string> -o <string>"
    echo "                     [-nt <int>] [-cr <string>] [-wcr <string>]"
    echo "                     [-sv <string>] [-sg <string>] [-mc <string>]"
    echo "                     [-egastr <int>] [-egacred <string>]"
    echo "                     [-asperausr <string>] [-asperapwd <string>]"
    echo "                     [-asperaserv <string>] [-egadecrpwd <string>]"
    echo "                     [-debug] [--help]"
    echo ""
    echo "-r <string>          File with reference genome"
    echo "-n <string>          Normal bam file"
    echo "-t <string>          Tumor bam file"
    echo "-extn <string>       External database id of normal bam file to download"
    echo "-extt <string>       External database id of tumor bam file to download"
    echo "-g <string>          Sample gender (XX|XY)"
    echo "-a <string>          File with analysis steps to be performed."
    echo "                     Expected format:"
    echo "                     <stepname> <account> <partition> <cpus> <mem> <time> <jobdeps=stepname1:...>"
    echo "-o <string>          Output directory"
    echo "-nt <int>            Number of download tries per file"
    echo "-cr <string>         bgzipped and tabixed bed file to specify regions to call for"
    echo "                     Manta and Strelka"
    echo "-wcr <string>        Reference file in npz format for WisecondorX"
    echo "-sv <string>         SNP vcf file required by Facets"
    echo "-sg <string>         SNP GC correction file required by AscatNGS"
    echo "-mc <string>         Name of male sex chromosome required by AscatNGS"
    echo "-egastr <int>        Number of streams used by the EGA download client"
    echo "                     (50 by default)"
    echo "-egacred <string>    File with EGA download client credentials"
    echo "-asperausr <string>  Username for Aspera server"
    echo "-asperapwd <string>  Password for Aspera server"
    echo "-asperaserv <string> Name of Aspera server"
    echo "-egadecrpwd <string> File with EGA decryptor password"
    echo "-debug               After ending, do not delete temporary files"
    echo "                     (for debugging purposes)"
    echo "--help               Display this help and exit"
}

########
read_pars()
{
    r_given=0
    n_given=0
    t_given=0
    extn_given=0
    extt_given=0
    a_given=0
    g_given=0
    gender="XX"
    o_given=0
    download_tries=5
    cr_given=0
    callregf="NONE"
    wcr_given=0
    wcref="NONE"
    sv_given=0
    snpvcf="NONE"
    sg_given=0
    snpgccorr="NONE"
    mc_given=0
    malesexchr="Y"
    egastr_given=0
    egastr=50
    egacred_given=0
    egacred="cred.json"
    asperausr_given=0
    asperausr="NONE"
    asperapwd_given=0
    asperapwd="NONE"
    asperaserv_given=0
    asperaserv="NONE"
    egadecrpwd_given=0
    egadecrpwd="NONE"
    debug=0
    while [ $# -ne 0 ]; do
        case $1 in
            "--help") usage
                      exit 1
                      ;;
            "--version") version
                         exit 1
                         ;;
            "-r") shift
                  if [ $# -ne 0 ]; then
                      ref=$1
                      r_given=1
                  fi
                  ;;
            "-n") shift
                  if [ $# -ne 0 ]; then
                      normalbam=$1
                      n_given=1
                  fi
                  ;;
            "-t") shift
                  if [ $# -ne 0 ]; then
                      tumorbam=$1
                      t_given=1
                  fi
                  ;;
            "-extn") shift
                  if [ $# -ne 0 ]; then
                      extid_normalbam=$1
                      extn_given=1
                  fi
                  ;;
            "-extt") shift
                  if [ $# -ne 0 ]; then
                      extid_tumorbam=$1
                      extt_given=1
                  fi
                  ;;
            "-g") shift
                  if [ $# -ne 0 ]; then
                      gender=$1
                      g_given=1
                  fi
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
            "-nt") shift
                  if [ $# -ne 0 ]; then
                      download_tries=$1
                      nt_given=1
                  fi
                  ;;
            "-cr") shift
                  if [ $# -ne 0 ]; then
                      callregf=$1
                      cr_given=1
                  fi
                  ;;
            "-wcr") shift
                  if [ $# -ne 0 ]; then
                      wcref=$1
                      wcr_given=1
                  fi
                  ;;
            "-sv") shift
                  if [ $# -ne 0 ]; then
                      snpvcf=$1
                      sv_given=1
                  fi
                  ;;
            "-sg") shift
                  if [ $# -ne 0 ]; then
                      snpgccorr=$1
                      sg_given=1
                  fi
                  ;;
            "-mc") shift
                  if [ $# -ne 0 ]; then
                      malesexchr=$1
                      mc_given=1
                  fi
                  ;;
            "-egastr") shift
                  if [ $# -ne 0 ]; then
                      egastr=$1
                      egastr_given=1
                  fi
                  ;;
            "-egacred") shift
                  if [ $# -ne 0 ]; then
                      egacred=$1
                      egacred_given=1
                  fi
                  ;;
            "-asperausr") shift
                  if [ $# -ne 0 ]; then
                      asperausr=$1
                      asperausr_given=1
                  fi
                  ;;
            "-asperapwd") shift
                  if [ $# -ne 0 ]; then
                      asperapwd=$1
                      asperapwd_given=1
                  fi
                  ;;
            "-asperaserv") shift
                  if [ $# -ne 0 ]; then
                      asperaserv=$1
                      asperaserv_given=1
                  fi
                  ;;
            "-egadecrpwd") shift
                  if [ $# -ne 0 ]; then
                      egadecrpwd=$1
                      egadecrpwd_given=1
                  fi
                  ;;
            "-debug") debug=1
                      debug_opt="-debug"
                      ;;
        esac
        shift
    done   
}

########
check_pars()
{
    if [ ${r_given} -eq 0 ]; then   
        echo "Error! -r parameter not given!" >&2
        exit 1
    else
        if [ ! -f ${ref} ]; then
            echo "Error! file ${ref} does not exist" >&2
            exit 1
        fi
    fi

    if [ ${n_given} -eq 0 -a ${extn_given} -eq 0 ]; then
        echo "Error, -n or -extn options should be given" >&2
        exit 1
    fi

    if [ ${n_given} -eq 1 -a ${extn_given} -eq 1 ]; then
        echo "Error, -n and -extn options cannot be given simultaneously" >&2
        exit 1
    fi

    if [ ${t_given} -eq 0 -a ${extt_given} -eq 0 ]; then
        echo "Error, -t or -extt options should be given" >&2
        exit 1
    fi

    if [ ${t_given} -eq 1 -a ${extt_given} -eq 1 ]; then
        echo "Error, -t and -extt options cannot be given simultaneously" >&2
        exit 1
    fi

    if [ ${n_given} -eq 1 ]; then
        if [ ! -f ${normalbam} ]; then
            echo "Error! file ${normalbam} does not exist" >&2
            exit 1
        fi
    fi

    if [ ${t_given} -eq 1 ]; then
        if [ ! -f ${tumorbam} ]; then
            echo "Error! file ${tumorbam} does not exist" >&2
            exit 1
        fi
    fi

    if [ ${g_given} -eq 0 ]; then   
        echo "Error! -g parameter not given!" >&2
        exit 1
    fi

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
print_pars()
{
    if [ ${r_given} -eq 1 ]; then
        echo "-r is ${ref}" >&2
    fi

    if [ ${n_given} -eq 1 ]; then
        echo "-n is ${normalbam}" >&2
    fi

    if [ ${t_given} -eq 1 ]; then
        echo "-t is ${tumorbam}" >&2
    fi

    if [ ${extn_given} -eq 1 ]; then
        echo "-extn is ${extid_normalbam}" >&2
    fi

    if [ ${extt_given} -eq 1 ]; then
        echo "-extt is ${extid_tumorbam}" >&2
    fi

    if [ ${a_given} -eq 1 ]; then
        echo "-a is ${afile}" >&2
    fi

    if [ ${g_given} -eq 1 ]; then
        echo "-g is ${gender}" >&2
    fi

    if [ ${o_given} -eq 1 ]; then
        echo "-o is ${outd}" >&2
    fi

    if [ ${sg_given} -eq 1 ]; then
        echo "-sg is ${snpgccorr}" >&2
    fi

    if [ ${mc_given} -eq 1 ]; then
        echo "-mc is ${malesexchr}" >&2
    fi

    if [ ${cr_given} -eq 1 ]; then
        echo "-cr is ${callregf}" >&2
    fi

    if [ ${wcr_given} -eq 1 ]; then
        echo "-wcr is ${wcref}" >&2
    fi

    if [ ${egastr_given} -eq 1 ]; then
        echo "-egastr is ${egastr}" >&2
    fi

    if [ ${egacred_given} -eq 1 ]; then
        echo "-egacred is ${egacred}" >&2
    fi

    if [ ${asperausr_given} -eq 1 ]; then
        echo "-asperausr is ${asperausr}" >&2
    fi

    if [ ${asperapwd_given} -eq 1 ]; then
        echo "-asperapwd is ${asperapwd}" >&2
    fi

    if [ ${asperaserv_given} -eq 1 ]; then
        echo "-asperaserv is ${asperaserv}" >&2
    fi

    if [ ${egadecrpwd_given} -eq 1 ]; then
        echo "-egadecrpwd is ${egadecrpwd}" >&2
    fi
}

########
create_dirs()
{
    mkdir -p ${outd} || { echo "Error! cannot create output directory" >&2; return 1; }

    mkdir -p ${outd}/scripts || { echo "Error! cannot create scripts directory" >&2; return 1; }

    mkdir -p ${outd}/data || { echo "Error! cannot create data directory" >&2; return 1; }
}

########
set_bam_filenames()
{
    if [ ${extn_given} -eq 1 ]; then
        normalbam=${outd}/data/normal.bam
    fi

    if [ ${extt_given} -eq 1 ]; then
        tumorbam=${outd}/data/tumor.bam
    fi
}

########
get_pars_manta_germline()
{
    echo "$ref $normalbam ${callregf} ${step_outd} $cpus"
}

########
get_pars_manta_somatic()
{
    echo "$ref $normalbam $tumorbam ${callregf} ${step_outd} $cpus"
}

########
get_pars_strelka_germline()
{
    echo "$ref $normalbam ${callregf} ${step_outd} $cpus"
}

########
get_pars_strelka_somatic()
{
    local manta_dep=`find_dependency_for_step ${jobdeps_spec} manta_somatic`
    local manta_outd=`get_outd_for_dep ${outd} "${manta_dep}"`
    echo "$ref $normalbam $tumorbam ${callregf} ${step_outd} "${manta_outd}" $cpus"
}

########
get_pars_msisensor()
{
    echo "$ref $normalbam $tumorbam ${step_outd} $cpus"
}

########
get_pars_platypus_germline()
{
    echo "$ref $normalbam ${step_outd}"
}

########
get_pars_cnvkit()
{
    echo "$ref $normalbam $tumorbam ${step_outd} $cpus"
}

########
get_pars_wisecondorx()
{
    echo "$wcref $tumorbam ${step_outd} $cpus"
}

########
get_pars_facets()
{
    echo "$normalbam $tumorbam $snpvcf ${step_outd}"
}

########
get_pars_ascatngs()
{
    echo "$ref $normalbam $tumorbam $gender $malesexchr $snpgccorr ${step_outd} $cpus"
}

########
get_pars_download_ega_norm_bam()
{
    echo "$normalbam ${extid_normalbam} $egastr $egacred ${download_tries} ${step_outd}"
}

########
get_pars_download_ega_tum_bam()
{
    echo "$tumorbam ${extid_tumorbam} $egastr $egacred ${download_tries} ${step_outd}"
}

########
get_pars_download_aws_norm_bam()
{
    echo "$normalbam ${extid_normalbam} ${download_tries} ${step_outd}"
}

########
get_pars_download_aws_tum_bam()
{
    echo "$tumorbam ${extid_tumorbam} ${download_tries} ${step_outd}"
}

########
get_pars_download_collab_norm_bam()
{
    echo "$normalbam ${extid_normalbam} ${download_tries} ${step_outd}"
}

########
get_pars_download_collab_tum_bam()
{
    echo "$tumorbam ${extid_tumorbam} ${download_tries} ${step_outd}"
}

########
get_pars_download_ega_asp_norm_bam()
{
    echo "$normalbam ${extid_normalbam} ${asperausr} ${asperapwd} ${asperaserv} ${egadecrpwd} ${download_tries} ${step_outd}"
}

########
get_pars_download_ega_asp_tum_bam()
{
    echo "$tumorbam ${extid_tumorbam} ${asperausr} ${asperapwd} ${asperaserv} ${egadecrpwd} ${download_tries} ${step_outd}"
}

########
get_pars_index_norm_bam()
{
    echo "$normalbam ${step_outd}"
}

########
get_pars_index_tum_bam()
{
    echo "$tumorbam ${step_outd}"
}

########
get_pars_sort_norm_bam()
{
    echo "$normalbam ${step_outd} $cpus"
}

########
get_pars_sort_tum_bam()
{
    echo "$tumorbam ${step_outd} $cpus"
}

########
get_pars_filter_norm_bam_contigs()
{
    echo "$ref $normalbam ${step_outd}"
}

########
get_pars_filter_tum_bam_contigs()
{
    echo "$ref $tumorbam ${step_outd}"
}

########
get_pars_delete_bam_files()
{
    echo "$normalbam $tumorbam ${step_outd}"
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

########
check_script_is_older_than_lib()
{
    local script_filename=$1
    # Check if script exists
    if [ -f ${script_filename} ]; then
        # script exists
        lib_timestamp=`get_lib_timestamp`
        script_timestamp=`get_file_timestamp ${script_filename}`
        if [ ${script_timestamp} -lt ${lib_timestamp} ]; then
            echo 1
        else
            echo 0
        fi
    else
        # script does not exist
        echo 0
    fi
}

########
execute_step()
{
    # Initialize variables
    local dirname=$1
    local stepname=$2
    local account=$3
    local partition=$4
    local cpus=$5
    local mem=$6
    local time=$7
    local jobdeps_spec=$8
    step_outd=`get_step_dirname ${outd} ${stepname}`
    
    # Execute step

    ## Obtain step status
    local status=`${bindir}/get_analysis_status -d ${dirname} -s "${stepname}"`
    echo "STEP: ${stepname} ; STATUS: ${status}" >&2

    ## Decide whether step should be executed
    if [ "${status}" != "FINISHED" ]; then
        # Initialize script variables
        local script_filename=`get_script_filename ${stepname}`
        local step_function=`get_step_function ${stepname}`
        local script_pars_funcname=`get_script_pars_funcname ${stepname}`
        local script_pars=`${script_pars_funcname}`
                
        # Create script
        create_script ${script_filename} ${step_function} "${script_pars}"

        # Archive script
        archive_script ${script_filename}

        # Execute script
        reset_outdir_for_step ${dirname} ${stepname} || return 1
        local jobdeps="`get_jobdeps ${jobdeps_spec}`"
        local stepname_jid=${stepname}_jid
        launch ${script_filename} ${account} ${partition} ${cpus} ${mem} ${time} "${jobdeps}" ${stepname_jid} || return 1
        
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
execute_steps_in_afile()
{
    # Read input parameters
    local dirname=$1
    local afile=$2

    # step_jids will store the job ids of the analysis steps
    step_jids=""
    
    # Read information about the steps to be executed
    while read entry; do
        entry_ok=`analysis_entry_is_ok "$entry"`
        if [ ${entry_ok} = "yes" ]; then
            # Extract entry information
            local stepname=`extract_stepname_from_entry "$entry"`
            local account=`extract_account_from_entry "$entry"`
            local partition=`extract_partition_from_entry "$entry"`
            cpus=`extract_cpus_from_entry "$entry"`
            mem=`extract_mem_from_entry "$entry"`
            local time=`extract_time_from_entry "$entry"`
            jobdeps_spec=`extract_jobdeps_spec_from_entry "$entry"`

            # Execute step
            execute_step ${dirname} ${stepname} ${account} ${partition} ${cpus} ${mem} ${time} ${jobdeps_spec} || return 1
        fi
    done < ${afile}
}

########

if [ $# -eq 0 ]; then
    print_desc
    exit 1
fi

read_pars $@ || exit 1

check_pars || exit 1

print_pars || exit 1

create_dirs || exit 1

set_bam_filenames || exit 1

execute_steps_in_afile ${outd} ${afile} || exit 1
