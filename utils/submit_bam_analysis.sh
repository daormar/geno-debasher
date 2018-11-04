# *- bash -*

# INCLUDE BASH LIBRARY
. ${bindir}/bam_utils_lib.sh

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
get_pars_manta_somatic()
{
    echo "$ref $normalbam $tumorbam ${callregf} ${step_outd} $cpus"
}

########
get_pars_strelka_somatic()
{
    local_manta_dep=`find_dependency_for_step ${jobdeps_spec} manta_somatic`
    local_manta_outd=`get_outd_for_dep ${outd} "${local_manta_dep}"`
    echo "$ref $normalbam $tumorbam ${callregf} ${step_outd} "${local_manta_outd}" $cpus"
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
get_pars_delete_bam_files()
{
    echo "$normalbam $tumorbam ${step_outd}"
}

########
get_jobdeps_from_detailed_spec()
{
    local_jobdeps_spec=$1

    local_jdeps=""
    # Iterate over the elements of the job specification: type1:stepname1,...,typen:stepnamen
    for dep_spec in `echo ${local_jobdeps_spec} | $SED 's/,/ /g'`; do
        local_deptype=`echo ${dep_spec} | $AWK -F ":" '{print $1}'`
        local_step=`echo ${dep_spec} | $AWK -F ":" '{print $2}'`
        
        # Check if there is a jid for the step
        local_step_jid=${local_step}_jid
        if [ ! -z "${!local_step_jid}" ]; then
            if [ -z "${local_jdeps}" ]; then
                local_jdeps=${local_deptype}":"${!local_step_jid}
            else
                local_jdeps=${local_jdeps}","${local_deptype}":"${!local_step_jid}
            fi
        fi
    done
    echo ${local_jdeps}
}

########
get_jobdeps()
{
    local_jobdeps_spec=$1
    case ${local_jobdeps_spec} in
            "afterok:all") apply_deptype_to_jobids ${step_jids} afterok
                    ;;
            "none") echo ""
                    ;;
            *) get_jobdeps_from_detailed_spec ${local_jobdeps_spec}
               ;;
    esac
}

########
execute_step()
{
    # Initialize variables
    local_dirname=$1
    local_stepname=$2
    local_account=$3
    local_partition=$4
    local_cpus=$5
    local_mem=$6
    local_time=$7
    local_jobdeps_spec=$8
    step_outd=`get_step_dirname ${outd} ${local_stepname}`
    
    # Execute step
    local_script_pars=`get_pars_${local_stepname}`
    script_filename=${local_dirname}/scripts/execute_${local_stepname}
    create_script ${script_filename} execute_${local_stepname} "${local_script_pars}"
    script_modified=`check_script_was_modified ${script_filename}`
    local_status=`${bindir}/get_analysis_status -d ${local_dirname} -s "${local_stepname}"`
    echo "STEP: ${local_stepname} ; STATUS: ${local_status}" >&2
    if [ "${local_status}" != "FINISHED" ]; then
        reset_outdir_for_step ${local_dirname} ${local_stepname} || return 1
        local_jobdeps="`get_jobdeps ${local_jobdeps_spec}`"
        local_stepname_jid=${local_stepname}_jid
        launch ${script_filename} ${local_account} ${local_partition} ${local_cpus} ${local_mem} ${local_time} "${local_jobdeps}" ${local_stepname_jid} || return 1
        
        # Update variables storing jids
        step_jids="${step_jids}:${!local_stepname_jid}"
    else
        if [ ${script_modified} -eq 1 ]; then
            echo "Warning: script was changed for this step with respect to last execution. See changes in file ${script_filename}.diff">&2
        fi
    fi
}

########
execute_steps_in_afile()
{
    # Read input parameters
    local_dirname=$1
    local_afile=$2

    # step_jids will store the job ids of the analysis steps
    step_jids=""
    
    # Read information about the steps to be executed
    while read entry; do
        entry_ok=`analysis_entry_is_ok "$entry"`
        if [ ${entry_ok} = "yes" ]; then
            # Extract entry information
            local_stepname=`extract_stepname_from_entry "$entry"`
            local_account=`extract_account_from_entry "$entry"`
            local_partition=`extract_partition_from_entry "$entry"`
            cpus=`extract_cpus_from_entry "$entry"`
            mem=`extract_mem_from_entry "$entry"`
            local_time=`extract_time_from_entry "$entry"`
            jobdeps_spec=`extract_jobdeps_spec_from_entry "$entry"`

            # Execute step
            execute_step ${local_dirname} ${local_stepname} ${local_account} ${local_partition} ${cpus} ${mem} ${local_time} ${jobdeps_spec} || return 1
        fi
    done < ${local_afile}
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
