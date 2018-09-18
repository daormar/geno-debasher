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
    echo "                     -n <string>|-egan <string> -t <string>|-egat <string>"
    echo "                     -g <string> -a <string> -o <string>"
    echo "                     [-wcr <string>] [-sv <string>]"
    echo "                     [-sg <string>] [-mc <string>]"
    echo "                     [-egastr <int>] [-egacred <string>]"
    echo "                     [-debug] [--help]"
    echo ""
    echo "-r <string>          File with reference genome"
    echo "-n <string>          Normal bam file"
    echo "-t <string>          Tumor bam file"
    echo "-egan <string>       EGA id of normal bam file to download"
    echo "-egat <string>       EGA id of tumor bam file to download"
    echo "-g <string>          Sample gender (XX|XY)"
    echo "-a <string>          File with analysis steps to be performed."
    echo "                     Expected format:"
    echo "                     <stepname> <partition> <cpus> <mem> <time> <jobdeps=stepname1:...>"
    echo "-o <string>          Output directory"
    echo "-wcr <string>        Reference file in npz format for WisecondorX"
    echo "-sv <string>         SNP vcf file required by Facets"
    echo "-sg <string>         SNP GC correction file required by AscatNGS"
    echo "-mc <string>         Name of male sex chromosome required by AscatNGS"
    echo "-egastr <int>        Number of streams used by the EGA download client"
    echo "                     (50 by default)"
    echo "-egacred <string>    File with EGA download client credentials"
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
    egan_given=0
    egat_given=0
    a_given=0
    g_given=0
    gender="XX"
    o_given=0
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
            "-egan") shift
                  if [ $# -ne 0 ]; then
                      egaid_normalbam=$1
                      egan_given=1
                  fi
                  ;;
            "-egat") shift
                  if [ $# -ne 0 ]; then
                      egaid_tumorbam=$1
                      egat_given=1
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

    if [ ${n_given} -eq 0 -a ${egan_given} -eq 0 ]; then
        echo "Error, -n or -egan options should be given" >&2
    fi

    if [ ${n_given} -eq 1 -a ${egan_given} -eq 1 ]; then
        echo "Error, -n and -egan options cannot be given simultaneously" >&2
    fi

    if [ ${t_given} -eq 0 -a ${egat_given} -eq 0 ]; then
        echo "Error, -t or -egat options should be given" >&2
    fi

    if [ ${t_given} -eq 1 -a ${egat_given} -eq 1 ]; then
        echo "Error, -t and -egat options cannot be given simultaneously" >&2
    fi

    if [ ${n_given} -eq 1 ]; then
        if [ ! -f ${normalbam} ]; then
            echo "Error! file ${normalbam} does not exist" >&2
            exit 1
        fi
    fi

    if [ ${t_given} -eq 1 ]; then
        if [ -a ! -f ${tumorbam} ]; then
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

    if [ ${egan_given} -eq 1 ]; then
        echo "-egan is ${egaid_normalbam}" >&2
    fi

    if [ ${egat_given} -eq 1 ]; then
        echo "-egat is ${egaid_tumorbam}" >&2
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

    if [ ${wcr_given} -eq 1 ]; then
        echo "-wcr is ${wcref}" >&2
    fi

    if [ ${egastr_given} -eq 1 ]; then
        echo "-egastr is ${egastr}" >&2
    fi

    if [ ${egacred_given} -eq 1 ]; then
        echo "-egacred is ${egacred}" >&2
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
    if [ ${egan_given} -eq 1 ]; then
        normalbam=${outd}/data/normal.bam
    fi

    if [ ${egat_given} -eq 1 ]; then
        tumorbam=${outd}/data/tumor.bam
    fi
}

########
get_pars_manta_somatic()
{
    echo "$ref $normalbam $tumorbam ${step_outd} $cpus"
}

########
get_pars_strelka_somatic()
{
#    manta_dep=`find_dependency_for_step ${local_jobdeps_spec} manta_somatic`
#    local_manta_outd=`get_outd_for_dep ${outd} "${manta_dep}"`
    echo "$ref $normalbam $tumorbam ${step_outd} "${local_manta_outd}" $cpus"
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
    echo "$normalbam $tumorbam $snpvcf outd$"
}

########
get_pars_ascatngs()
{
    echo "$ref $normalbam $tumorbam $gender $malesexchr $snpgccorr ${step_outd} $cpus"
}

########
get_pars_download_ega_norm_bam()
{
    echo "$normalbam ${egaid_normalbam} $egastr $egacred ${step_outd}"
}

########
get_pars_download_ega_tum_bam()
{
    echo "$tumorbam ${egaid_tumorbam} $egastr $egacred ${step_outd}"
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
    create_script ${local_dirname}/scripts/execute_${local_stepname} execute_${local_stepname} "${local_script_pars}"
    local_status=`${bindir}/get_analysis_status -d ${local_dirname} -s "${local_stepname}"`
    echo "STEP: ${local_stepname} ; STATUS: ${local_status}" >&2
    if [ "${local_status}" != "FINISHED" ]; then
        reset_outdir_for_step ${local_dirname} ${local_stepname} || exit 1
        local_jobdeps="`get_jobdeps ${local_jobdeps_spec}`"
        local_stepname_jid=${local_stepname}_jid
        launch ${local_dirname}/scripts/execute_${local_stepname} ${local_account} ${local_partition} ${local_cpus} ${local_mem} ${local_time} "${local_jobdeps}" ${local_stepname_jid}
        
        # Update variables storing jids
        step_jids="${step_jids}:${!local_stepname_jid}"
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
            execute_step ${local_dirname} ${local_stepname} ${local_account} ${local_partition} ${cpus} ${mem} ${local_time} ${jobdeps_spec} || exit 1
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
