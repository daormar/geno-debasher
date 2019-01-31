# *- bash -*

# INCLUDE BASH LIBRARY
. ${PANPIPE_HOME_DIR}/bin/panpipe_lib || exit 1

########
print_desc()
{
    echo "analyze_dataset analyses samples of a given dataset"
    echo "type \"analyze_dataset --help\" to get usage information"
}

########
usage()
{
    echo "analyze_dataset      -r <string> -m <string>"
    echo "                     -p <string> -o <string>"
    echo "                     [-wcr <string>] [-sv <string>]"
    echo "                     [-sg <string>] [-mc <string>]"
    echo "                     [-egastr <int>] [-egacred <string>]"
    echo "                     [-asperausr <string>] [-asperapwd <string>]"
    echo "                     [-asperaserv <string>] [-egadecrpwd <string>]"
    echo "                     [-p] [--help]"
    echo ""
    echo "-r <string>          File with reference genome"
    echo "-m <string>          File with metadata, one entry per line."
    echo "                     Format: ID PHENOTYPE GENDER ; ID PHENOTYPE GENDER"
    echo "-p <string>          File with pipeline steps to be performed"
    echo "-o <string>          Output directory"
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
    echo "--help               Display this help and exit"
}

########
read_pars()
{
    r_given=0
    m_given=0
    p_given=0
    o_given=0
    cr_given=0
    callregf=${NOFILE}
    wcr_given=0
    wcref=${NOFILE}
    sv_given=0
    snpvcf=${NOFILE}
    sg_given=0
    snpgccorr=${NOFILE}
    mc_given=0
    malesexchr="Y"
    egastr_given=0
    egastr=50
    egacred_given=0
    egacred=${NOFILE}
    asperausr_given=0
    asperausr=${NOFILE}
    asperapwd_given=0
    asperapwd=${NOFILE}
    asperaserv_given=0
    asperaserv=${NOFILE}
    egadecrpwd_given=0
    egadecrpwd=${NOFILE}
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
            "-m") shift
                  if [ $# -ne 0 ]; then
                      metadata=$1
                      m_given=1
                  fi
                  ;;
            "-p") shift
                  if [ $# -ne 0 ]; then
                      pfile=$1
                      p_given=1
                  fi
                  ;;
            "-o") shift
                  if [ $# -ne 0 ]; then
                      outd=$1
                      o_given=1
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
    
    if [ ${m_given} -eq 0 ]; then
        echo "Error, -m option should be given" >&2
    fi

    if [ ${m_given} -eq 1 ]; then
        if [ ! -f ${metadata} ]; then
            echo "Error! file ${metadata} does not exist" >&2
            exit 1
        fi
    fi

    if [ ${p_given} -eq 0 ]; then   
        echo "Error! -a parameter not given!" >&2
        exit 1
    else
        if [ ! -f ${pfile} ]; then
            echo "Error! file ${pfile} does not exist" >&2
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

    if [ ${cr_given} -eq 1 ]; then
        if [ "${callregf}" != ${NOFILE} -a ! -f ${callregf} ]; then
            echo "Error! file ${callregf} does not exist" >&2
            exit 1
        fi
    fi

    if [ ${wcr_given} -eq 1 ]; then
        if [ "${wcref}" != ${NOFILE} -a ! -f ${wcref} ]; then
            echo "Error! file ${wcref} does not exist" >&2
            exit 1
        fi
    fi

    if [ ${sv_given} -eq 1 ]; then
        if [ "${snpvcf}" != ${NOFILE} -a ! -f ${snpvcf} ]; then
            echo "Error! file ${snpvcf} does not exist" >&2
            exit 1
        fi
    fi

    if [ ${sg_given} -eq 1 ]; then
        if [ "${snpgccorr}" != ${NOFILE} -a ! -f ${snpgccorr} ]; then
            echo "Error! file ${snpgccorr} does not exist" >&2
            exit 1
        fi
    fi

    if [ ${egacred_given} -eq 1 ]; then
        if [ "${egacred}" != ${NOFILE} -a ! -f ${egacred} ]; then
            echo "Error! file ${egacred} does not exist" >&2
            exit 1
        fi
    fi

    if [ ${egadecrpwd_given} -eq 1 ]; then
        if [ "${egadecrpwd}" != ${NOFILE} -a ! -f ${egadecrpwd} ]; then
            echo "Error! file ${egadecrpwd} does not exist" >&2
            exit 1
        fi
    fi
}

########
absolutize_file_paths()
{
    if [ ${r_given} -eq 1 ]; then   
        ref=`get_absolute_path ${ref}`
    fi
    
    if [ ${m_given} -eq 1 ]; then
        metadata=`get_absolute_path ${metadata}`
    fi

    if [ ${p_given} -eq 1 ]; then   
        pfile=`get_absolute_path ${pfile}`
    fi

    if [ ${o_given} -eq 1 ]; then
        outd=`get_absolute_path ${outd}`
    fi

    if [ ${cr_given} -eq 1 -a "${callrefg}" != ${NOFILE} ]; then
        callregf=`get_absolute_path ${callregf}`
    fi

    if [ ${wcr_given} -eq 1 -a "${wcref}" != ${NOFILE} ]; then
        wcref=`get_absolute_path ${wcref}`
    fi

    if [ ${sv_given} -eq 1 -a "${snpvcf}" != ${NOFILE} ]; then
        snpvcf=`get_absolute_path ${snpvcf}`
    fi

    if [ ${sg_given} -eq 1 -a "${snpgccorr}" != ${NOFILE} ]; then
        snpgccorr=`get_absolute_path ${snpgccorr}`
    fi

    if [ ${egacred_given} -eq 1 -a "${egacred}" != ${NOFILE} ]; then
        egacred=`get_absolute_path ${egacred}`
    fi

    if [ ${egadecrpwd_given} -eq 1 -a "${egadecrpwd}" != ${NOFILE} ]; then
        egadecrpwd=`get_absolute_path ${egadecrpwd}`
    fi
}

########
print_pars()
{
    if [ ${r_given} -eq 1 ]; then
        echo "-r is ${ref}" >&2
    fi

    if [ ${m_given} -eq 1 ]; then
        echo "-m is ${metadata}" >&2
    fi

    if [ ${p_given} -eq 1 ]; then
        echo "-p is ${pfile}" >&2
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
extract_normal_sample_info()
{
    local entry=$1
    local sample1=`echo ${entry} | $AWK -F ";" '{print $1}' | $GREP 'Normal\|normal'`
    local sample2=`echo ${entry} | $AWK -F ";" '{print $2}' | $GREP 'Normal\|normal'`

    if [ ! -z "${sample1}" ]; then
        echo ${sample1}
    else
        if [ ! -z "${sample2}" ]; then
            echo ${sample2}
        else
            echo ""
        fi
    fi
}

########
extract_tumor_sample_info()
{
    local entry=$1
    local sample1=`echo ${entry} | $AWK -F ";" '{print $1}' | $GREP 'Tumour\|tumour\|Tumor\|tumor'`
    local sample2=`echo ${entry} | $AWK -F ";" '{print $2}' | $GREP 'Tumour\|tumour\|Tumor\|tumor'`

    if [ ! -z "${sample1}" ]; then
        echo ${sample1}
    else
        if [ ! -z "${sample2}" ]; then
            echo ${sample2}
        else
            echo ""
        fi
    fi
}

########
entry_is_ok()
{
    local entry=$1
    local nsample=`extract_normal_sample_info "${entry}"`
    local tsample=`extract_tumor_sample_info "${entry}"`

    if [ ! -z "${nsample}" -a ! -z "${tsample}" ]; then
        echo "yes"
    else
        echo "no"
    fi
}

########
extract_id_from_sample_info()
{
    local sample_info=$1
    echo ${sample_info} | $AWK '{print $1}'
}

########
extract_gender_from_sample_info()
{
    local sample_info=$1
    local tmp=`echo ${sample_info} | $GREP 'Female\|female'`
    if [ ! -z "${tmp}" ]; then
        echo "female"
    else
        echo "male"
    fi
}

########
get_outd_name()
{
    local norm_id=$1
    local tum_id=$2

    # If id contains a file path, retain file name only
    local norm_id_wo_pathinfo=`$BASENAME ${norm_id}`
    local tum_id_wo_pathinfo=`$BASENAME ${tum_id}`
    
    echo ${norm_id_wo_pathinfo}"_"${tum_id_wo_pathinfo}
}

########
process_pars()
{
    # Read metadata file
    entry_num=1
    while read entry; do
        entry_ok=`entry_is_ok "$entry"`
        if [ ${entry_ok} = "yes" ]; then
            # Extract sample info
            normal_sample_info=`extract_normal_sample_info "$entry"`
            normal_id=`extract_id_from_sample_info "${normal_sample_info}"`
            
            tumor_sample_info=`extract_tumor_sample_info "$entry"`
            tumor_id=`extract_id_from_sample_info "${tumor_sample_info}"`

            gender=`extract_gender_from_sample_info "${normal_sample_info}"`

            # Obtain value for -g option
            if [ ${gender} = "male" ]; then
                gender_opt="XY"
            else                
                gender_opt="XX"
            fi
            
            # Set name of output directory for analysis
            analysis_outd=`get_outd_name ${normal_id} ${tumor_id}`
            
            # Print command to execute pipeline
            echo ${bindir}/pipe_exec -r ${ref} -extn ${normal_id} -extt ${tumor_id} -p ${pfile} -g ${gender_opt} -o ${outd}/${analysis_outd} -cr ${callregf} -wcr ${wcref} -sv ${snpvcf} -sg ${snpgccorr} -mc ${malesexchr} -egastr ${egastr} -egacred ${egacred} -asperausr ${asperausr} -asperapwd ${asperapwd} -asperaserv ${asperaserv} -egadecrpwd ${egadecrpwd}
        else
            echo "Error in entry number ${entry_num}"
        fi

        entry_num=`expr ${entry_num} + 1`
        
    done < ${metadata}
}

########

if [ $# -eq 0 ]; then
    print_desc
    exit 1
fi

read_pars $@ || exit 1

check_pars || exit 1

absolutize_file_paths || exit 1

print_pars || exit 1

process_pars || exit 1
