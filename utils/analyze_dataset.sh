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
    echo "analyze_dataset       --pfile <string> -outdir <string>"
    echo "                      --sched <string> --metadata <string>"
    echo "                      [--dflt-nodes <string>] --ppl-opts <string>"
    echo "                      [--help]"
    echo ""
    echo "--pfile <string>      File with pipeline steps to be performed"
    echo "--outdir <string>     Output directory"
    echo "--sched <string>      Scheduler used to execute the pipelines"
    echo "--metadata <string>   File with metadata, one entry per line."
    echo "                      Format: ID PHENOTYPE GENDER ; ID PHENOTYPE GENDER"
    echo "--dflt-nodes <string> Default set of nodes used to execute the pipeline"
    echo "--ppl-opts <string>   File containing a string with pipeline options"
    echo "--help                Display this help and exit"
}

########
read_pars()
{
    pfile_given=0
    outdir_given=0
    sched_given=0
    metadata_given=0
    dflt_nodes_given=0
    ppl_opts_given=0
    while [ $# -ne 0 ]; do
        case $1 in
            "--help") usage
                      exit 1
                      ;;
            "--pfile") shift
                  if [ $# -ne 0 ]; then
                      pfile=$1
                      pfile_given=1
                  fi
                  ;;
            "--outdir") shift
                  if [ $# -ne 0 ]; then
                      outd=$1
                      outdir_given=1
                  fi
                  ;;
            "--sched") shift
                  if [ $# -ne 0 ]; then
                      sched=$1
                      sched_given=1
                  fi
                  ;;
            "--metadata") shift
                  if [ $# -ne 0 ]; then
                      metadata=$1
                      metadata_given=1
                  fi
                  ;;
            "--dflt-nodes") shift
                  if [ $# -ne 0 ]; then
                      dflt_nodes=$1
                      dflt_nodes_given=1
                  fi
                  ;;
            "--ppl-opts") shift
                  if [ $# -ne 0 ]; then
                      ppl_opts=$1
                      ppl_opts_given=1
                  fi
                  ;;
        esac
        shift
    done   
}

########
check_pars()
{
    if [ ${pfile_given} -eq 0 ]; then   
        echo "Error! --pfile parameter not given!" >&2
        exit 1
    else
        if [ ! -f ${pfile} ]; then
            echo "Error! file ${pfile} does not exist" >&2
            exit 1
        fi
    fi

    if [ ${outdir_given} -eq 0 ]; then
        echo "Error! --outdir parameter not given!" >&2
        exit 1
    else
        if [ -d ${outd} ]; then
            echo "Warning! output directory does exist" >&2 
        fi
    fi

    if [ ${sched_given} -eq 0 ]; then
        echo "Error, --sched option should be given" >&2
        exit 1
    fi
    
    if [ ${metadata_given} -eq 0 ]; then
        echo "Error, --metadata option should be given" >&2
        exit 1
    fi

    if [ ${metadata_given} -eq 1 ]; then
        if [ ! -f ${metadata} ]; then
            echo "Error! file ${metadata} does not exist" >&2
            exit 1
        fi
    fi

    if [ ${ppl_opts_given} -eq 0 ]; then
        echo "Error, --ppl-opts option should be given" >&2
        exit 1
    fi

    if [ ${ppl_opts_given} -eq 1 ]; then
        if [ ! -f ${ppl_opts} ]; then
            echo "Error! file ${ppl_opts} does not exist" >&2
            exit 1
        fi
    fi
}

########
absolutize_file_paths()
{
    if [ ${pfile_given} -eq 1 ]; then   
        pfile=`get_absolute_path ${pfile}`
    fi

    if [ ${outdir_given} -eq 1 ]; then
        outd=`get_absolute_path ${outd}`
    fi

    if [ ${metadata_given} -eq 1 ]; then
        metadata=`get_absolute_path ${metadata}`
    fi

    if [ ${ppl_opts_given} -eq 1 ]; then
        ppl_opts=`get_absolute_path ${ppl_opts}`
    fi
}

########
print_pars()
{
    if [ ${pfile_given} -eq 1 ]; then
        echo "--pfile is ${pfile}" >&2
    fi

    if [ ${outdir_given} -eq 1 ]; then
        echo "--outdir is ${outd}" >&2
    fi

    if [ ${sched_given} -eq 1 ]; then
        echo "--sched is ${sched}" >&2
    fi

    if [ ${metadata_given} -eq 1 ]; then
        echo "--metadata is ${metadata}" >&2
    fi

    if [ ${dflt_nodes_given} -eq 1 ]; then
        echo "--dflt-nodes is ${dflt_nodes}" >&2
    fi

    if [ ${ppl_opts_given} -eq 1 ]; then
        echo "--ppl-opts is ${ppl_opts}" >&2
    fi
}

########
check_sample_is_normal()
{
    $GREP 'Normal\|normal\|Non-tumor\|non-tumor'
}

########
extract_normal_sample_info()
{
    local entry=$1
    local sample1=`echo ${entry} | $AWK -F ";" '{print $1}' | check_sample_is_normal`
    local sample2=`echo ${entry} | $AWK -F ";" '{print $2}' | check_sample_is_normal`

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
check_sample_is_tumor()
{
    $GREP 'Tumour\|tumour\|Tumor\|tumor' | $GREP -v 'Non-tumor\|non-tumor'
}

########
extract_tumor_sample_info()
{
    local entry=$1
    local sample1=`echo ${entry} | $AWK -F ";" '{print $1}' | check_sample_is_tumor`
    local sample2=`echo ${entry} | $AWK -F ";" '{print $2}' | check_sample_is_tumor`

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
get_ppl_opts_str()
{
    cat ${ppl_opts}
}

########
get_dflt_nodes_opt()
{
    if [ ${dflt_nodes_given} -eq 1 ]; then
        echo "--dflt-nodes ${dflt_nodes}"
    else
        echo ""
    fi
}

########
process_pars()
{
    # Set options
    ppl_opts_str=`get_ppl_opts_str`
    
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

            # Obtain --dflt-nodes option
            dflt_nodes_opt=`get_dflt_nodes_opt`
            
            # Print command to execute pipeline
            echo ${PANPIPE_HOME_DIR}/bin/pipe_exec --pfile ${pfile} --outdir ${outd}/${analysis_outd} --sched ${sched} ${dflt_nodes_opt} -extn ${normal_id} -extt ${tumor_id} -g ${gender_opt} ${ppl_opts_str}
        else
            echo "Error in entry number ${entry_num}"
        fi

        entry_num=$((entry_num + 1))
        
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
