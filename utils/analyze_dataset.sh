# Geno-PanPipe package
# Copyright (C) 2019,2020 Daniel Ortiz-Mart\'inez
#
# This library is free software; you can redistribute it and/or
# modify it under the terms of the GNU Lesser General Public License
# as published by the Free Software Foundation; either version 3
# of the License, or (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with this program; If not, see <http://www.gnu.org/licenses/>.

# *- bash -*

# INCLUDE BASH LIBRARY
. "${PANPIPE_HOME_DIR}"/bin/panpipe_lib || exit 1

########
print_desc()
{
    echo "analyze_dataset analyses samples of a given dataset"
    echo "type \"analyze_dataset --help\" to get usage information"
}

########
usage()
{
    echo "analyze_dataset       --pfile <string> --outdir <string>"
    echo "                      --sched <string> --metadata <string>"
    echo "                      [--dflt-nodes <string>] --ppl-opts <string>"
    echo "                      [--lcxx <string> --lcxy <string>]"
    echo "                      [--local-bams] [--help]"
    echo ""
    echo "--pfile <string>      File with pipeline processes to be performed"
    echo "--outdir <string>     Output directory"
    echo "--sched <string>      Scheduler used to execute the pipelines"
    echo "--metadata <string>   File with metadata, one entry per line."
    echo "                      Format: ID PHENOTYPE GENDER ; ID PHENOTYPE GENDER"
    echo "--dflt-nodes <string> Default set of nodes used to execute the pipeline"
    echo "--ppl-opts <string>   File containing a string with pipeline options"
    echo "--lcxx <string>       File containing list of contigs of interest for XX samples"
    echo "--lcxy <string>       File containing list of contigs of interest for XY samples"
    echo "--local-bams          ID's for bam files contained in metadata are considered"
    echo "                      as local file names"
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
    lcxx_given=0
    lcxy_given=0
    local_bams_given=0
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
            "--lcxx") shift
                  if [ $# -ne 0 ]; then
                      lcxx=$1
                      lcxx_given=1
                  fi
                  ;;
            "--lcxy") shift
                  if [ $# -ne 0 ]; then
                      lcxy=$1
                      lcxy_given=1
                  fi
                  ;;
            "--local-bams")
                  if [ $# -ne 0 ]; then
                      local_bams_given=1
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
        if [ ! -f "${pfile}" ]; then
            echo "Warning! file ${pfile} does not exist" >&2
        fi
    fi

    if [ ${outdir_given} -eq 0 ]; then
        echo "Error! --outdir parameter not given!" >&2
        exit 1
    else
        if [ -d "${outd}" ]; then
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
        if [ ! -f "${ppl_opts}" ]; then
            echo "Error! file ${ppl_opts} does not exist" >&2
            exit 1
        fi
    fi

    if [ ${lcxx_given} -eq 1 -a  ${lcxy_given} -eq 0 ]; then
        echo "Error, --lcxx and --lcxy should be given simultaneously" >&2
        exit 1
    fi

    if [ ${lcxx_given} -eq 0 -a  ${lcxy_given} -eq 1 ]; then
        echo "Error, --lcxx and --lcxy should be given simultaneously" >&2
        exit 1
    fi

    if [ ${lcxx_given} -eq 1 ]; then
        if [ ! -f ${lcxx} ]; then
            echo "Error! file ${lcxx} does not exist" >&2
            exit 1
        fi
    fi

    if [ ${lcxy_given} -eq 1 ]; then
        if [ ! -f ${lcxy} ]; then
            echo "Error! file ${lcxy} does not exist" >&2
            exit 1
        fi
    fi
}

########
absolutize_file_paths()
{
    if [ ${pfile_given} -eq 1 ]; then
        pfile=`get_absolute_path "${pfile}"`
    fi

    if [ ${outdir_given} -eq 1 ]; then
        outd=`get_absolute_path "${outd}"`
    fi

    if [ ${metadata_given} -eq 1 ]; then
        metadata=`get_absolute_path "${metadata}"`
    fi

    if [ ${ppl_opts_given} -eq 1 ]; then
        ppl_opts=`get_absolute_path "${ppl_opts}"`
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
filter_normal_sample()
{
    "$GREP" 'Normal\|normal\|Non-tumor\|non-tumor'
}

########
extract_normal_sample_info()
{
    local entry=$1
    local sample1=`echo ${entry} | "$AWK" -F ";" '{print $1}' | filter_normal_sample`
    local sample2=`echo ${entry} | "$AWK" -F ";" '{print $2}' | filter_normal_sample`

    if [ ! -z "${sample1}" ]; then
        echo "${sample1}"
    else
        if [ ! -z "${sample2}" ]; then
            echo "${sample2}"
        else
            echo ""
        fi
    fi
}

########
filter_tumor_sample()
{
    "$GREP" 'Tumour\|tumour\|Tumor\|tumor' | "$GREP" -v 'Non-tumor\|non-tumor'
}

########
extract_tumor_sample_info()
{
    local entry=$1
    local sample1=`echo ${entry} | "$AWK" -F ";" '{print $1}' | filter_tumor_sample`
    local sample2=`echo ${entry} | "$AWK" -F ";" '{print $2}' | filter_tumor_sample`

    if [ ! -z "${sample1}" ]; then
        echo "${sample1}"
    else
        if [ ! -z "${sample2}" ]; then
            echo "${sample2}"
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
    echo "${sample_info}" | "$AWK" '{print $1}'
}

########
extract_gender_from_sample_info()
{
    local sample_info=$1
    local tmp=`echo ${sample_info} | "$GREP" 'Female\|female'`
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
    local norm_id_wo_pathinfo=`"$BASENAME" "${norm_id}"`
    local tum_id_wo_pathinfo=`"$BASENAME" "${tum_id}"`

    echo "${norm_id_wo_pathinfo}_${tum_id_wo_pathinfo}"
}

########
get_ppl_opts_str()
{
    cat "${ppl_opts}"
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
get_lc_opt()
{
    gender_opt=$1
    lcxx=$2
    lcxy=$3

    if [ "${lcxx}" = "" -o "${lcxy}" = "" ]; then
        echo ""
    else
        if [ "${gender_opt}" = "XX" ]; then
            echo "-lc ${lcxx}"
        else
            echo "-lc ${lcxy}"
        fi
    fi
}

########
esc_dq()
{
    "$SED" 's/"/\\\"/g' <<< "$1"
}

########
process_pars()
{
    # Set options
    ppl_opts_str=`get_ppl_opts_str`

    # Get pipe_exec path
    local pipe_exec_path
    panpipe_exec_path=`get_panpipe_exec_path`

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

            # Obtain value for -lc option
            lc_opt=`get_lc_opt ${gender_opt} ${lcxx} ${lcxy}`

            # Set name of output directory for analysis
            analysis_outd=`get_outd_name "${normal_id}" "${tumor_id}"`

            # Obtain --dflt-nodes option
            dflt_nodes_opt=`get_dflt_nodes_opt`

            # Determine whether the normal and tumor ids correspond to
            # locally stored file names or not
            if [ ${local_bams_given} -eq 1 ]; then
                nopt="-n"
                topt="-t"
            else
                nopt="-extn"
                topt="-extt"
            fi

            # Print command to execute pipeline
            normalize_cmd "\"$(esc_dq "${panpipe_exec_path}")\" --pfile \"$(esc_dq "${pfile}")\" --outdir \"$(esc_dq "${outd}/${analysis_outd}")\" --sched ${sched} ${dflt_nodes_opt} ${nopt} \"$(esc_dq "${normal_id}")\" ${topt} \"$(esc_dq "${tumor_id}")\" -g ${gender_opt} ${lc_opt} ${ppl_opts_str}"
        else
            echo "Error in entry number ${entry_num}"
        fi

        entry_num=$((entry_num + 1))

    done < "${metadata}"
}

########

if [ $# -eq 0 ]; then
    print_desc
    exit 1
fi

read_pars "$@" || exit 1

check_pars || exit 1

absolutize_file_paths || exit 1

print_pars || exit 1

process_pars || exit 1
