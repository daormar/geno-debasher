# *- bash -*

# INCLUDE BASH LIBRARY
. ${bindir}/bam_utils_lib.sh

########
print_desc()
{
    echo "analyze_ega_study analyses samples of ega study"
    echo "type \"analyze_ega_study --help\" to get usage information"
}

########
usage()
{
    echo "analyze_ega_study    -r <string> -e <string>"
    echo "                     -a <string> -o <string>"
    echo "                     [-wcr <string>] [-sv <string>]"
    echo "                     [-sg <string>] [-mc <string>]"
    echo "                     [-egastr <int>] [-egacred <string>]"
    echo "                     [-debug] [--help]"
    echo ""
    echo "-r <string>          File with reference genome"
    echo "-e <string>          File with processed EGA metadata using the query_ega_metadata"
    echo "                     tool (option -f 3)"
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
    e_given=0
    a_given=0
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
            "-e") shift
                  if [ $# -ne 0 ]; then
                      egadata=$1
                      e_given=1
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

    if [ ${e_given} -eq 0 ]; then   
        echo "Error! -e parameter not given!" >&2
        exit 1
    else
        if [ ! -f ${egadata} ]; then
            echo "Error! file ${egadata} does not exist" >&2
            exit 1
        fi
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

    if [ ${e_given} -eq 1 ]; then
        echo "-e is ${egadata}" >&2
    fi

    if [ ${a_given} -eq 1 ]; then
        echo "-a is ${afile}" >&2
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
}

########
extract_normal_sample_info()
{
    local_entry=$1
    sample1=`echo ${local_entry} | $AWK -F ";" '{print $1}'` | $GREP 'Normal\|normal'
    sample2=`echo ${local_entry} | $AWK -F ";" '{print $2}'` | $GREP 'Normal\|normal'

    if [ ! -z "${sample1}" ]; then
        echo ${sample1}
    fi

    if [ ! -z "${sample2}" ]; then
        echo ${sample2}
    fi

    echo ""
}

########
extract_tumor_sample_info()
{
    local_entry=$1
    sample1=`echo ${local_entry} | $AWK -F ";" '{print $1}'` | $GREP 'Tumour\|tumour'
    sample2=`echo ${local_entry} | $AWK -F ";" '{print $2}'` | $GREP 'Tumour\|tumour'

    if [ ! -z "${sample1}" ]; then
        echo ${sample1}
    fi

    if [ ! -z "${sample2}" ]; then
        echo ${sample2}
    fi

    echo ""
}

########
egadata_entry_is_ok()
{
    local_entry=$1
    nsample=`extract_normal_sample_info ${local_entry}`
    tsample=`extract_tumor_sample_info ${local_entry}`

    if [ ! -z "${nsample}" -a ! -z "${tsample}" ]; then
        echo "yes"
    else
        echo "no"
    fi
}

########
extract_egafid_from_sample_info()
{
    local_sample_info=$1
    echo ${local_sample_info} | $AWK '{print $3}'
}

########
extract_gender_from_sample_info()
{
    local_sample_info=$1
    echo ${local_sample_info} | $AWK '{print $NF}'
}

########
process_pars()
{
    # Read EGA data file
    entry_num=1
    while read entry; do
        entry_ok=`egadata_entry_is_ok "$entry"`
        if [ ${entry_ok} = "yes" ]; then
            # Extract sample info
            normal_sample_info=`extract_normal_sample_info $entry`
            egan_id=`extract_egafid_from_sample_info ${normal_sample_info}`
            
            tumor_sample_info=`extract_tumor_sample_info $entry`
            egat_id=`extract_egafid_from_sample_info ${normal_sample_info}`

            gender=`extract_gender_from_sample_info ${normal_sample_info}`

            # Obtain value for -g option
            if [ ${gender} = "gender=male" ]; then
                gender_opt="XY"
            else                
                gender_opt="XX"
            fi
            
            # Set name of output directory for analysis
            outd=${egan_id}"_"${egat_id}
            
            # Submit bam analysis for normal and tumor samples
            echo ${bindir}/submit_bam_analysis -r ${ref} -egan ${egan_id} -egat ${egat_id} -a ${afile} -g ${gender_opt} -o ${outd} -wcr ${wcref} -sv ${snpvcf} -sg ${snpgccorr} -mc ${malesexchr} -egastr ${egastr} -egacred ${egacred}
        else
            echo "Error in entry number ${entry_num}"
        fi

        entry_num=`expr ${entry_num} + 1`
        
    done < ${egadata}
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

process_pars || exit 1
