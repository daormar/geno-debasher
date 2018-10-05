# *- bash -*

# INCLUDE BASH LIBRARY
. ${bindir}/bam_utils_lib.sh

########
print_desc()
{
    echo "analyze_study analyses samples of a given study"
    echo "type \"analyze_ega_study --help\" to get usage information"
}

########
usage()
{
    echo "analyze_study        -r <string> -e <string>|-i <string>"
    echo "                     -a <string> -o <string>"
    echo "                     [-wcr <string>] [-sv <string>]"
    echo "                     [-sg <string>] [-mc <string>]"
    echo "                     [-egastr <int>] [-egacred <string>]"
    echo "                     [-p] [-debug] [--help]"
    echo ""
    echo "-r <string>          File with reference genome"
    echo "-e <string>          File with processed EGA metadata using the"
    echo "                     query_ega_metadata tool (option -f 3)"
    echo "-i <string>          File with processed ICGC metadata using the"
    echo "                     query_icgc_metadata tool (option -f 4)"
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
    echo "-p                   Only print the commands executing the analysis"
    echo "-debug               After ending, do not delete temporary files"
    echo "                     (for debugging purposes)"
    echo "--help               Display this help and exit"
}

########
read_pars()
{
    r_given=0
    e_given=0
    i_given=0
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
    p_given=0
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
            "-i") shift
                  if [ $# -ne 0 ]; then
                      icgcdata=$1
                      i_given=1
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
            "-p") p_given=1
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
    
    if [ ${e_given} -eq 0 -a ${i_given} -eq 0 ]; then
        echo "Error, -e or -i options should be given" >&2
    fi

    if [ ${e_given} -eq 1 -a ${i_given} -eq 1 ]; then
        echo "Error, -e and -i options cannot be given simultaneously" >&2
    fi

    if [ ${e_given} -eq 1 -a ! -f ${egadata} ]; then
            echo "Error! file ${egadata} does not exist" >&2
            exit 1
    fi

    if [ ${i_given} -eq 1 -a ! -f ${icgcdata} ]; then
            echo "Error! file ${icgcdata} does not exist" >&2
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

    if [ ${e_given} -eq 1 ]; then
        echo "-e is ${egadata}" >&2
    fi

    if [ ${i_given} -eq 1 ]; then
        echo "-e is ${icgcdata}" >&2
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
extract_normal_sample_info_ega()
{
    local_entry=$1
    sample1=`echo ${local_entry} | $AWK -F ";" '{print $1}' | $GREP 'Normal\|normal'`
    sample2=`echo ${local_entry} | $AWK -F ";" '{print $2}' | $GREP 'Normal\|normal'`

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
extract_tumor_sample_info_ega()
{
    local_entry=$1
    sample1=`echo ${local_entry} | $AWK -F ";" '{print $1}' | $GREP 'Tumour\|tumour\|Tumor\|tumor'`
    sample2=`echo ${local_entry} | $AWK -F ";" '{print $2}' | $GREP 'Tumour\|tumour\|Tumor\|tumor'`

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
egadata_entry_is_ok()
{
    local_entry=$1
    nsample=`extract_normal_sample_info_ega "${local_entry}"`
    tsample=`extract_tumor_sample_info_ega "${local_entry}"`

    if [ ! -z "${nsample}" -a ! -z "${tsample}" ]; then
        echo "yes"
    else
        echo "no"
    fi
}

########
extract_fid_from_sample_info_ega()
{
    local_sample_info_ega=$1
    echo ${local_sample_info_ega} | $AWK '{print $3}'
}

########
extract_gender_from_sample_info_ega()
{
    local_sample_info_ega=$1
    echo ${local_sample_info_ega} | $AWK '{print $NF}'
}

########
analyze_ega_study()
{
    # Read EGA data file
    entry_num=1
    while read entry; do
        entry_ok=`egadata_entry_is_ok "$entry"`
        if [ ${entry_ok} = "yes" ]; then
            # Extract sample info
            normal_sample_info_ega=`extract_normal_sample_info_ega "$entry"`
            egan_id=`extract_fid_from_sample_info_ega "${normal_sample_info_ega}"`
            
            tumor_sample_info_ega=`extract_tumor_sample_info_ega "$entry"`
            egat_id=`extract_fid_from_sample_info_ega "${tumor_sample_info_ega}"`

            gender=`extract_gender_from_sample_info_ega "${normal_sample_info_ega}"`

            # Obtain value for -g option
            if [ ${gender} = "gender=male" ]; then
                gender_opt="XY"
            else                
                gender_opt="XX"
            fi
            
            # Set name of output directory for analysis
            analysis_outd=${egan_id}"_"${egat_id}
            
            # Submit bam analysis for normal and tumor samples
            if [ ${p_given} -eq 0 ]; then
                ${bindir}/submit_bam_analysis -r ${ref} -extn ${egan_id} -extt ${egat_id} -a ${afile} -g ${gender_opt} -o ${outd}/${analysis_outd} -wcr ${wcref} -sv ${snpvcf} -sg ${snpgccorr} -mc ${malesexchr} -egastr ${egastr} -egacred ${egacred}
            else
                echo ${bindir}/submit_bam_analysis -r ${ref} -extn ${egan_id} -extt ${egat_id} -a ${afile} -g ${gender_opt} -o ${outd}/${analysis_outd} -wcr ${wcref} -sv ${snpvcf} -sg ${snpgccorr} -mc ${malesexchr} -egastr ${egastr} -egacred ${egacred}
            fi
        else
            echo "Error in entry number ${entry_num}"
        fi

        entry_num=`expr ${entry_num} + 1`
        
    done < ${egadata}
}

########
extract_normal_sample_info_icgc()
{
    local_entry=$1
    sample1=`echo ${local_entry} | $AWK -F ";" '{print $1}' | $GREP 'Normal\|normal'`
    sample2=`echo ${local_entry} | $AWK -F ";" '{print $2}' | $GREP 'Normal\|normal'`

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
extract_tumor_sample_info_icgc()
{
    local_entry=$1
    sample1=`echo ${local_entry} | $AWK -F ";" '{print $1}' | $GREP 'Tumour\|tumour\|Tumor\|tumor'`
    sample2=`echo ${local_entry} | $AWK -F ";" '{print $2}' | $GREP 'Tumour\|tumour\|Tumor\|tumor'`

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
icgcdata_entry_is_ok()
{
    local_entry=$1
    nsample=`extract_normal_sample_info_icgc "${local_entry}"`
    tsample=`extract_tumor_sample_info_icgc "${local_entry}"`

    if [ ! -z "${nsample}" -a ! -z "${tsample}" ]; then
        echo "yes"
    else
        echo "no"
    fi
}

########
extract_oid_from_sample_info_icgc()
{
    local_sample_info_icgc=$1
    echo ${local_sample_info_icgc} | $AWK '{print $2}'
}

########
extract_gender_from_sample_info_icgc()
{
    local_sample_info_icgc=$1
    echo ${local_sample_info_icgc} | $AWK '{print $NF}'
}

########
analyze_icgc_study()
{
    # Read ICGC data file
    entry_num=1
    while read entry; do
        entry_ok=`icgcdata_entry_is_ok "$entry"`
        if [ ${entry_ok} = "yes" ]; then
            # Extract sample info
            normal_sample_info_icgc=`extract_normal_sample_info_icgc "$entry"`
            icgcn_id=`extract_oid_from_sample_info_icgc "${normal_sample_info_icgc}"`
            
            tumor_sample_info_icgc=`extract_tumor_sample_info_icgc "$entry"`
            icgct_id=`extract_oid_from_sample_info_icgc "${tumor_sample_info_icgc}"`

            gender=`extract_gender_from_sample_info_icgc "${normal_sample_info_icgc}"`

            # Obtain value for -g option
            if [ ${gender} = "male" ]; then
                gender_opt="XY"
            else                
                gender_opt="XX"
            fi
            
            # Set name of output directory for analysis
            analysis_outd=${icgcn_id}"_"${icgct_id}
            
            # Submit bam analysis for normal and tumor samples
            if [ ${p_given} -eq 0 ]; then
                ${bindir}/submit_bam_analysis -r ${ref} -extn ${icgcn_id} -extt ${icgct_id} -a ${afile} -g ${gender_opt} -o ${outd}/${analysis_outd} -wcr ${wcref} -sv ${snpvcf} -sg ${snpgccorr} -mc ${malesexchr}
            else
                echo ${bindir}/submit_bam_analysis -r ${ref} -extn ${icgcn_id} -extt ${icgct_id} -a ${afile} -g ${gender_opt} -o ${outd}/${analysis_outd} -wcr ${wcref} -sv ${snpvcf} -sg ${snpgccorr} -mc ${malesexchr}
            fi
        else
            echo "Error in entry number ${entry_num}"
        fi

        entry_num=`expr ${entry_num} + 1`
        
    done < ${icgcdata}
}

########
process_pars()
{
    if [ ${e_given} -eq 1 ]; then
        analyze_ega_study
    fi

    if [ ${i_given} -eq 1 ]; then
        analyze_icgc_study
    fi
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
