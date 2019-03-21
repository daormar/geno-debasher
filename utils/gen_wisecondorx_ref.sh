# *- bash -*

# INCLUDE BASH LIBRARY
. ${PANPIPE_HOME_DIR}/bin/panpipe_lib || exit 1

########
print_desc()
{
    echo "gen_wisecondorx_ref generates the ref file required by Wisecondorx"
    echo "type \"gen_wisecondorx_ref --help\" to get usage information"
}

########
usage()
{
    echo "gen_wisecondorx_ref  -egalist <string> -o <string> -T <string>"
    echo "                     -i <string>"
    echo "                     [-egastr <int>] [-egacred <string>]"
    echo "                     [--debug] [--help]"
    echo ""
    echo "-egalist <string>    File with list of EGA ids of normal bam files (one"
    echo "                     per line)"
    echo "-o <string>          Output file"
    echo "-T <string>          Directory for temporary files"
    echo "-i <string>          File with information about the steps to execute"
    echo "-egastr <int>        Number of streams used by the EGA download client"
    echo "                     (50 by default)"
    echo "-egacred <string>    File with EGA download client credentials"
    echo "--debug              After ending, do not delete temporary files"
    echo "                     (for debugging purposes)"
    echo "--help               Display this help and exit"
}

########
read_pars()
{
    egalist_given=0
    o_given=0
    T_given=0
    i_given=0
    egastr_given=0
    egastr=50
    egacred_given=0
    egacred=${NOFILE}
    debug=0
    while [ $# -ne 0 ]; do
        case $1 in
            "--help") usage
                      exit 1
                      ;;
            "-egalist") shift
                  if [ $# -ne 0 ]; then
                      egalist=$1
                      egalist_given=1
                  fi
                  ;;
            "-o") shift
                  if [ $# -ne 0 ]; then
                      outf=$1
                      o_given=1
                  fi
                  ;;
            "-T") shift
                  if [ $# -ne 0 ]; then
                      tdir=$1
                      T_given=1
                  fi
                  ;;
            "-i") shift
                  if [ $# -ne 0 ]; then
                      infofile=$1
                      i_given=1
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
            "--debug") debug=1
                      ;;
        esac
        shift
    done   
}

########
check_pars()
{
    if [ ${egalist_given} -eq 0 ]; then   
        echo "Error! -egalist parameter not given!" >&2
        exit 1
    else
        if [ ! -f ${egalist} ]; then
            echo "Error! file ${egalist} does not exist" >&2
            exit 1
        fi
    fi

    if [ ${T_given} -eq 0 ]; then
        echo "Error! -T parameter not given!" >&2
        exit 1
    else
        if [ ! -d ${tdir} ]; then
            echo "Error! directory for temporaries does not exist" >&2
            exit 1
        fi
    fi

    if [ ${i_given} -eq 0 ]; then
        echo "Error! -i parameter not given!" >&2
        exit 1
    else
        if [ ! -f ${infofile} ]; then
            echo "Error! file with step information does not exist" >&2
            exit 1
        fi
    fi
}

########
create_dirs()
{
    tmpdir=`mktemp -d $tdir/gen_wisecondorx_ref_XXXX`
    if [ ! -d "${tmpdir}" ]; then
        echo "Error! cannot create directory for temporary files" >&2; return 1;
    fi
    mkdir -p ${tmpdir}/scripts || { echo "Error! cannot create directory for scripts" >&2; return 1; }
    mkdir -p ${tmpdir}/data || { echo "Error! cannot create directory to store data files" >&2; return 1; }
}

########
bam_download_and_npz_conv()
{
    # Initialize variables
    local datadir=$1
    local egaid=$2
    local egastr=$3
    local egacred=$4

    ### Download bam file

    # Activate conda environment
    conda activate pyega3 || exit 1

    # Download file
    pyega3 -c ${egastr} -cf ${egacred} fetch ${egaid} ${datadir}/${egaid}.bam || exit 1

    # Deactivate conda environment
    conda deactivate

    ### Index bam file
    
    # Activate conda environment
    conda activate samtools || exit 1

    # Index file
    samtools index ${datadir}/${egaid}.bam || exit 1

    # Deactivate conda environment
    conda deactivate
    
    ### Convert bam file into npz file

    # Activate conda environment
    conda activate wisecondorx || exit 1

    # Execute conversion
    BINSIZE=5000
    WisecondorX convert ${datadir}/${egaid}.bam ${datadir}/${egaid}.npz --binsize $BINSIZE || exit 1

    # Deactivate conda environment
    conda deactivate
}

########
gen_reffile_wisecondorx()
{
    # Initialize variables
    local datadir=$1
    local outf=$2

    # Activate conda environment
    conda activate wisecondorx || exit 1

    # Convert tumor bam file into npz
    BINSIZE=5000
    WisecondorX newref ${datadir}/*.npz ${outf} --binsize $BINSIZE || exit 1

    # Deactivate conda environment
    conda deactivate
}

########
remove_dir()
{
    # Initialize variables
    local dir=$1

    # Remove directory
    rm -rf ${dir} || exit 1
}

########
extract_egaid_from_entry()
{
    local entry=$1
    echo ${entry} | $AWK '{print $1}'
}

########
get_pars_bam_download_and_npz_conv()
{
    echo "${tmpdir}/data $egaid $egastr $egacred"
}

########
get_pars_gen_reffile_wisecondorx()
{
    echo "${tmpdir}/data $outf"
}

########
get_pars_remove_dir()
{
    echo "${tmpdir}"
}

########
process_pars()
{
    # Read file with list of EGA ids
    local ids=""
    while read entry; do
        # Extract EGA id
        egaid=`extract_egaid_from_entry $entry`
        
        # Process EGA id
        local stepname=bam_download_and_npz_conv
        local job_array_list="1"
        local stepinfo=`get_step_info ${stepname} ${infofile}`
        local stepdeps=""
        local script_pars=`get_pars_${stepname}`
        launch_step ${tmpdir} ${stepname} ${job_array_list} "${stepinfo}" "${stepdeps}" "${script_pars}" id

        # Update variables storing ids
        if [ -z "${ids}" ]; then
            ids=${id}
        else
            ids="${ids},${id}"
        fi

    done < ${egalist}

    # Generate reference file
    local stepname=gen_reffile_wisecondorx
    local job_array_list="1"
    local stepinfo=`get_step_info ${stepname} ${infofile}`
    local stepdeps=`apply_deptype_to_stepids ${ids} afterok`
    local script_pars=`get_pars_${stepname}`
    launch_step ${tmpdir} ${stepname} ${job_array_list} "${stepinfo}" "${stepdeps}" "${script_pars}" id

    # Update variables storing ids
    ids="${ids},${id}"

    if [ ${debug} -eq 0 ]; then
        # Remove directory with temporary files
        local stepname=remove_dir
        local job_array_list="1"
        local stepinfo=`get_step_info ${stepname} ${infofile}`
        local stepdeps=`apply_deptype_to_stepids ${ids} afterok`
        local script_pars=`get_pars_${stepname}`
        launch_step ${tmpdir} ${stepname} ${job_array_list} "${stepinfo}" "${stepdeps}" "${script_pars}" id
    fi
}

########

if [ $# -eq 0 ]; then
    print_desc
    exit 1
fi

read_pars $@ || exit 1

check_pars || exit 1

create_dirs || exit 1

process_pars || exit 1
