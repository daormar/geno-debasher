# *- bash -*

# INCLUDE BASH LIBRARY
. ${bindir}/bam_utils_lib.sh

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
    echo "                     -a <string> -p <string>"
    echo "                     [-egastr <int>] [-egacred <string>]"
    echo "                     [--debug] [--help]"
    echo ""
    echo "-egalist <string>    File with list of EGA ids of normal bam files (one"
    echo "                     per line)"
    echo "-o <string>          Output file"
    echo "-T <string>          Directory for temporary files"
    echo "-a <string>          Slurm account name used to submit jobs"
    echo "-p <string>          Partition where the jobs should be executed"
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
    egalist_given=0
    o_given=0
    T_given=0
    a_given=0
    p_given=0
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
            "-a") shift
                  if [ $# -ne 0 ]; then
                      account=$1
                      a_given=1
                  fi
                  ;;
            "-p") shift
                  if [ $# -ne 0 ]; then
                      partition=$1
                      p_given=1
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

    if [ ${a_given} -eq 0 ]; then
        echo "Error! -a parameter not given!" >&2
        exit 1
    fi

    if [ ${p_given} -eq 0 ]; then
        echo "Error! -p parameter not given!" >&2
        exit 1
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
    local_datadir=$1
    local_egaid=$2
    local_egastr=$3
    local_egacred=$4

    ### Download bam file

    # Activate conda environment
    conda activate pyega3 || exit 1

    # Download file
    pyega3 -c ${local_egastr} -cf ${local_egacred} fetch ${local_egaid} ${local_datadir}/${local_egaid}.bam || exit 1

    # Deactivate conda environment
    conda deactivate
    
    ### Convert bam file into npz file

    # Activate conda environment
    conda activate wisecondorx || exit 1

    # Execute conversion
    BINSIZE=5000
    WisecondorX convert ${local_datadir}/${local_egaid}.bam ${local_datadir}/${local_egaid}.npz --binsize $BINSIZE || exit 1

    # Deactivate conda environment
    conda deactivate
}

########
gen_reffile_wisecondorx()
{
    # Initialize variables
    local_datadir=$1
    local_outf=$2

    # Activate conda environment
    conda activate wisecondorx || exit 1

    # Convert tumor bam file into npz
    BINSIZE=5000
    WisecondorX newref ${local_datadir}/*.npz ${local_outf} --binsize $BINSIZE || exit 1

    # Deactivate conda environment
    conda deactivate
}

########
remove_dir()
{
    # Initialize variables
    local_dir=$1

    # Remove directory
    rm -rf ${local_dir} || exit 1
}

########
extract_egaid_from_entry()
{
    local_entry=$1
    echo ${local_entry} | $AWK '{print $1}'
}

########
process_pars()
{
    # Initialize variables
    local_cpus=1
    local_mem=1024
    local_time=12:00:00
    local_jids=""
    
    # Read file with list of EGA ids
    while read entry; do
        # Extract EGA id
        egaid=`extract_egaid_from_entry $entry`
        
        # Process EGA id
        local_jobdeps=""
        create_script ${tmpdir}/scripts/bam_download_and_npz_conv bam_download_and_npz_conv "${tmpdir}/data $egaid $egastr $egacred"
        launch ${tmpdir}/scripts/bam_download_and_npz_conv ${account} ${partition} ${local_cpus} ${local_mem} ${local_time} "${local_jobdeps}" local_jid

        # Update variables storing jids
        local_jids="${local_jids},${local_jid}"

    done < ${egalist}

    # Generate reference file
    create_script ${tmpdir}/scripts/gen_reffile_wisecondorx gen_reffile_wisecondorx "${tmpdir}/data $outf"
    local_job_deps=`apply_deptype_to_jobids ${local_jids} afterok`
    launch ${tmpdir}/scripts/gen_reffile_wisecondorx ${account} ${partition} ${local_cpus} ${local_mem} ${local_time} "${local_job_deps}" local_jid
    local_jids="${local_jids},${local_jid}"

    if [ ${debug} -eq 0 ]; then
        # Remove directory with temporary files
        create_script ${tmpdir}/scripts/remove_dir remove_dir "${tmpdir}"
        local_job_deps=`apply_deptype_to_jobids ${local_jids} afterok`
        launch ${tmpdir}/scripts/remove_dir ${account} ${partition} ${local_cpus} ${local_mem} ${local_time} "${local_job_deps}" local_jid
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
