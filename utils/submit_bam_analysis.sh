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
    echo "submit_bam_analysis       -r <string> -n <string> -t <string> -o <string>"
    echo "                          -a <string>"
    echo "                          [-debug] [--help]"
    echo ""
    echo "-r <string>               File with reference genome."
    echo "-n <string>               File with normal bam."
    echo "-t <string>               File with tumor bam."
    echo "-o <string>               Output directory."
    echo "-a <string>               File with analysis steps to be performed."
    echo "                          Expected format:"
    echo "                          <stepname> <cpus> <mem> <time>"
    echo "-debug                    After ending, do not delete temporary files"
    echo "                          (for debugging purposes)."
    echo "--help                    Display this help and exit."
}

########
read_pars()
{
    r_given=0
    n_given=0
    t_given=0
    o_given=0
    a_given=0
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
            "-o") shift
                  if [ $# -ne 0 ]; then
                      outd=$1
                      o_given=1
                  fi
                  ;;
            "-a") shift
                  if [ $# -ne 0 ]; then
                      afile=$1
                      a_given=1
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

    if [ ${n_given} -eq 0 ]; then   
        echo "Error! -n parameter not given!" >&2
        exit 1
    else
        if [ ! -f ${normalbam} ]; then
            echo "Error! file ${normalbam} does not exist" >&2
            exit 1
        fi
    fi

    if [ ${t_given} -eq 0 ]; then
        echo "Error! -t parameter not given!" >&2
        exit 1
    else
        if [ ! -f ${tumorbam} ]; then
            echo "Error! file ${tumorbam} does not exist" >&2
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

    if [ ${a_given} -eq 0 ]; then   
        echo "Error! -s parameter not given!" >&2
        exit 1
    else
        if [ ! -f ${afile} ]; then
            echo "Error! file ${afile} does not exist" >&2
            exit 1
        fi
    fi
}

########
create_dirs()
{
    mkdir -p ${outd} || { echo "Error! cannot create output directory" >&2; return 1; }

    mkdir -p ${outd}/scripts || { echo "Error! cannot create scripts directory" >&2; return 1; }
}

########
execute_manta_somatic()
{
    # Initialize variables
    MANTA_OUTD=`get_step_dirname ${outd} ${stepname}`

    # Activate conda environment
    conda activate manta
    
    # Configure Manta
    configManta.py --normalBam ${normalbam} --tumorBam ${tumorbam} --referenceFasta ${ref} --runDir ${MANTA_OUTD} > ${MANTA_OUTD}/configManta.log 2>&1 || exit 1

    # Execute Manta
    ${MANTA_OUTD}/runWorkflow.py -m local -j ${cpus} > ${MANTA_OUTD}/runWorkflow.log 2>&1 || exit 1

    # Create file indicating that execution was finished
    touch ${MANTA_OUTD}/finished
}

########
execute_strelka_somatic()
{
    # Initialize variables
    STRELKA_OUTD=`get_step_dirname ${outd} ${stepname}`

    # Activate conda environment
    conda activate strelka

    # Configure Strelka
    configureStrelkaSomaticWorkflow.py --normalBam ${normalbam} --tumorBam ${tumorbam} --referenceFasta ${ref} --runDir ${STRELKA_OUTD} > ${STRELKA_OUTD}/configureStrelkaSomaticWorkflow.log 2>&1 || exit 1

    # Execute Strelka
    ${STRELKA_OUTD}/runWorkflow.py -m local -j ${cpus} > ${STRELKA_OUTD}/runWorkflow.log 2>&1 || exit 1
    
    # Create file indicating that execution was finished
    touch ${STRELKA_OUTD}/finished
}

########
execute_msisensor()
{
    # Initialize variables
    MSISENSOR_OUTD=`get_step_dirname ${outd} ${stepname}`
    
    # Activate conda environment
    conda activate msisensor

    # Create homopolymer and microsatellites file
    msisensor scan -d ${ref} -o ${MSISENSOR_OUTD}/msisensor.list > ${MSISENSOR_OUTD}/msisensor_scan.log 2>&1 || exit 1

    # Run MSIsensor analysis
    msisensor msi -d ${MSISENSOR_OUTD}/msisensor.list -n ${normalbam} -t ${tumorbam} -o ${MSISENSOR_OUTD}/output -l 1 -q 1 -b ${cpus} > ${MSISENSOR_OUTD}/msisensor_msi.log 2>&1 || exit 1
    
    # Create file indicating that execution was finished
    touch ${MSISENSOR_OUTD}/finished
}

########
execute_platypus_germline()
{
    # Initialize variables
    PLATYPUS_OUTD=`get_step_dirname ${outd} ${stepname}`
    
    # # Activate conda environment
    # conda activate platypus

    # Run Platypus
    python ${PLATYPUS_BUILD_DIR}/bin/Platypus.py callVariants --bamFiles=${normalbam} --refFile=${ref} --output=$3 --verbosity=1 > ${PLATYPUS_OUTD}/platypus.log 2>&1 || exit 1

    # Create file indicating that execution was finished
    touch ${PLATYPUS_OUTD}/finished
}

########
execute_cnvkit()
{
    # Initialize variables
    CNVKIT_OUTD=`get_step_dirname ${outd} ${stepname}`
    
    # Activate conda environment
    conda activate cnvkit

    # Run cnvkit
    cnvkit.py batch ${tumorbam} -n ${normalbam} -m wgs -f ${ref}  -d ${CNVKIT_OUTD} -p ${cpus} > ${CNVKIT_OUTD}/cnvkit.log 2>&1 || exit 1

    # Create file indicating that execution was finished
    touch ${CNVKIT_OUTD}/finished
}

########
execute_steps_in_afile()
{
    local_dirname=$1
    local_afile=$2
    
    # Read information about the steps to be executed
    while read entry; do
        entry_ok=`entry_is_ok "$entry"`
        if [ ${entry_ok} = "yes" ]; then
            # Extract entry information
            stepname=`extract_stepname_from_entry "$entry"`
            cpus=`extract_cpus_from_entry "$entry"`
            mem=`extract_mem_from_entry "$entry"`
            time=`extract_time_from_entry "$entry"`

            # Execute step
            execute_step ${local_dirname} ${stepname} ${cpus} ${mem} ${time}
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

create_dirs || exit 1

execute_steps_in_afile ${outd} ${afile}
