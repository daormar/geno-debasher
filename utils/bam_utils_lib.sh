# *- bash -*


#############
# CONSTANTS #
#############

NOFILE="_NONE_"

#####################
# GENERAL FUNCTIONS #
#####################

init_bash_shebang_var()
{
    echo "#!${BASH}"
}

########
is_absolute_path()
{
    case $1 in
        /*) echo 1 ;;
        *) echo 0 ;;
    esac
}

########
get_absolute_path()
{
    local file=$1
    # Check if an absolute path was given
    local absolute=`is_absolute_path $file`
    if [ $absolute -eq 1 ]; then
        echo $file
    else
        local oldpwd=$PWD
        local basetmp=`$BASENAME $PWD/$file`
        local dirtmp=`$DIRNAME $PWD/$file`
        cd $dirtmp
        local result=${PWD}/${basetmp}
        cd $oldpwd
        echo $result
    fi
}

########
exclude_readonly_vars()
{
    $AWK -F "=" 'BEGIN{
                         readonlyvars["BASHOPTS"]=1
                         readonlyvars["BASH_VERSINFO"]=1
                         readonlyvars["EUID"]=1
                         readonlyvars["PPID"]=1
                         readonlyvars["SHELLOPTS"]=1
                         readonlyvars["UID"]=1
                        }
                        {
                         if(!($1 in readonlyvars)) printf"%s\n",$0
                        }'
}

########
exclude_bashisms()
{
    $AWK '{if(index($1,"=(")==0) printf"%s\n",$0}'
}

########
create_script()
{
    # Init variables
    local name=$1
    local command=$2
    local script_pars=$3
    
    # Write bash shebang
    local BASH_SHEBANG=`init_bash_shebang_var`
    echo ${BASH_SHEBANG} > ${name} || return 1

    # Write SLURM commands
    echo "#SBATCH --job-name=${command}" >> ${name} || return 1
    echo "#SBATCH --output=${name}.slurm_out" >> ${name} || return 1

    # Write environment variables
    set | exclude_readonly_vars | exclude_bashisms >> ${name} || return 1
    
    # Write command to be executed
    echo "${command} ${script_pars}" >> ${name} || return 1

    # Give execution permission
    chmod u+x ${name} || return 1
}

########
get_lib_timestamp()
{
    $GREP "Lib time stamp:" ${bindir}/bam_utils_lib | $AWK '{if(NF==5) printf"%s",$NF}'
}

########
get_file_timestamp()
{
    local file=$1
    ${PYTHON} -c "import os; t=os.path.getmtime('${file}'); print int(t)"
}

########
get_account_opt()
{
    local account=$1

    if [ -z "${account}" ]; then
        echo ""
    else
        echo "-A ${account}"
    fi
}

########
get_partition_opt()
{
    local partition=$1

    if [ -z "${partition}" ]; then
        echo ""
    else
        echo "--partition=${partition}"
    fi
}

########
get_script_filename() 
{
    local stepname=$1
    
    echo ${dirname}/scripts/${stepname}
}

########
remove_suffix_from_stepname()
{
    local stepname=$1
    
    echo ${stepname} | $AWK '{if(index($1,"__")==0){print $1} else{printf "%s\n",substr($1,1,index($1,"__")-1)}}'
}

########
get_step_function()
{
    local stepname=$1

    local stepname_wo_suffix=`remove_suffix_from_stepname ${stepname}`
    
    echo "${stepname_wo_suffix}"
}

########
get_script_pars_funcname()
{
    local stepname=$1

    local stepname_wo_suffix=`remove_suffix_from_stepname ${stepname}`

    echo get_pars_${stepname_wo_suffix}
}

########
find_dependency_for_step()
{
    local jobdeps_spec=$1
    local stepname_part=$2

    for local_dep in `echo ${jobdeps_spec} | $SED 's/,/ /g'`; do
        local stepname_part_in_dep=`echo ${local_dep} | $AWK -F ":" '{print $2}'`
        if [ ${stepname_part_in_dep} = ${stepname_part} ]; then
            echo ${local_dep}
        fi
    done
}

########
get_outd_for_dep()
{
    local outd=$1
    local dep=$2

    if [ -z "${dep}" ]; then
        echo ""
    else
        local stepname_part=`echo ${dep} | $AWK -F ":" '{print $2}'`
        get_step_dirname ${outd} ${stepname_part}
    fi
}

########
apply_deptype_to_jobids()
{
    # Initialize variables
    local jids=$1
    local deptype=$2

    # Apply deptype
    result=""
    for local_jid in `echo ${jids} | $SED 's/,/ /g'`; do
        if [ -z "" ]; then
            result=${deptype}:${local_jid}
        else
            result=${result}","${deptype}:${local_jid}
        fi
    done

    echo $result
}

########
get_slurm_dependency_opt()
{
    local jobdeps=$1

    # Create dependency option
    if [ -z "${jobdeps}" ]; then
        echo ""
    else
        echo "--dependency=${jobdeps}"
    fi
}

########
launch()
{
    # Initialize variables
    local file=$1
    local account=$2
    local partition=$3
    local cpus=$4
    local mem=$5
    local time=$6
    local jobdeps=$7
    local outvar=$8

    # Launch file
    if [ -z "${SBATCH}" ]; then
        ${file} > ${file}.log 2>&1 || return 1
        eval "${outvar}=\"\""
    else
        account_opt=`get_account_opt ${account}`
        partition_opt=`get_partition_opt ${partition}`
        dependency_opt=`get_slurm_dependency_opt "${jobdeps}"`
        local jid=$($SBATCH --cpus-per-task=${cpus} --mem=${mem} --time ${time} --parsable ${account_opt} ${partition_opt} ${dependency_opt} ${file})
        eval "${outvar}='${jid}'"
    fi
}

########
launch_step()
{
    # Initialize variables
    local stepname=$1
    local stepinfo=$2
    local jobdeps=$3
    local script_pars=$4
    local jid=$5

    # Create script
    create_script ${tmpdir}/scripts/${stepname} ${stepname} "${script_pars}" || return 1

    # Retrieve requirements
    local account=`extract_account_from_entry "$stepinfo"`
    local partition=`extract_partition_from_entry "$stepinfo"`
    local cpus=`extract_cpus_from_entry "$stepinfo"`
    local mem=`extract_mem_from_entry "$stepinfo"`
    local time=`extract_time_from_entry "$stepinfo"`

    # Launch script
    launch ${tmpdir}/scripts/${stepname} ${account} ${partition} ${cpus} ${mem} ${time} "${jobdeps}" ${jid} || return 1
}

########
get_step_info()
{
    local stepname=$1
    local infofile=$2

    $AWK -v stepname=${stepname} '{if($1==stepname) print $0}' ${infofile}
}

########
analysis_entry_is_ok()
{
    local entry=$1
    echo "${entry}" | $AWK '{if(NF>=4) print"yes\n"; else print"no\n"}'
}

########
extract_stepname_from_entry()
{
    local entry=$1
    echo "${entry}" | $AWK '{print $1}'
}

########
extract_account_from_entry()
{
    local entry=$1
    echo "${entry}" | $AWK '{print $2}'
}

########
extract_partition_from_entry()
{
    local entry=$1
    echo "${entry}" | $AWK '{print $3}'
}

########
extract_cpus_from_entry()
{
    local entry=$1
    echo "${entry}" | $AWK '{print $4}'
}

########
extract_mem_from_entry()
{
    local entry=$1
    echo "${entry}" | $AWK '{print $5}'
}

########
extract_time_from_entry()
{
    local entry=$1
    echo "${entry}" | $AWK '{print $6}'
}

########
extract_jobdeps_spec_from_entry()
{
    local entry=$1
    echo "${entry}" | $AWK '{print substr($7,9)}'
}

########
get_step_dirname()
{
    local dirname=$1
    local stepname=$2
    echo ${dirname}/${stepname}
}

########
reset_outdir_for_step() 
{
    local dirname=$1
    local stepname=$2
    local outd=`get_step_dirname ${dirname} ${stepname}`

    if [ -d ${outd} ]; then
        echo "Warning: ${stepname} output directory already exists but analysis was not finished, directory content will be removed">&2
        rm -rf ${outd}/* || { echo "Error! could not clear output directory" >&2; return 1; }
    else
        mkdir ${outd} || { echo "Error! cannot create output directory" >&2; return 1; }
    fi
}

########
get_step_status()
{
    local dirname=$1
    local stepname=$2
    local stepdirname=`get_step_dirname ${dirname} ${stepname}`
    
    if [ -d ${stepdirname} ]; then
        if [ -f ${stepdirname}/finished ]; then
            echo "FINISHED"
        else
            echo "UNFINISHED"
        fi
    else
        echo "TO-DO"
    fi
}

########
display_begin_step_message()
{
    if [ -z "${SLURM_JOB_ID}" ]; then
        echo "Step started at `date`" >&2
    else
        echo "Step started at `date` (SLURM_JOB_ID= ${SLURM_JOB_ID})" >&2
    fi
}

########
display_end_step_message()
{
    echo "Step finished at `date`" >&2
}

##################
# ANALYSIS STEPS #
##################

########
get_callreg_opt()
{
    local callregf=$1

    if [ ${callregf} = ${NOFILE} ]; then
        echo ""
    else
        echo "--callRegions ${callregf}"
    fi
}

########
manta_germline()
{
    display_begin_step_message

    # Initialize variables
    local ref=$1
    local normalbam=$2
    local callregf=$3
    local step_outd=$4
    local cpus=$5

    # Define --callRegions option
    call_reg_opt=`get_callreg_opt "${callregf}"`

    # Activate conda environment
    conda activate manta > ${step_outd}/conda_activate.log 2>&1 || exit 1

    # Configure Manta
    configManta.py --bam ${normalbam} --referenceFasta ${ref} ${call_reg_opt} --runDir ${step_outd} > ${step_outd}/configManta.log 2>&1 || exit 1

    # Execute Manta
    ${step_outd}/runWorkflow.py -m local -j ${cpus} > ${step_outd}/runWorkflow.log 2>&1 || exit 1

    # Deactivate conda environment
    conda deactivate > ${step_outd}/conda_deactivate.log 2>&1

    # Create file indicating that execution was finished
    touch ${step_outd}/finished

    display_end_step_message
}

########
manta_somatic()
{
    display_begin_step_message

    # Initialize variables
    local ref=$1
    local normalbam=$2
    local tumorbam=$3
    local callregf=$4
    local step_outd=$5
    local cpus=$6

    # Define --callRegions option
    call_reg_opt=`get_callreg_opt "${callregf}"`

    # Activate conda environment
    conda activate manta > ${step_outd}/conda_activate.log 2>&1 || exit 1
    
    # Configure Manta
    configManta.py --normalBam ${normalbam} --tumorBam ${tumorbam} --referenceFasta ${ref} ${call_reg_opt} --runDir ${step_outd} > ${step_outd}/configManta.log 2>&1 || exit 1

    # Execute Manta
    ${step_outd}/runWorkflow.py -m local -j ${cpus} > ${step_outd}/runWorkflow.log 2>&1 || exit 1

    # Deactivate conda environment
    conda deactivate > ${step_outd}/conda_deactivate.log 2>&1

    # Create file indicating that execution was finished
    touch ${step_outd}/finished

    display_end_step_message
}

########
get_indel_cand_opt()
{
    local manta_outd=$1

    if [ -z "${manta_outd}" ]; then
        echo ""
    else
        manta_indel_file="${manta_outd}/results/variants/candidateSmallIndels.vcf.gz"
        if [ -f ${manta_indel_file} ]; then
            echo "--indelCandidates ${manta_indel_file}"
        else
            echo "WARNING: Manta indel file for Strelka not found! (${manta_indel_file})" >&2
            echo ""
        fi
    fi
}

########
strelka_germline()
{
    display_begin_step_message

    # Initialize variables
    local ref=$1
    local normalbam=$2
    local callregf=$3
    local step_outd=$4
    local cpus=$5

    # Define --callRegions option
    call_reg_opt=`get_callreg_opt "${callregf}"`

    # Activate conda environment
    conda activate strelka > ${step_outd}/conda_activate.log 2>&1 || exit 1

    # Configure Strelka
    configureStrelkaGermlineWorkflow.py --bam ${normalbam} --referenceFasta ${ref} ${call_reg_opt} --runDir ${step_outd} > ${step_outd}/configureStrelkaGermlineWorkflow.log 2>&1 || exit 1

    # Execute Strelka
    ${step_outd}/runWorkflow.py -m local -j ${cpus} > ${step_outd}/runWorkflow.log 2>&1 || exit 1

    # Deactivate conda environment
    conda deactivate > ${step_outd}/conda_deactivate.log 2>&1

    # Create file indicating that execution was finished
    touch ${step_outd}/finished

    display_end_step_message
}

########
strelka_somatic()
{
    display_begin_step_message

    # Initialize variables
    local ref=$1
    local normalbam=$2
    local tumorbam=$3
    local callregf=$4
    local step_outd=$5
    local manta_outd=$6
    local cpus=$7

    # Define --indelCandidates option if output from Manta is available
    indel_cand_opt=`get_indel_cand_opt "${manta_outd}"`

    # Define --callRegions option
    call_reg_opt=`get_callreg_opt "${callregf}"`

    # Activate conda environment
    conda activate strelka > ${step_outd}/conda_activate.log 2>&1 || exit 1

    # Configure Strelka
    configureStrelkaSomaticWorkflow.py --normalBam ${normalbam} --tumorBam ${tumorbam} --referenceFasta ${ref} ${indel_cand_opt} ${call_reg_opt} --runDir ${step_outd} > ${step_outd}/configureStrelkaSomaticWorkflow.log 2>&1 || exit 1

    # Execute Strelka
    ${step_outd}/runWorkflow.py -m local -j ${cpus} > ${step_outd}/runWorkflow.log 2>&1 || exit 1

    # Deactivate conda environment
    conda deactivate > ${step_outd}/conda_deactivate.log 2>&1

    # Create file indicating that execution was finished
    touch ${step_outd}/finished

    display_end_step_message
}

########
msisensor()
{
    display_begin_step_message

    # Initialize variables
    local ref=$1
    local normalbam=$2
    local tumorbam=$3
    local step_outd=$4
    local cpus=$5

    # Activate conda environment
    conda activate msisensor > ${step_outd}/conda_activate.log 2>&1 || exit 1

    # Create homopolymer and microsatellites file
    msisensor scan -d ${ref} -o ${step_outd}/msisensor.list > ${step_outd}/msisensor_scan.log 2>&1 || exit 1

    # Run MSIsensor analysis
    msisensor msi -d ${step_outd}/msisensor.list -n ${normalbam} -t ${tumorbam} -o ${step_outd}/output -l 1 -q 1 -b ${cpus} > ${step_outd}/msisensor_msi.log 2>&1 || exit 1

    # Deactivate conda environment
    conda deactivate > ${step_outd}/conda_deactivate.log 2>&1
    
    # Create file indicating that execution was finished
    touch ${step_outd}/finished

    display_end_step_message
}

########
platypus_germline_conda()
{
    display_begin_step_message
    
    # Initialize variables
    local ref=$1
    local normalbam=$2
    local step_outd=$3

    # Activate conda environment
    conda activate platypus > ${step_outd}/conda_activate.log 2>&1 || exit 1

    # Run Platypus
    Platypus.py callVariants --bamFiles=${normalbam} --refFile=${ref} --output=${step_outd}/output.vcf --verbosity=1 > ${step_outd}/platypus.log 2>&1 || exit 1

    # Deactivate conda environment
    conda deactivate > ${step_outd}/conda_deactivate.log 2>&1

    # Create file indicating that execution was finished
    touch ${step_outd}/finished        

    display_end_step_message
}

########
platypus_germline_local()
{
    display_begin_step_message

    # Initialize variables
    local ref=$1
    local normalbam=$2
    local step_outd=$3

    # Run Platypus
    python ${PLATYPUS_HOME_DIR}/bin/Platypus.py callVariants --bamFiles=${normalbam} --refFile=${ref} --output=${step_outd}/output.vcf --verbosity=1 > ${step_outd}/platypus.log 2>&1 || exit 1

    # Create file indicating that execution was finished
    touch ${step_outd}/finished    

    display_end_step_message
}

########
platypus_germline()
{
    # Initialize variables
    local ref=$1
    local normalbam=$2
    local step_outd=$3

    if [ -z "${PLATYPUS_HOME_DIR}" ]; then
        platypus_germline_conda ${ref} ${normalbam} ${step_outd}
    else
        platypus_germline_local ${ref} ${normalbam} ${step_outd}
    fi
}

########
cnvkit()
{
    display_begin_step_message

    # Initialize variables
    local ref=$1
    local normalbam=$2
    local tumorbam=$3
    local step_outd=$4
    local cpus=$5
    
    # Activate conda environment
    conda activate cnvkit > ${step_outd}/conda_activate.log 2>&1 || exit 1

    # Run cnvkit
    cnvkit.py batch ${tumorbam} -n ${normalbam} -m wgs -f ${ref}  -d ${step_outd} -p ${cpus} > ${step_outd}/cnvkit.log 2>&1 || exit 1

    # Deactivate conda environment
    conda deactivate > ${step_outd}/conda_deactivate.log 2>&1

    # Create file indicating that execution was finished
    touch ${step_outd}/finished

    display_end_step_message
}

########
wisecondorx()
{
    display_begin_step_message

    # Initialize variables
    local wcref=$1
    local tumorbam=$2
    local step_outd=$3
    local cpus=$4
    
    # Activate conda environment
    conda activate wisecondorx > ${step_outd}/conda_activate.log 2>&1 || exit 1

    # Convert tumor bam file into npz
    BINSIZE=5000
    WisecondorX convert ${tumorbam} ${step_outd}/tumor.npz --binsize $BINSIZE > ${step_outd}/wisecondorx_convert.log 2>&1 || exit 1
    
    # Use WisecondorX for prediction
    WisecondorX predict ${step_outd}/tumor.npz ${wcref} ${step_outd}/out > ${step_outd}/wisecondorx_predict.log 2>&1 || exit 1

    # Deactivate conda environment
    conda deactivate > ${step_outd}/conda_deactivate.log 2>&1

    # Create file indicating that execution was finished
    touch ${step_outd}/finished

    display_end_step_message
}

########
facets()
{
    display_begin_step_message

    # Initialize variables
    local normalbam=$1
    local tumorbam=$2
    local snpvcf=$3
    local step_outd=$4

    # Activate conda environment if needed
    if [ -z "${FACETS_HOME_DIR}" ]; then
        conda activate facets > ${step_outd}/conda_activate.log 2>&1 || exit 1
    fi
        
    # Execute snp-pileup
    if [ -z "${FACETS_HOME_DIR}" ]; then
        snp-pileup ${snpvcf} ${step_outd}/snp-pileup-counts.csv ${normalbam} ${tumorbam} > ${step_outd}/snp-pileup.log 2>&1 || exit 1
    else
        ${FACETS_HOME_DIR}/inst/extcode/snp-pileup ${snpvcf} ${step_outd}/snp-pileup-counts.csv ${normalbam} ${tumorbam} > ${step_outd}/snp-pileup.log 2>&1 || exit 1
    fi
    
    # Execute facets
    ${bindir}/run_facets -c ${step_outd}/snp-pileup-counts.csv > ${step_outd}/facets.out 2> ${step_outd}/run_facets.log || exit 1

    # Deactivate conda environment if needed
    if [ -z "${FACETS_HOME_DIR}" ]; then
        conda deactivate > ${step_outd}/conda_deactivate.log 2>&1
    fi

    # Create file indicating that execution was finished
    touch ${step_outd}/finished

    display_end_step_message
}

########
ascatngs()
{
    display_begin_step_message

    # Initialize variables
    local ref=$1
    local normalbam=$2
    local tumorbam=$3
    local gender=$4
    local malesexchr=$5
    local snpgccorr=$6
    local step_outd=$7
    local cpus=$8
    
    # Activate conda environment
    conda activate ascatngs > ${step_outd}/conda_activate.log 2>&1 || exit 1

    # Run cnvkit
    ascat.pl -n ${normalbam} -t ${tumorbam} -r ${ref} -sg ${snpgccorr} -pr WGS -g ${gender} -gc ${malesexchr} -cpus ${cpus} -o ${step_outd} > ${step_outd}/ascat.log 2>&1 || exit 1

    # Deactivate conda environment
    conda deactivate > ${step_outd}/conda_deactivate.log 2>&1

    # Create file indicating that execution was finished
    touch ${step_outd}/finished

    display_end_step_message
}

########
ega_download_retry()
{
    # Initialize variables
    local egastr=$1
    local egacred=$2
    local egaid=$3
    local outf=$4
    local download_tries=$5
    local step_outd=`${DIRNAME} ${outf}`
    
    # Start download with multiple tries
    local ntry=1
    while [ ${ntry} -le ${download_tries} ]; do
        echo "Starting download try number ${ntry}..." >&2

        # Remove previously downloaded file (if any)
        if [ -f ${outf} ]; then
            rm ${outf}
        fi

        # Download file
        pyega3 -c ${egastr} -cf ${egacred} fetch ${egaid} ${outf} > ${step_outd}/pyega3.log 2>&1
        
        # Check if download was successful
        if [ $? -eq 0 -a -f ${outf} ]; then
            return 0
        fi

        # Save log file
        cp ${step_outd}/pyega3.log ${step_outd}/pyega3.log.attempt${ntry}

        ntry=`expr ${ntry} + 1`
    done

    echo "All download attempts failed!" >&2

    return 1
}

########
download_ega_norm_bam()
{
    display_begin_step_message

    # Initialize variables
    local normalbam=$1
    local egaid_normalbam=$2
    local egastr=$3
    local egacred=$4
    local download_tries=$5
    local step_outd=$6

    # Activate conda environment
    conda activate pyega3 > ${step_outd}/conda_activate.log 2>&1 || exit 1

    # Download file (with multiple tries)
    ega_download_retry ${egastr} ${egacred} ${egaid_normalbam} ${step_outd}/normal.bam ${download_tries} || exit 1

    # Move file
    mv ${step_outd}/normal.bam ${normalbam} || exit 1
    
    # Deactivate conda environment
    conda deactivate > ${step_outd}/conda_deactivate.log 2>&1

    # Create file indicating that execution was finished
    touch ${step_outd}/finished

    display_end_step_message
}

########
download_ega_tum_bam()
{
    display_begin_step_message

    # Initialize variables
    local tumorbam=$1
    local egaid_tumorbam=$2
    local egastr=$3
    local egacred=$4
    local download_tries=$5
    local step_outd=$6

    # Activate conda environment
    conda activate pyega3 > ${step_outd}/conda_activate.log 2>&1 || exit 1

    # Download file (with multiple tries)
    ega_download_retry ${egastr} ${egacred} ${egaid_tumorbam} ${step_outd}/tumor.bam ${download_tries} || exit 1

    # Move file
    mv ${step_outd}/tumor.bam ${tumorbam} || exit 1

    # Deactivate conda environment
    conda deactivate > ${step_outd}/conda_deactivate.log 2>&1

    # Create file indicating that execution was finished
    touch ${step_outd}/finished

    display_end_step_message
}

########
find_bam_filename()
{
    local step_outd=$1
    local result=""
    
    for f in ${step_outd}/*.bam; do
        if [ -f $f ]; then
            result=$f
        fi
    done

    echo ${result}
}

########
download_aws_norm_bam()
{
    display_begin_step_message

    # Initialize variables
    local normalbam=$1
    local icgcid_normalbam=$2
    local download_tries=$3
    local step_outd=$4

    # Download file
    ${ICGCSTOR_HOME_DIR}/bin/icgc-storage-client --profile aws download --object-id ${icgcid_normalbam} --output-dir ${step_outd} > ${step_outd}/icgc-storage-client.log 2>&1 || exit 1

    # Find bam file name
    local bam_file_name=`find_bam_filename ${step_outd}`
    
    if [ -z "${bam_file_name}" ]; then
        echo "Error: bam file not found after download process was completed" >&2
        exit 1
    fi

    # Move file
    mv ${bam_file_name} ${normalbam} || exit 1

    # Create file indicating that execution was finished
    touch ${step_outd}/finished

    display_end_step_message
}

########
download_aws_tum_bam()
{
    display_begin_step_message

    # Initialize variables
    local tumorbam=$1
    local icgcid_tumorbam=$2
    local download_tries=$3
    local step_outd=$4

    # Download file
    ${ICGCSTOR_HOME_DIR}/bin/icgc-storage-client --profile aws download --object-id ${icgcid_tumorbam} --output-dir ${step_outd} > ${step_outd}/icgc-storage-client.log 2>&1 || exit 1

    # Find bam file name
    local bam_file_name=`find_bam_filename ${step_outd}`
    
    if [ -z "${bam_file_name}" ]; then
        echo "Error: bam file not found after download process was completed" >&2
        exit 1
    fi

    # Move file
    mv ${bam_file_name} ${tumorbam} || exit 1

    # Create file indicating that execution was finished
    touch ${step_outd}/finished

    display_end_step_message
}

########
download_collab_norm_bam()
{
    display_begin_step_message

    # Initialize variables
    local normalbam=$1
    local icgcid_normalbam=$2
    local download_tries=$3
    local step_outd=$4

    # Download file
    ${ICGCSTOR_HOME_DIR}/bin/icgc-storage-client --profile collab download --object-id ${icgcid_normalbam} --output-dir ${step_outd} > ${step_outd}/icgc-storage-client.log 2>&1 || exit 1

    # Find bam file name
    local bam_file_name=`find_bam_filename ${step_outd}`
    
    if [ -z "${bam_file_name}" ]; then
        echo "Error: bam file not found after download process was completed" >&2
        exit 1
    fi

    # Move file
    mv ${bam_file_name} ${normalbam} || exit 1

    # Create file indicating that execution was finished
    touch ${step_outd}/finished

    display_end_step_message
}

########
download_collab_tum_bam()
{
    display_begin_step_message

    # Initialize variables
    local tumorbam=$1
    local icgcid_tumorbam=$2
    local download_tries=$3
    local step_outd=$4

    # Download file
    ${ICGCSTOR_HOME_DIR}/bin/icgc-storage-client --profile collab download --object-id ${icgcid_tumorbam} --output-dir ${step_outd} > ${step_outd}/icgc-storage-client.log 2>&1 || exit 1

    # Find bam file name
    local bam_file_name=`find_bam_filename ${step_outd}`
    
    if [ -z "${bam_file_name}" ]; then
        echo "Error: bam file not found after download process was completed" >&2
        exit 1
    fi

    # Move file
    mv ${bam_file_name} ${tumorbam} || exit 1

    # Create file indicating that execution was finished
    touch ${step_outd}/finished

    display_end_step_message
}

########
download_ega_asp_norm_bam()
{
    display_begin_step_message

    # Initialize variables
    local normalbam=$1
    local normalbam_file=$2
    local aspera_user=$3
    local aspera_passwd=$4
    local aspera_server=$5
    local egadecrypt_pwd=$6
    local download_tries=$7
    local step_outd=$8
    local max_trans_rate=100m
    
    # Download file
    ASPERA_SCP_PASS=${aspera_passwd} ${ASPERA_HOME_DIR}/bin/ascp --ignore-host-key -QTl ${max_trans_rate} ${aspera_user}@${aspera_server}:${normalbam_file} ${step_outd}/normal.bam.crypt > ${step_outd}/ascp.log 2>&1 || exit 1

    # Decrypt file
    $JAVA -jar ${EGADECRYPT_HOME_DIR}/decryptor.jar ${egadecrypt_pwd} ${step_outd}/normal.bam.crypt > ${step_outd}/decryptor.log 2>&1 || exit 1
    
    # Obtain file name
    local bam_file_name=`find_bam_filename ${step_outd}`
    
    if [ -z "${bam_file_name}" ]; then
        echo "Error: bam file not found after download process was completed" >&2
        exit 1
    fi

    # Move file
    mv ${bam_file_name} ${normalbam} || exit 1

    # Remove encrypted file
    rm ${step_outd}/normal.bam.crypt || exit 1
    
    # Create file indicating that execution was finished
    touch ${step_outd}/finished

    display_end_step_message
}

########
download_ega_asp_tum_bam()
{
    display_begin_step_message

    # Initialize variables
    local tumorbam=$1
    local tumorbam_file=$2
    local aspera_user=$3
    local aspera_passwd=$4
    local aspera_server=$5
    local egadecrypt_pwd=$6
    local download_tries=$7
    local step_outd=$8
    local max_trans_rate=100m

    # Download file
    ASPERA_SCP_PASS=${aspera_passwd} ${ASPERA_HOME_DIR}/bin/ascp --ignore-host-key -QTl ${max_trans_rate} ${aspera_user}@${aspera_server}:${tumorbam_file} ${step_outd}/tumor.bam.crypt > ${step_outd}/ascp.log 2>&1 || exit 1

    # Decrypt file
    $JAVA -jar ${EGADECRYPT_HOME_DIR}/decryptor.jar ${egadecrypt_pwd} ${step_outd}/tumor.bam.crypt > ${step_outd}/decryptor.log 2>&1 || exit 1

    # Obtain file name
    local bam_file_name=`find_bam_filename ${step_outd}`
    
    if [ -z "${bam_file_name}" ]; then
        echo "Error: bam file not found after download process was completed" >&2
        exit 1
    fi

    # Move file
    mv ${bam_file_name} ${tumorbam} || exit 1

    # Remove encrypted file
    rm ${step_outd}/tumor.bam.crypt || exit 1

    # Create file indicating that execution was finished
    touch ${step_outd}/finished

    display_end_step_message
}

########
index_norm_bam()
{
    display_begin_step_message

    # Initialize variables
    local normalbam=$1
    local step_outd=$2

    # Remove previous index if one was created
    if [ -f ${normalbam}.bai ]; then
        rm ${normalbam}.bai || exit 1
    fi
        
    # Index normal bam file

    # Activate conda environment
    conda activate base > ${step_outd}/conda_activate.log 2>&1 || exit 1

    # Execute samtools
    samtools index ${normalbam} > ${step_outd}/samtools.log 2>&1 || exit 1

    # Deactivate conda environment
    conda deactivate > ${step_outd}/conda_deactivate.log 2>&1
    
    # Create file indicating that execution was finished
    touch ${step_outd}/finished

    display_end_step_message
}

########
sort_norm_bam()
{
    display_begin_step_message

    # Initialize variables
    local normalbam=$1
    local step_outd=$2
    local cpus=$3

    # Activate conda environment
    conda activate base > ${step_outd}/conda_activate.log 2>&1 || exit 1

    # Verify if bam file is already sorted
    local bam_is_sorted=`samtools view -H ${normalbam} | $GREP SO:coordinate | wc -l` || exit 1
    if [ ${bam_is_sorted} -eq 1 ]; then
        echo "Warning: bam file is already sorted"
    else
        # Execute samtools
        samtools sort -T ${step_outd} -o ${step_outd}/sorted.bam -m 2G -@ ${cpus} ${normalbam} >  ${step_outd}/samtools.log 2>&1 || exit 1
        # NOTE: -m option is used here to increase the maximum memory per
        # thread. One lateral efect of this is that the number of tmp files
        # generated is decreased. This constitutes one possible way to avoid
        # the "Too many open files" error reported by samtools

        # Replace initial bam file by the sorted one
        mv ${step_outd}/sorted.bam ${normalbam} 2> ${step_outd}/mv.log || exit 1
    fi
    
    # Deactivate conda environment
    conda deactivate > ${step_outd}/conda_deactivate.log 2>&1

    # Create file indicating that execution was finished
    touch ${step_outd}/finished

    display_end_step_message
}

########
filter_norm_bam_contigs()
{
    display_begin_step_message

    # Initialize variables
    local ref=$1
    local normalbam=$2
    local step_outd=$3

    # Activate conda environment
    conda activate base > ${step_outd}/conda_activate.log 2>&1 || exit 1

    # Generate bed file for genome reference
    ${bindir}/gen_bed_for_genome -r ${ref} -o ${step_outd}/genref
    
    # Filter normal bam file
    samtools view -b -L ${step_outd}/genref.bed ${normalbam} > ${step_outd}/filtered.bam 2> ${step_outd}/samtools_view.log || exit 1

    # Replace initial bam file by the filtered one
    mv ${step_outd}/filtered.bam ${normalbam} 2> ${step_outd}/mv.log || exit 1

    # Deactivate conda environment
    conda deactivate > ${step_outd}/conda_deactivate.log 2>&1

    # Create file indicating that execution was finished
    touch ${step_outd}/finished

    display_end_step_message
}

########
filter_tum_bam_contigs()
{
    display_begin_step_message

    # Initialize variables
    local ref=$1
    local tumorbam=$2
    local step_outd=$3

    # Activate conda environment
    conda activate base > ${step_outd}/conda_activate.log 2>&1 || exit 1

    # Generate bed file for genome reference
    ${bindir}/gen_bed_for_genome -r ${ref} -o ${step_outd}/genref
    
    # Filter tumor bam file
    samtools view -b -L ${step_outd}/genref.bed ${tumorbam} > ${step_outd}/filtered.bam 2> ${step_outd}/samtools_view.log || exit 1

    # Replace initial bam file by the filtered one
    mv ${step_outd}/filtered.bam ${tumorbam} 2> ${step_outd}/mv.log || exit 1

    # Deactivate conda environment
    conda deactivate > ${step_outd}/conda_deactivate.log 2>&1

    # Create file indicating that execution was finished
    touch ${step_outd}/finished

    display_end_step_message
}

########
sort_tum_bam()
{
    display_begin_step_message

    # Initialize variables
    local tumorbam=$1
    local step_outd=$2
    local cpus=$3
    
    # Activate conda environment
    conda activate base > ${step_outd}/conda_activate.log 2>&1 || exit 1

    # Verify if bam file is already sorted
    local bam_is_sorted=`samtools view -H ${tumorbam} | $GREP SO:coordinate | wc -l` || exit 1
    if [ ${bam_is_sorted} -eq 1 ]; then
        echo "Warning: bam file is already sorted"
    else
        # Execute samtools
        samtools sort -T ${step_outd} -o ${step_outd}/sorted.bam -m 2G -@ ${cpus} ${tumorbam} >  ${step_outd}/samtools.log 2>&1 || exit 1
        # NOTE: -m option is used here to increase the maximum memory per
        # thread. One lateral efect of this is that the number of tmp files
        # generated is decreased. This constitutes one possible way to avoid
        # the "Too many open files" error reported by samtools

        # Replace initial bam file by the sorted one
        mv ${step_outd}/sorted.bam ${tumorbam} 2> ${step_outd}/mv.log || exit 1
    fi
    
    # Deactivate conda environment
    conda deactivate > ${step_outd}/conda_deactivate.log 2>&1

    # Create file indicating that execution was finished
    touch ${step_outd}/finished

    display_end_step_message
}

########
index_tum_bam()
{
    display_begin_step_message

    # Initialize variables
    local tumorbam=$1
    local step_outd=$2

    # Remove previous index if one was created
    if [ -f ${tumorbam}.bai ]; then
        rm ${tumorbam}.bai || exit 1
    fi

    # Index tumor bam file

    # Activate conda environment
    conda activate base > ${step_outd}/conda_activate.log 2>&1 || exit 1
    
    # Execute samtools
    samtools index ${tumorbam} > ${step_outd}/samtools.log 2>&1 || exit 1

    # Deactivate conda environment
    conda deactivate > ${step_outd}/conda_deactivate.log 2>&1
    
    # Create file indicating that execution was finished
    touch ${step_outd}/finished

    display_end_step_message
}

########
delete_bam_files()
{
    display_begin_step_message

    # Initialize variables
    local normalbam=$1
    local tumorbam=$2
    local step_outd=$3

    # Delete normal bam file
    rm ${normalbam} > ${step_outd}/rm_norm.log 2>&1 || exit 1
    
    # Delete tumor bam file
    rm ${tumorbam} > ${step_outd}/rm_tum.log 2>&1 || exit 1

    # Create file indicating that execution was finished
    touch ${step_outd}/finished

    display_end_step_message
}
