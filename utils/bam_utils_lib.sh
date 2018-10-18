# *- bash -*

########
init_bash_shebang_var()
{
    echo "#!${BASH}"
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
    local_name=$1
    local_command=$2
    local_script_pars=$3

    # Save previous file (if any)
    if [ -f ${local_name} ]; then
        cp ${local_name} ${local_name}.previous
    fi
    
    # Write bash shebang
    local_BASH_SHEBANG=`init_bash_shebang_var`
    echo ${local_BASH_SHEBANG} > ${local_name} || return 1

    # Write SLURM commands
    echo "#SBATCH --job-name=${local_command}" >> ${local_name} || return 1
    echo "#SBATCH --output=${local_name}.slurm_out" >> ${local_name} || return 1

    # Write environment variables
    set | exclude_readonly_vars | exclude_bashisms >> ${local_name} || return 1
    
    # Write command to be executed
    echo "${local_command} ${local_script_pars}" >> ${local_name} || return 1

    # Give execution permission
    chmod u+x ${local_name} || return 1

    # Archive script with date info
    curr_date=`date '+%Y_%m_%d'`
    cp ${local_name} ${local_name}.${curr_date}
}

########
get_account_opt()
{
    local_account=$1

    if [ -z "${local_account}" ]; then
        echo ""
    else
        echo "-A ${local_account}"
    fi
}

########
get_partition_opt()
{
    local_partition=$1

    if [ -z "${local_partition}" ]; then
        echo ""
    else
        echo "--partition=${local_partition}"
    fi
}

########
find_dependency_for_step()
{
    local_jobdeps_spec=$1
    local_stepname_part=$2

    for local_dep in `echo ${local_jobdeps_spec} | $SED 's/,/ /g'`; do
        local_stepname_part_in_dep=`echo ${local_dep} | $AWK -F ":" '{print $2}'`
        if [ ${local_stepname_part_in_dep} = ${local_stepname_part} ]; then
            echo ${local_dep}
        fi
    done
}

########
get_outd_for_dep()
{
    local_outd=$1
    local_dep=$2

    if [ -z "${local_dep}" ]; then
        echo ""
    else
        local_stepname_part=`echo ${local_dep} | $AWK -F ":" '{print $2}'`
        get_step_dirname ${local_outd} ${local_stepname_part}
    fi
}

########
apply_deptype_to_jobids()
{
    # Initialize variables
    local_jids=$1
    local_deptype=$2

    # Apply deptype
    result=""
    for jid in `echo ${local_jids} | $SED 's/,/ /g'`; do
        if [ -z "" ]; then
            result=${local_deptype}:${jid}
        else
            result=${result}","${local_deptype}:${jid}
        fi
    done

    echo $result
}

########
get_slurm_dependency_opt()
{
    local_jobdeps=$1

    # Create dependency option
    if [ -z "${local_jobdeps}" ]; then
        echo ""
    else
        echo "--dependency=${local_jobdeps}"
    fi
}

########
launch()
{
    local_file=$1
    local_account=$2
    local_partition=$3
    local_cpus=$4
    local_mem=$5
    local_time=$6
    local_jobdeps=$7
    local_outvar=$8
    
    if [ -z "${SBATCH}" ]; then
        ${local_file} > ${local_file}.log 2>&1 || return 1
        eval "${local_outvar}=\"\""
    else
        account_opt=`get_account_opt ${local_account}`
        partition_opt=`get_partition_opt ${local_partition}`
        dependency_opt=`get_slurm_dependency_opt "${local_jobdeps}"`
        local_jid=$($SBATCH --cpus-per-task=${local_cpus} --mem=${local_mem} --time ${local_time} --parsable ${account_opt} ${partition_opt} ${dependency_opt} ${local_file})
        eval "${local_outvar}='${local_jid}'"
    fi
}

########
launch_step()
{
    # Initialize variables
    local_stepname=$1
    local_stepinfo=$2
    local_jobdeps=$3
    local_script_pars=$4
    local_jid=$5

    # Create script
    create_script ${tmpdir}/scripts/${local_stepname} ${local_stepname} "${local_script_pars}" || return 1

    # Retrieve requirements
    local_account=`extract_account_from_entry "$local_stepinfo"`
    local_partition=`extract_partition_from_entry "$local_stepinfo"`
    local_cpus=`extract_cpus_from_entry "$local_stepinfo"`
    local_mem=`extract_mem_from_entry "$local_stepinfo"`
    local_time=`extract_time_from_entry "$local_stepinfo"`

    # Launch script
    launch ${tmpdir}/scripts/${local_stepname} ${local_account} ${local_partition} ${local_cpus} ${local_mem} ${local_time} "${local_jobdeps}" ${local_jid} || return 1
}

########
get_step_info()
{
    local_stepname=$1
    local_infofile=$2

    $AWK -v stepname=${local_stepname} '{if($1==stepname) print $0}' ${local_infofile}
}

########
analysis_entry_is_ok()
{
    local_entry=$1
    echo "${local_entry}" | $AWK '{if(NF>=4) print"yes\n"; else print"no\n"}'
}

########
extract_stepname_from_entry()
{
    local_entry=$1
    echo "${local_entry}" | $AWK '{print $1}'
}

########
extract_account_from_entry()
{
    local_entry=$1
    echo "${local_entry}" | $AWK '{print $2}'
}

########
extract_partition_from_entry()
{
    local_entry=$1
    echo "${local_entry}" | $AWK '{print $3}'
}

########
extract_cpus_from_entry()
{
    local_entry=$1
    echo "${local_entry}" | $AWK '{print $4}'
}

########
extract_mem_from_entry()
{
    local_entry=$1
    echo "${local_entry}" | $AWK '{print $5}'
}

########
extract_time_from_entry()
{
    local_entry=$1
    echo "${local_entry}" | $AWK '{print $6}'
}

########
extract_jobdeps_spec_from_entry()
{
    local_entry=$1
    echo "${local_entry}" | $AWK '{print substr($7,9)}'
}

########
get_step_dirname()
{
    local_dirname=$1
    local_stepname=$2
    echo ${local_dirname}/${local_stepname}
}

########
reset_outdir_for_step() 
{
    local_dirname=$1
    local_stepname=$2
    local_outd=`get_step_dirname ${local_dirname} ${local_stepname}`

    if [ -d ${local_outd} ]; then
        echo "Warning: ${local_stepname} output directory already exists but analysis was not finished, directory content will be removed">&2
        rm -rf ${local_outd}/* || { echo "Error! could not clear output directory" >&2; return 1; }
    else
        mkdir ${local_outd} || { echo "Error! cannot create output directory" >&2; return 1; }
    fi
}

########
get_step_status()
{
    local_dirname=$1
    local_stepname=$2
    local_stepdirname=`get_step_dirname ${local_dirname} ${local_stepname}`
    
    if [ -d ${local_stepdirname} ]; then
        if [ -f ${local_stepdirname}/finished ]; then
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
execute_manta_somatic()
{
    display_begin_step_message

    # Initialize variables
    local_ref=$1
    local_normalbam=$2
    local_tumorbam=$3
    local_step_outd=$4
    local_cpus=$5
    
    # Activate conda environment
    conda activate manta 2> ${local_step_outd}/conda_activate.log || exit 1
    
    # Configure Manta
    configManta.py --normalBam ${local_normalbam} --tumorBam ${local_tumorbam} --referenceFasta ${local_ref} --runDir ${local_step_outd} > ${local_step_outd}/configManta.log 2>&1 || exit 1

    # Execute Manta
    ${local_step_outd}/runWorkflow.py -m local -j ${local_cpus} > ${local_step_outd}/runWorkflow.log 2>&1 || exit 1

    # Deactivate conda environment
    conda deactivate

    # Create file indicating that execution was finished
    touch ${local_step_outd}/finished

    display_end_step_message
}

########
get_indel_cand_opt()
{
    local_manta_outd=$1

    if [ -z "${local_manta_outd}" ]; then
        echo ""
    else
        manta_indel_file="${local_manta_outd}/results/variants/candidateSmallIndels.vcf.gz"
        if [ -f ${manta_indel_file} ]; then
            echo "--indelCandidates ${manta_indel_file}"
        else
            echo "WARNING: Manta indel file for Strelka not found! (${manta_indel_file})" >&2
            echo ""
        fi
    fi
}

########
execute_strelka_somatic()
{
    display_begin_step_message

    # Initialize variables
    local_ref=$1
    local_normalbam=$2
    local_tumorbam=$3
    local_step_outd=$4
    local_manta_outd=$5
    local_cpus=$6

    # Define --indelCandidates option if output from Manta is available
    indel_cand_opt=`get_indel_cand_opt "${local_manta_outd}"`

    # Activate conda environment
    conda activate strelka 2> ${local_step_outd}/conda_activate.log || exit 1

    # Configure Strelka
    configureStrelkaSomaticWorkflow.py --normalBam ${local_normalbam} --tumorBam ${local_tumorbam} --referenceFasta ${local_ref} ${indel_cand_opt} --runDir ${local_step_outd} > ${local_step_outd}/configureStrelkaSomaticWorkflow.log 2>&1 || exit 1

    # Execute Strelka
    ${local_step_outd}/runWorkflow.py -m local -j ${local_cpus} > ${local_step_outd}/runWorkflow.log 2>&1 || exit 1

    # Deactivate conda environment
    conda deactivate

    # Create file indicating that execution was finished
    touch ${local_step_outd}/finished

    display_end_step_message
}

########
execute_msisensor()
{
    display_begin_step_message

    # Initialize variables
    local_ref=$1
    local_normalbam=$2
    local_tumorbam=$3
    local_step_outd=$4
    local_cpus=$5

    # Activate conda environment
    conda activate msisensor 2> ${local_step_outd}/conda_activate.log || exit 1

    # Create homopolymer and microsatellites file
    msisensor scan -d ${local_ref} -o ${local_step_outd}/msisensor.list > ${local_step_outd}/msisensor_scan.log 2>&1 || exit 1

    # Run MSIsensor analysis
    msisensor msi -d ${local_step_outd}/msisensor.list -n ${local_normalbam} -t ${local_tumorbam} -o ${local_step_outd}/output -l 1 -q 1 -b ${local_cpus} > ${local_step_outd}/msisensor_msi.log 2>&1 || exit 1

    # Deactivate conda environment
    conda deactivate
    
    # Create file indicating that execution was finished
    touch ${local_step_outd}/finished

    display_end_step_message
}

########
execute_platypus_germline_conda()
{
    display_begin_step_message
    
    # Initialize variables
    local_ref=$1
    local_normalbam=$2
    local_step_outd=$3

    # Activate conda environment
    conda activate platypus 2> ${local_step_outd}/conda_activate.log || exit 1

    # Run Platypus
    Platypus.py callVariants --bamFiles=${local_normalbam} --refFile=${local_ref} --output=${local_step_outd}/output.vcf --verbosity=1 > ${local_step_outd}/platypus.log 2>&1 || exit 1

    # Deactivate conda environment
    conda deactivate

    # Create file indicating that execution was finished
    touch ${local_step_outd}/finished        

    display_end_step_message
}

########
execute_platypus_germline_local()
{
    display_begin_step_message

    # Initialize variables
    local_ref=$1
    local_normalbam=$2
    local_step_outd=$3

    # Run Platypus
    python ${PLATYPUS_HOME_DIR}/bin/Platypus.py callVariants --bamFiles=${local_normalbam} --refFile=${local_ref} --output=${local_step_outd}/output.vcf --verbosity=1 > ${local_step_outd}/platypus.log 2>&1 || exit 1

    # Create file indicating that execution was finished
    touch ${local_step_outd}/finished    

    display_end_step_message
}

########
execute_platypus_germline()
{
    # Initialize variables
    local_ref=$1
    local_normalbam=$2
    local_step_outd=$3

    if [ -z "${PLATYPUS_HOME_DIR}" ]; then
        execute_platypus_germline_conda ${local_ref} ${local_normalbam} ${local_step_outd}
    else
        execute_platypus_germline_local ${local_ref} ${local_normalbam} ${local_step_outd}
    fi
}

########
execute_cnvkit()
{
    display_begin_step_message

    # Initialize variables
    local_ref=$1
    local_normalbam=$2
    local_tumorbam=$3
    local_step_outd=$4
    local_cpus=$5
    
    # Activate conda environment
    conda activate cnvkit 2> ${local_step_outd}/conda_activate.log || exit 1

    # Run cnvkit
    cnvkit.py batch ${local_tumorbam} -n ${local_normalbam} -m wgs -f ${local_ref}  -d ${local_step_outd} -p ${local_cpus} > ${local_step_outd}/cnvkit.log 2>&1 || exit 1

    # Deactivate conda environment
    conda deactivate

    # Create file indicating that execution was finished
    touch ${local_step_outd}/finished

    display_end_step_message
}

########
execute_wisecondorx()
{
    display_begin_step_message

    # Initialize variables
    local_wcref=$1
    local_tumorbam=$2
    local_step_outd=$3
    local_cpus=$4
    
    # Activate conda environment
    conda activate wisecondorx 2> ${local_step_outd}/conda_activate.log || exit 1

    # Convert tumor bam file into npz
    BINSIZE=5000
    WisecondorX convert ${local_tumorbam} ${local_step_outd}/tumor.npz --binsize $BINSIZE > ${local_step_outd}/wisecondorx_convert.log 2>&1 || exit 1
    
    # Use WisecondorX for prediction
    WisecondorX predict ${local_step_outd}/tumor.npz ${local_wcref} ${local_step_outd}/out > ${local_step_outd}/wisecondorx_predict.log 2>&1 || exit 1

    # Deactivate conda environment
    conda deactivate

    # Create file indicating that execution was finished
    touch ${local_step_outd}/finished

    display_end_step_message
}

########
execute_facets()
{
    display_begin_step_message

    # Initialize variables
    local_normalbam=$1
    local_tumorbam=$2
    local_snpvcf=$3
    local_step_outd=$4

    # Execute snp-pileup
    ${FACETS_HOME_DIR}/inst/extcode/snp-pileup ${snpvcf} ${local_step_outd}/snp-pileup-counts.csv ${local_normalbam} ${local_tumorbam} > ${local_step_outd}/snp-pileup.log 2>&1 || exit 1
    
    # Execute facets
    ${bindir}/run_facets -c ${local_step_outd}/snp-pileup-counts.csv > ${local_step_outd}/facets.out 2> ${local_step_outd}/run_facets.log || exit 1

    # Create file indicating that execution was finished
    touch ${local_step_outd}/finished

    display_end_step_message
}

########
execute_ascatngs()
{
    display_begin_step_message

    # Initialize variables
    local_ref=$1
    local_normalbam=$2
    local_tumorbam=$3
    local_gender=$4
    local_malesexchr=$5
    local_snpgccorr=$6
    local_step_outd=$7
    local_cpus=$8
    
    # Activate conda environment
    conda activate ascatngs 2> ${local_step_outd}/conda_activate.log || exit 1

    # Run cnvkit
    ascat.pl -n ${local_normalbam} -t ${local_tumbam} -r ${local_ref} -sg ${local_snpgccorr} -pr WGS -g ${local_gender} -gc ${local_malesexchr} -cpus ${local_cpus} > ${local_step_outd}/ascat.log 2>&1 || exit 1

    # Deactivate conda environment
    conda deactivate

    # Create file indicating that execution was finished
    touch ${local_step_outd}/finished

    display_end_step_message
}

########
execute_download_ega_norm_bam()
{
    display_begin_step_message

    # Initialize variables
    local_normalbam=$1
    local_egaid_normalbam=$2
    local_egastr=$3
    local_egacred=$4
    local_step_outd=$5

    # Activate conda environment
    conda activate pyega3 2> ${local_step_outd}/conda_activate.log || exit 1

    # Download file
    pyega3 -c ${local_egastr} -cf ${local_egacred} fetch ${local_egaid_normalbam} ${local_step_outd}/normal.bam > ${local_step_outd}/pyega3.log 2>&1 || exit 1

    # Move file
    mv ${local_step_outd}/normal.bam ${local_normalbam}
    
    # Deactivate conda environment
    conda deactivate

    # Create file indicating that execution was finished
    touch ${local_step_outd}/finished

    display_end_step_message
}

########
execute_download_ega_tum_bam()
{
    display_begin_step_message

    # Initialize variables
    local_tumorbam=$1
    local_egaid_tumorbam=$2
    local_egastr=$3
    local_egacred=$4
    local_step_outd=$5

    # Activate conda environment
    conda activate pyega3 2> ${local_step_outd}/conda_activate.log || exit 1

    # Download file
    pyega3 -c ${local_egastr} -cf ${local_egacred} fetch ${local_egaid_tumorbam} ${local_step_outd}/tumor.bam > ${local_step_outd}/pyega3.log 2>&1 || exit 1

    # Move file
    mv ${local_step_outd}/tumor.bam ${local_tumorbam}

    # Deactivate conda environment
    conda deactivate

    # Create file indicating that execution was finished
    touch ${local_step_outd}/finished

    display_end_step_message
}

########
find_bam_filename()
{
    local_step_outd=$1
    local_result=""
    
    for f in ${local_step_outd}/*.bam; do
        if [ -f $f ]; then
            local_result=$f
        fi
    done

    echo ${local_result}
}

########
execute_download_aws_norm_bam()
{
    display_begin_step_message

    # Initialize variables
    local_normalbam=$1
    local_icgcid_normalbam=$2
    local_step_outd=$3

    # Download file
    ${ICGCSTOR_HOME_DIR}/bin/icgc-storage-client download --object-id ${local_icgcid_normalbam} --output-dir ${local_step_outd} || exit 1

    # Find bam file name
    local_bam_file_name=`find_bam_filename ${local_step_outd}`
    
    if [ -z "${local_bam_file_name}" ]; then
        echo "Error: bam file not found after download process was completed" >&2
        exit 1
    fi

    # Move file
    mv ${local_bam_file_name} ${local_normalbam} || exit 1

    # Create file indicating that execution was finished
    touch ${local_step_outd}/finished

    display_end_step_message
}

########
execute_download_aws_tum_bam()
{
    display_begin_step_message

    # Initialize variables
    local_tumorbam=$1
    local_icgcid_tumorbam=$2
    local_step_outd=$3

    # Download file
    ${ICGCSTOR_HOME_DIR}/bin/icgc-storage-client download --object-id ${local_icgcid_tumorbam} --output-dir ${local_step_outd} || exit 1

    # Find bam file name
    local_bam_file_name=`find_bam_filename ${local_step_outd}`
    
    if [ -z "${local_bam_file_name}" ]; then
        echo "Error: bam file not found after download process was completed" >&2
        exit 1
    fi

    # Move file
    mv ${local_bam_file_name} ${local_tumorbam} || exit 1

    # Create file indicating that execution was finished
    touch ${local_step_outd}/finished

    display_end_step_message
}

########
execute_index_norm_bam()
{
    display_begin_step_message

    # Initialize variables
    local_normalbam=$1
    local_step_outd=$2

    # Index normal bam file if necessary
    if [ ! -f ${local_normalbam}.bai ]; then

        # Activate conda environment
        conda activate base 2> ${local_step_outd}/conda_activate.log || exit 1

        # Execute samtools
        samtools index ${local_normalbam} > ${local_step_outd}/samtools.log 2>&1 || exit 1

        # Deactivate conda environment
        conda deactivate
    fi
    
    # Create file indicating that execution was finished
    touch ${local_step_outd}/finished

    display_end_step_message
}

########
execute_sort_norm_bam()
{
    display_begin_step_message

    # Initialize variables
    local_normalbam=$1
    local_step_outd=$2
    local_cpus=$3

    # Activate conda environment
    conda activate base 2> ${local_step_outd}/conda_activate.log || exit 1

    # Execute samtools
    samtools sort -T ${local_step_outd} -o ${local_step_outd}/sorted.bam -m 2G -@ ${local_cpus} ${local_normalbam} >  ${local_step_outd}/samtools.log 2>&1 || exit 1
    # NOTE: -m option is used here to increase the maximum memory per
    # thread. One lateral efect of this is that the number of tmp files
    # generated is decreased. This constitutes one possible way to avoid
    # the "Too many open files" error reported by samtools
    
    # Deactivate conda environment
    conda deactivate

    # Replace initial bam file by the sorted one
    mv ${local_step_outd}/sorted.bam ${local_normalbam} 2> ${local_step_outd}/mv.log || exit 1

    # Create file indicating that execution was finished
    touch ${local_step_outd}/finished

    display_end_step_message
}

########
execute_sort_tum_bam()
{
    display_begin_step_message

    # Initialize variables
    local_tumorbam=$1
    local_step_outd=$2
    local_cpus=$3

    # Activate conda environment
    conda activate base 2> ${local_step_outd}/conda_activate.log || exit 1

    # Execute samtools
    samtools sort -T ${local_step_outd} -o ${local_step_outd}/sorted.bam -m 2G -@ ${local_cpus} ${local_tumorbam} >  ${local_step_outd}/samtools.log 2>&1 || exit 1
    # NOTE: -m option is used here to increase the maximum memory per
    # thread. One lateral efect of this is that the number of tmp files
    # generated is decreased. This constitutes one possible way to avoid
    # the "Too many open files" error reported by samtools
    
    # Deactivate conda environment
    conda deactivate

    # Replace initial bam file by the sorted one
    mv ${local_step_outd}/sorted.bam ${local_tumorbam} 2> ${local_step_outd}/mv.log || exit 1

    # Create file indicating that execution was finished
    touch ${local_step_outd}/finished

    display_end_step_message
}

########
execute_index_tum_bam()
{
    display_begin_step_message

    # Initialize variables
    local_tumorbam=$1
    local_step_outd=$2

    # Index tumor bam file if necessary
    if [ ! -f ${local_tumorbam}.bai ]; then

        # Activate conda environment
        conda activate base 2> ${local_step_outd}/conda_activate.log || exit 1

        # Execute samtools
        samtools index ${local_tumorbam} > ${local_step_outd}/samtools.log 2>&1 || exit 1

        # Deactivate conda environment
        conda deactivate
    fi
    
    # Create file indicating that execution was finished
    touch ${local_step_outd}/finished

    display_end_step_message
}

########
execute_delete_bam_files()
{
    display_begin_step_message

    # Initialize variables
    local_normalbam=$1
    local_tumorbam=$2
    local_step_outd=$3

    # Delete normal bam file
    rm ${local_normalbam} > ${local_step_outd}/rm_norm.log 2>&1 || exit 1
    
    # Delete tumor bam file
    rm ${local_tumorbam} > ${local_step_outd}/rm_tum.log 2>&1 || exit 1

    # Create file indicating that execution was finished
    touch ${local_step_outd}/finished

    display_end_step_message
}
