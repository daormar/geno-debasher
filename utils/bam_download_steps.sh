# *- bash -*

#############
# CONSTANTS #
#############

DEFAULT_NUMBER_OF_DOWNLOAD_TRIES=5
DEFAULT_NUMBER_OF_EGA_DOWNLOAD_STREAMS=50
DEFAULT_NUMBER_OF_GDC_DOWNLOAD_PROCS=50
DEFAULT_ASP_MAX_TRANS_RATE=100m

######################
# BAM DOWNLOAD STEPS #
######################

########
copy_norm_bam_explain_cmdline_opts()
{
    # -extn option
    description="Path to local normal bam file to be copied"
    explain_cmdline_opt "-extn" "<string>" "$description"
}


########
copy_norm_bam_define_opts()
{
    # Initialize variables
    local cmdline=$1
    local stepspec=$2
    local optlist=""

    # Define the -step-outd option, the output directory for the step
    local step_outd=`get_step_outdir_given_stepspec "$stepspec"`
    define_opt "-step-outd" ${step_outd} optlist || exit 1

    # -extn option
    define_cmdline_opt "$cmdline" "-extn" optlist || exit 1

    # -normalbam option
    local abs_datadir=`get_absolute_shdirname ${DATADIR_BASENAME}`
    local normalbam=${abs_datadir}/normal.bam
    define_opt "-normalbam" $normalbam optlist || exit 1

    # Save option list
    save_opt_list optlist    
}

########
copy_norm_bam()
{
    # Initialize variables
    local normalbam=`read_opt_value_from_line "$*" "-normalbam"`
    local local_normalbam=`read_opt_value_from_line "$*" "-extn"`
    local step_outd=`read_opt_value_from_line "$*" "-step-outd"`

    # Copy file
    logmsg "* Copying file..."
    cp ${local_normalbam} ${normalbam} || exit 1
}

########
copy_tum_bam_explain_cmdline_opts()
{
    # -extt option
    description="Path to local tumor bam file to be copied"
    explain_cmdline_req_opt "-extt" "<string>" "$description"
}

########
copy_tum_bam_define_opts()
{
    # Initialize variables
    local cmdline=$1
    local stepspec=$2
    local optlist=""

    # Define the -step-outd option, the output directory for the step
    local step_outd=`get_step_outdir_given_stepspec "$stepspec"`
    define_opt "-step-outd" ${step_outd} optlist || exit 1

    # -extt option
    define_cmdline_opt "$cmdline" "-extt" optlist || exit 1

    # -tumorbam option
    local abs_datadir=`get_absolute_shdirname ${DATADIR_BASENAME}`
    local tumorbam=${abs_datadir}/tumor.bam
    define_opt "-tumorbam" $tumorbam optlist || exit 1

    # Save option list
    save_opt_list optlist    
}

########
copy_tum_bam()
{
    # Initialize variables
    local tumorbam=`read_opt_value_from_line "$*" "-tumorbam"`
    local local_tumorbam=`read_opt_value_from_line "$*" "-extt"`
    local step_outd=`read_opt_value_from_line "$*" "-step-outd"`

    # Copy file
    logmsg "* Copying file..."
    cp ${local_tumorbam} ${tumorbam} || exit 1
}

########
download_ega_norm_bam_explain_cmdline_opts()
{
    # -extn option
    description="External database id of normal bam file to download"
    explain_cmdline_opt "-extn" "<string>" "$description"

    # -egastr option
    description="Number of streams used by the EGA download client (${DEFAULT_NUMBER_OF_EGA_DOWNLOAD_STREAMS} by default)"
    explain_cmdline_opt "-egastr" "<int>" "$description"

    # -egacred option
    description="File with EGA download client credentials"
    explain_cmdline_req_opt "-egacred" "<string>" "$description"

    # -nt option
    description="Number of download tries per file (${DEFAULT_NUMBER_OF_DOWNLOAD_TRIES} by default)"
    explain_cmdline_opt "-nt" "<int>" "$description"
}


########
download_ega_norm_bam_define_opts()
{
    # Initialize variables
    local cmdline=$1
    local stepspec=$2
    local optlist=""

    # Define the -step-outd option, the output directory for the step
    local step_outd=`get_step_outdir_given_stepspec "$stepspec"`
    define_opt "-step-outd" ${step_outd} optlist || exit 1

    # -extn option
    define_cmdline_opt "$cmdline" "-extn" optlist || exit 1

    # -egastr option
    define_cmdline_nonmandatory_opt "$cmdline" "-egastr" ${DEFAULT_NUMBER_OF_EGA_DOWNLOAD_STREAMS} optlist || exit 1

    # -egacred option
    define_cmdline_opt "$cmdline" "-egacred" optlist || exit 1

    # -nt option
    define_cmdline_nonmandatory_opt "$cmdline" "-nt" ${DEFAULT_NUMBER_OF_DOWNLOAD_TRIES} optlist || exit 1

    # -normalbam option
    local abs_datadir=`get_absolute_shdirname ${DATADIR_BASENAME}`
    local normalbam=${abs_datadir}/normal.bam
    define_opt "-normalbam" $normalbam optlist || exit 1

    # Save option list
    save_opt_list optlist    
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
        logmsg "Starting download try number ${ntry}..."

        # Remove previously downloaded file (if any)
        if [ -f ${outf} ]; then
            rm ${outf}
        fi

        # Download file
        pyega3 -c ${egastr} -cf ${egacred} fetch ${egaid} ${outf} 2>&1
        
        # Check if download was successful
        if [ $? -eq 0 -a -f ${outf} ]; then
            return 0
        fi

        ntry=$((ntry+1))
    done

    logmsg "All download attempts failed!"

    return 1
}

########
download_ega_norm_bam()
{
    # Initialize variables
    local normalbam=`read_opt_value_from_line "$*" "-normalbam"`
    local egaid_normalbam=`read_opt_value_from_line "$*" "-extn"`
    local egastr=`read_opt_value_from_line "$*" "-egastr"`
    local egacred=`read_opt_value_from_line "$*" "-egacred"`
    local download_tries=`read_opt_value_from_line "$*" "-nt"`
    local step_outd=`read_opt_value_from_line "$*" "-step-outd"`

    # Activate conda environment
    logmsg "* Activating conda environment..."
    conda activate pyega3 2>&1 || exit 1

    # Download file (with multiple tries)
    ega_download_retry ${egastr} ${egacred} ${egaid_normalbam} ${step_outd}/normal.bam ${download_tries} || exit 1

    # Move file
    mv ${step_outd}/normal.bam ${normalbam} || exit 1
    
    # Deactivate conda environment
    logmsg "* Deactivating conda environment..."
    conda deactivate 2>&1
}

########
download_ega_norm_bam_conda_envs()
{
    define_conda_env pyega3 pyega3.yml
}

########
download_ega_tum_bam_explain_cmdline_opts()
{
    # -extt option
    description="External database id of tumor bam file to download"
    explain_cmdline_req_opt "-extt" "<string>" "$description"

    # -egastr option
    description="Number of streams used by the EGA download client (${DEFAULT_NUMBER_OF_EGA_DOWNLOAD_STREAMS} by default)"
    explain_cmdline_opt "-egastr" "<int>" "$description"

    # -egacred option
    description="File with EGA download client credentials"
    explain_cmdline_req_opt "-egacred" "<string>" "$description"

    # -nt option
    description="Number of download tries per file (${DEFAULT_NUMBER_OF_DOWNLOAD_TRIES} by default)"
    explain_cmdline_opt "-nt" "<int>" "$description"
}

########
download_ega_tum_bam_define_opts()
{
    # Initialize variables
    local cmdline=$1
    local stepspec=$2
    local optlist=""

    # Define the -step-outd option, the output directory for the step
    local step_outd=`get_step_outdir_given_stepspec "$stepspec"`
    define_opt "-step-outd" ${step_outd} optlist || exit 1

    # -extt option
    define_cmdline_opt "$cmdline" "-extt" optlist || exit 1

    # -egastr option
    define_cmdline_nonmandatory_opt "$cmdline" "-egastr" ${DEFAULT_NUMBER_OF_EGA_DOWNLOAD_STREAMS} optlist || exit 1

    # -egacred option
    define_cmdline_opt "$cmdline" "-egacred" optlist || exit 1

    # -nt option
    define_cmdline_nonmandatory_opt "$cmdline" "-nt" ${DEFAULT_NUMBER_OF_DOWNLOAD_TRIES} optlist || exit 1

    # -tumorbam option
    local abs_datadir=`get_absolute_shdirname ${DATADIR_BASENAME}`
    local tumorbam=${abs_datadir}/tumor.bam
    define_opt "-tumorbam" $tumorbam optlist || exit 1

    # Save option list
    save_opt_list optlist    
}

########
download_ega_tum_bam()
{
    # Initialize variables
    local tumorbam=`read_opt_value_from_line "$*" "-tumorbam"`
    local egaid_tumorbam=`read_opt_value_from_line "$*" "-extt"`
    local egastr=`read_opt_value_from_line "$*" "-egastr"`
    local egacred=`read_opt_value_from_line "$*" "-egacred"`
    local download_tries=`read_opt_value_from_line "$*" "-nt"`
    local step_outd=`read_opt_value_from_line "$*" "-step-outd"`

    # Activate conda environment
    logmsg "* Activating conda environment..."
    conda activate pyega3 2>&1 || exit 1

    # Download file (with multiple tries)
    ega_download_retry ${egastr} ${egacred} ${egaid_tumorbam} ${step_outd}/tumor.bam ${download_tries} || exit 1

    # Move file
    mv ${step_outd}/tumor.bam ${tumorbam} || exit 1

    # Deactivate conda environment
    logmsg "* Deactivating conda environment..."
    conda deactivate 2>&1
}

########
download_ega_tum_bam_conda_envs()
{
    define_conda_env pyega3 pyega3.yml
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
download_ega_asp_norm_bam_explain_cmdline_opts()
{
    # -extn option
    description="External database id of normal bam file to download"
    explain_cmdline_req_opt "-extn" "<string>" "$description"

    # -asperausr option
    description="Username for Aspera server"
    explain_cmdline_req_opt "-asperausr" "<string>" "$description"

    # -asperapwd option
    description="Password for Aspera server"
    explain_cmdline_req_opt "-asperapwd" "<string>" "$description"

    # -asperaserv option
    description="Name of Aspera server"
    explain_cmdline_req_opt "-asperaserv" "<string>" "$description"

    # -egadecrpwd option
    description="File with EGA decryptor password"
    explain_cmdline_req_opt "-egadecrpwd" "<string>" "$description"
    
    # -nt option
    description="Number of download tries per file (${DEFAULT_NUMBER_OF_DOWNLOAD_TRIES} by default)"
    explain_cmdline_opt "-nt" "<int>" "$description"
}

########
download_ega_asp_norm_bam_define_opts()
{
    # Initialize variables
    local cmdline=$1
    local stepspec=$2
    local optlist=""

    # Define the -step-outd option, the output directory for the step
    local step_outd=`get_step_outdir_given_stepspec "$stepspec"`
    define_opt "-step-outd" ${step_outd} optlist || exit 1
    
    # -extn option
    define_cmdline_opt "$cmdline" "-extn" optlist || exit 1

    # -asperausr option
    define_cmdline_opt "$cmdline" "-asperausr" optlist || exit 1

    # -asperapwd option
    define_cmdline_opt "$cmdline" "-asperapwd" optlist || exit 1

    # -asperaserv option
    define_cmdline_opt "$cmdline" "-asperaserv" optlist || exit 1

    # -egadecrpwd option
    define_cmdline_opt "$cmdline" "-egadecrpwd" optlist || exit 1

    # -nt option
    define_cmdline_nonmandatory_opt "$cmdline" "-nt" ${DEFAULT_NUMBER_OF_DOWNLOAD_TRIES} optlist || exit 1

    # -normalbam option
    local abs_datadir=`get_absolute_shdirname ${DATADIR_BASENAME}`
    local normalbam=${abs_datadir}/normal.bam
    define_opt "-normalbam" $normalbam optlist || exit 1

    # Save option list
    save_opt_list optlist
}

########
download_ega_asp_norm_bam()
{
    # Initialize variables
    local normalbam=`read_opt_value_from_line "$*" "-normalbam"`
    local normalbam_file=`read_opt_value_from_line "$*" "-extn"`
    local aspera_user=`read_opt_value_from_line "$*" "-asperausr"`
    local aspera_passwd=`read_opt_value_from_line "$*" "-asperapwd"`
    local aspera_server=`read_opt_value_from_line "$*" "-asperaserv"`
    local egadecrypt_pwd=`read_opt_value_from_line "$*" "-egadecrpwd"`
    local download_tries=`read_opt_value_from_line "$*" "-nt"`
    local step_outd=`read_opt_value_from_line "$*" "-step-outd"`
    local max_trans_rate=${DEFAULT_ASP_MAX_TRANS_RATE}
    
    # Download file
    logmsg "* Executing ascp (${normalbam_file})..."
    ASPERA_SCP_PASS=${aspera_passwd} ${ASPERA_HOME_DIR}/bin/ascp --ignore-host-key -QTl ${max_trans_rate} ${aspera_user}@${aspera_server}:${normalbam_file} ${step_outd}/normal.bam.crypt 2>&1 || exit 1

    # Decrypt file
    logmsg "* Executing decryptor.jar..."
    $JAVA -jar ${EGADECRYPT_HOME_DIR}/decryptor.jar ${egadecrypt_pwd} ${step_outd}/normal.bam.crypt 2>&1 || exit 1
    
    # Obtain file name
    local bam_file_name=`find_bam_filename ${step_outd}`
    
    if [ -z "${bam_file_name}" ]; then
        logmsg "Error: bam file not found after download process was completed"
        exit 1
    fi

    # Move file
    mv ${bam_file_name} ${normalbam} || exit 1

    # Remove encrypted file
    rm ${step_outd}/normal.bam.crypt || exit 1
}

########
download_ega_asp_tum_bam_explain_cmdline_opts()
{
    # -extt option
    description="External database id of tumor bam file to download"
    explain_cmdline_req_opt "-extt" "<string>" "$description"

    # -asperausr option
    description="Username for Aspera server"
    explain_cmdline_req_opt "-asperausr" "<string>" "$description"

    # -asperapwd option
    description="Password for Aspera server"
    explain_cmdline_req_opt "-asperapwd" "<string>" "$description"

    # -asperaserv option
    description="Name of Aspera server"
    explain_cmdline_req_opt "-asperaserv" "<string>" "$description"

    # -egadecrpwd option
    description="File with EGA decryptor password"
    explain_cmdline_req_opt "-egadecrpwd" "<string>" "$description"
    
    # -nt option
    description="Number of download tries per file (${DEFAULT_NUMBER_OF_DOWNLOAD_TRIES} by default)"
    explain_cmdline_opt "-nt" "<int>" "$description"
}

########
download_ega_asp_tum_bam_define_opts()
{
    # Initialize variables
    local cmdline=$1
    local stepspec=$2
    local optlist=""

    # Define the -step-outd option, the output directory for the step
    local step_outd=`get_step_outdir_given_stepspec "$stepspec"`
    define_opt "-step-outd" ${step_outd} optlist || exit 1
    
    # -extt option
    define_cmdline_opt "$cmdline" "-extt" optlist || exit 1

    # -asperausr option
    define_cmdline_opt "$cmdline" "-asperausr" optlist || exit 1

    # -asperapwd option
    define_cmdline_opt "$cmdline" "-asperapwd" optlist || exit 1

    # -asperaserv option
    define_cmdline_opt "$cmdline" "-asperaserv" optlist || exit 1

    # -egadecrpwd option
    define_cmdline_infile_opt "$cmdline" "-egadecrpwd" optlist || exit 1

    # -nt option
    define_cmdline_nonmandatory_opt "$cmdline" "-nt" ${DEFAULT_NUMBER_OF_DOWNLOAD_TRIES} optlist || exit 1

    # -tumorbam option
    local abs_datadir=`get_absolute_shdirname ${DATADIR_BASENAME}`
    local tumorbam=${abs_datadir}/tumor.bam
    define_opt "-tumorbam" $tumorbam optlist || exit 1

    # Save option list
    save_opt_list optlist
}

########
download_ega_asp_tum_bam()
{
    # Initialize variables
    local tumorbam=`read_opt_value_from_line "$*" "-tumorbam"`
    local tumorbam_file=`read_opt_value_from_line "$*" "-extt"`
    local aspera_user=`read_opt_value_from_line "$*" "-asperausr"`
    local aspera_passwd=`read_opt_value_from_line "$*" "-asperapwd"`
    local aspera_server=`read_opt_value_from_line "$*" "-asperaserv"`
    local egadecrypt_pwd=`read_opt_value_from_line "$*" "-egadecrpwd"`
    local download_tries=`read_opt_value_from_line "$*" "-nt"`
    local step_outd=`read_opt_value_from_line "$*" "-step-outd"`
    local max_trans_rate=${DEFAULT_ASP_MAX_TRANS_RATE}

    # Download file
    logmsg "* Executing ascp (${tumorbam_file})..."
    ASPERA_SCP_PASS=${aspera_passwd} ${ASPERA_HOME_DIR}/bin/ascp --ignore-host-key -QTl ${max_trans_rate} ${aspera_user}@${aspera_server}:${tumorbam_file} ${step_outd}/tumor.bam.crypt 2>&1 || exit 1

    # Decrypt file
    logmsg "* Executing decryptor.jar..."
    $JAVA -jar ${EGADECRYPT_HOME_DIR}/decryptor.jar ${egadecrypt_pwd} ${step_outd}/tumor.bam.crypt 2>&1 || exit 1

    # Obtain file name
    local bam_file_name=`find_bam_filename ${step_outd}`
    
    if [ -z "${bam_file_name}" ]; then
        logmsg "Error: bam file not found after download process was completed"
        exit 1
    fi

    # Move file
    mv ${bam_file_name} ${tumorbam} || exit 1

    # Remove encrypted file
    rm ${step_outd}/tumor.bam.crypt || exit 1
}

########
decrypt_ega_norm_bam_explain_cmdline_opts()
{
    # -extn option
    description="File name of encrypted normal bam file to process"
    explain_cmdline_req_opt "-extn" "<string>" "$description"

    # -egadecrpwd option
    description="File with EGA decryptor password"
    explain_cmdline_req_opt "-egadecrpwd" "<string>" "$description"
}

########
decrypt_ega_norm_bam_define_opts()
{
    # Initialize variables
    local cmdline=$1
    local stepspec=$2
    local optlist=""

    # Define the -step-outd option, the output directory for the step
    local step_outd=`get_step_outdir_given_stepspec "$stepspec"`
    define_opt "-step-outd" ${step_outd} optlist || exit 1
    
    # -extn option
    define_cmdline_opt "$cmdline" "-extn" optlist || exit 1

    # -egadecrpwd option
    define_cmdline_opt "$cmdline" "-egadecrpwd" optlist || exit 1

    # -normalbam option
    local abs_datadir=`get_absolute_shdirname ${DATADIR_BASENAME}`
    local normalbam=${abs_datadir}/normal.bam
    define_opt "-normalbam" $normalbam optlist || exit 1

    # Save option list
    save_opt_list optlist
}

########
decrypt_ega_norm_bam()
{
    # Initialize variables
    local normalbam=`read_opt_value_from_line "$*" "-normalbam"`
    local normalbam_file=`read_opt_value_from_line "$*" "-extn"`
    local egadecrypt_pwd=`read_opt_value_from_line "$*" "-egadecrpwd"`
    local step_outd=`read_opt_value_from_line "$*" "-step-outd"`
    
    # Decrypt file
    logmsg "* Executing decryptor.jar..."
    $JAVA -jar ${EGADECRYPT_HOME_DIR}/decryptor.jar ${egadecrypt_pwd} --output-folder ${step_outd} ${normalbam_file} 2>&1 || exit 1
    
    # Obtain file name
    local bam_file_name=`find_bam_filename ${step_outd}`
    
    if [ -z "${bam_file_name}" ]; then
        logmsg "Error: bam file not found after decryption process was completed"
        exit 1
    fi

    # Move file
    mv ${bam_file_name} ${normalbam} || exit 1
}

########
decrypt_ega_tum_bam_explain_cmdline_opts()
{
    # -extt option
    description="File name of encrypted tumor bam file to process"
    explain_cmdline_req_opt "-extt" "<string>" "$description"

    # -egadecrpwd option
    description="File with EGA decryptor password"
    explain_cmdline_req_opt "-egadecrpwd" "<string>" "$description"
}

########
decrypt_ega_tum_bam_define_opts()
{
    # Initialize variables
    local cmdline=$1
    local stepspec=$2
    local optlist=""

    # Define the -step-outd option, the output directory for the step
    local step_outd=`get_step_outdir_given_stepspec "$stepspec"`
    define_opt "-step-outd" ${step_outd} optlist || exit 1
    
    # -extt option
    define_cmdline_opt "$cmdline" "-extt" optlist || exit 1

    # -egadecrpwd option
    define_cmdline_infile_opt "$cmdline" "-egadecrpwd" optlist || exit 1

    # -tumorbam option
    local abs_datadir=`get_absolute_shdirname ${DATADIR_BASENAME}`
    local tumorbam=${abs_datadir}/tumor.bam
    define_opt "-tumorbam" $tumorbam optlist || exit 1

    # Save option list
    save_opt_list optlist
}

########
decrypt_ega_tum_bam()
{
    # Initialize variables
    local tumorbam=`read_opt_value_from_line "$*" "-tumorbam"`
    local tumorbam_file=`read_opt_value_from_line "$*" "-extt"`
    local egadecrypt_pwd=`read_opt_value_from_line "$*" "-egadecrpwd"`
    local step_outd=`read_opt_value_from_line "$*" "-step-outd"`
    
    # Decrypt file
    logmsg "* Executing decryptor.jar..."
    $JAVA -jar ${EGADECRYPT_HOME_DIR}/decryptor.jar ${egadecrypt_pwd} --output-folder ${step_outd} ${tumorbam_file} 2>&1 || exit 1

    # Obtain file name
    local bam_file_name=`find_bam_filename ${step_outd}`
    
    if [ -z "${bam_file_name}" ]; then
        logmsg "Error: bam file not found after decryption process was completed"
        exit 1
    fi

    # Move file
    mv ${bam_file_name} ${tumorbam} || exit 1
}

########
decsingle_ega_norm_bam_explain_cmdline_opts()
{
    # -extn option
    description="File name of encrypted normal bam file to process"
    explain_cmdline_req_opt "-extn" "<string>" "$description"

    # -decsinglepwd option
    description="Password for bam file to be processed with decSINGLE tool"
    explain_cmdline_req_opt "-decsinglepwd" "<string>" "$description"
}

########
decsingle_ega_norm_bam_define_opts()
{
    # Initialize variables
    local cmdline=$1
    local stepspec=$2
    local optlist=""

    # Define the -step-outd option, the output directory for the step
    local step_outd=`get_step_outdir_given_stepspec "$stepspec"`
    define_opt "-step-outd" ${step_outd} optlist || exit 1
    
    # -extn option
    define_cmdline_opt "$cmdline" "-extn" optlist || exit 1

    # -egadecrpwd option
    define_cmdline_opt "$cmdline" "-decsinglepwd" optlist || exit 1

    # -normalbam option
    local abs_datadir=`get_absolute_shdirname ${DATADIR_BASENAME}`
    local normalbam=${abs_datadir}/normal.bam
    define_opt "-normalbam" $normalbam optlist || exit 1

    # Save option list
    save_opt_list optlist
}

########
decsingle_ega_norm_bam()
{
    # Initialize variables
    local normalbam=`read_opt_value_from_line "$*" "-normalbam"`
    local normalbam_file=`read_opt_value_from_line "$*" "-extn"`
    local decsingle_pwd=`read_opt_value_from_line "$*" "-decsinglepwd"`
    local step_outd=`read_opt_value_from_line "$*" "-step-outd"`

    # Activate conda environment
    logmsg "* Activating conda environment..."
    conda activate java-jre 2>&1 || exit 1

    # Decrypt file
    logmsg "* Executing decryptor.jar..."
    cat ${normalbam_file} | $JAVA -cp ${DECSINGLE_HOME_DIR} decSINGLE ${egadecrypt_pwd} --output-folder ${step_outd} ${normalbam_file} <(echo ${decsingle_pwd}) > ${normalbam} || exit 1

    # Deactivate conda environment
    logmsg "* Deactivating conda environment..."
    conda deactivate 2>&1
}

########
decsingle_ega_norm_bam_conda_envs()
{
    define_conda_env java-jre java-jre.yml
}

########
decsingle_ega_tum_bam_explain_cmdline_opts()
{
    # -extn option
    description="File name of encrypted tumor bam file to process"
    explain_cmdline_req_opt "-extt" "<string>" "$description"

    # -decsinglepwd option
    description="Password for bam file to be processed with decSINGLE tool"
    explain_cmdline_req_opt "-decsinglepwd" "<string>" "$description"
}

########
decsingle_ega_tum_bam_define_opts()
{
    # Initialize variables
    local cmdline=$1
    local stepspec=$2
    local optlist=""

    # Define the -step-outd option, the output directory for the step
    local step_outd=`get_step_outdir_given_stepspec "$stepspec"`
    define_opt "-step-outd" ${step_outd} optlist || exit 1
    
    # -extt option
    define_cmdline_opt "$cmdline" "-extt" optlist || exit 1

    # -egadecrpwd option
    define_cmdline_opt "$cmdline" "-decsinglepwd" optlist || exit 1

    # -tumorbam option
    local abs_datadir=`get_absolute_shdirname ${DATADIR_BASENAME}`
    local tumorbam=${abs_datadir}/tumor.bam
    define_opt "-tumorbam" $tumorbam optlist || exit 1

    # Save option list
    save_opt_list optlist
}

########
decsingle_ega_tum_bam()
{
    # Initialize variables
    local tumorbam=`read_opt_value_from_line "$*" "-tumorbam"`
    local tumorbam_file=`read_opt_value_from_line "$*" "-extt"`
    local decsingle_pwd=`read_opt_value_from_line "$*" "-decsinglepwd"`
    local step_outd=`read_opt_value_from_line "$*" "-step-outd"`

    # Activate conda environment
    logmsg "* Activating conda environment..."
    conda activate java-jre 2>&1 || exit 1

    # Decrypt file
    logmsg "* Executing decryptor.jar..."
    cat ${tumorbam_file} | $JAVA -cp ${DECSINGLE_HOME_DIR} decSINGLE ${egadecrypt_pwd} --output-folder ${step_outd} ${tumorbam_file} <(echo ${decsingle_pwd}) > ${tumorbam} || exit 1

    # Deactivate conda environment
    logmsg "* Deactivating conda environment..."
    conda deactivate 2>&1
}

########
decsingle_ega_tum_bam_conda_envs()
{
    define_conda_env java-jre java-jre.yml
}

########
download_aws_norm_bam_explain_cmdline_opts()
{
    # -extn option
    description="External database id of normal bam file to download"
    explain_cmdline_req_opt "-extn" "<string>" "$description"

    # -nt option
    description="Number of download tries per file (${DEFAULT_NUMBER_OF_DOWNLOAD_TRIES} by default)"
    explain_cmdline_opt "-nt" "<int>" "$description"
}

########
download_aws_norm_bam_define_opts()
{
    # Initialize variables
    local cmdline=$1
    local stepspec=$2
    local optlist=""

    # Define the -step-outd option, the output directory for the step
    local step_outd=`get_step_outdir_given_stepspec "$stepspec"`
    define_opt "-step-outd" ${step_outd} optlist || exit 1

    # -extn option
    define_cmdline_opt "$cmdline" "-extn" optlist || exit 1

    # -nt option
    define_cmdline_nonmandatory_opt "$cmdline" "-nt" ${DEFAULT_NUMBER_OF_DOWNLOAD_TRIES} optlist || exit 1

    # -normalbam option
    local abs_datadir=`get_absolute_shdirname ${DATADIR_BASENAME}`
    local normalbam=${abs_datadir}/normal.bam
    define_opt "-normalbam" $normalbam optlist || exit 1

    # Save option list
    save_opt_list optlist    
}

########
download_aws_norm_bam()
{
    # Initialize variables
    local normalbam=`read_opt_value_from_line "$*" "-normalbam"`
    local icgcid_normalbam=`read_opt_value_from_line "$*" "-extn"`
    local download_tries=`read_opt_value_from_line "$*" "-nt"`
    local step_outd=`read_opt_value_from_line "$*" "-step-outd"`

    # Download file
    logmsg "* Executing score-client..."
    ${ICGCSTOR_HOME_DIR}/bin/score-client --profile aws download --object-id ${icgcid_normalbam} --output-dir ${step_outd} 2>&1 || exit 1

    # Find bam file name
    local bam_file_name=`find_bam_filename ${step_outd}`
    
    if [ -z "${bam_file_name}" ]; then
        logmsg "Error: bam file not found after download process was completed"
        exit 1
    fi

    # Move file
    mv ${bam_file_name} ${normalbam} || exit 1
}

########
download_aws_tum_bam_explain_cmdline_opts()
{
    # -extt option
    description="External database id of tumor bam file to download"
    explain_cmdline_req_opt "-extt" "<string>" "$description"

    # -nt option
    description="Number of download tries per file (${DEFAULT_NUMBER_OF_DOWNLOAD_TRIES} by default)"
    explain_cmdline_opt "-nt" "<int>" "$description"
}

########
download_aws_tum_bam_define_opts()
{
    # Initialize variables
    local cmdline=$1
    local stepspec=$2
    local optlist=""

    # Define the -step-outd option, the output directory for the step
    local step_outd=`get_step_outdir_given_stepspec "$stepspec"`
    define_opt "-step-outd" ${step_outd} optlist || exit 1

    # -extt option
    define_cmdline_opt "$cmdline" "-extt" optlist || exit 1

    # -nt option
    define_cmdline_nonmandatory_opt "$cmdline" "-nt" ${DEFAULT_NUMBER_OF_DOWNLOAD_TRIES} optlist || exit 1

    # -tumorbam option
    local abs_datadir=`get_absolute_shdirname ${DATADIR_BASENAME}`
    local tumorbam=${abs_datadir}/tumor.bam
    define_opt "-tumorbam" $tumorbam optlist || exit 1

    # Save option list
    save_opt_list optlist    
}

########
download_aws_tum_bam()
{
    # Initialize variables
    local tumorbam=`read_opt_value_from_line "$*" "-tumorbam"`
    local icgcid_tumorbam=`read_opt_value_from_line "$*" "-extt"`
    local download_tries=`read_opt_value_from_line "$*" "-nt"`
    local step_outd=`read_opt_value_from_line "$*" "-step-outd"`

    # Download file
    logmsg "* Executing score-client..."
    ${ICGCSTOR_HOME_DIR}/bin/score-client --profile aws download --object-id ${icgcid_tumorbam} --output-dir ${step_outd} 2>&1 || exit 1

    # Find bam file name
    local bam_file_name=`find_bam_filename ${step_outd}`
    
    if [ -z "${bam_file_name}" ]; then
        logmsg "Error: bam file not found after download process was completed"
        exit 1
    fi

    # Move file
    mv ${bam_file_name} ${tumorbam} || exit 1
}

########
download_collab_norm_bam_explain_cmdline_opts()
{
    # -extn option
    description="External database id of normal bam file to download"
    explain_cmdline_req_opt "-extn" "<string>" "$description"

    # -nt option
    description="Number of download tries per file (${DEFAULT_NUMBER_OF_DOWNLOAD_TRIES} by default)"
    explain_cmdline_opt "-nt" "<int>" "$description"
}

########
download_collab_norm_bam_define_opts()
{
    # Initialize variables
    local cmdline=$1
    local stepspec=$2
    local optlist=""

    # Define the -step-outd option, the output directory for the step
    local step_outd=`get_step_outdir_given_stepspec "$stepspec"`
    define_opt "-step-outd" ${step_outd} optlist || exit 1

    # -extn option
    define_cmdline_opt "$cmdline" "-extn" optlist || exit 1

    # -nt option
    define_cmdline_nonmandatory_opt "$cmdline" "-nt" ${DEFAULT_NUMBER_OF_DOWNLOAD_TRIES} optlist || exit 1

    # -normalbam option
    local abs_datadir=`get_absolute_shdirname ${DATADIR_BASENAME}`
    local normalbam=${abs_datadir}/normal.bam
    define_opt "-normalbam" $normalbam optlist || exit 1

    # Save option list
    save_opt_list optlist    
}

########
download_collab_norm_bam()
{
    # Initialize variables
    local normalbam=`read_opt_value_from_line "$*" "-normalbam"`
    local icgcid_normalbam=`read_opt_value_from_line "$*" "-extn"`
    local download_tries=`read_opt_value_from_line "$*" "-nt"`
    local step_outd=`read_opt_value_from_line "$*" "-step-outd"`

    # Download file
    logmsg "* Executing score-client..."
    ${ICGCSTOR_HOME_DIR}/bin/score-client --profile collab download --object-id ${icgcid_normalbam} --output-dir ${step_outd} 2>&1 || exit 1

    # Find bam file name
    local bam_file_name=`find_bam_filename ${step_outd}`
    
    if [ -z "${bam_file_name}" ]; then
        logmsg "Error: bam file not found after download process was completed"
        exit 1
    fi

    # Move file
    mv ${bam_file_name} ${normalbam} || exit 1
}

########
download_collab_tum_bam_explain_cmdline_opts()
{
    # -extt option
    description="External database id of tumor bam file to download"
    explain_cmdline_req_opt "-extt" "<string>" "$description"

    # -nt option
    description="Number of download tries per file (${DEFAULT_NUMBER_OF_DOWNLOAD_TRIES} by default)"
    explain_cmdline_opt "-nt" "<int>" "$description"    
}

########
download_collab_tum_bam_define_opts()
{
    # Initialize variables
    local cmdline=$1
    local stepspec=$2
    local optlist=""

    # Define the -step-outd option, the output directory for the step
    local step_outd=`get_step_outdir_given_stepspec "$stepspec"`
    define_opt "-step-outd" ${step_outd} optlist || exit 1
    
    # -extt option
    define_cmdline_opt "$cmdline" "-extt" optlist || exit 1

    # -nt option
    define_cmdline_nonmandatory_opt "$cmdline" "-nt" ${DEFAULT_NUMBER_OF_DOWNLOAD_TRIES} optlist || exit 1

    # -tumorbam option
    local abs_datadir=`get_absolute_shdirname ${DATADIR_BASENAME}`
    local tumorbam=${abs_datadir}/tumor.bam
    define_opt "-tumorbam" $tumorbam optlist || exit 1

    # Save option list
    save_opt_list optlist    
}

########
download_collab_tum_bam()
{
    # Initialize variables
    local tumorbam=`read_opt_value_from_line "$*" "-tumorbam"`
    local icgcid_tumorbam=`read_opt_value_from_line "$*" "-extt"`
    local download_tries=`read_opt_value_from_line "$*" "-nt"`
    local step_outd=`read_opt_value_from_line "$*" "-step-outd"`

    # Download file
    logmsg "* Executing score-client..."
    ${ICGCSTOR_HOME_DIR}/bin/score-client --profile collab download --object-id ${icgcid_tumorbam} --output-dir ${step_outd} 2>&1 || exit 1

    # Find bam file name
    local bam_file_name=`find_bam_filename ${step_outd}`
    
    if [ -z "${bam_file_name}" ]; then
        logmsg "Error: bam file not found after download process was completed"
        exit 1
    fi

    # Move file
    mv ${bam_file_name} ${tumorbam} || exit 1
}

########
download_gdc_norm_bam_explain_cmdline_opts()
{
    # -extn option
    description="External database id of normal bam file to download"
    explain_cmdline_opt "-extn" "<string>" "$description"

    # -gdprocs option
    description="Number of processes used by the gdc download client (${DEFAULT_NUMBER_OF_GDC_DOWNLOAD_PROCS} by default)"
    explain_cmdline_opt "-gdprocs" "<int>" "$description"

    # -gdctok option
    description="GDC API auth token file"
    explain_cmdline_req_opt "-gdctok" "<string>" "$description"

    # -nt option
    description="Number of download tries per file (${DEFAULT_NUMBER_OF_DOWNLOAD_TRIES} by default)"
    explain_cmdline_opt "-nt" "<int>" "$description"
}

########
download_gdc_norm_bam_define_opts()
{
    # Initialize variables
    local cmdline=$1
    local stepspec=$2
    local optlist=""

    # Define the -step-outd option, the output directory for the step
    local step_outd=`get_step_outdir_given_stepspec "$stepspec"`
    define_opt "-step-outd" ${step_outd} optlist || exit 1

    # -extn option
    define_cmdline_opt "$cmdline" "-extn" optlist || exit 1

    # -gdprocs option
    define_cmdline_nonmandatory_opt "$cmdline" "-gdprocs" ${DEFAULT_NUMBER_OF_GDC_DOWNLOAD_PROCS} optlist || exit 1

    # -gdctok option
    define_cmdline_opt "$cmdline" "-gdctok" optlist || exit 1

    # -nt option
    define_cmdline_nonmandatory_opt "$cmdline" "-nt" ${DEFAULT_NUMBER_OF_DOWNLOAD_TRIES} optlist || exit 1

    # -normalbam option
    local abs_datadir=`get_absolute_shdirname ${DATADIR_BASENAME}`
    local normalbam=${abs_datadir}/normal.bam
    define_opt "-normalbam" $normalbam optlist || exit 1

    # Save option list
    save_opt_list optlist    
}

########
get_gdc_bamfname()
{
    local gdcid=$1
    local outd=$2
    local outf

    outf=`find ${outd} -name "*.bam"` || return 1

    echo ${outf}
}

########
download_gdc_norm_bam()
{
    # Initialize variables
    local normalbam=`read_opt_value_from_line "$*" "-normalbam"`
    local gdcid_normalbam=`read_opt_value_from_line "$*" "-extn"`
    local gdprocs=`read_opt_value_from_line "$*" "-gdprocs"`
    local gdctok=`read_opt_value_from_line "$*" "-gdctok"`
    local download_tries=`read_opt_value_from_line "$*" "-nt"`
    local step_outd=`read_opt_value_from_line "$*" "-step-outd"`

    # Activate conda environment
    logmsg "* Activating conda environment..."
    conda activate gdc-client 2>&1 || exit 1

    # Download file (with multiple tries)
    gdc-client download -n ${gdprocs} -t ${gdctok} -d ${step_outd} --retry-amount ${download_tries} ${gdcid_normalbam} 2>/dev/null || exit 1
    
    # Deactivate conda environment
    logmsg "* Deactivating conda environment..."
    conda deactivate 2>&1

    # Move file
    local gdc_bamfname
    gdc_bamfname=`get_gdc_bamfname ${gdcid_normalbam} ${step_outd}` || exit 1
    mv ${gdc_bamfname} ${normalbam} || exit 1
}

########
download_gdc_norm_bam_conda_envs()
{
    define_conda_env gdc-client gdc-client.yml
}

########
download_gdc_tum_bam_explain_cmdline_opts()
{
    # -extt option
    description="External database id of tumor bam file to download"
    explain_cmdline_opt "-extt" "<string>" "$description"

    # -gdprocs option
    description="Number of processes used by the gdc download client (${DEFAULT_NUMBER_OF_GDC_DOWNLOAD_PROCS} by default)"
    explain_cmdline_opt "-gdprocs" "<int>" "$description"

    # -gdctok option
    description="GDC API auth token file"
    explain_cmdline_req_opt "-gdctok" "<string>" "$description"

    # -nt option
    description="Number of download tries per file (${DEFAULT_NUMBER_OF_DOWNLOAD_TRIES} by default)"
    explain_cmdline_opt "-nt" "<int>" "$description"
}

########
download_gdc_tum_bam_define_opts()
{
    # Initialize variables
    local cmdline=$1
    local stepspec=$2
    local optlist=""

    # Define the -step-outd option, the output directory for the step
    local step_outd=`get_step_outdir_given_stepspec "$stepspec"`
    define_opt "-step-outd" ${step_outd} optlist || exit 1

    # -extt option
    define_cmdline_opt "$cmdline" "-extt" optlist || exit 1

    # -gdprocs option
    define_cmdline_nonmandatory_opt "$cmdline" "-gdprocs" ${DEFAULT_NUMBER_OF_GDC_DOWNLOAD_PROCS} optlist || exit 1

    # -gdctok option
    define_cmdline_opt "$cmdline" "-gdctok" optlist || exit 1

    # -nt option
    define_cmdline_nonmandatory_opt "$cmdline" "-nt" ${DEFAULT_NUMBER_OF_DOWNLOAD_TRIES} optlist || exit 1

    # -tumorbam option
    local abs_datadir=`get_absolute_shdirname ${DATADIR_BASENAME}`
    local tumorbam=${abs_datadir}/tumor.bam
    define_opt "-tumorbam" $tumorbam optlist || exit 1

    # Save option list
    save_opt_list optlist    
}

########
download_gdc_tum_bam()
{
    # Initialize variables
    local tumorbam=`read_opt_value_from_line "$*" "-tumorbam"`
    local gdcid_tumorbam=`read_opt_value_from_line "$*" "-extt"`
    local gdprocs=`read_opt_value_from_line "$*" "-gdprocs"`
    local gdctok=`read_opt_value_from_line "$*" "-gdctok"`
    local download_tries=`read_opt_value_from_line "$*" "-nt"`
    local step_outd=`read_opt_value_from_line "$*" "-step-outd"`

    # Activate conda environment
    logmsg "* Activating conda environment..."
    conda activate gdc-client 2>&1 || exit 1

    # Download file (with multiple tries)
    gdc-client download -n ${gdprocs} -t ${gdctok} -d ${step_outd} --retry-amount ${download_tries} ${gdcid_tumorbam} 2>/dev/null || exit 1
    
    # Deactivate conda environment
    logmsg "* Deactivating conda environment..."
    conda deactivate 2>&1

    # Move file
    local gdc_bamfname
    gdc_bamfname=`get_gdc_bamfname ${gdcid_tumorbam} ${step_outd}` || exit 1
    mv ${gdc_bamfname} ${tumorbam} || exit 1
}

########
download_gdc_tum_bam_conda_envs()
{
    define_conda_env gdc-client gdc-client.yml
}
