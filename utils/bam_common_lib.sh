# *- bash -*

#############
# CONSTANTS #
#############

DATADIR_BASENAME="data"

#############
# FUNCTIONS #
#############

########
get_normal_bam_filename()
{
    local cmdline=$1
    local given=0
    local normalbam
    normalbam=`read_opt_value_from_line "$cmdline" "-n"` && given=1
    if [ $given -eq 1 ]; then
        # -n option was given
        file_exists $normalbam || { errmsg "file $normalbam does not exist" ; return 1; }
        echo $normalbam
    else
        # Check -extn option
        check_opt_given "$cmdline" "-extn" || { errmsg "-n or -extn option should be given" ; return 1; }
        local abs_datadir=`get_absolute_shdirname ${DATADIR_BASENAME}`
        normalbam=${abs_datadir}/normal.bam
        echo $normalbam
    fi
}

########
get_tumor_bam_filename()
{
    local cmdline=$1
    local given=0
    local tumorbam
    tumorbam=`read_opt_value_from_line "$cmdline" "-t"` && given=1
    if [ $given -eq 1 ]; then
        # -t option was given
        file_exists $tumorbam || { errmsg "file $tumorbam does not exist" ; return 1; }
        echo $tumorbam
    else
        # Check -extt option
        check_opt_given "$cmdline" "-extt" || { errmsg "-t or -extt option should be given" ; return 1; }
        local abs_datadir=`get_absolute_shdirname ${DATADIR_BASENAME}`
        tumorbam=${abs_datadir}/tumor.bam
        echo $tumorbam
    fi
}
