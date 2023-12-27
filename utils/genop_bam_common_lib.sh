# Geno-DeBasher package
# Copyright (C) 2019-2024 Daniel Ortiz-Mart\'inez
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

#############
# CONSTANTS #
#############

DATADIR_BASENAME="data"
SPLITDIR_BASENAME="split"
SUMMARYDIR_BASENAME="summary"
GERM_SNVS_SUM_DIR_BASENAME="summary/germline_snvs"
SUMMARY_FILE_EXT="sum"

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
        file_exists "$normalbam" || { errmsg "file $normalbam does not exist" ; return 1; }
        echo "$normalbam"
    else
        # Check -extn option
        check_opt_given "$cmdline" "-extn" || { errmsg "-n or -extn option should be given" ; return 1; }
        local abs_datadir=`get_absolute_shdirname "${DATADIR_BASENAME}"`
        normalbam="${abs_datadir}"/normal.bam
        echo "$normalbam"
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
        file_exists "$tumorbam" || { errmsg "file $tumorbam does not exist" ; return 1; }
        echo "$tumorbam"
    else
        # Check -extt option
        check_opt_given "$cmdline" "-extt" || { errmsg "-t or -extt option should be given" ; return 1; }
        local abs_datadir=`get_absolute_shdirname "${DATADIR_BASENAME}"`
        tumorbam="${abs_datadir}"/tumor.bam
        echo "$tumorbam"
    fi
}

########
create_summary_file()
{
    # Initialize variables
    summarydir=$1
    label=$2
    vcf=$3

    # Create file
    echo "$vcf" > "${summarydir}/${label}.${SUMMARY_FILE_EXT}"
}

########
slurm_to_java_mem_spec()
{
    local mem=$1
    echo "${mem}" | "${AWK}" '{if(substr($1,length($1),1) ~ /^[0-9]/) printf"%sM",$1; else printf"%s",$1}'
}
