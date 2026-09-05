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

######################
# CLEANING PROCESSES #
######################

########
delete_bam_files_explain_opts()
{
    # -process-outd option
    local description="output directory"
    explain_opt "-process-outd" "<file>" "$description"

    # -datadir option
    local description="data directory"
    explain_opt "-datadir" "<file>" "$description"
}

########
delete_bam_files_identify_cmdline_opts()
{
    :
}

########
delete_bam_files_define_opts()
{
    # Initialize variables
    local cmdline=$1
    local process_spec=$2
    local process_name=$3
    local process_outdir=$4
    local optlist=""

    # Define the -process-outd option, the output directory for the process
    define_opt "-process-outd" "${process_outdir}" optlist || return 1

    # -datadir option
    local abs_datadir=`get_absolute_shdirname "${GENODB_BAM_COMMON_DATADIR_BASENAME}"`
    define_opt "-datadir" "${abs_datadir}" optlist || return 1

    # Save option list
    save_opt_list optlist
}

########
delete_bam_files()
{
    # Initialize variables
    local process_outd=`read_opt_value_from_func_args "-process-outd" "$@"`
    local abs_datadir=`read_opt_value_from_func_args "-datadir" "$@"`

    # Delete bam files
    "${RM}" -f "${abs_datadir}"/*.bam || return 1
}

########
clear_datadir_explain_opts()
{
    # -process-outd option
    local description="output directory"
    explain_opt "-process-outd" "<file>" "$description"

    # -datadir option
    local description="data directory"
    explain_opt "-datadir" "<file>" "$description"
}

########
clear_datadir_identify_cmdline_opts()
{
    :
}

########
clear_datadir_define_opts()
{
    # Initialize variables
    local cmdline=$1
    local process_spec=$2
    local process_name=$3
    local process_outdir=$4
    local optlist=""

    # Define the -process-outd option, the output directory for the process
    define_opt "-process-outd" "${process_outdir}" optlist || return 1

    # -datadir option
    local abs_datadir=`get_absolute_shdirname ${GENODB_BAM_COMMON_DATADIR_BASENAME}`
    define_opt "-datadir" "${abs_datadir}" optlist || return 1

    # Save option list
    save_opt_list optlist
}

########
clear_datadir()
{
    # Initialize variables
    local process_outd=`read_opt_value_from_func_args "-process-outd" "$@"`
    local abs_datadir=`read_opt_value_from_func_args "-datadir" "$@"`

    # Delete bam files
    "${RM}" -rf "${abs_datadir}"/* || return 1

    # Print README.txt file
    echo "NOTE: This directory was cleared by means of the 'clear_datadir' process" > "${abs_datadir}"/README.txt || return 1
}
