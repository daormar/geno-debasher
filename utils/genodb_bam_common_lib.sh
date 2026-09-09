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

GENODB_BAM_COMMON_SUMMARY_FILE_EXT="sum"

#############
# FUNCTIONS #
#############

########
genodb_bam_common::create_summary_file()
{
    # Initialize variables
    summarydir=$1
    label=$2
    vcf=$3

    # Create file
    echo "$vcf" > "${summarydir}/${label}.${GENODB_BAM_COMMON_SUMMARY_FILE_EXT}"
}

########
genodb_bam_common::slurm_to_java_mem_spec()
{
    local mem=$1
    echo "${mem}" | "${AWK}" '{if(substr($1,length($1),1) ~ /^[0-9]/) printf"%sM",$1; else printf"%s",$1}'
}
