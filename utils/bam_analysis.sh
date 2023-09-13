# Bio-PanPipe package
# Copyright (C) 2019,2020 Daniel Ortiz-Mart\'inez
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

# INCLUDE BASH LIBRARY
. "${biopanpipe_bindir}"/bam_common_lib || exit 1

#################
# CFG FUNCTIONS #
#################

########
bam_analysis_shared_dirs()
{
    define_shared_dir "${DATADIR_BASENAME}"
    define_shared_dir "${SUMMARYDIR_BASENAME}"
    define_shared_dir "${GERM_SNVS_SUM_DIR_BASENAME}"
}

########
bam_analysis_fifos()
{
    :
}

# INCLUDE BASH FILES IMPLEMENTING PROCESSES
. "${biopanpipe_bindir}"/genref_processes || exit 1
. "${biopanpipe_bindir}"/bam_download_processes || exit 1
. "${biopanpipe_bindir}"/bam_manip_processes || exit 1
. "${biopanpipe_bindir}"/bam_analysis_processes || exit 1
. "${biopanpipe_bindir}"/bam_summariz_processes || exit 1
. "${biopanpipe_bindir}"/cleaning_processes || exit 1
