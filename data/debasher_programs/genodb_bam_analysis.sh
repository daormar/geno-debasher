genodebasher_pkglibdir="/home/dortiz/bio/software/geno-debasher/lib/geno-debasher"
genodebasher_debasher_programs_dir="/home/dortiz/bio/software/geno-debasher/share/geno-debasher/debasher_programs"
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

# INCLUDE BASH LIBRARY
. "${genodebasher_pkglibdir}"/genodb_bam_common_lib || exit 1

#################
# CFG FUNCTIONS #
#################

########
genodb_bam_analysis_shared_dirs()
{
    define_shared_dir "data"
    define_shared_dir "split"
    define_shared_dir "summary/germline_snvs"
}

# INCLUDE BASH FILES IMPLEMENTING PROCESSES
. "${genodebasher_debasher_programs_dir}"/genodb_genref_processes.sh || exit 1
. "${genodebasher_debasher_programs_dir}"/genodb_bam_download_processes.sh || exit 1
. "${genodebasher_debasher_programs_dir}"/genodb_bam_manip_processes.sh || exit 1
. "${genodebasher_debasher_programs_dir}"/genodb_bam_filter_processes.sh || exit 1
. "${genodebasher_debasher_programs_dir}"/genodb_bam_analysis_processes.sh || exit 1
. "${genodebasher_debasher_programs_dir}"/genodb_bam_summariz_processes.sh || exit 1
. "${genodebasher_debasher_programs_dir}"/genodb_cleaning_processes.sh || exit 1

# EXTRA FILES
. "${genodebasher_debasher_programs_dir}"/genodb_bam_ascat_processes.sh || exit 1
. "${genodebasher_debasher_programs_dir}"/genodb_bam_facets_processes.sh || exit 1
