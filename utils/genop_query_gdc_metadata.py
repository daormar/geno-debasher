"""
Geno-PanPipe package
Copyright 2019,2020 Daniel Ortiz-Mart\'inez

This library is free software; you can redistribute it and/or
modify it under the terms of the GNU Lesser General Public License
as published by the Free Software Foundation; either version 3
of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License
along with this program; If not, see <http://www.gnu.org/licenses/>.
"""

# *- python -*

# import modules
import io
import sys
import os
import getopt
import operator
import csv

##################################################
def take_pars():
    flags = {}
    values = {}
    flags["m_given"] = False
    flags["f_given"] = False
    flags["verbose"] = False

    try:
        opts, args = getopt.getopt(sys.argv[1:], "m:f:v", [
                                   "metadatafile=", "format="])
    except getopt.GetoptError:
        print_help()
        sys.exit(2)
    if(len(opts) == 0):
        print_help()
        sys.exit()
    else:
        for opt, arg in opts:
            if opt in ("-m", "--metadatafile"):
                values["metadatafile"] = arg
                flags["m_given"] = True
            elif opt in ("-f", "--format"):
                values["format"] = int(arg)
                flags["f_given"] = True
            elif opt in ("-v", "--verbose"):
                flags["verbose"] = True
    return (flags, values)

##################################################
def check_pars(flags, values):
    if(flags["m_given"] == False):
        print("Error! -m parameter not given", file=sys.stderr)
        sys.exit(2)

    if(flags["f_given"] == False):
        print("Error! -f parameter not given", file=sys.stderr)
        sys.exit(2)

##################################################
def print_help():
    print("genop_query_gdc_metadata -m <string> -f <int>", file=sys.stderr)
    print("               [-v]", file=sys.stderr)
    print("", file=sys.stderr)
    print("-m <string>    File with metadata information", file=sys.stderr)
    print("-f <int>       Output format:", file=sys.stderr)
    print("                1: case_id,sample_id,uid,filename,tissue,phenotype,gender", file=sys.stderr)
    print("                2: Same as 1 but sorted by case_id", file=sys.stderr)
    print("                3: Same as 2 but entries for same case_id appear in same line", file=sys.stderr)
    print("                4: Same as 3 but only uid, phenotype and gender is listed", file=sys.stderr)
    print("-v             Verbose mode", file=sys.stderr)

##################################################
def get_pheno_info(field):
    if "Normal" in field:
        return "phenotype=non-tumor"
    elif "Tumor" in field:
        return "phenotype=tumor"
    else:
        return None

##################################################
def get_gender_info(field):
    if "male" in field:
        return "gender=male"
    elif "female" in field:
        return "gender=female"
    else:
        return None

##################################################
def get_tissue_info(field):
    return "tissue="+field.replace(" ", "_")

##################################################
def extract_metadata_info(filename):
    file = open(filename, 'r')
    metadata_info = csv.DictReader(file, dialect='excel-tab')
    return metadata_info

##################################################
def get_info_in_basic_format(metadata_info):
    formatted_info = []

    # Populate formatted_info structure
    for record in metadata_info:
        formatted_info.append((record["cases.0.case_id"], record["cases.0.samples.0.sample_id"], record["file_name"], record["file_id"], get_tissue_info(
            record["cases.0.diagnoses.0.tissue_or_organ_of_origin"]), get_pheno_info(record["cases.0.samples.0.tissue_type"]), get_gender_info(record["cases.0.demographic.gender"])))

    return formatted_info

##################################################
def filter_fields(entry, field_list):
    if(field_list):
        filtered_entry = []
        for field in field_list:
            filtered_entry.append(entry[field])
        return filtered_entry
    else:
        return entry

##################################################
def group_formatted_info_by_donor(formatted_info, field_list):
    # Create and populate map to make grouping easier
    group_map = {}
    for elem in formatted_info:
        case_id_idx = 0
        if(elem[case_id_idx] in group_map):
            group_map[elem[case_id_idx]].append(
                filter_fields(elem, field_list))
        else:
            group_map[elem[case_id_idx]] = []
            group_map[elem[case_id_idx]].append(
                filter_fields(elem, field_list))
    # Create grouped info
    formatted_info_grouped = []
    for key in group_map:
        tmplist = []
        for elem in group_map[key]:
            if(tmplist):
                tmplist.append((";"))
            tmplist.append(elem)
        flattmplist = [item for sublist in tmplist for item in sublist]
        formatted_info_grouped.append(flattmplist)

    return formatted_info_grouped

##################################################
def format_info(format, metadata_info):
    if(format == 1):
        return get_info_in_basic_format(metadata_info)
    elif(format == 2):
        formatted_info = get_info_in_basic_format(metadata_info)
        return sorted(formatted_info, key=operator.itemgetter(0))
    elif(format == 3):
        formatted_info = get_info_in_basic_format(metadata_info)
        return group_formatted_info_by_donor(formatted_info, [])
    elif(format == 4):
        formatted_info = get_info_in_basic_format(metadata_info)
        return group_formatted_info_by_donor(formatted_info, [3, 5, 6])

##################################################
def print_info(formatted_info):
    for elem in formatted_info:
        row = ""
        for i in range(len(elem)):
            if(i == 0):
                row = elem[i]
            else:
                row = row+" "+elem[i]
        print(row)

##################################################
def process_pars(flags, values):
    # Extract info from files
    metadata_info = extract_metadata_info(values["metadatafile"])

    # Format information
    formatted_info = format_info(values["format"], metadata_info)

    # Print information
    print_info(formatted_info)

##################################################
def main(argv):
    # take parameters
    (flags, values) = take_pars()

    # check parameters
    check_pars(flags, values)

    # process parameters
    process_pars(flags, values)


if __name__ == "__main__":
    main(sys.argv)
