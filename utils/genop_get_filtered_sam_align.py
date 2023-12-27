"""
Geno-DeBasher package
Copyright 2019-2024 Daniel Ortiz-Mart\'inez

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
import getopt
import operator

##################################################
def take_pars():
    flags = {}
    values = {}
    flags["l_given"] = False

    try:
        opts, args = getopt.getopt(sys.argv[1:], "l:", ["listc="])
    except getopt.GetoptError:
        print_help()
        sys.exit(2)
    if(len(opts) == 0):
        print_help()
        sys.exit()
    else:
        for opt, arg in opts:
            if opt in ("-l", "--listc"):
                values["listc"] = arg
                flags["l_given"] = True
    return (flags, values)

##################################################
def check_pars(flags, values):
    if(flags["l_given"] == False):
        print("Error! -l parameter not given", file=sys.stderr)
        sys.exit(2)

##################################################
def print_help():
    print("genop_get_filtered_sam_align -l <string>", file=sys.stderr)
    print("", file=sys.stderr)
    print("-l <string>    List of contigs to keep (one contig name per line)", file=sys.stderr)

##################################################
def getContigsToKeep(listc):
    file = open(listc, 'r')
    contigs_to_keep = {}
    for line in file:
        line = line.strip("\n")
        fields = line.split()
        contigs_to_keep[fields[0]] = 1
    return contigs_to_keep

##################################################
def extract_safield(fields):
    for i in range(10, len(fields)):
        if fields[i].startswith("SA:Z:"):
            return i, fields[i]
    return -1, ""

##################################################
def extract_xafield(fields):
    for i in range(10, len(fields)):
        if fields[i].startswith("XA:Z:"):
            return i, fields[i]
    return -1, ""

##################################################
def alig_contains_contig_to_keep(alig, contigs_to_keep):
    alig_fields = alig.split(",")
    if alig_fields >= 1:
        if alig_fields[0] in contigs_to_keep:
            return True
        else:
            return False
    else:
        return False

##################################################
def filter_chim_aligs(chim_aligs, contigs_to_keep):
    filtered_chim_aligs = ""
    chim_aligs_array = chim_aligs.split(";")
    for alig in chim_aligs_array:
        if alig_contains_contig_to_keep(alig, contigs_to_keep):
            filtered_chim_aligs = filtered_chim_aligs+alig+";"

    return filtered_chim_aligs

##################################################
def filter_sa_xa_field(field, contigs_to_keep):
    field_parts = field.split(":")
    if len(field_parts) == 3:
        chim_aligs = field_parts[2]
        filtered_chim_aligs = filter_chim_aligs(chim_aligs, contigs_to_keep)
        if filtered_chim_aligs == "":
            return ""
        else:
            return field_parts[0]+":"+field_parts[1]+":"+filtered_chim_aligs
    else:
        return field

##################################################
def replace_field(fields, field_idx, new_field):
    filtered_entry = ""
    for i in range(len(fields)):
        if i == field_idx:
            field_to_add = new_field
        else:
            field_to_add = fields[i]
        if field_to_add != "":
            if filtered_entry == "":
                filtered_entry = fields[i]
            else:
                filtered_entry = filtered_entry+"\t"+field_to_add
    return filtered_entry

##################################################
def obtain_filtered_entry_sa(entry, contigs_to_keep):
    fields = entry.split()
    sa_idx, safield = extract_safield(fields)
    if sa_idx < 0:
        return False, entry
    else:
        filtered_safield = filter_sa_xa_field(safield, contigs_to_keep)
        if safield == filtered_safield:
            return False, entry
        else:
            filtered_entry = replace_field(fields, sa_idx, filtered_safield)
            return True, filtered_entry

##################################################
def obtain_filtered_entry_xa(entry, contigs_to_keep):
    fields = entry.split()
    xa_idx, xafield = extract_xafield(fields)
    if xa_idx < 0:
        return False, entry
    else:
        filtered_xafield = filter_sa_xa_field(xafield, contigs_to_keep)
        if xafield == filtered_xafield:
            return False, entry
        else:
            filtered_entry = replace_field(fields, xa_idx, filtered_xafield)
            return True, filtered_entry

##################################################
def obtain_filtered_entry(entry, contigs_to_keep):
    filter_sa_required, filtered_entry = obtain_filtered_entry_sa(
        entry, contigs_to_keep)
    filter_xa_required, filtered_entry = obtain_filtered_entry_xa(
        filtered_entry, contigs_to_keep)
    return filter_sa_required or filter_xa_required, filtered_entry

##################################################
def process_pars(flags, values):
    contigs_to_keep = getContigsToKeep(values["listc"])

    # read stdin line by line
    lineno = 1
    for line in sys.stdin:
        line = line.strip("\n")
        filter_required, filtered_entry = obtain_filtered_entry(
            line, contigs_to_keep)
        if filter_required:
            print(filtered_entry)
        else:
            print(line)
        lineno += 1

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
