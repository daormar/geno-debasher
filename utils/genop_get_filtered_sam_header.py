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
import getopt
import operator

##################################################
def take_pars():
    flags = {}
    values = {}
    flags["h_given"] = False
    flags["l_given"] = False

    try:
        opts, args = getopt.getopt(sys.argv[1:], "h:l:", [
                                   "samhdrf=", "listc="])
    except getopt.GetoptError:
        print_help()
        sys.exit(2)
    if(len(opts) == 0):
        print_help()
        sys.exit()
    else:
        for opt, arg in opts:
            if opt in ("-h", "--samhdrf"):
                values["samhdrf"] = arg
                flags["h_given"] = True
            elif opt in ("-l", "--listc"):
                values["listc"] = arg
                flags["l_given"] = True
    return (flags, values)

##################################################
def check_pars(flags, values):
    if(flags["h_given"] == False):
        print("Error! -h parameter not given", file=sys.stderr)
        sys.exit(2)

    if(flags["l_given"] == False):
        print("Error! -l parameter not given", file=sys.stderr)
        sys.exit(2)

##################################################
def print_help():
    print("genop_get_filtered_sam_header -h <string> -l <string>", file=sys.stderr)
    print("", file=sys.stderr)
    print("-h <string>    File with SAM header", file=sys.stderr)
    print("-l <string>    List of contigs to keep (one contig name per line)", file=sys.stderr)

##################################################
def get_contigs_to_keep(listc):
    file = open(listc, 'r')
    contigs_to_keep = {}
    for line in file:
        line = line.strip("\n")
        fields = line.split()
        contigs_to_keep[fields[0]] = 1
    return contigs_to_keep

##################################################
def process_pars(flags, values):
    contigs_to_keep = get_contigs_to_keep(values["listc"])

    # Filter header
    file = open(values["samhdrf"], 'r')
    # read file line by line
    for line in file:
        line = line.strip("\n")
        fields = line.split()
        if(fields[0] != "@SQ"):
            print(line)
        else:
            seqname = fields[1]
            seqname_fields = seqname.split(":")
            if(seqname_fields[1] in contigs_to_keep):
                print(line)

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
