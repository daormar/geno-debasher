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
import operator
import getopt

##################################################
def take_pars():
    flags = {}
    values = {}
    flags["f_given"] = False
    values["file"] = ""

    try:
        opts, args = getopt.getopt(sys.argv[1:], "f:", ["file="])
    except getopt.GetoptError:
        print_help()
        sys.exit(2)
    if(len(opts) > 1):
        print_help()
        sys.exit()
    else:
        for opt, arg in opts:
            if opt in ("-f", "--file"):
                values["file"] = arg
                flags["f_given"] = True
    return (flags, values)

##################################################
def check_pars(flags, values):
    None

##################################################
def print_help():
    print("genop_get_two_elem_combs [-f <string>]", file=sys.stderr)
    print("", file=sys.stderr)
    print("-f <string>        File containing output of query metadata tools (if not", file=sys.stderr)
    print("                   given, input is read from stdin)", file=sys.stderr)

##################################################
def reverse(string):
    string = "".join(reversed(string))
    return string

##################################################
def get_num_ones(binary_num):
    num_ones = 0
    for number in binary_num:
        if number == "1":
            num_ones = num_ones+1
    return num_ones

##################################################
def get_positions_of_ones(binary_num):
    positions = []
    reversed_binary_num = reverse(binary_num)
    for i in range(len(binary_num)-2):
        if reversed_binary_num[i] == "1":
            positions.append(i)
    return positions

##################################################
def print_combination(positions, fields):
    sample1 = fields[positions[0]].strip()
    sample2 = fields[positions[1]].strip()
    print((sample1+" ; "+sample2))

##################################################
def process_entry_with_more_than_two_samples(fields):
    for i in range(2**len(fields)):
        bin_i = bin(i)
        num_ones = get_num_ones(bin_i)
        if(num_ones == 2):
            positions = get_positions_of_ones(bin_i)
            print_combination(positions, fields)

##################################################
def process_pars(flags, values):
    # Determine stream to be processed
    if values["file"] == "":
        stream = sys.stdin
    else:
        stream = open(values["file"], 'r')

    # Process output of tools to query metadata. For those entries with
    # more than two samples, generate all possible combinations of two
    # elements without repetitions
    for line in stream:
        line = line.strip("\n")
        fields = line.split(";")
        if len(fields) == 1:
            print("ERROR")
        elif len(fields) == 2:
            print(line)
        else:
            process_entry_with_more_than_two_samples(fields)

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
