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
import pysam

##################################################
def take_pars():
    flags = {}
    values = {}
    flags["f_given"] = False
    flags["c_given"] = False
    values["altaft"] = -1

    try:
        opts, args = getopt.getopt(sys.argv[1:], "f:c:", ["vcf=", "altaft="])
    except getopt.GetoptError:
        print_help()
        sys.exit(2)
    if(len(opts) == 0):
        print_help()
        sys.exit()
    else:
        for opt, arg in opts:
            if opt in ("-f", "--vcf"):
                values["vcf"] = arg
                flags["f_given"] = True
            elif opt in ("-c", "--altaft"):
                values["altaft"] = float(arg)
                flags["c_given"] = True
    return (flags, values)

##################################################
def check_pars(flags, values):
    if(flags["f_given"] == False):
        print("Error! -f parameter not given", file=sys.stderr)
        sys.exit(2)

##################################################
def print_help():
    print("filter_vcf -f <string> -c <int>", file=sys.stderr)
    print("", file=sys.stderr)
    print("-f <string>  name of VCF. To avoid possible problems, it is better that the", file=sys.stderr)
    print("             file is bgzipped and tabix indexed", file=sys.stderr)
    print("-c <float>   Alternate allele frequencies threshold", file=sys.stderr)

##################################################
def check_altaft(line, altaft):
    if altaft < 0:
        return True
    else:
        if "CAF" in line.info:
            # CAF info available
            if len(line.info['CAF']) <= 1:
                # Discard malformed entries
                return True
            else:
                # Find at least one alternate allele frequency equal or
                # above threshold
                for i in range(1, len(line.info['CAF'])):
                    af_str = line.info['CAF'][i]
                    if af_str == ".":
                        af_str = "0.0"
                    if float(af_str) >= altaft:
                        return True
                return False
        else:
            # No CAF info available
            return True

##################################################
def filt_vcf(file, altaft):
    # Open files
    invcf = pysam.VariantFile(file)
    outvcf = pysam.VariantFile("-", mode='w', header=invcf.header)

    # Process input vcf
    for line in invcf:
        # Check alternate allele frequency threshold
        altaft_passed = check_altaft(line, altaft)
        if altaft_passed:
            outvcf.write(line)

    # Close files
    outvcf.close()
    invcf.close()


##################################################
def process_pars(flags, values):
    filt_vcf(values["vcf"], values["altaft"])

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
