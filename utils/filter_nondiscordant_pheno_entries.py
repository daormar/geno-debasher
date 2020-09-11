"""
Bio-PanPipe package
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
import io, sys, operator, getopt

##################################################
def take_pars():
    flags={}
    values={}
    flags["f_given"]=False
    values["file"] = ""
    flags["h_given"]=False
    values["healthy_label"] = "non-tumor"
    flags["t_given"]=False
    values["tumor_label"] = "tumor"

    try:
        opts, args = getopt.getopt(sys.argv[1:],"f:h:t:",["file=","healthy-label=","tumor-label="])
    except getopt.GetoptError:
        print_help()
        sys.exit(2)
    if(len(opts)>3):
        print_help()
        sys.exit()
    else:
        for opt, arg in opts:
            if opt in ("-f", "--file"):
                values["file"] = arg
                flags["f_given"]=True
            elif opt in ("-h", "--healthy-label"):
                values["healthy_label"] = arg
                flags["h_given"]=True
            elif opt in ("-t", "--tumor-label"):
                values["tumor_label"] = arg
                flags["t_given"]=True
    return (flags,values)

##################################################
def check_pars(flags,values):
    None
    
##################################################
def print_help():
    print("filter_nondiscordant_pheno_entries [-f <string>] [-h <string>] [-t <string>]", file=sys.stderr)
    print("", file=sys.stderr)
    print("-f <string>        File containing output of query metadata tools (if not", file=sys.stderr)
    print("                   given, input is read from stdin)", file=sys.stderr)
    print("-h <string>        Label to identify healthy samples (\"non-tumor\" by default)", file=sys.stderr)
    print("-t <string>        Label to identify tumor samples (\"tumor\" by default)", file=sys.stderr)

##################################################
def process_pars(flags,values):
    # Determine stream to be processed
    if values["file"]=="":
        stream = sys.stdin
    else:
        stream = open(values["file"], 'r')

    # Filter entries
    for line in stream:
        line=line.strip("\n")
        if line.find("="+values["tumor_label"])!=-1 and line.find("="+values["healthy_label"])!=-1:
            print(line)

##################################################
def main(argv):
    # take parameters
    (flags,values)=take_pars()

    # check parameters
    check_pars(flags,values)

    # process parameters
    process_pars(flags,values)

if __name__ == "__main__":
    main(sys.argv)
