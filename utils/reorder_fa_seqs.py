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
import io, sys, getopt, operator

##################################################
def take_pars():
    flags={}
    values={}
    flags["f_given"]=False
    flags["l_given"]=False

    try:
        opts, args = getopt.getopt(sys.argv[1:],"f:l:",["fasta=","listseq="])
    except getopt.GetoptError:
        print_help()
        sys.exit(2)
    if(len(opts)==0):
        print_help()
        sys.exit()
    else:
        for opt, arg in opts:
            if opt in ("-f", "--fasta"):
                values["fasta"] = arg
                flags["f_given"]=True
            elif opt in ("-l", "--listseq"):
                values["listseq"] = arg
                flags["l_given"]=True
    return (flags,values)

##################################################
def check_pars(flags,values):
    if(flags["f_given"]==False):
        print("Error! -f parameter not given", file=sys.stderr)
        sys.exit(2)

    if(flags["l_given"]==False):
        print("Error! -l parameter not given", file=sys.stderr)
        sys.exit(2)

##################################################
def print_help():
    print("reorder_fa_seqs -f <string> -l <string>", file=sys.stderr)
    print("", file=sys.stderr)
    print("-f <string>    File with genome reference", file=sys.stderr)
    print("-l <string>    File with ordered list of sequence names", file=sys.stderr)

##################################################
def read_seqs(filename):
    file = open(filename, 'r')
    seqs=[]
    for line in file:
        line=line.strip("\n")
        fields=line.split()
        seqname=fields[0]
        seqs.append(seqname)
    return seqs
    
##################################################
def extract_seq_name(line):
    fields=line.split()
    return fields[0][1:]
    
##################################################
def load_fasta(filename):
    fasta= open(filename)
    seq_dict= {}
    while True:
        line= fasta.readline()
        if line == '':
            break
        if line.strip().startswith('>'):
            seq_name= line.strip()[1:]
            seq_name= extract_seq_name(line.strip("\n"))
            seq_dict[seq_name]= []
        else:
            seq_dict[seq_name].append(line.strip())
    fasta.close()
    return seq_dict

##################################################
def print_reordered_fasta(seqs, fasta):
    for seq_name in seqs:
        seq_name = seq_name.strip()
        print('>' + seq_name)
        print('\n'.join(fasta[seq_name]))

##################################################
def process_pars(flags,values):
    # Initialize variables
    seqs=read_seqs(values["listseq"])

    # Load fasta
    fasta = load_fasta(values["fasta"])

    # Print reordered fasta
    print_reordered_fasta(seqs, fasta)
    
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
