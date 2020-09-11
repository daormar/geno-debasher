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
import sys, gzip, getopt

##################################################
def take_pars():
    flags={}
    values={}
    flags["m_given"]=False
    flags["g_given"]=False
    flags["v_given"]=False
    flags["l_given"]=False

    try:
        opts, args = getopt.getopt(sys.argv[1:],"m:g:v:l:",["maf=","gap=","vcf=","listc="])
    except getopt.GetoptError:
        print_help()
        sys.exit(2)
    if(len(opts)==0):
        print_help()
        sys.exit()
    else:
        for opt, arg in opts:
            if opt in ("-m", "--maf"):
                values["maf"] = float(arg)
                flags["m_given"]=True
            elif opt in ("-g", "--gap"):
                values["gap"] = int(arg)
                flags["g_given"]=True
            elif opt in ("-v", "--vcf"):
                values["vcf"] = arg
                flags["v_given"]=True
            elif opt in ("-l", "--listc"):
                values["listc"] = arg
                flags["l_given"]=True
    return (flags,values)

##################################################
def check_pars(flags,values):
    if(flags["m_given"]==False):
        print("Error! -m parameter not given", file=sys.stderr)
        sys.exit(2)
    if(flags["g_given"]==False):
        print("Error! -g parameter not given", file=sys.stderr)
        sys.exit(2)
    if(flags["v_given"]==False):
        print("Error! -v parameter not given", file=sys.stderr)
        sys.exit(2)
    if(flags["l_given"]==False):
        print("Error! -l parameter not given", file=sys.stderr)
        sys.exit(2)

##################################################
def print_help():
    print("create_snv_pos_ascat -m <float> -g <int> -v <string> -l <string>", file=sys.stderr)
    print("", file=sys.stderr)
    print("-m <float>     MAF value", file=sys.stderr)
    print("-g <int>       GAP value", file=sys.stderr)
    print("-v <string>    VCF file name", file=sys.stderr)
    print("-l <string>    List of contigs (one contig name per line)", file=sys.stderr)

##################################################
def get_contigs(listc):
    file = open(listc, 'r')
    contigs={}
    for line in file:
        line=line.strip("\n")
        fields=line.split()
        contigs[fields[0]]=1
    return contigs

##################################################
def process_vcf(maf,gap,vcf,contigs):
    with gzip.open(vcf, "r") as fi:
        prev_event_contig=""
        # Process vcf lines
        for lin in fi :
            cols = lin.split("\t")
            if len(cols) == 8 and cols[0] in contigs:
                # Initialize entry variables
                entry_contig=cols[0]
                entry_pos=int(cols[1])
                entry_id=cols[2]
                entry_maf=-1
                entry_is_snv=False

                # Extract info
                info = cols[7].split(";")
                for i in info:
                    # Check SNV status
                    if(i=="TSA=SNV"):
                        entry_is_snv=True
                    # Check MAF
                    aux = i.split("=")
                    if(aux[0] == "MAF"):
                        entry_maf=float(aux[1])

                # Print data when appropriate
                if(entry_is_snv and entry_maf >= maf and (prev_event_contig!=entry_contig or (prev_event_contig==entry_contig and entry_pos-gap > prev_event_pos))):
                    print("{}\t{}\t{}".format(entry_id, entry_contig, entry_pos))
                    prev_event_contig=entry_contig
                    prev_event_pos=entry_pos

##################################################
def process_pars(flags,values):
    contigs=get_contigs(values["listc"])
    process_vcf(values["maf"],values["gap"],values["vcf"],contigs)

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
