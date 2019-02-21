# *- python -*

# import modules
import io, sys, getopt, operator

##################################################
def take_pars():
    flags={}
    values={}
    flags["c_given"]=False
    flags["g_given"]=False
    flags["l_given"]=False

    try:
        opts, args = getopt.getopt(sys.argv[1:],"c:g:l:",["contig=","genref=","listc="])
    except getopt.GetoptError:
        print_help()
        sys.exit(2)
    if(len(opts)==0):
        print_help()
        sys.exit()
    else:
        for opt, arg in opts:
            if opt in ("-c", "--contig"):
                values["contig"] = arg
                flags["c_given"]=True
            if opt in ("-g", "--genref"):
                values["genref"] = arg
                flags["g_given"]=True
            if opt in ("-l", "--listc"):
                values["listc"] = arg
                flags["l_given"]=True
    return (flags,values)

##################################################
def check_pars(flags,values):
    if(flags["g_given"]==False):
        print >> sys.stderr, "Error! -g parameter not given"
        sys.exit(2)
    if(flags["c_given"]==False and flags["l_given"]==False):
        print >> sys.stderr, "Error! -c or -l parameter not given"
        sys.exit(2)
    if(flags["c_given"]==True and flags["l_given"]==True):
        print >> sys.stderr, "Error! -c and -l parameters cannot be given simultaneously"
        sys.exit(2)

##################################################
def print_help():
    print >> sys.stderr, "filter_contig_from_genref -c <string> {-l <string> | -g <string>}"
    print >> sys.stderr, ""
    print >> sys.stderr, "-c <string>    Name of contig to remove"
    print >> sys.stderr, "-l <string>    List of contigs to keep (one contig name per line)"
    print >> sys.stderr, "-g <string>    File with genome reference"

##################################################
def getContigsToKeep(listc):
    file = open(listc, 'r')
    contigsToKeep={}
    for line in file:
        line=line.strip("\n")
        fields=line.split()
        contigsToKeep[fields[0]]=1
    return contigsToKeep
    
##################################################
def process_pars(flags,values):
    # Initialize variables
    if(flags["c_given"]):
        contigToFilter=values["contig"]
    else:
        contigToFilter=""
        
    if(flags["l_given"]):
        contigsToKeep=getContigsToKeep(values["listc"])
    else:
        contigsToKeep={}

    # Filter genome
    file = open(values["genref"], 'r')
    skip=False
    # read file line by line
    for line in file:
        line=line.strip("\n")
        fields=line.split()
        if(">" in fields[0]):
            contigname=fields[0][1:]
            if(contigname==contigToFilter or (not contigname in contigsToKeep)):
                skip=True
            else:
                skip=False
            
        if(not skip):
            print line
                                   
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
