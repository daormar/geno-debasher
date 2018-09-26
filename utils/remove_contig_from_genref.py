# *- python -*

# import modules
import io, sys, getopt, operator

##################################################
def take_pars():
    flags={}
    values={}
    flags["c_given"]=False
    flags["g_given"]=False

    try:
        opts, args = getopt.getopt(sys.argv[1:],"c:g:",["contig=","genref="])
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
    return (flags,values)

##################################################
def check_pars(flags,values):
    if(flags["c_given"]==False):
        print >> sys.stderr, "Error! -c parameter not given"
        sys.exit(2)

    if(flags["g_given"]==False):
        print >> sys.stderr, "Error! -g parameter not given"
        sys.exit(2)

##################################################
def print_help():
    print >> sys.stderr, "remove_contig_from_genref -c <string> -g <string>"
    print >> sys.stderr, ""
    print >> sys.stderr, "-c <string>    Contig name"
    print >> sys.stderr, "-g <string>    File with genome reference"

##################################################
def process_pars(flags,values):
    file = open(values["genref"], 'r')
    contigstr=">"+values["contig"]
    skip=False
    # read file line by line
    for line in file:
        line=line.strip("\n")
        fields=line.split()
        if(contigstr==fields[0]):
            skip=True
        elif(">" in fields[0]):
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
