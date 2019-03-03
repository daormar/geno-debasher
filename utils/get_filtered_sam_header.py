# *- python -*

# import modules
import io, sys, getopt, operator

##################################################
def take_pars():
    flags={}
    values={}
    flags["h_given"]=False
    flags["l_given"]=False

    try:
        opts, args = getopt.getopt(sys.argv[1:],"h:l:",["samhdrf=","listc="])
    except getopt.GetoptError:
        print_help()
        sys.exit(2)
    if(len(opts)==0):
        print_help()
        sys.exit()
    else:
        for opt, arg in opts:
            if opt in ("-h", "--samhdrf"):
                values["samhdrf"] = arg
                flags["h_given"]=True
            if opt in ("-l", "--listc"):
                values["listc"] = arg
                flags["l_given"]=True
    return (flags,values)

##################################################
def check_pars(flags,values):
    if(flags["h_given"]==False):
        print >> sys.stderr, "Error! -h parameter not given"
        sys.exit(2)

    if(flags["l_given"]==False):
        print >> sys.stderr, "Error! -l parameter not given"
        sys.exit(2)

##################################################
def print_help():
    print >> sys.stderr, "get_filtered_sam_header -h <string> -l <string>"
    print >> sys.stderr, ""
    print >> sys.stderr, "-h <string>    File with SAM header"
    print >> sys.stderr, "-l <string>    List of contigs to keep (one contig name per line)"

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
    contigsToKeep=getContigsToKeep(values["listc"])

    # Filter header
    file = open(values["samhdrf"], 'r')
    # read file line by line
    for line in file:
        line=line.strip("\n")
        fields=line.split()
        if(fields[0] != "@SQ"):
            print line
        else:
            seqname=fields[1]
            seqname_fields=seqname.split(":")
            if(seqname_fields[1] in contigsToKeep):
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
