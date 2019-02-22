# *- python -*

# import modules
import io, sys, getopt, operator

##################################################
def take_pars():
    flags={}
    values={}
    flags["m_given"]=False
    flags["f_given"]=False

    try:
        opts, args = getopt.getopt(sys.argv[1:],"m:f:",["mapf=","vcf="])
    except getopt.GetoptError:
        print_help()
        sys.exit(2)
    if(len(opts)==0):
        print_help()
        sys.exit()
    else:
        for opt, arg in opts:
            if opt in ("-m", "--mapf"):
                values["mapf"] = arg
                flags["m_given"]=True
            if opt in ("-f", "--vcf"):
                values["vcf"] = arg
                flags["f_given"]=True
    return (flags,values)

##################################################
def check_pars(flags,values):
    if(flags["m_given"]==False):
        print >> sys.stderr, "Error! -m parameter not given"
        sys.exit(2)

    if(flags["f_given"]==False):
        print >> sys.stderr, "Error! -f parameter not given"
        sys.exit(2)

##################################################
def print_help():
    print >> sys.stderr, "map_vcf_contnames -m <string> -f <string>"
    print >> sys.stderr, ""
    print >> sys.stderr, "-m <string>    File with contig to contig mapping"
    print >> sys.stderr, "-f <string>    VCF file)"

##################################################
def getContigMap(mapf):
    file = open(mapf, 'r')
    contigMap={}
    for line in file:
        line=line.strip("\n")
        fields=line.split()
        if(len(fields)==2):
            contigMap[fields[0]]=fields[1]
    return contigMap

##################################################
def remove_prefix(prefix,text):
    if text.startswith(prefix):
        return text[len(prefix):]
    return text

##################################################
def process_pars(flags,values):
    contigMap=getContigMap(values["mapf"])

    # Filter genome
    file = open(values["vcf"], 'r')
    # read file line by line
    for line in file:
        line=line.strip("\n")
        if line.startswith("#"):
            print line
        else:
            fields=line.split()
            contig=fields[0]
            if(contig in contigMap):
                line_wo_contig=remove_prefix(contig,line)
                print contigMap[contig]+line_wo_contig
            else:
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
