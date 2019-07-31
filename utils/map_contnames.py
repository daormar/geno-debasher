# *- python -*

# import modules
import io, sys, getopt, operator, csv

##################################################
def take_pars():
    flags={}
    values={}
    flags["m_given"]=False
    flags["f_given"]=False
    flags["c_given"]=False
    flags["invert"]=False

    try:
        opts, args = getopt.getopt(sys.argv[1:],"m:f:c:i",["mapf=","file=","col=","invert"])
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
            elif opt in ("-f", "--file"):
                values["file"] = arg
                flags["f_given"]=True
            elif opt in ("-c", "--col"):
                values["col"] = int(arg)
                flags["f_given"]=True
            elif opt in ("-i", "--invert"):
                flags["invert"]=True
    return (flags,values)

##################################################
def check_pars(flags,values):
    if(flags["m_given"]==False):
        print >> sys.stderr, "Error! -m parameter not given"
        sys.exit(2)

    if(flags["f_given"]==False):
        print >> sys.stderr, "Error! -f parameter not given"
        sys.exit(2)

    if(flags["c_given"]==False):
        print >> sys.stderr, "Error! -c parameter not given"
        sys.exit(2)

##################################################
def print_help():
    print >> sys.stderr, "map_contnames -m <string> -f <string> -c <int> [--invert]"
    print >> sys.stderr, ""
    print >> sys.stderr, "-m <string>    File with contig to contig mapping"
    print >> sys.stderr, "-f <string>    File where mapping needs to be done"
    print >> sys.stderr, "-c <int>       Column to be mapped (numbered from 0)"
    print >> sys.stderr, "--invert       Invert mapping"   

##################################################
def get_contig_map(mapf):
    file = open(mapf, 'r')
    contig_map={}
    inv_contig_map={}
    for line in file:
        line=line.strip("\n")
        fields=line.split()
        if(len(fields)==2):
            contig_map[fields[0]]=fields[1]
            inv_contig_map[fields[1]]=fields[0]
    return contig_map,inv_contig_map

##################################################
def map_contigs(filename,contig_map,inv_contig_map,col,invert):
    # Select contig map to be used
    if invert:
        selected_contig_map=inv_contig_map
    else:
        selected_contig_map=contig_map

    # read file line by line
    file = open(filename, 'r')
    for line in file:
        line=line.strip("\n")
        sniffer = csv.Sniffer()
        dialect = sniffer.sniff(line)
        fields=line.split(dialect.delimiter)
        if col<len(fields):
            contig=fields[col]
            if(contig in selected_contig_map):
                modified_line=fields[0]
                for i in range(1,len(fields)):
                    if i==col:
                        modified_line=modified_line+dialect.delimiter+selected_contig_map[contig]
                    else:
                        modified_line=modified_line+dialect.delimiter+fields[i]
                print modified_line
            else:
                print line

##################################################
def process_pars(flags,values):
    contig_map,inv_contig_map=get_contig_map(values["mapf"])
    map_contigs(values["file"],contig_map,inv_contig_map,values["col"],flags["invert"])
                                   
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
