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
    print >> sys.stderr, "filter_nondiscordant_pheno_entries [-f <string>] [-h <string>] [-t <string>]"
    print >> sys.stderr, ""
    print >> sys.stderr, "-f <string>        File containing output of query metadata tools (if not"
    print >> sys.stderr, "                   given, input is read from stdin)"
    print >> sys.stderr, "-h <string>        Label to identify healthy samples (\"non-tumor\" by default)"
    print >> sys.stderr, "-t <string>        Label to identify tumor samples (\"tumor\" by default)"

##################################################
def process_pars(flags,values):
    # Determine stream to be processed
    if values["file"]=="":
        stream = sys.stdin
    else:
        stream = open(values["file"], 'r')

    # Process output of tools to query metadata. For those entries with
    # more than two samples, generate all possible combinations of two
    # elements without repetitions
    for line in stream:
        line=line.strip("\n")
        if(line.find("="+values["tumor_label"])!=-1 and line.find("="+values["healthy_label"])!=-1):
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
