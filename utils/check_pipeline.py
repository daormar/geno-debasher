# *- python -*

# import modules
import io, sys, getopt, operator

##################################################
def take_pars():
    flags={}
    values={}
    flags["a_given"]=False
    values["verbose"]=False

    try:
        opts, args = getopt.getopt(sys.argv[1:],"a:v",["afile="])
    except getopt.GetoptError:
        print_help()
        sys.exit(2)
    if(len(opts)==0):
        print_help()
        sys.exit()
    else:
        for opt, arg in opts:
            if opt in ("-a", "--afile"):
                values["afile"] = arg
                flags["a_given"]=True
            elif opt in ("-v", "--verbose"):
                flags["verbose"]=True
    return (flags,values)

##################################################
def check_pars(flags,values):
    if(flags["a_given"]==False):
        print >> sys.stderr, "Error! -a parameter not given"
        sys.exit(2)

##################################################
def print_help():
    print >> sys.stderr, "check_pipeline -a <string>"
    print >> sys.stderr, "               [-v]"
    print >> sys.stderr, ""
    print >> sys.stderr, "-a <string>    Analysis file"
    print >> sys.stderr, "-v             Verbose mode"

##################################################
def process_pars(flags,values):
    # TO-BE-DONE
    print "TO-BE-DONE"

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
