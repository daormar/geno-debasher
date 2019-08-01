# *- python -*

# import modules
import io, sys, getopt, operator, pysam

##################################################
def take_pars():
    flags={}
    values={}
    flags["f_given"]=False
    flags["c_given"]=False
    values["altaft"]=-1
    
    try:
        opts, args = getopt.getopt(sys.argv[1:],"f:c:",["vcf=","altaft="])
    except getopt.GetoptError:
        print_help()
        sys.exit(2)
    if(len(opts)==0):
        print_help()
        sys.exit()
    else:
        for opt, arg in opts:
            if opt in ("-f", "--vcf"):
                values["vcf"] = arg
                flags["f_given"]=True
            elif opt in ("-c", "--altaft"):
                values["altaft"] = float(arg)
                flags["c_given"]=True
    return (flags,values)

##################################################
def check_pars(flags,values):
    if(flags["f_given"]==False):
        print >> sys.stderr, "Error! -f parameter not given"
        sys.exit(2)

##################################################
def print_help():
    print >> sys.stderr, "filter_vcf -f <string> -c <int>"
    print >> sys.stderr, ""
    print >> sys.stderr, "-f <string>  name of VCF. To avoid possible problems, it is better that the"
    print >> sys.stderr, "             file is bgzipped and tabix indexed"
    print >> sys.stderr, "-c <float>   Alternate allele frequencies threshold"

##################################################
def check_altaft(line,altaft):
    if altaft<0:
        return True
    else:
        if "CAF" in line.info:
            # CAF info available
            if len(line.info['CAF']) <= 1:
                # Discard malformed entries
                return True
            else:
                # Find at least one alternate allele frequency equal or
                # above threshold
                for i in range(1,len(line.info['CAF'])):
                    af_str=line.info['CAF'][i]
                    if af_str==".":
                        af_str="0.0"
                    if float(af_str) >= altaft:
                        return True
                return False
        else:
            # No CAF info available
            return True

##################################################
def filt_vcf(file,altaft):
    # Open files
    invcf= pysam.VariantFile(file)
    outvcf= pysam.VariantFile("-", mode= 'w', header= invcf.header)

    # Process input vcf
    for line in invcf:
        # Check alternate allele frequency threshold
       altaft_passed=check_altaft(line,altaft)
       if altaft_passed:
           outvcf.write(line)
        
    # Close files
    outvcf.close()
    invcf.close()

    
##################################################
def process_pars(flags,values):
    filt_vcf(values["vcf"],values["altaft"])
    
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
