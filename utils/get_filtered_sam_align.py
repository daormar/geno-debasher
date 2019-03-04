# *- python -*

# import modules
import io, sys, getopt, operator

##################################################
def take_pars():
    flags={}
    values={}
    flags["l_given"]=False

    try:
        opts, args = getopt.getopt(sys.argv[1:],"l:",["listc="])
    except getopt.GetoptError:
        print_help()
        sys.exit(2)
    if(len(opts)==0):
        print_help()
        sys.exit()
    else:
        for opt, arg in opts:
            if opt in ("-l", "--listc"):
                values["listc"] = arg
                flags["l_given"]=True
    return (flags,values)

##################################################
def check_pars(flags,values):
    if(flags["l_given"]==False):
        print >> sys.stderr, "Error! -l parameter not given"
        sys.exit(2)

##################################################
def print_help():
    print >> sys.stderr, "get_filtered_sam_header -l <string>"
    print >> sys.stderr, ""
    print >> sys.stderr, "-l <string>    List of contigs to keep (one contig name per line)"

##################################################
def getContigsToKeep(listc):
    file = open(listc, 'r')
    contigs_to_keep={}
    for line in file:
        line=line.strip("\n")
        fields=line.split()
        contigs_to_keep[fields[0]]=1
    return contigs_to_keep

##################################################
def extract_safield(fields):    
    for i in range(10,len(fields)):
        if fields[i].startswith("SA:Z:"):
            return i,fields[i]
    return -1,""

##################################################
def alig_contains_contig_to_keep(alig,contigs_to_keep):
    alig_fields=alig.split(",")
    if alig_fields>=1:
        if alig_fields[0] in contigs_to_keep:
            return True
        else:
            return False
    else:
        return False

##################################################
def filter_chim_aligs(chim_aligs,contigs_to_keep):
    filtered_chim_aligs=""
    chim_aligs_array=chim_aligs.split(";")
    for alig in chim_aligs_array:
        if alig_contains_contig_to_keep(alig,contigs_to_keep):
            filtered_chim_aligs=filtered_chim_aligs+alig+";"
        
    return filtered_chim_aligs

##################################################
def filter_safield(safield,contigs_to_keep):
    safield_parts=safield.split(":")
    if len(safield_parts)==3:
        chim_aligs=safield_parts[2]
        filtered_chim_aligs=filter_chim_aligs(chim_aligs,contigs_to_keep)
        if filtered_chim_aligs=="":
            return ""
        else:
            return safield_parts[0]+":"+safield_parts[1]+":"+filtered_chim_aligs
    else:
        return safield

##################################################
def replace_filtered_sa_field(fields,sa_idx,filtered_safield):
    filtered_entry=""
    for i in range(len(fields)):
        if i==sa_idx:
            field_to_add=filtered_safield
        else:
            field_to_add=fields[i]
        if field_to_add!="":
            if filtered_entry=="":
                filtered_entry=fields[i]
            else:
                filtered_entry=filtered_entry+"\t"+field_to_add
    return filtered_entry
    
##################################################
def obtain_filtered_entry(entry,contigs_to_keep):
    fields=entry.split()
    sa_idx,safield=extract_safield(fields)
    if sa_idx<0:
        return False,""
    else:
        filtered_safield=filter_safield(safield,contigs_to_keep)
        if safield==filtered_safield:
            return False,""
        else:
            filtered_entry=replace_filtered_sa_field(fields,sa_idx,filtered_safield)
            return True,filtered_entry

##################################################
def process_pars(flags,values):
    contigs_to_keep=getContigsToKeep(values["listc"])

    # read stdin line by line
    lineno=1
    for line in sys.stdin:
        line=line.strip("\n")
        filter_required,filtered_entry=obtain_filtered_entry(line,contigs_to_keep)
        if filter_required:
            print filtered_entry
        else:
            print line
        lineno+=1
                
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
