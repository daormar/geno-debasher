# *- python -*

# import modules
import io, sys, getopt, operator

##################################################
class donor_data:
    def __init__(self):
        self.donor_sex=None

##################################################
class awsmanif_data:
    def __init__(self):
        self.object_id=None
        self.filename=None
        self.donor_id=None

##################################################
def take_pars():
    flags={}
    values={}
    flags["d_given"]=False
    flags["a_given"]=False
    flags["f_given"]=False
    values["verbose"]=True

    try:
        opts, args = getopt.getopt(sys.argv[1:],"d:a:f:v",["donorinfofile=","awsmanif=","format="])
    except getopt.GetoptError:
        print_help()
        sys.exit(2)
    if(len(opts)==0):
        print_help()
        sys.exit()
    else:
        for opt, arg in opts:
            if opt in ("-d", "--donorinfofile"):
                values["donorinfo"] = arg
                flags["d_given"]=True
            if opt in ("-a", "--awsmanif"):
                values["awsmanif"] = arg
                flags["a_given"]=True
            if opt in ("-f", "--format"):
                values["format"] = int(arg)
                flags["f_given"]=True
            elif opt in ("-v", "--verbose"):
                verbose=1
    return (flags,values)

##################################################
def check_pars(flags,values):
    if(flags["d_given"]==False):
        print >> sys.stderr, "Error! -d parameter not given"
        sys.exit(2)

    if(flags["a_given"]==False):
        print >> sys.stderr, "Error! -a parameter not given"
        sys.exit(2)

    if(flags["f_given"]==False):
        print >> sys.stderr, "Error! -f parameter not given"
        sys.exit(2)

##################################################
def print_help():
    print >> sys.stderr, "query_icgc_metadata -d <string> -a <string> -f <int> [-v]"
    print >> sys.stderr, ""
    print >> sys.stderr, "-d <string>    File with donor information"
    print >> sys.stderr, "-a <string>    File with aws manifest"
    print >> sys.stderr, "-f <int>       Output format:"
    print >> sys.stderr, "                1: "
    print >> sys.stderr, "-v             Verbose mode"

##################################################
def extract_donor_info(filename):
    #TBD
    # sample_info_map={}
    # file = open(filename, 'r')
    # # read file line by line
    # for line in file:
    #     line=line.strip("\n")
    #     fields=line.split()
    #     sample_accession=fields[1]
    #     sd=sample_data()
    #     sd.sample_alias=fields[0]
    #     sd.filename=fields[2]
    #     sd.fileaccession=fields[3]
    #     sample_info_map[sample_accession]=sd
    # return sample_info_map

##################################################
# def extract_attribute_info(line):
#     fields=line.split()
#     # Remove blanks from attributes
#     attribute_str=""
#     for i in range(1,len(fields)):
#         attribute_str=attribute_str+"_"+fields[i]
#     attr_fields=attribute_str.split(";")
#     donor_id=attr_fields[1]
#     phenotype=attr_fields[2]
#     gender=attr_fields[4]

#     return (donor_id,phenotype,gender)
        
##################################################
def extract_awsmanif(filename):
    #TBD
    # analysis_info_map={}
    # file = open(filename, 'r')
    # # read file line by line
    # for line in file:
    #     line=line.strip("\n")
    #     fields=line.split()
    #     ega_sample_id=fields[0]
    #     ad=analysis_data()
    #     (ad.donor_id,ad.phenotype,ad.gender)=extract_attribute_info(line)
    #     analysis_info_map[ega_sample_id]=ad
    # return analysis_info_map

##################################################
# def get_info_in_basic_format(sample_info_map,analysis_info_map,study_info_map):
#     formatted_info=[]

#     # Populate formatted_info structure
#     for sample_accession in study_info_map:
#         ega_sample_id=study_info_map[sample_accession].ega_sample_id
#         filename=study_info_map[sample_accession].filename
#         fileaccession=sample_info_map[sample_accession].fileaccession
#         donor_id=analysis_info_map[ega_sample_id].donor_id
#         phenotype=analysis_info_map[ega_sample_id].phenotype
#         gender=analysis_info_map[ega_sample_id].gender
#         formatted_info.append((sample_accession,ega_sample_id,fileaccession,filename,donor_id,phenotype,gender))

#     return formatted_info

##################################################
# def group_formatted_info_by_donor(formatted_info):
#     # Create and populate map to make grouping easier
#     group_map={}
#     for elem in formatted_info:
#         if(elem[4] in group_map):
#             group_map[elem[4]].append(elem)
#         else:
#             group_map[elem[4]]=[]
#             group_map[elem[4]].append(elem)
#     # Created grouped info
#     formatted_info_grouped=[]
#     for key in group_map:
#         tmplist=[]
#         for elem in group_map[key]:
#             if(tmplist):
#                 tmplist.append((";"))
#             tmplist.append(elem)
#         flattmplist=[item for sublist in tmplist for item in sublist]
#         formatted_info_grouped.append(flattmplist)

#     return formatted_info_grouped

##################################################
def format_info(format,donor_info_map,awsmanif_map):
    # if(format==1):
    #     return get_info_in_basic_format(sample_info_map,analysis_info_map,study_info_map)
    # elif(format==2):
    #     formatted_info=get_info_in_basic_format(sample_info_map,analysis_info_map,study_info_map)
    #     return sorted(formatted_info, key=operator.itemgetter(4))
    # elif(format==3):
    #     formatted_info=get_info_in_basic_format(sample_info_map,analysis_info_map,study_info_map)
    #     return group_formatted_info_by_donor(formatted_info)
    
##################################################
def print_info(formatted_info):
    # for elem in formatted_info:
    #     row=""
    #     for i in range(len(elem)):
    #         if(i==0):
    #             row=elem[i]
    #         else:
    #             row=row+" "+elem[i]
    #     print row

##################################################
def process_pars(flags,values):
    # Extract info from files
    donor_info_map=extract_donor_info(values["donorinfofile"])
    awsmanif_map=extract_awsmanif(values["awsmanif"])

    # Format information
    formatted_info=format_info(values["format"],donor_info_map,awsmanif_map)
    
    # Print information
    print_info(formatted_info)

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
