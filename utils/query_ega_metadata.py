"""
Geno-PanPipe package
Copyright 2019,2020 Daniel Ortiz-Mart\'inez
 
This library is free software; you can redistribute it and/or
modify it under the terms of the GNU Lesser General Public License
as published by the Free Software Foundation; either version 3
of the License, or (at your option) any later version.
 
This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU Lesser General Public License for more details.
 
You should have received a copy of the GNU Lesser General Public License
along with this program; If not, see <http://www.gnu.org/licenses/>.
"""
 
# *- python -*

# import modules
import io, sys, os, getopt, operator

# Global variables
ATTR_NOT_FOUND="ATTR_NOT_FOUND"
NO_ASPERA_FILENAME_AVAILABLE="NO_ASPERA_FILENAME_AVAILABLE"

##################################################
class sample_data:
    def __init__(self):
        self.filename=None
        self.fileaccession=None
        # NOTE: information indexed by sample_accession

##################################################
class analysis_data:
    def __init__(self):
        self.donor_id=None
        self.phenotype=None
        self.gender=None
        # NOTE: information indexed by ega_sample_id

##################################################
class study_data:
    def __init__(self):
        self.ega_sample_id=None
        # NOTE: information indexed by sample_accession
        
##################################################
def take_pars():
    flags={}
    values={}
    flags["s_given"]=False
    flags["a_given"]=False
    flags["t_given"]=False
    flags["f_given"]=False
    flags["p_given"]=False
    flags["verbose"]=False

    try:
        opts, args = getopt.getopt(sys.argv[1:],"s:a:t:f:p:v",["sampleinfofile=","analysisinfofile=","studyinfofile=","format=","asperacontent="])
    except getopt.GetoptError:
        print_help()
        sys.exit(2)
    if(len(opts)==0):
        print_help()
        sys.exit()
    else:
        for opt, arg in opts:
            if opt in ("-s", "--sampleinfofile"):
                values["sampleinfofile"] = arg
                flags["s_given"]=True
            elif opt in ("-a", "--analysisinfofile"):
                values["analysisinfofile"] = arg
                flags["a_given"]=True
            elif opt in ("-t", "--studyinfofile"):
                values["studyinfofile"] = arg
                flags["t_given"]=True
            elif opt in ("-f", "--format"):
                values["format"] = int(arg)
                flags["f_given"]=True
            elif opt in ("-p", "--asperacontent"):
                values["asperacontent"] = arg
                flags["p_given"]=True
            elif opt in ("-v", "--verbose"):
                flags["verbose"]=True
    return (flags,values)

##################################################
def check_pars(flags,values):
    if(flags["s_given"]==False):
        print("Error! -s parameter not given", file=sys.stderr)
        sys.exit(2)

    if(flags["a_given"]==False):
        print("Error! -a parameter not given", file=sys.stderr)
        sys.exit(2)

    if(flags["t_given"]==False):
        print("Error! -t parameter not given", file=sys.stderr)
        sys.exit(2)

    if(flags["f_given"]==False):
        print("Error! -f parameter not given", file=sys.stderr)
        sys.exit(2)

##################################################
def print_help():
    print("query_ega_metadata -s <string> -a <string> -t <string> -f <int>", file=sys.stderr)
    print("                   [-p <string>] [-v]", file=sys.stderr)
    print("", file=sys.stderr)
    print("-s <string>    File with sample information", file=sys.stderr)
    print("-a <string>    File with analysis information", file=sys.stderr)
    print("-t <string>    File with study information", file=sys.stderr)
    print("-f <int>       Output format:", file=sys.stderr)
    print("                1: SAMPLE_ACCESSION EGA_SAMPLE_ID FILE_ACCESSION FILENAME ASPERA_BOX_FILENAME DONOR_ID PHENOTYPE GENDER", file=sys.stderr)
    print("                2: Same as 1 but sorted by donor_id", file=sys.stderr)
    print("                3: Same as 2 but entries for same donor_id appear in same line", file=sys.stderr)
    print("                4: Same as 3 but only EGA_SAMPLE_ID, PHENOTYPE and GENDER are listed", file=sys.stderr)
    print("                5: Same as 3 but only ASPERA_BOX_FILENAME, PHENOTYPE and GENDER are listed (-p option is required)", file=sys.stderr)
    print("-p <string>    File listing Aspera box content", file=sys.stderr)
    print("-v             Verbose mode", file=sys.stderr)

##################################################
def standardize_filename(filename):
    # Remove extension
    if os.path.splitext(filename)[1]==".cip":
        result=os.path.splitext(filename)[0]
    else:
        result=filename
        
    # Retain base file name
    result=os.path.basename(result)
    
    return result

##################################################
def extract_sample_info(filename):
    sample_info_map={}
    file = open(filename, 'r')
    # read file line by line
    for line in file:
        line=line.strip("\n")
        fields=line.split()
        sample_accession=fields[1]
        sd=sample_data()
        sd.filename=standardize_filename(fields[2])
        sd.fileaccession=fields[3]
        sample_info_map[sample_accession]=sd
    return sample_info_map

##################################################
def extract_attribute_info(line):
    fields=line.split()
    # Remove blanks from attributes
    attribute_str=""
    for i in range(1,len(fields)):
        if i>1:
            attribute_str=attribute_str+"_"
        attribute_str=attribute_str+fields[i]

    # Initialize output values
    donor_id=ATTR_NOT_FOUND
    phenotype=ATTR_NOT_FOUND
    gender=ATTR_NOT_FOUND
    
    # Iterate over attributes
    attr_fields=attribute_str.split(";")
    for attr_field in attr_fields:
        # Extract key-value information
        keyvalue=attr_field.split("=")
        if len(keyvalue)==2:
            key=keyvalue[0]
            value=keyvalue[1]
            if(key=="donor_id"):
                donor_id=attr_field
            elif(key=="phenotype"):
                phenotype=attr_field
            elif(key=="gender"):
                gender=attr_field

    return (donor_id,phenotype,gender)
        
##################################################
def extract_analysis_info(filename):
    analysis_info_map={}
    file = open(filename, 'r')
    # read file line by line
    for line in file:
        line=line.strip("\n")
        fields=line.split()
        ega_sample_id=fields[0]
        ad=analysis_data()
        (ad.donor_id,ad.phenotype,ad.gender)=extract_attribute_info(line)
        analysis_info_map[ega_sample_id]=ad
    return analysis_info_map

##################################################
def extract_study_info(filename):
    study_info_map={}
    file = open(filename, 'r')
    # read file line by line
    for line in file:
        line=line.strip("\n")
        fields=line.split()
        sample_accession=fields[len(fields)-3]
        sd=study_data()
        sd.ega_sample_id=fields[len(fields)-1]
        study_info_map[sample_accession]=sd
    return study_info_map

##################################################
def load_aspera_box_info(filename):
    aspera_box_info=[]
    file = open(filename, 'r')
    # read file line by line
    for line in file:
        line=line.strip("\n")
        fields=line.split()
        aspera_box_info.append(fields[0])
    # Output is a list of aspera box filenames. Each one contains as
    # substring one of the filenames stored in sample data
    return aspera_box_info

##################################################
def match_aspera_filename(egafilename,aspera_box_info):
    matching = [s for s in aspera_box_info if egafilename in s]
    if(matching):
        return matching[0]
    else:
        return NO_ASPERA_FILENAME_AVAILABLE

##################################################
def get_info_in_basic_format(sample_info_map,analysis_info_map,study_info_map,aspera_box_info):
    formatted_info=[]

    # Populate formatted_info structure
    for sample_accession in study_info_map:
        ega_sample_id=study_info_map[sample_accession].ega_sample_id
        filename=sample_info_map[sample_accession].filename
        fileaccession=sample_info_map[sample_accession].fileaccession
        donor_id=analysis_info_map[ega_sample_id].donor_id
        phenotype=analysis_info_map[ega_sample_id].phenotype
        gender=analysis_info_map[ega_sample_id].gender
        aspera_filename=match_aspera_filename(filename,aspera_box_info)
        formatted_info.append((sample_accession,ega_sample_id,fileaccession,filename,aspera_filename,donor_id,phenotype,gender))

    return formatted_info

##################################################
def filter_fields(entry,field_list):
    if(field_list):
        filtered_entry=[]
        for field in field_list:
            filtered_entry.append(entry[field])
        return filtered_entry
    else:
        return entry    

##################################################
def group_formatted_info_by_donor(formatted_info,field_list):
    # Create and populate map to make grouping easier
    group_map={}
    for elem in formatted_info:
        donor_id_idx=5
        if(elem[donor_id_idx] in group_map):
            group_map[elem[donor_id_idx]].append(filter_fields(elem,field_list))
        else:
            group_map[elem[donor_id_idx]]=[]
            group_map[elem[donor_id_idx]].append(filter_fields(elem,field_list))
    # Create grouped info
    formatted_info_grouped=[]
    for key in group_map:
        tmplist=[]
        for elem in group_map[key]:
            if(tmplist):
                tmplist.append((";"))
            tmplist.append(elem)
        flattmplist=[item for sublist in tmplist for item in sublist]
        formatted_info_grouped.append(flattmplist)

    return formatted_info_grouped

##################################################
def format_info(format,sample_info_map,analysis_info_map,study_info_map,aspera_box_info):
    if(format==1):
        return get_info_in_basic_format(sample_info_map,analysis_info_map,study_info_map,aspera_box_info)
    elif(format==2):
        formatted_info=get_info_in_basic_format(sample_info_map,analysis_info_map,study_info_map,aspera_box_info)
        return sorted(formatted_info, key=operator.itemgetter(4))
    elif(format==3):
        formatted_info=get_info_in_basic_format(sample_info_map,analysis_info_map,study_info_map,aspera_box_info)
        return group_formatted_info_by_donor(formatted_info,[])
    elif(format==4):
        formatted_info=get_info_in_basic_format(sample_info_map,analysis_info_map,study_info_map,aspera_box_info)
        return group_formatted_info_by_donor(formatted_info,[2,6,7])
    elif(format==5):
        formatted_info=get_info_in_basic_format(sample_info_map,analysis_info_map,study_info_map,aspera_box_info)
        return group_formatted_info_by_donor(formatted_info,[4,6,7])

##################################################
def print_info(formatted_info):
    for elem in formatted_info:
        row=""
        for i in range(len(elem)):
            if(i==0):
                row=elem[i]
            else:
                row=row+" "+elem[i]
        print(row)

##################################################
def process_pars(flags,values):
    # Extract info from files
    sample_info_map=extract_sample_info(values["sampleinfofile"])
    analysis_info_map=extract_analysis_info(values["analysisinfofile"])
    study_info_map=extract_study_info(values["studyinfofile"])
    if(flags["p_given"]):
        aspera_box_info=load_aspera_box_info(values["asperacontent"])
    else:
        aspera_box_info=[]
        
    # Format information
    formatted_info=format_info(values["format"],sample_info_map,analysis_info_map,study_info_map,aspera_box_info)
    
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
