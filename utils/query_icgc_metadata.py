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
import io, sys, getopt, operator, json

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
class table_data:
    def __init__(self):
        self.phenotype=None

##################################################
def take_pars():
    flags={}
    values={}
    flags["d_given"]=False
    flags["a_given"]=False
    flags["t_given"]=False
    flags["f_given"]=False
    flags["verbose"]=False

    try:
        opts, args = getopt.getopt(sys.argv[1:],"d:a:t:f:v",["donorinfofile=","awsmanif=","table=","format="])
    except getopt.GetoptError:
        print_help()
        sys.exit(2)
    if(len(opts)==0):
        print_help()
        sys.exit()
    else:
        for opt, arg in opts:
            if opt in ("-d", "--donorinfo"):
                values["donorinfo"] = arg
                flags["d_given"]=True
            elif opt in ("-a", "--awsmanif"):
                values["awsmanif"] = arg
                flags["a_given"]=True
            elif opt in ("-t", "--table"):
                values["table"] = arg
                flags["t_given"]=True
            elif opt in ("-f", "--format"):
                values["format"] = int(arg)
                flags["f_given"]=True
            elif opt in ("-v", "--verbose"):
                flags["verbose"]=True
    return (flags,values)

##################################################
def check_pars(flags,values):
    if(flags["d_given"]==False):
        print("Error! -d parameter not given", file=sys.stderr)
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
    print("query_icgc_metadata -d <string> -a <string> -t <string> -f <int> [-v]", file=sys.stderr)
    print("", file=sys.stderr)
    print("-d <string>    File with donor information", file=sys.stderr)
    print("-a <string>    File with aws manifest", file=sys.stderr)
    print("-t <string>    Table file in json format", file=sys.stderr)
    print("-f <int>       Output format:", file=sys.stderr)
    print("                1: FILE_ID OBJECT_ID FILENAME DONOR_ID PHENOTYPE GENDER", file=sys.stderr)
    print("                2: Same as 1 but sorted by donor_id", file=sys.stderr)
    print("                3: Same as 2 but entries for same donor_id appear in same line", file=sys.stderr)
    print("                4: Same as 3 but mini-bam files are excluded", file=sys.stderr)
    print("                5: Same as 3 but only OBJECT_ID, PHENOTYPE and GENDER are listed", file=sys.stderr)
    print("-v             Verbose mode", file=sys.stderr)

##################################################
def extract_donor_info(filename):
    donor_info_map={}
    file = open(filename, 'r')
    lineno=1
    # read file line by line
    for line in file:
        if(lineno>1):
            line=line.strip("\n")
            fields=line.split()
            donor_id=fields[0]
            dd=donor_data()
            dd.donor_sex=fields[4]
            donor_info_map[donor_id]=dd
        lineno+=1
    return donor_info_map
        
##################################################
def extract_awsmanif(filename):
    awsmanif_map={}
    file = open(filename, 'r')
    lineno=1
    # read file line by line
    for line in file:
        if(lineno>1):
            line=line.strip("\n")
            fields=line.split()
            file_id=fields[1]
            ad=awsmanif_data()
            ad.object_id=fields[2]
            ad.filename=fields[4]
            ad.donor_id=fields[8]
            awsmanif_map[file_id]=ad
        lineno+=1
    return awsmanif_map

##################################################
def extract_table(filename):
    table_map={}
    file=open(filename,"read")
    table=json.load(file)
    for row in table:
        file_id=row["id"]
        td=table_data()
        td.donor_id=row["donors"][0]["donorId"]
        td.phenotype=row["donors"][0]["specimenType"][0]
        table_map[file_id]=td
        
    return table_map

##################################################
def get_phenotype_code(phenotype):
    if("Normal" in phenotype or "normal" in phenotype):
        return "normal"
    elif("Tumour" in phenotype or "tumour" in phenotype or "Tumor" in phenotype or "tumor" in phenotype):
        return "tumor"
    
##################################################
def get_info_in_basic_format(donor_info_map,awsmanif_map,table_map,exclude_mini_bam_files):
    formatted_info=[]

    # Populate formatted_info structure
    for file_id in awsmanif_map:
        object_id=awsmanif_map[file_id].object_id
        filename=awsmanif_map[file_id].filename
        donor_id=awsmanif_map[file_id].donor_id
        donor_sex=donor_info_map[donor_id].donor_sex
        phenotype=table_map[file_id].phenotype
        phenotype_code=get_phenotype_code(phenotype)
        if(not exclude_mini_bam_files or ("mini" not in filename)):
            formatted_info.append((file_id,object_id,filename,donor_id,phenotype_code,donor_sex))

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
        if(elem[3] in group_map):
            group_map[elem[3]].append(filter_fields(elem,field_list))
        else:
            group_map[elem[3]]=[]
            group_map[elem[3]].append(filter_fields(elem,field_list))
    # Created grouped info
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
def format_info(format,donor_info_map,awsmanif_map,table_map):
    if(format==1):
        exclude_mini_bam_files=False
        return get_info_in_basic_format(donor_info_map,awsmanif_map,table_map,exclude_mini_bam_files)
    elif(format==2):
        exclude_mini_bam_files=False
        formatted_info=get_info_in_basic_format(donor_info_map,awsmanif_map,table_map,exclude_mini_bam_files)
        return sorted(formatted_info, key=operator.itemgetter(3))
    elif(format==3):
        exclude_mini_bam_files=False
        formatted_info=get_info_in_basic_format(donor_info_map,awsmanif_map,table_map,exclude_mini_bam_files)
        return group_formatted_info_by_donor(formatted_info,[])
    elif(format==4):
        exclude_mini_bam_files=True
        formatted_info=get_info_in_basic_format(donor_info_map,awsmanif_map,table_map,exclude_mini_bam_files)
        return group_formatted_info_by_donor(formatted_info,[])
    elif(format==5):
        exclude_mini_bam_files=True
        formatted_info=get_info_in_basic_format(donor_info_map,awsmanif_map,table_map,exclude_mini_bam_files)
        return group_formatted_info_by_donor(formatted_info,[1,4,5])

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
    donor_info_map=extract_donor_info(values["donorinfo"])
    awsmanif_map=extract_awsmanif(values["awsmanif"])
    table_map=extract_table(values["table"])

    # Format information
    formatted_info=format_info(values["format"],donor_info_map,awsmanif_map,table_map)
    
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
