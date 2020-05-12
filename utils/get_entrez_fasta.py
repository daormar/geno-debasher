"""
Bio-PanPipe package
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
import io, sys, getopt, operator, requests, time
import xml.etree.ElementTree as ET

##################################################
def take_pars():
    flags={}
    values={}
    flags["a_given"]=False

    try:
        opts, args = getopt.getopt(sys.argv[1:],"a:",["accession="])
    except getopt.GetoptError:
        print_help()
        sys.exit(2)
    if(len(opts)==0):
        print_help()
        sys.exit()
    else:
        for opt, arg in opts:
            if opt in ("-a", "--accession"):
                values["accession"] = arg
                flags["a_given"]=True
    return (flags,values)

##################################################
def check_pars(flags,values):
    if(flags["a_given"]==False):
        print >> sys.stderr, "Error! -a parameter not given"
        sys.exit(2)

##################################################
def print_help():
    print >> sys.stderr, "get_entrez_fasta -a <string>"
    print >> sys.stderr, ""
    print >> sys.stderr, "-a <string>    Accession number of the FASTA data to download"

##################################################
def get_info(url,num_retries,time_between_retries=1):
    success=False
    for i in range(num_retries):
        print >> sys.stderr, "Getting data from",url,"( Attempt:",i+1,")"
        try:
            req=requests.get(url)
            req.raise_for_status()
            success=True
            break
        except requests.exceptions.RequestException as e:
            print >> sys.stderr, e
            time.sleep(time_between_retries)
    if success:
        return req
    else:
        print >> sys.stderr, "Maximum number of attempts exceeded, aborting"
        sys.exit(1)
        
##################################################
def extract_esearch_info(accession):
    # Get information
    url="https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"+"esearch.fcgi?db=nuccore&term="+accession+"&usehistory=y"
    num_retries=5
    time_between_retries=1
    req=get_info(url,num_retries,time_between_retries)

    # Process information
    root = ET.fromstring(req.content)

    # Check for errors
    for child in root:
        if child.tag=="ErrorList":
            print >> sys.stderr, "Error while extracting esearch information, aborting (see request result below)"
            print >> sys.stderr,req.content
            sys.exit(1)
            
    # Extract QueryKey and WebEnv fields
    key = None
    web = None
    for child in root:
        if child.tag=="QueryKey":
            key=child.text
        elif child.tag=="WebEnv":
            web=child.text

    # Check for possible errors
    if key is None or web is None:
        print >> sys.stderr, "Something went wrong while extracting esearch information, aborting"
        sys.exit(1)

    return key,web
    
##################################################
def post_efetch_info(key,web):
    url="https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"+"efetch.fcgi?db=nuccore&query_key="+key+"&WebEnv="+web+"&rettype=fasta&retmode=text";

    # Get information
    num_retries=5
    req=get_info(url,num_retries)

    # If request content is in xml format, then retrieval has gone wrong
    try:
        root = ET.fromstring(req.content)
        print >> sys.stderr, "Error while extracting efetch information, aborting (see request result below)"
        print >> sys.stderr,req.content
        sys.exit(1)
    except ET.ParseError as e:
        pass
        
    return req
    
##################################################
def process_pars(flags,values):
    key,web=extract_esearch_info(values["accession"])
        
    req=post_efetch_info(key,web)
    print req.content
    
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
