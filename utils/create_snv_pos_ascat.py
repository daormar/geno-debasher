import gzip
import sys

# Read options
arguments= len(sys.argv) -1
if(arguments != 3):
  print "create_snv_pos_ascat <maf_value> <gap_value> <file.vcf>"
  sys.exit(1)
  
maf=float(sys.argv[1])
gap=int(sys.argv[2])
vcf=sys.argv[3]

# Initialize variables
chroms = ["1", "2", "3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","X","Y"]
prev_event_chrom=""

# Process parameters
with gzip.open(vcf, "r") as fi :
  # Process vcf line
  for lin in fi :
    cols = lin.split("\t")
    if len(cols) == 8 and cols[0] in chroms :
      # Initialize entry variables
      entry_chrom=cols[0]
      entry_pos=int(cols[1])
      entry_id=cols[2]
      entry_maf=-1
      entry_is_snv=False

      # Extract info
      info = cols[7].split(";")
      for i in info:
        # Check SNV status
        if(i=="TSA=SNV"):
          entry_is_snv=True
        # Check MAF
        aux = i.split("=")
        if(aux[0] == "MAF"):
          entry_maf=float(aux[1])

      # Print data when appropriate
      if(entry_is_snv and entry_maf >= maf and (prev_event_chrom!=entry_chrom or (prev_event_chrom==entry_chrom and entry_pos-gap > prev_event_pos))):
        print "{}\t{}\t{}".format(entry_id, entry_chrom, entry_pos)
        prev_event_chrom=entry_chrom
        prev_event_pos=entry_pos
