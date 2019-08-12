# *- python -*

# import modules
import io, sys, operator

##################################################
def reverse(string): 
    string = "".join(reversed(string)) 
    return string

##################################################
def get_num_ones(binary_num):
    num_ones=0
    for number in binary_num:
        if number=="1":
            num_ones=num_ones+1
    return num_ones

##################################################
def get_positions_of_ones(binary_num):
    positions=[]
    reversed_binary_num=reverse(binary_num)
    for i in range(len(binary_num)-2):
        if reversed_binary_num[i]=="1":
            positions.append(i)
    return positions

##################################################
def print_combination(positions,fields):
    sample1=fields[positions[0]].strip()
    sample2=fields[positions[1]].strip()
    print(sample1+" ; "+sample2)
    
##################################################
def process_entry_with_more_than_two_samples(fields):
    for i in range(2**len(fields)):
        bin_i=bin(i)
        num_ones=get_num_ones(bin_i)
        if(num_ones==2):
            positions=get_positions_of_ones(bin_i)
            print_combination(positions,fields)

##################################################
def main(argv):
    # Process output of tools to query metadata. For those entries with
    # more than two samples, generate all possible combinations of two
    # elements without repetitions
    for line in sys.stdin:
        line=line.strip("\n")
        fields=line.split(";")
        if len(fields)==1:
            print("ERROR")
        elif len(fields)==2:
            print(line)
        else:
            process_entry_with_more_than_two_samples(fields)

if __name__ == "__main__":
    main(sys.argv)
