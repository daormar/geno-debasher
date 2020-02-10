# *- python -*

# import modules
import io, sys

##################################################
def main(argv):
    for line in sys.stdin:
        line=line.strip("\n")
        if(line.find("=tumor")!=-1 and line.find("=non-tumor")!=-1):
            print(line)

if __name__ == "__main__":
    main(sys.argv)
