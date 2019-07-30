# *- bash -*

########
print_desc()
{
    echo "postproc_snpgcc post-processes snp gc correction file so as to make it compatible with run_ascat tool"
    echo "type \"postproc_snpgcc --help\" to get usage information"
}

########
usage()
{
    echo "postproc_snpgcc   -s <string> [--help]"
    echo ""
    echo "-s <string>       SNP gc correction file"
    echo "--help            Display this help and exit"
}

########
read_pars()
{
    s_given=0
    while [ $# -ne 0 ]; do
        case $1 in
            "--help") usage
                      exit 1
                      ;;
            "-s") shift
                  if [ $# -ne 0 ]; then
                      snpgccfile=$1
                      s_given=1
                  fi
                  ;;
        esac
        shift
    done   
}

########
check_pars()
{
    if [ ${s_given} -eq 0 ]; then   
        echo "Error! -s parameter not given!" >&2
        exit 1
    else
        if [ ! -f ${snpgccfile} ]; then
            echo "Error! file ${snpgccfile} does not exist" >&2
            exit 1
        fi
    fi
}

########
replace_first_column()
{
    local file=$1
    ${AWK} '{
              if(FNR==1) print $0
              else
              {
               firstfield=sprintf("snp%d",FNR-1) 
               $1=firstfield
               printf "%s",firstfield
               for(i=2;i<=NF;++i)
                printf "\t%s",$i
               printf"\n"
              }
            }' ${file}
}

########
process_pars()
{
    replace_first_column ${snpgccfile}
}

########

if [ $# -eq 0 ]; then
    print_desc
    exit 1
fi

read_pars $@ || exit 1

check_pars || exit 1

process_pars || exit 1
