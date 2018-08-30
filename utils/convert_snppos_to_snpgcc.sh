# *- bash -*

if [ $# -ne 3 ]; then
    echo "convert_snppos_to_snpgcc <ref.fa> <snpposfile.tsv> <outfile>" >&2
else
    # Initialize variables
    ref=$1
    snpposfile=$2
    outfile=$3
    TMPDIR=/tmp

    # Create directories
    mkdir ${TMPDIR}/splitPos ${TMPDIR}/splitGc ${TMPDIR}/splitGcLogs

    # Split file
    $SPLIT --number=l/10 -d ${snpposfile} ${TMPDIR}/splitPos/snpPos.
    
    # Process fragments
    conda activate ascatngs
    export PERL5LIB=/opt/anaconda3/envs/ascatngs/lib/perl5:/home/dortiz/bio/software/ascatngs/build/lib/perl5
    for file in `ls ${TMPDIR}/splitPos/`; do
        ascatSnpPanelGcCorrections.pl ${ref} ${TMPDIR}/splitPos/${file} > ${TMPDIR}/splitGc/${file} 2> ${TMPDIR}/splitGcLogs/${file}.log &
    done

    # Wait until all processes finish
    wait
    
    # Merge solutions
    $HEAD -n 1 ${TMPDIR}/splitGc/snpPos.00 > ${outfile}
    cat ${TMPDIR}/splitGc/snpPos.* | $GREP -vP 'Chr\tPosition' >> ${outfile}
    
    # Remove temporary directories
    rm -rf ${TMPDIR}/splitPos ${TMPDIR}/splitGc ${TMPDIR}/splitGcLogs
fi
