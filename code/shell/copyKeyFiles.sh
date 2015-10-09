# Copy key files in output directory $1 to $2

if [[ -z $1 || -z $2 ]] ; then
    echo "Missing argument 1 and/or 2"
    echo "Usage: source copyKeyFiles.sh <srcDir> <destDir>"
    echo "Example: source code/shell/copyKeyFiles.sh outTestHP11 goodTest"
    return
fi

mkdir -p $2

function checkIfFilesExistAndCopy {
    SRCDIR=$1
    DSTDIR=$2
    PATTERN=$3
    FILES=$SRCDIR/$PATTERN
    for f in $FILES
    do
        if [ ! -f $f ]; then
            echo "File $f does not exist, not copied"
        else
            echo "Copying $f"
            cp $f $DSTDIR/
        fi
    done
    }

checkIfFilesExistAndCopy $1 $2 'LCRs_*.tsv'
checkIfFilesExistAndCopy $1 $2 'BadKmers_*.tsv'
checkIfFilesExistAndCopy $1 $2 'IndelGroupsOverlapping_*.tsv'
checkIfFilesExistAndCopy $1 $2 'IndelGroupsNonoverlapping_*.tsv'
checkIfFilesExistAndCopy $1 $2 'NonvalidatedMarkers_*.tsv'
checkIfFilesExistAndCopy $1 $2 'MarkerErrors_*.tsv'
checkIfFilesExistAndCopy $1 $2 'MarkersOverlapping_*.tsv'
checkIfFilesExistAndCopy $1 $2 'MarkersNonoverlapping_*.tsv'
checkIfFilesExistAndCopy $1 $2 'MarkersNonoverlapping_*.pdf'
checkIfFilesExistAndCopy $1 $2 'MarkerCounts_*.pdf'
checkIfFilesExistAndCopy $1 $2 'MarkerDensity_*.png'
checkIfFilesExistAndCopy $1 $2 'LCRs_*.dotplot.png'
checkIfFilesExistAndCopy $1 $2 'MarkersAnnotated_*'
