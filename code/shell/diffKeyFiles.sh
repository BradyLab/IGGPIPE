# Run diff on key files in directory $1 to $1/$2

if [[ -z $1 || -z $2 ]] ; then
    echo "Missing argument 1 and/or 2"
    echo "Usage: source diffKeyFiles.sh <srcDestDir> <diffSubdir>"
    echo "Example: source code/shell/diffKeyFiles.sh outTestHP11 goodTest"
    return
fi

mkdir -p $1/$2

function checkIfFilesExistAndDiff {
    MAINDIR=$1
    SUBDIR=$2
    PATTERN=$3
    DEST=$MAINDIR/$SUBDIR/
    FILES=$MAINDIR/$PATTERN
    for f in $FILES
    do
        if [ ! -f $f ]; then
            echo "File $f does not exist, not diffed"
        else
            echo "Diffing $f"
            diff $f $DEST
        fi
    done
    }

checkIfFilesExistAndDiff $1 $2 'LCRs_*.tsv'
checkIfFilesExistAndDiff $1 $2 'BadKmers_*.tsv'
checkIfFilesExistAndDiff $1 $2 'IndelGroupsOverlapping_*.tsv'
checkIfFilesExistAndDiff $1 $2 'IndelGroupsNonoverlapping_*.tsv'
checkIfFilesExistAndDiff $1 $2 'NonvalidatedMarkers_*.tsv'
checkIfFilesExistAndDiff $1 $2 'MarkerErrors_*.tsv'
checkIfFilesExistAndDiff $1 $2 'MarkersOverlapping_*.tsv'
checkIfFilesExistAndDiff $1 $2 'MarkersNonoverlapping_*.tsv'
checkIfFilesExistAndDiff $1 $2 'MarkersNonoverlapping_*.pdf'
checkIfFilesExistAndDiff $1 $2 'MarkerCounts_*.pdf'
checkIfFilesExistAndDiff $1 $2 'MarkerDensity_*.png'
checkIfFilesExistAndDiff $1 $2 'MarkersAnnotated_*'
