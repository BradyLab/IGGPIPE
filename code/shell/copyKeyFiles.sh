# Copy key files in test directory $1 to $1/$2

if [[ -z $1 || -z $2 ]] ; then
    echo Missing argument 1 and/or 2
    echo Usage: source copyKeyFiles.sh "<srcDestDir>" "<destSubdir>"
    return
fi

mkdir -p $1/$2

function checkIfFilesExistAndCopy {
    MAINDIR=$1
    SUBDIR=$2
    PATTERN=$3
    DEST=$MAINDIR/$SUBDIR/
    FILES=$MAINDIR/$PATTERN
    for f in $FILES
    do
        if [ ! -f $f ]; then
            echo "File $f does not exist, not copied"
        else
            echo "Copying $f"
            cp $f $DEST
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
checkIfFilesExistAndCopy $1 $2 'MarkersAnnotated_*'
