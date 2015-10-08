# Copy key files in test directory $1 to $1/$2

if [[ -z $1 || -z $2 ]] ; then
    echo Missing argument 1 and/or 2
    echo Usage: source copyKeyFiles.sh "<srcDestDir>" "<destSubdir>"
    return
fi

mkdir -p $1/$2

function checkIfFileExistsAndCopy {
    if [ ! -f $1 ];
       then
        echo "File $1 does not exist, not copied"
        return
    fi
    cp $1 $2
    }

checkIfFileExistsAndCopy $1/LCRs_*.tsv $1/$2/
checkIfFileExistsAndCopy $1/BadKmers_*.tsv $1/$2/
checkIfFileExistsAndCopy $1/IndelGroupsOverlapping_*.tsv $1/$2/
checkIfFileExistsAndCopy $1/IndelGroupsNonoverlapping_*.tsv $1/$2/
checkIfFileExistsAndCopy $1/NonvalidatedMarkers_*.tsv $1/$2/
checkIfFileExistsAndCopy $1/MarkerErrors_*.tsv $1/$2/
checkIfFileExistsAndCopy $1/MarkersOverlapping_*.tsv $1/$2/
checkIfFileExistsAndCopy $1/MarkersNonoverlapping_*.tsv $1/$2/
checkIfFileExistsAndCopy $1/MarkerCounts_*.pdf $1/$2/
checkIfFileExistsAndCopy $1/MarkerDensity_*.png $1/$2/
