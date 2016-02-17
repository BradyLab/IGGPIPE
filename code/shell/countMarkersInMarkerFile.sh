# Count the number of markers present a marker file given by $1, with the value
# of NDAMIN given by $2.

if [[ -z $1 || -z $2 ]] ; then
    echo "Missing argument 1 and/or 2"
    echo "Usage: source countMarkersInMarkerFile.sh <markerFile> <NDAMIN>"
    echo "Example: source code/shell/countMarkersInMarkerFile.sh outTestHP11/MarkersOverlapping_K11k2L100D10_2000A100_2000d10_100N2F0X20V3000W8M3G1.tsv 2"
    return
fi

if [[ ! -f $1 ]] ; then
    echo "File $1 does not exist"
    return
fi

# The second column of the marker file gives the value of NDA.  If $2 is not 2,
# loop from 2 to $2 and count the number of lines whose second column equals the
# loop value.  Then do a global count too.

if [[ $2 -ne 2 ]] ; then
    for i in `seq 2 $2` ; do
        echo "Number of markers in $1 with NDA = $i:"
        tail +2 $1 | cut -f 2 | grep -E "^$i$" | wc -l
    done
fi

echo "Total number of markers in $1:"
tail +2 $1 | cut -f 2 | wc -l
