# Count the number of distinct LCRs present in the LCR file given by $1, with
# the number of genomes given by $2.

if [[ -z $1 || -z $2 ]] ; then
    echo "Missing argument 1 and/or 2"
    echo "Usage: source countLCRsInLCRfile.sh <LCRfile> <N_GENOMES>"
    echo "Example: source code/shell/countLCRsInLCRfile.sh outTestHP11/LCRs_K11k2L100D10_2000.tsv 2"
    return
fi

if [[ ! -f $1 ]] ; then
    echo "File $1 does not exist"
    return
fi

# This is not as simple as a line count, because each line of an LCR file
# contains one k-mer, which is a member of one LCR, whose ID is included
# as the last value on each line.  Thus, the number of unique LCR IDs must
# be counted.  The number of columns in the LCRs file is 2+5*N_GENOMES,
# since the file contains five columns of data for each genome that is
# processed, plus two additional columns (k-mer is first column, LCR ID
# is last column).

echo 'LAST_COLUMN=`echo "2+5*$2" | bc`'
LAST_COLUMN=`echo "2+5*$2" | bc`
# Use tail to discard first line, which is the column name header line.
echo "tail +2 $1 | cut -f $LAST_COLUMN | uniq | wc -l"
tail +2 $1 | cut -f $LAST_COLUMN | uniq | wc -l
