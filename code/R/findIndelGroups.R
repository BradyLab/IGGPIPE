#######################################################################################
# See usage below for description.
# Author: Ted Toal
# Date: 2013-2015
# Brady Lab, UC Davis
#######################################################################################

# Enclose everything in braces so stop statements will work correctly.
{

# Pathname separator.
PATHSEP = ifelse(grepl("/", Sys.getenv("HOME")), "/", "\\")

# cat() that immediately flushes to console.
catnow = function(...)
    {
    cat(...)
    flush.console()
    return(invisible(0))
    }

# Print an R object's value.
objPrint = function(x, title="")
    {
    sink("temp.txt")
    print(x)
    sink()
    x.S = readLines("temp.txt")
    x.S = sub(" $", "", x.S)
    multiline = (length(x.S) > 1)
    if (!multiline)
        x.S = sub("[1] ", "", x.S, fixed=TRUE)
    x.S = paste(x.S, collapse="\n")
    if (title != "")
        {
        if (!multiline)
            catnow(title, ": ", sep="")
        else
            catnow(title, ":\n", sep="")
        }
    catnow(x.S, "\n", sep="")
    }

# Define function to display info when investigating.
inv = function(a, title="") { if (investigate) objPrint(a, title) }

# Get arguments.
testing = 0
#testing = 1 # For testing only.
#testing = 2 # For testing only.
{
if (testing == 0)
    args = commandArgs(TRUE)
else if (testing == 1)
    {
    args = c("~/Documents/UCDavis/BradyLab/Genomes/kmers/IGGPIPE",
        "outTestHP11/LCRs_K11k2L100D10_2000.tsv", 100, 2000, 10, 100, 2, 0,
        "MAX",
        "outTestHP11/IndelGroupsOverlapping_K11k2L100D10_2000A100_2000d10_100N2F0.tsv",
        "outTestHP11/IndelGroupsNonoverlapping_K11k2L100D10_2000A100_2000d10_100N2F0.tsv",
        "HP", TRUE,
        "outTestHP11/GenomeData/Genome_1.idlens", "outTestHP11/GenomeData/Genome_2.idlens")
    }
else if (testing == 2)
    {
    args = c("~/Documents/UCDavis/BradyLab/Genomes/kmers/IGGPIPE",
        "outHP14/LCRs_K14k2L400D10_1500.tsv", 100, 2000, 10, 100, 2, 0,
        "MAX",
        "outHP14/IndelGroupsOverlapping_K14k2L400D10_1500A400_1500d50_300N2F0.tsv",
        "outHP14/IndelGroupsNonoverlapping_K14k2L400D10_1500A400_1500d50_300N2F0.tsv",
        "HP", TRUE,
        "outHP14/GenomeData/Genome_1.idlens", "outHP14/GenomeData/Genome_2.idlens")
    }
else stop("Unknown value for 'testing'")
}

NexpectedMin = 15
if (length(args) < NexpectedMin)
    {
    usage = c(
        "Read a data frame of k-mers that have been identified as defining 'LCRs' or",
        "'locally contiguous blocks' between two or more genomes.  Analyze them to find",
        "ones that contain Indel Groups within the appropriate size range for PCR.",
        "Amplicon size limit arguments are approximate, not exact.",
        "Write the resulting LCRs (with additional info) to a file.",
        "",
        "Usage: Rscript findIndelGroups.R <wd> <inFile> <Amin> <Amax> <ADmin> <ADmax> \\",
        "       <NDAmin> <minFlank> <minMax> <outOverlappingIndelGroups> \\",
        "       <outNonoverlappingIndelGroups> <ltrs> <investigate> <idlensFile1> ...",
        "",
        "Arguments:",
        "   <wd>        : Path of R working directory, specify other file paths relative to this.",
        "   <inFile>    : Input file containing the LCR k-mers.",
        "   <Amin>      : Minimum IGG marker amplicon size in any genome.",
        "   <Amax>      : Maximum IGG marker amplicon size in any genome.  Size will always be less than <Dmax> (findLCRs.R)",
        "   <ADmin>     : When smallest amplicon size is <Amin>, this is minimum ADDITIONAL bp's of next larger amplicon.",
        "   <ADmax>     : When largest amplicon size is <Amax>, this is minimum FEWER bp's of next smaller amplicon.",
        "                 Interpolation is linear between <ADmin> and <ADmax>.",
        "   <NDAmin>    : Minimum number of distinct amplicon sizes for each marker.  Set this to the number of genomes if you",
        "                 want distinct amplicon sizes for each genome.  Set it to a smaller value (>= 2) if fewer distinct",
        "                 sizes are acceptable.",
        "   <minFlank>  : Minimum number of common unique LCR k-mers to the left of left-side marker anchor k-mer and to the",
        "                 right of the right-side marker anchor k-mer.  May be 0 for none.",
        "   <minMax>    : Either MIN or MAX to indicate the method to be used to remove overlapping Indel Groups.  MIN means",
        "                 that the smallest Indel Groups are retained over larger ones; MAX means the larger are retained.",
        "   <outOverlappingIndelGroups>    : Name of output file to which overlapping Indel Groups are to be written.",
        "   <outNonoverlappingIndelGroups> : Name of output file to which non-overlapping Indel Groups are to be written.",
        "   <ltrs>        : String of genome designator letters, one letter each in same order as genomes were specified for input file.",
        "   <investigate> : FALSE for normal operation, TRUE for more verbose debugging output",
        "   <idlensFile1> : File containing sequence IDs/lengths for genome 1.",
        "   ...           : Additional sequence IDs/lengths files for genomes 2..Ngenomes."
        )
    for (S in usage)
        catnow(S, "\n", sep="")
    stop("Try again with correct number of arguments")
    }

catnow("findIndelGroups.R arguments:\n")
workingDirectory = args[1]
catnow("  workingDirectory: ", workingDirectory, "\n")
if (!dir.exists(workingDirectory))
    stop("Directory doesn't exist: ", workingDirectory)
setwd(workingDirectory)

inFile = args[2]
catnow("  inFile: ", inFile, "\n")
if (!file.exists(inFile))
    stop("File doesn't exist: ", inFile)

Amin = as.integer(args[3])
catnow("  Amin: ", Amin, "\n")
if (is.na(Amin) || Amin < 1)
    stop("Amin must be >= 1")

Amax = as.integer(args[4])
catnow("  Amax: ", Amax, "\n")
if (is.na(Amax) || Amax < Amin+10)
    stop("Amax must be >= ", (Amin+10))

ADmin = as.integer(args[5])
catnow("  ADmin: ", ADmin, "\n")
if (is.na(ADmin) || ADmin < 1)
    stop("ADmin must be >= 1")

ADmax = as.integer(args[6])
catnow("  ADmax: ", ADmax, "\n")
if (is.na(ADmax) || ADmax < 1)
    stop("ADmax must be >= 1")

NDAmin = as.integer(args[7])
catnow("  NDAmin: ", NDAmin, "\n")
if (is.na(NDAmin) || NDAmin < 2)
    stop("NDAmin must be >= 0")

minFlank = as.integer(args[8])
catnow("  minFlank: ", minFlank, "\n")
if (is.na(minFlank) || minFlank < 0)
    stop("minFlank must be >= 0")

minMax = args[9]
catnow("  minMax: ", minMax, "\n")
if (!minMax %in% c("MIN", "MAX"))
    stop("minMax must be MIN or MAX")

outOverlappingIndelGroups = args[10]
catnow("  outOverlappingIndelGroups: ", outOverlappingIndelGroups, "\n")
if (is.na(outOverlappingIndelGroups))
    stop("outOverlappingIndelGroups must be specified")

outNonoverlappingIndelGroups = args[11]
catnow("  outNonoverlappingIndelGroups: ", outNonoverlappingIndelGroups, "\n")
if (is.na(outNonoverlappingIndelGroups))
    stop("outNonoverlappingIndelGroups must be specified")

genomeLtrs = args[12]
catnow("  genomeLtrs: ", genomeLtrs, "\n")
if (is.na(genomeLtrs))
    stop("genomeLtrs must be specified")

investigate = as.logical(args[13])
catnow("  investigate: ", investigate, "\n")
if (is.na(investigate))
    stop("investigate must be TRUE or FALSE")

idlensFiles = args[14:length(args)]
for (idlensFile in idlensFiles)
    {
    if (!file.exists(idlensFile))
        stop("File doesn't exist: ", idlensFile)
    catnow("  idlensFile: ", idlensFile, "\n")
    }

########################################
# Initialize.
########################################

# The required separation in amplicon size between a genome with an amplicon
# size of A1 and the next-larger amplicon size A2 is a linear function of A1.
# Compute the offset and scale factor for this linear relationship.
# Relationship: amplicon difference D = m * amplicon size A + b
par.m = (ADmax-ADmin)/(Amax-Amin)
par.b = ADmin - par.m*Amin
minAmplDiff = function(Asize) par.m*Asize+par.b # Asize can be a vector.
#minAmplDiff(400)
#minAmplDiff(1500)

# Read LCR k-mer data.
{
if (testing == 0)
    df = read.table(inFile, header=TRUE, row.names=1, sep="\t", stringsAsFactors=FALSE)
else
    df = read.table(inFile, header=TRUE, nrows=100000, row.names=1, sep="\t", stringsAsFactors=FALSE)
}
inv(dim(df), "input data dim")
inv(colnames(df), "input data columns")
inv(head(df), "input data head")

# The k-mer length we are working with.
kmerLen = nchar(rownames(df)[1])
inv(kmerLen, "k")

# Get the genome names from the columns names.
genomes = sub("\\.seqID", "", colnames(df)[grepl("\\.seqID", colnames(df))])
if (length(genomes) < 2)
    stop("Expected to recognize at least two genomes in the data column names")
refGenome = genomes[1]
otherGenomes = genomes[-1]
Ngenomes = length(genomes)
if (NDAmin > Ngenomes)
    stop("NDAmin must be <= number of genomes")
genomeLtrs = unlist(strsplit(genomeLtrs, "", fixed=TRUE))
if (length(genomeLtrs) != Ngenomes)
    stop("Number of genome letters must be equal to the number of genomes")
names(genomeLtrs) = genomes
inv(genomeLtrs, "genomeLtrs")
if (length(idlensFiles) != Ngenomes)
    stop("Number of <idlensFile>'s must be equal to the number of genomes")
names(idlensFiles) = genomes

# Create vectors of column names in df, indexed by genome, for the .seqID, .pos,
# .strand, .contig, and .contigPos columns.
makeColVec = function(S)
    {
    V = paste(genomes, S, sep="")
    names(V) = genomes
    return(V)
    }
idCol = makeColVec(".seqID")
posCol = makeColVec(".pos")
strandCol = makeColVec(".strand")
contigCol = makeColVec(".contig")
contigPosCol = makeColVec(".contigPos")
refIdCol = idCol[refGenome]
refPosCol = posCol[refGenome]
refStrandCol = strandCol[refGenome]
refContigCol = contigCol[refGenome]
refContigPosCol = contigPosCol[refGenome]

########################################
# Find Indel Groups.
########################################

# Assign a name that will be different for each candidate LCR.  The LCRs are
# assigned numbers (LCR column), but the reference genome ID must be included
# because the same set of numbers is reused for each different reference ID.
# Remove the LCR column because it is now useless and even dangerous.
df$LCRname = paste(df[,refIdCol], df$LCR, sep="_")
df = df[, colnames(df) != "LCR"]

# Sort df by the name column so that all k-mers of a candidate LCR will be
# adjacent to one another.
df = df[order(df$LCRname),]

# Add column NkmerFromLCRstart containing the number of this k-mer from the START
# of the LCR, the first one being 1.
numKmersInLCR = tapply(1:nrow(df), df$LCRname, length)
# Following needs unlist() and c() depending on whether all elements of numKmersInLCR are equal
df$NkmerFromLCRstart = unlist(c(sapply(numKmersInLCR, function(N) 1:N)))

# Add column NkmerFromLCRend containing the number of this k-mer from the END of
# the LCR, the last one being 1.
# Following needs unlist() and c() depending on whether all elements of numKmersInLCR are equal
df$NkmerFromLCRend = unlist(c(sapply(numKmersInLCR, function(N) N:1)))

# Add columns "(genome).difPos" giving absolute difference in "(genome).pos"
# from one df row to the next, i.e. is the spacing between one k-mer and the
# next.  All positions are either increasing or decreasing within one genome
# of one LCR.  The first k-mer of each LCR has a value of 0 for these columns.
difPosCol = makeColVec(".difPos")
refDifPosCol = difPosCol[refGenome]
for (genome in genomes)
    {
    df[,difPosCol[genome]] = abs(df[,posCol[genome]] - c(0, df[-nrow(df), posCol[genome]]))
    df[df$NkmerFromLCRstart == 1, difPosCol[genome]] = 0
    }

# Add columns "(genome).difDif" giving absolute difference in "(genome).difPos"
# and (refGenome).difPos.  This is the difference in size of the stretch between
# this k-mer and the next one, in each genome compared to reference genome.  But
# note that either genome may be the larger distance, since we are taking absolute
# value.  The first k-mer of each LCR has a value of 0 for these columns.
difDifCol = makeColVec(".difDif")
for (genome in otherGenomes)
    {
    df[,difDifCol[genome]] = abs(df[,difPosCol[genome]] - df[,refDifPosCol])
    df[df$NkmerFromLCRstart == 1, difPosCol[genome]] = 0
    }

# Remove candidate LCRs that don't have at least 2*minFlank+2 k-mers in them.
inv(length(unique(df$LCRname)), "Number of LCRs before removing those with too few k-mers")
numKmersInLCR = tapply(1:nrow(df), df$LCRname, length)
keepLCRs = names(numKmersInLCR)[numKmersInLCR >= 2*minFlank+2]
df = df[df$LCRname %in% keepLCRs,]
inv(length(unique(df$LCRname)), "Number of LCRs after removing those with too few k-mers")

# Also get rid of candidate LCRs that have insufficient total distance polymorphism
# between the k-mers, excluding the minFlank k-mers on each side.
inv(length(unique(df$LCRname)), "Number of LCRs before removing those with too little polymorphism")
for (genome in otherGenomes)
    {
    # We will avoid slow tapply() and instead use cumsum() on the difDif columns.
    # If we subtract the cumulative sum at the minFlank-th k-mer from the END of
    # each LCR FROM the cumulative sum at the minFlank-th k-mer from the START of
    # each LCR and then add on the difDif value of that same start k-mer, we will
    # have the maximum POSSIBLE size difference between any amplicon of that genome
    # and the amplicon from the reference genome.  If that size difference is < ADmin,
    # then there is no way any amplicon from that LCR can be >= ADmin in that genome,
    # so discard the LCR.
    cumSum = cumsum(df[,difDifCol[genome]])
    difDifAtStartKmer = df[,difDifCol[genome]][df$NkmerFromLCRstart == minFlank+1]
    cumSumThruStartKmer = cumSum[df$NkmerFromLCRstart == minFlank+1]
    cumSumThruEndKmer = cumSum[df$NkmerFromLCRend == minFlank+1]
    cumSumLCRs = cumSumThruEndKmer - cumSumThruStartKmer + difDifAtStartKmer
    tooLittleChangeLCRs = df$LCRname[df$NkmerFromLCRend == 1][cumSumLCRs < ADmin]
    #cat("Number of LCRs of genome", genome, "with too little polymorphism:", length(tooLittleChangeLCRs), "\n")
    df = df[!df$LCRname %in% tooLittleChangeLCRs,,drop=FALSE]
    }
inv(length(unique(df$LCRname)), "Number of LCRs after removing those with too little polymorphism")

# Locate ALL pairs of k-mers that satisfy the requirements for an IGG marker
# (Amin, Amax, ADmin, ADmax).
#
# This is challenging to do efficiently, but I've found a way (not as efficient
# when Ngenomes > 3).  Pay careful attention.
#
# Each LCR is processed as follows (and all LCRs are processed in parallel):
#
# 1. Each row of df is a common unique k-mer of some LCR.  Get the df row indexes
# of ALL df k-mers (for all LCRs) EXCEPT THE FIRST minFlank+1 k-mers and the last
# minFlank k-mers of each LCR into vector rightSideKmers.  These are the candidate
# right-side k-mers for each Indel Group.  Each one can be a right-side k-mer for
# more than one Indel Group, because several different left-side k-mers (from the
# same LCR) might work with one right-side k-mer to give valid Indel Groups.
# Note: the "minFlank" requirements serve to ensure that all Indel Groups have
# the required minimum number of flanking k-mers within the LCR.  Also note that
# the df row index is sufficient because the df row itself contains all needed
# info about the common unique k-mer, including which LCR it belongs to and which
# k-mer it is relative to the first k-mer of the LCR.
#
# 2. Initialize vector leftSideKmers with one lower k-mer index number than the
# indexes in rightSideKmers, i.e. simply subtract 1 from rightSideKmers to form
# leftSideKmers, so that leftSideKmers initially has the index of the k-mer
# immediately preceding the right-side k-mer.  The length of the vectors
# rightSideKmers and leftSideKmers is the same, and corresponding elements form
# a pair.  The two vectors initially hold every possible ADJACENT k-mer pair
# within each LCR that may be adjacent-k-mer marker candidates.  However, there
# will be many more pairs of candidates besides these.  The next steps describe
# how the analysis proceeds using these two vectors.
#
# 3. For each candidate k-mer pair in the two vectors, test to see if the pair
# satisfies ADmin/ADmax/NDAmin considering all genomes, and satisfies Amin and
# Amax.  If so, the pair become an Indel Group, and so their indexes are placed
# in data frame goodPairs, which accumulates the Indel Groups in its columns
# "left" and "right", which simply hold the row indexes into df of the k-mers
# for the left and right sides of the Indel Groups.
#
# 4. After testing the pairs in the two vectors, additional pairs need to be
# tested.  For each right-side k-mer, k-mers that PRECEDE the adjacent k-mer
# in the left-side vector and are in the same LCR are also candidate left-side
# k-mers for that right-side k-mer.  In order to test all these possibilities,
# we will back up the left-side k-mers one k-mer at a time (i.e. one df row at
# a time, i.e. one df row index at a time, i.e. subtract 1 from each left-side
# index) and re-test.  Eventually, the left-side k-mer indexes reach the index
# of the minFlank+1 k-mer in the LCR, and after that is tested, there is no
# longer any need for that PAIR of k-mers to remain in the left-side and
# right-side vectors, because the right-side k-mer has now been tested with
# every possible left-side k-mer.  Since only the left-side indexes are changed,
# the right-side indexes always represent the same right-side k-mer being tested
# against every candidate left-side k-mer.  Once it has been tested against all
# of them, removing that pair of indexes from both left-side and right-side
# vectors removes that right-side k-mer from further testing, but removing the
# left-side index does NOT mean the left-side k-mer is finished being tested,
# the other left-side vector entries for that LCR will eventually be decremented
# to be equal to the same left-side k-mer as is being removed.  And thus, every
# possible k-mer pair within each LCR is tested.  Each test involves the entire
# left-side and right-side vectors, so many tests are done "simultaneously"
# using R's vector capability.  Given this explanation, pay attention to the
# following steps.
#
# 5. If at any time Amax is exceeded for a k-mer pair, that pair is removed from
# the left-side and right-side vectors, and that right-k-mer is no longer
# considered (because every other left-side k-mer before that pair will exceed
# Amax by even more), nor is that left-side k-mer considered any longer (because
# every other right-side k-mer after that pair will exceed Amax by even more).
# As a consequence, the number of tested pairs in an LCR containing N k-mers can
# be far less than N^2/2, if it contains a lot of k-mers separated by more than
# Amax.  Every right-side k-mer beyond those in that LCR, in the right-side
# vector, will eventually also have its left-side k-mer pass over this same k-mer
# boundary (but passing to the next left-side k-mer BEFORE this one, since this
# one has been removed), at which time it too will be eliminated.  The right-side
# vector therefore holds onto right-side k-mers until they have been tested with
# every possible left-side k-mer that could form a valid candidate, and likewise
# the left-side vector holds onto left-side k-mers until they have been decremented
# through all possible left-side k-mers that could form a valid candidate.
#
# 6. After testing the pairs in the left-side and right-side vectors, pairs are
# removed from the vectors if the left-side k-mer was the minFlank+1 k-mer of the
# LCR, because in that case the corresponding right-side k-mer has now been tested
# with every possible valid left-side k-mer, and the left-side k-mer has been
# decremented through all those possible valid left-side k-mers.
#
# 7. ALL left-side k-mer indexes in the left-side vector are decremented by ONE
# so they now index a k-mer one row earlier or one k-mer more towards the start
# of the LCR.  Now return to step 3 to re-test the left-side/right-side k-mer
# pairs to see if they are good candidates.
#
# This algorithm allows all LCRs to be processed in parallel, and the left-side
# and right-side vectors gradually gets smaller and smaller, as potential
# right-side k-mers are eliminated from it along with the corresponding left-side
# k-mer index that reached an index incompatible with that fixed right-side k-mer.
# When the vector sizes reach 0 we are finished.
#

rightSideKmers = (1:nrow(df))[df$NkmerFromLCRstart > minFlank+1 & df$NkmerFromLCRend > minFlank]
leftSideKmers = rightSideKmers-1
if (testing != 0 && FALSE)
    {
    x = sample(1:length(rightSideKmers), 200)
    rightSideKmers = rightSideKmers[x]
    leftSideKmers = leftSideKmers[x]
    }
goodPairs = NULL
progressEverySecs = 5 # Show progress every this many seconds.
lastProgressTime = Sys.time()
loopCount = 0
catnow("Seeking good IGG marker candidate Indel Group regions:\n")
while (length(rightSideKmers) > 0)
    {
    # Calculate distance between left-side k-mer and right-side k-mer in each genome,
    # as an absolute value, unsigned.  Accumulate the distances in matrix difPos,
    # which has one column for each genome and one row for each element of the
    # left-side/right-side vectors.
    # At the same time, check too see if each amplicon resulting from the k-mer pair
    # exceeds the limits of Amin/Amax and record the results in vectors lessThanAmin
    # and moreThanAmax.
    lessThanAmin = rep(FALSE, length(rightSideKmers))
    moreThanAmax = rep(FALSE, length(rightSideKmers))
    difPos = NULL
    for (genome in genomes)
        {
        x = abs(df[rightSideKmers, posCol[genome]] - df[leftSideKmers, posCol[genome]])
        difPos = cbind(difPos, x)
        lessThanAmin[x < Amin] = TRUE
        moreThanAmax[x > Amax] = TRUE
        }

    # Remove k-mer pairs whose distance exceeds Amax from the left-side and right-
    # side vectors.  That right-side k-mer cannot form an amplicon of size < Amax
    # with any other k-mer that lies earlier in the LCR from the one it is paired
    # with now in the left-side vector.
    if (any(moreThanAmax))
        {
        rightSideKmers = rightSideKmers[!moreThanAmax]
        leftSideKmers = leftSideKmers[!moreThanAmax]
        lessThanAmin = lessThanAmin[!moreThanAmax]
        difPos = difPos[!moreThanAmax,,drop=FALSE]
        }
    rm(moreThanAmax)

    # Get k-mer pairs that satisfy Amin (i.e. distance is > Amin), if any, and
    # check to see if they satisfy ADmin/ADmax/NDAmin.
    if (any(!lessThanAmin))
        {
        tmp.rightSideKmers = rightSideKmers[!lessThanAmin]
        tmp.leftSideKmers = leftSideKmers[!lessThanAmin]
        tmp.difPos = difPos[!lessThanAmin,,drop=FALSE]

        # For k-mer pairs that satisfy Amin, we need to see if they satisfy
        # ADmin/ADmax/NDAmin.  An easy first step that will eliminate many pairs
        # from consideration is to see which ones have at least NDAmin distinct
        # difPos values.  If they don't, we can eliminate those from consideration.
        numDistinct = apply(tmp.difPos, 1, length)
        atLeastNDAmin = (numDistinct >= NDAmin)
        # Remove k-mer pairs that have fewer than NDAmin distinct values.
        tmp.rightSideKmers = tmp.rightSideKmers[atLeastNDAmin]
        tmp.leftSideKmers = tmp.leftSideKmers[atLeastNDAmin]
        tmp.difPos = tmp.difPos[atLeastNDAmin,,drop=FALSE]
        numDistinct = numDistinct[atLeastNDAmin]
        rm(atLeastNDAmin)

        # If there are still candidates left, we must test them further.
        if (nrow(tmp.difPos) > 0)
            {
            # Now check the remaining pairs to see if at least NDAmin of them
            # satisfy ADmin/ADmax.  If NDAmin < Ngenomes, it is ok if some don't.
            # This requires sorting the amplicon lengths of each pair, because
            # we need to compare sizes of next-larger amplicon to each amplicon.
            # But sorting takes way too long.  We need a faster way.  When
            # Ngenomes is small, as we expect it to usually be, we can just use
            # min() to do this.  Let's divide this up into three algorithms,
            # for when Ngenomes is 2, 3, or > 3, using pure sorting only for > 3.
            sortCols = function(df, col1, col2)
                {
                swapIt = (df[,col1] > df[,col2])
                t = df[swapIt,col1]
                df[swapIt,col1] = df[swapIt,col2]
                df[swapIt,col2] = t
                return(df)
                }
            if (Ngenomes > 3)
                tmp.sortedLens = t(apply(tmp.difPos, 1, sort, method="quick"))
            else if (Ngenomes == 2)
                tmp.sortedLens = sortCols(tmp.difPos, 1, 2)
            else # (Ngenomes == 3)
                {
                tmp.sortedLens = sortCols(tmp.difPos, 1, 2)
                tmp.sortedLens = sortCols(tmp.sortedLens, 2, 3)
                tmp.sortedLens = sortCols(tmp.sortedLens, 1, 2)
                }

            # Now that they are sorted, generate the differences of successive ones.
            # Then make sure those differences are larger than minAmplDiff(Asize),
            # where Asize is the size of the smaller of the two from which the
            # difference was made.  The number that are at least the minimum size
            # different from the following one must be at least NDAmin-1 (since
            # the last one is not considered).  Record the number in tmp.NDA,
            # where NDA = number distinct amplicons.
            diffs = t(apply(tmp.sortedLens, 1, diff))
            if (ncol(tmp.sortedLens) == 2)
                diffs = t(diffs) # Note: doing it 0 times leaves it as vector, we want matrix.
            minDiff = t(apply(tmp.sortedLens[,-ncol(tmp.sortedLens),drop=FALSE], 1, minAmplDiff))
            if (ncol(tmp.sortedLens) == 2)
                minDiff = t(minDiff) # Note: doing it 0 times leaves it as vector, we want matrix.
            tmp.NDA = 1 + apply(diffs >= minDiff, 1, sum)
            enoughDistinctAmplicons = (tmp.NDA >= NDAmin)

            # Reduce the vectors to contain those k-mer pairs that have enough distinct amplicon sizes.
            tmp.rightSideKmers = tmp.rightSideKmers[enoughDistinctAmplicons]
            tmp.leftSideKmers = tmp.leftSideKmers[enoughDistinctAmplicons]
            tmp.difPos = tmp.difPos[enoughDistinctAmplicons,,drop=FALSE]
            tmp.NDA = tmp.NDA[enoughDistinctAmplicons]

            # Add the successful Indel Groups to goodPairs.
            inv(nrow(tmp.difPos), "    New Indel Groups found")
            if (nrow(tmp.difPos) > 0)
                goodPairs = rbind(goodPairs,
                    data.frame(left=tmp.leftSideKmers, right=tmp.rightSideKmers, NDA=tmp.NDA))
            }
        }

    # Remove k-mer pairs where the left-side k-mer is at k-mer position minFlank+1
    # from the start of the LCR.
    keepKmer = (df$NkmerFromLCRstart[leftSideKmers] > minFlank+1)
    leftSideKmers = leftSideKmers[keepKmer]
    rightSideKmers = rightSideKmers[keepKmer]

    # Advance left-side k-mers to the next k-mer BACK towards START of the LCR.
    leftSideKmers = leftSideKmers - 1

    # Report progress if it is time.
    loopCount = loopCount + 1
    curTime = Sys.time()
    if (difftime(curTime, lastProgressTime, units="secs") >= progressEverySecs)
        {
        catnow("Loop ", format(loopCount, width=4), "  # right-side k-mers remaining to test:",
            format(length(rightSideKmers), width=9), "  # candidates found: ", nrow(goodPairs), "\n")
        lastProgressTime = curTime
        }
    }
inv(sum(!is.na(goodPairs)))

# Retrieve the df rows for the good left and right side k-mers into dfEL and dfER.
dfEL = df[goodPairs$left,]
dfER = df[goodPairs$right,]
rownames(dfEL) = NULL
rownames(dfER) = NULL
dfEL$kmer = rownames(df)[goodPairs$left]
dfER$kmer = rownames(df)[goodPairs$right]

# Using dfEL/dfER as Indel Groups, make data frame dfIGs containing the needed paired
# info from each of those two data frames.  Add column NDA from goodPairs.  Include
# columns for expected amplicon length between exact k-mer start positions for
# each genome.  Add "pct" columns for each genome: percent of pos1 position from
# start of that sequence ID.  We must read the idlens files for sequence lengths.
# Note: .I means Indel Groups, i.e. these column vectors apply to dfIGs, not df.
makeColVec.I = function(S)
    {
    V = paste(genomeLtrs, S, sep="")
    names(V) = genomes
    return(V)
    }
idCol.I = makeColVec.I("id")
pos1Col.I = makeColVec.I("pos1")
pos2Col.I = makeColVec.I("pos2")
str1Col.I = makeColVec.I("s1")
str2Col.I = makeColVec.I("s2")
ctg1Col.I = makeColVec.I("ctg1")
ctg2Col.I = makeColVec.I("ctg2")
kkLenCol.I = makeColVec.I("kkLen")
pctCol.I = makeColVec.I("pct")

dfIGs = data.frame(kmer1=dfEL$kmer, kmer2=dfER$kmer, NDA=goodPairs$NDA, stringsAsFactors=FALSE)
for (genome in genomes)
    {
    dfIdLens = read.table(idlensFiles[genome], header=TRUE, row.names=1, sep="\t", stringsAsFactors=FALSE)

    dfIGs[,idCol.I[genome]] = dfEL[,idCol[genome]]
    dfIGs[,pos1Col.I[genome]] = dfEL[,posCol[genome]]
    dfIGs[,pos2Col.I[genome]] = dfER[,posCol[genome]]
    dfIGs[,str1Col.I[genome]] = dfEL[,strandCol[genome]]
    dfIGs[,str2Col.I[genome]] = dfER[,strandCol[genome]]
    dfIGs[,ctg1Col.I[genome]] = dfEL[,contigCol[genome]]
    dfIGs[,ctg2Col.I[genome]] = dfER[,contigCol[genome]]
    dfIGs[,kkLenCol.I[genome]] = abs(dfIGs[,pos2Col.I[genome]] - dfIGs[,pos1Col.I[genome]] + 1)
    dfIGs[,pctCol.I[genome]] = signif(100*dfIGs[, pos1Col.I[genome]]/dfIdLens[dfIGs[,idCol.I[genome]], "len"], 2)
    }
rownames(dfIGs) = NULL

# Put the data frame in order by reference genome position.
dfIGs = dfIGs[order(dfIGs[, idCol.I[refGenome]], dfIGs[, pos1Col.I[refGenome]]),]
rownames(dfIGs) = NULL
if (nrow(dfIGs) == 0)
    stop("There are no Indel Groups.")

########################################
# Write the overlapping Indel Group data to a file.
########################################

write.table(dfIGs, outOverlappingIndelGroups, row.names=FALSE, quote=FALSE, sep="\t")
# dfIGs = read.table(outOverlappingIndelGroups, header=TRUE, row.names=NULL, sep="\t", stringsAsFactors=FALSE)

catnow("\n")
catnow("Found", nrow(dfIGs), "Indel Group region k-mer pairs with distance between", Amin, "and", Amax, "and\n")
catnow("delta distance at least", ADmin, "at distance", Amin, "and at least", ADmax, "at distance", Amax, "\n")
catnow("in genomes", genomes, " (Indel Group regions can overlap, with each perhaps containing multiple actual Indel Groups)\n")
catnow("\n")

numPerRefSeqId = tapply(1:nrow(dfIGs), dfIGs[,idCol.I[refGenome]], length)
inv(numPerRefSeqId, "Number of Indel Groups per reference genome sequence ID")

########################################
# Remove overlapping Indel Groups, guided by the value of minMax, which is either
# MIN or MAX.  If MIN, retain smallest Indel Groups of each group of overlapping
# ones, and if MAX, retain largest.
########################################

inv(nrow(dfIGs), "Number of Indel Groups including overlapping ones")
df = dfIGs
# Test each genome, one by one, for Indel Group overlaps, and remove Indel Groups
# to get rid of them.
for (genome in genomes)
    {
    inv(genome, "Remove overlaps")
    id.Col = idCol.I[genome]
    pos1.Col = pos1Col.I[genome]
    pos2.Col = pos2Col.I[genome]

    # The method will be:
    #   1. Find, for each Indel Group X with end position X.E, the row index of
    #       the Indel Group Y with starting position Y.S such that X.E > Y.S.
    #       (We count X.E == Y.S as non-overlapping, i.e. Indel Group AB does
    #       not overlap Indel Group BC).  Then, set new column "overlap" TRUE
    #       if that row index is not equal to the row index of Indel Group X.
    #   2. Find the start and end index of each continuous set of TRUE overlap
    #       values, those being the start and end of a group of Indel Groups
    #       that all overlap one another in some manner.
    #   3. For each such start and end index, find the index of the Indel Group
    #       in that row index range which has the smallest (MIN) or largest (MAX)
    #       length.
    #   4. Again for each such start and end index, remove all Indel Groups in
    #       that index range which overlap with the one with the smallest or
    #       largest length.
    #   5. Repeat steps 1-4 until there are no more overlaps.

    # Loop until no more Indel Groups are found to overlap in this genome.
    while (TRUE)
        {
        inv(nrow(df), "Loop with # Indel Groups remaining")

        # Copy start or end position, whichever is smaller (it varies depending on strand)
        # into new column "start", and vice-versa, copy the larger position into new
        # column "end".
        df$start = df[, pos1.Col]
        df$end = df[, pos2.Col]
        pos2IsSmaller = (df[,pos1.Col] > df[,pos2.Col])
        df$start[pos2IsSmaller] = df[pos2IsSmaller, pos2.Col]
        df$end[pos2IsSmaller] = df[pos2IsSmaller, pos1.Col]

        # Add new column "len" equal to length of the Indel Group segments.
        df$len = df$end - df$start + 1

        # Sort by ID and start position.
        df = df[order(df[, id.Col], df$start),]
        N = nrow(df)

        # Find index of Indel Group that each Indel Group overlaps through and
        # put it in column thruIdx.
        df$thruIdx = NA
        for (id in unique(df[, id.Col]))
            {
            thisId = (df[, id.Col] == id)
            df$thruIdx[thisId] = match(TRUE, thisId) - 1 + findInterval(df$end[thisId], df$start[thisId]+1)
            }

        # Set thruIdx of Indel Groups that do not overlap even the next Indel
        # Group to 0.
        df$thruIdx[df$thruIdx == 1:nrow(df)] = 0

        # Get the set of indexes of all Indel Groups which overlap at least one
        # other Indel Group.
        overlapIdxs = sapply(1:nrow(df), function(i)
            {
            if (df$thruIdx[i] == 0)
                return(0)
            return(i:df$thruIdx[i])
            })
        overlapIdxs = unique(unlist(overlapIdxs))
        # Remove index 0, which comes from non-overlapping Indel Groups.
        overlapIdxs = overlapIdxs[overlapIdxs != 0]

        # If no Indel Groups overlap, break out of loop.
        inv(length(overlapIdxs), "Number of overlapping Indel Groups")
        if (length(overlapIdxs) == 0)
            break

        # Set column "overlap" TRUE for each of those overlapping Indel Groups.
        df$overlap = FALSE
        df$overlap[overlapIdxs] = TRUE

        # Get the start and end index of each group of mutually overlapping Indel Groups.
        startOverlap = which(!c(FALSE, df$overlap[-N]) & df$overlap)
        endOverlap = which(df$overlap & !c(df$overlap[-1], FALSE))
        if (length(startOverlap) != length(endOverlap)) stop("Expected equal start/end overlap vectors")

        # Find the index within each group of the Indel Group with the smallest
        # or largest length.  Then get the indexes within the group of the Indel
        # Groups that overlap that Indel Group with the smallest or largest length.
        idxsToRemove = sapply(1:length(startOverlap), function(i)
            {
            idxs = startOverlap[i]:endOverlap[i]
            len.SL = ifelse(minMax == "MIN", min(df$len[idxs]), max(df$len[idxs]))
            idx.SL = idxs[df$len[idxs] == len.SL][1] # If more than one, pick the first.
            return(idxs[!(df$end[idxs] <= df$start[idx.SL]) & !(df$start[idxs] >= df$end[idx.SL]) & idxs != idx.SL])
            })
        idxsToRemove = unlist(idxsToRemove)

        # Remove those Indel Groups.
        df = df[-idxsToRemove,]
        inv(length(idxsToRemove), "Number of Indel Groups deleted")
        inv(nrow(df), "Number of Indel Groups remaining")
        }
    }
dfNoOverlaps = df
dfNoOverlaps = dfNoOverlaps[, !colnames(dfNoOverlaps) %in% c("start", "end", "len", "thruIdx", "overlap")]
inv(nrow(dfNoOverlaps), "Number of Indel Groups with overlapping Indel Groups removed")

########################################
# Write the non-overlapping Indel Group data to a file.
########################################

write.table(dfNoOverlaps, outNonoverlappingIndelGroups, row.names=FALSE, quote=FALSE, sep="\t")
# dfNoOverlaps = read.table(outNonoverlappingIndelGroups, header=TRUE, row.names=NULL, sep="\t", stringsAsFactors=FALSE)

catnow("Number of non-overlapping Indel Groups:", nrow(dfNoOverlaps), "\n")

}
