#######################################################################################
# See usage below for description.
# Author: Ted Toal
# Date: 2013-2015
# Brady Lab, UC Davis
#
# Possible later improvements to pipeline:
# 1. Enhance to look for  1-2 1-3 2-3 3-3 etc. mappings.
# 2. It would make sense for k-mers (and LCB's) to only need to be common between 2
#   genomes rather than all.  But this is a big change to the code. Each k-mer might
#   be NA in one or more genomes. Each LCB would essentially be between only two
#   genomes, always.
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
#testing = 3 # For testing only.
{
if (testing == 0)
    args = commandArgs(TRUE)
else if (testing == 1)
    {
    args = c("~/Documents/UCDavis/BradyLab/Genomes/kmers/IGGPIPE",
        "outTestHP11/Kmers/common.unique.kmers", "HP", 2, 100, 10, 2000,
        "outTestHP11/LCRs_K11k2L100D10_2000.tsv",
        "outTestHP11/BadKmers_K11k2L100D10_2000.tsv", TRUE)
    }
else if (testing == 2)
    {
    args = c("~/Documents/UCDavis/BradyLab/Genomes/kmers/IGGPIPE",
        "outHP14/Kmers/common.unique.kmers", "HP", 2, 400, 10, 1500,
        "outHP14/LCRs_K14k2L400D10_1500.tsv",
        "outHP14/BadKmers_K14k2L400D10_1500.tsv", TRUE)
    }
else if (testing == 3)
    {
    args = c("~/Documents/UCDavis/BradyLab/Genomes/kmers/IGGPIPE",
        "outHPT14/Kmers/common.unique.kmers", "HPT", 2, 300, 5, 1500,
        "outHPT14/LCRs_K14k2L300D5_1500.tsv",
        "outHPT14/BadKmers_K14k2L300D5_1500.tsv", TRUE)
    }
else stop("Unknown value for 'testing'")
}

Nexpected = 10
if (length(args) != Nexpected)
    {
    usage = c(
        "Analyze the positions of a set of k-mers that occur once only in each of N genomes,",
        "with genome1 being the 'reference genome'.  Find blocks of k-mers that are located",
        "on the same contig and near and contiguous to one another in each genome, and call",
        "these 'LCRs' or 'locally contiguous regions'.  Write them out to a file.",
        "",
        "Usage: Rscript findLCRs.R <wd> <inKmerFile> <ltrs> <kMin> <Lmin> <Dmin> <Dmax> \\",
        "       <outLcbFile> <outBadKmers> <investigate>",
        "",
        "Arguments:",
        "   <wd>    : Path of R working directory, specify other file paths relative to this.",
        "   <inKmerFile> : Path of common unique k-mers file with positions of k-mers within",
        "                  each genome, sorted by reference genome position.",
        "   <ltrs>  : String of genome designator letters, one letter each, first letter is",
        "             for genome 1, the reference genome, second letter for genome 2, etc.",
        "   <kMin>  : Minimum number of sequential k-mers to create an LCR.",
        "   <Lmin>  : Minimum length of an LCR in base-pairs.",
        "   <Dmin>  : Minimum distance in bp between two adjacent k-mers in an LCR.  When one",
        "             k-mer is included in an LCR, all k-mers within <Dmin> of it are ignored.",
        "   <Dmax>  : Maximum distance in bp between two adjacent k-mers in an LCR.",
        "   <outLcbFile>  : Name of output file to which LCRs are to be written.",
        "   <outBadKmers> : Name of output file to which rejected k-mers are to be written (for inspection).",
        "   <investigate> : FALSE for normal operation, TRUE for more verbose debugging output"
        )
    for (S in usage)
        catnow(S, "\n", sep="")
    stop("Try again with correct number of arguments")
    }

catnow("findLCRs.R arguments:\n")
workingDirectory = args[1]
catnow("  workingDirectory: ", workingDirectory, "\n")
if (!dir.exists(workingDirectory))
    stop("Directory doesn't exist: ", workingDirectory)
setwd(workingDirectory)

inKmerFile = args[2]
catnow("  inKmerFile: ", inKmerFile, "\n")
if (!file.exists(inKmerFile))
    stop("File doesn't exist: ", inKmerFile)

genomeLtrs = args[3]
catnow("  genomeLtrs: ", genomeLtrs, "\n")
if (is.na(genomeLtrs))
    stop("genomeLtrs must be specified")

kMin = as.integer(args[4])
catnow("  kMin: ", kMin, "\n")
if (is.na(kMin) || kMin < 2)
    stop("kMin must be > 1")

Lmin = as.integer(args[5])
catnow("  Lmin: ", Lmin, "\n")
if (is.na(Lmin) || Lmin < 100)
    stop("Lmin must be >= 100")

Dmin = as.integer(args[6])
catnow("  Dmin: ", Dmin, "\n")
if (is.na(Dmin) || Dmin >= 100)
    stop("Dmin must be < 100")

Dmax = as.integer(args[7])
catnow("  Dmax: ", Dmax, "\n")
if (is.na(Dmax) || Dmax < 100)
    stop("Dmax must be >= 100")

outLcbFile = args[8]
catnow("  outLcbFile: ", outLcbFile, "\n")

outBadKmers = args[9]
catnow("  outBadKmers: ", outBadKmers, "\n")

investigate = as.logical(args[10])
catnow("  investigate: ", investigate, "\n")
if (is.na(investigate))
    stop("investigate must be TRUE or FALSE")

########################################
# Initialization.
########################################

# Split out the genome letters.
genomeLtrs = unlist(strsplit(genomeLtrs, "", fixed=TRUE))
inv(genomeLtrs, "genomeLtrs")
Ngenomes = length(genomeLtrs)
if (Ngenomes < 2)
    stop("There must be at least two letters specified for the genome letters")
refGenome = genomeLtrs[1]
otherGenomes = genomeLtrs[-1]

# Make column names to use for the common unique k-mer file.  The column names
# have the genome name prepended, with a "." separator.
colNames = c()
for (genome in genomeLtrs)
    colNames = c(colNames, paste(genome, c(".seqID", ".pos", ".strand", ".contig", ".contigPos"), sep=""))

# Create vectors of column names in candLCRdf, indexed by genome letter, for the
# .seqID, .pos, .strand, .contig, and .contigPos columns.
makeColVec = function(S)
    {
    V = paste(genomeLtrs, S, sep="")
    names(V) = genomeLtrs
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
# Find LCRs.
########################################

# Smallest number of common unique k-mers to process at one time.  We do not
# process all of them at once because the memory requirements may be too intensive.
# However, we must process all k-mers of a given contig together.  We will do even
# more than that, and process all k-mers of a given reference genome sequence ID
# together.  The input common unique k-mer file is sorted by reference genome
# sequence ID and position.  Therefore, we will read in MIN_KMERS_AT_ONCE
# common unique k-mers from this file, then continue reading in MIN_KMERS_AT_ONCE
# additional k-mers at a time, until we encounter the next sequence ID.  Then, we
# will process all k-mers prededing that next sequence ID, then move on to repeat
# the process.
MIN_KMERS_AT_ONCE = 500000

# Open the input file and read the first MIN_KMERS_AT_ONCE k-mers.
inFile = file(inKmerFile, open="r")
kmerBufferDf = read.table(inFile, header=FALSE, row.names=1, sep="\t", nrows=MIN_KMERS_AT_ONCE, stringsAsFactors=FALSE)
if (nrow(kmerBufferDf) == 0)
    stop("No k-mer data in file ", inKmerFile)
if (ncol(kmerBufferDf) != length(colNames))
    stop("findLCRs expected ", length(colNames), " columns but found ", ncol(kmerBufferDf), ": ",
        paste(colNames, collapse=","), " vs. ", paste(colnames(kmerBufferDf), collapse=",")) 
colnames(kmerBufferDf) = colNames
colClasses = c("character", sapply(1:ncol(kmerBufferDf), function(i) class(kmerBufferDf[,i]))) # first "character" for row names.
NkmersRead = nrow(kmerBufferDf)
catnow("Initial read of", NkmersRead, "k-mers\n")
#catnow("colNames:", paste(colNames, collapse=","), "\n")
#catnow("colClasses:", paste(colClasses, collapse=","), "\n")

# Find out just what k is, i.e. how big are the k-mers?
k = nchar(rownames(kmerBufferDf)[1])
inv(k, "k-mer size k")

# Loop until there are no more k-mers, read additional k-mers from the input file
# until we have one or more complete sequence IDs and at least MIN_KMERS_AT_ONCE
# k-mers to process, or until the end of the input file is reach.  Find LCRs using
# those k-mers.
# Accumulate results in data frames cumLcbDf and cumDiscardDf.
cumLcbDf = NULL
cumDiscardDf = NULL
isEOF = FALSE
nextLCRnum = 1
while (nrow(kmerBufferDf) > 0)
    {
    catnow("Reading data...\n")

    ############################################################################
    # Gather into 'candLCRdf' a set of common unique k-mers to analyze, that
    # includes ALL k-mers on the reference sequence IDs that are included, thus
    # ensuring that k-mers that will be part of an LCR are not analyzed separately
    # from one another.  Use 'kmerBufferDf' as a buffer to hold k-mers read from
    # the input file that are still awaiting transfer to 'candLCRdf'.
    ############################################################################

    # Read more k-mers, we need at least one in kmerBufferDf after we move at
    # least MIN_KMERS_AT_ONCE to candLCRdf.  After this, kmerBufferDf has MORE
    # THAN MIN_KMERS_AT_ONCE rows, unless end of file.
    if (!isEOF)
        {
        cat("  Read...")
        tdf = read.table(inFile, header=FALSE, row.names=1, sep="\t", nrows=MIN_KMERS_AT_ONCE,
            col.names=c("kmers", colNames), colClasses=colClasses, stringsAsFactors=FALSE)
        if (nrow(tdf) == 0)
            isEOF = TRUE
        else
            {
            colnames(tdf) = colNames
            NkmersRead = NkmersRead + nrow(tdf)
            kmerBufferDf = rbind(kmerBufferDf, tdf)
            rm(tdf)
            }
        cat(" NkmersRead =", NkmersRead, " isEOF =", isEOF, "\n")
        }

    # Move at least MIN_KMERS_AT_ONCE k-mers from kmerBufferDf to candLCRdf, but
    # be sure to move ALL k-mers of a sequence ID.  The only time kmerBufferDf
    # becomes empty here is when we reached the end of the file
    N = nrow(kmerBufferDf)
    {
    if (N <= MIN_KMERS_AT_ONCE)
        {
        if (!isEOF) stop("findLCRs expected to be at end of file")
        candLCRdf = kmerBufferDf
        kmerBufferDf = kmerBufferDf[c(),]
        }
    else
        {
        # Transfer at least MIN_KMERS_AT_ONCE, but as many additional k-mers as
        # necessary to reach the end of the sequence ID that occurs at the
        # MIN_KMERS_AT_ONCE position.  If there are not enough sequence IDs,
        # read more lines until end of file or we find the next one.
        lastID = kmerBufferDf[MIN_KMERS_AT_ONCE, refIdCol]
        nextIDidx = MIN_KMERS_AT_ONCE + match(FALSE, lastID == kmerBufferDf[(MIN_KMERS_AT_ONCE+1):N, refIdCol])
        cat("  Read...")
        cnt = 0
        while (!isEOF && is.na(nextIDidx))
            {
            tdf = read.table(inFile, header=FALSE, row.names=1, sep="\t", nrows=MIN_KMERS_AT_ONCE,
                col.names=c("kmers", colNames), colClasses=colClasses, stringsAsFactors=FALSE)
            if (nrow(tdf) == 0)
                isEOF = TRUE
            else
                {
                colnames(tdf) = colNames
                NkmersRead = NkmersRead + nrow(tdf)
                kmerBufferDf = rbind(kmerBufferDf, tdf)
                nextIDidx = match(FALSE, lastID == kmerBufferDf[(N+1):nrow(tdf), refIdCol])
                rm(tdf)
                }
            cnt = cnt+1
            catnow(cnt, " ", sep="")
            }
        cat(" NkmersRead =", NkmersRead, " isEOF =", isEOF, "\n")
        # If nextIDidx is still NA, we reached end of file, so process all remaining k-mers.
        if (is.na(nextIDidx))
            {
            candLCRdf = kmerBufferDf
            kmerBufferDf = kmerBufferDf[c(),]
            }
        # Otherwise, take all k-mers preceding nextIDidx.
        else
            {
            candLCRdf = kmerBufferDf[1:(nextIDidx-1),]
            kmerBufferDf = kmerBufferDf[nextIDidx:nrow(kmerBufferDf),]
            }
        }
    }
    catnow(NkmersRead, "Common unique k-mers read so far\n")
    inv(dim(candLCRdf), "Common unique k-mers dim")
    inv(head(candLCRdf), "Common unique k-mers head")
    catnow("  Processing", nrow(candLCRdf), "common unique k-mers\n")
    IDs = unique(candLCRdf[, refIdCol])
    {
    if (length(IDs) > 5)
        catnow("   Number of reference sequence IDs:", length(IDs), "\n")
    else
        catnow("   Reference sequence IDs:", paste(IDs, collapse=","), "\n")
    }

    ############################################################################
    # Preparations for analysis.
    ############################################################################

    # findMers.cpp reports a k-mer position that is the position of the 5' end of
    # the k-mer on the strand where the k-mer is found, and it reports which strand
    # that is.  To simplify processing, we will convert all positions to be the
    # position of the 5' end of the + strand side of the k-mer, even if the k-mer
    # itself was on the - strand.
    catnow("  Canonicalizing positions\n")
    for (genome in genomeLtrs)
        {
        minusStrand = (candLCRdf[,strandCol[genome]] == "-")
        candLCRdf[minusStrand, posCol[genome]] = candLCRdf[minusStrand, posCol[genome]] - (k-1)
        }

    # Calculate the "strand group" of each k-mer and add column "strandGroup".
    # A strand group is the strand direction relative to the reference genome
    # strand direction, for all genomes.
    # If there are N genomes, there are 2^(N-1) possible strand groups, since
    # we are comparing each genome's strand to the reference genome strand, so
    # the reference genome will always match itself, while the other genomes will
    # either match or mismatch the reference genome.
    # Remember that in actuality each k-mer is on BOTH strands.  A k-mer's
    # canonical sequence (the one lower in alphabetical order) is the one in the
    # data frame row name, and its strand is in the data frame "strand" column.
    # We compare the "strand" column of each genome to the "strand" column of the
    # reference genome to obtain the strand group of that genome.
    # We use the R function xor() (exclusive-or function) to invert the strand
    # group when ref strand is '-' = 1.  We form a number from the strands of
    # each kmer using 0 for '-' strand and 1 for '+' strand, with the number
    # being represented as a binary number.  It is that binary number that is
    # stored in column "strandGroup".  The reference genome strand group number
    # is always 0.
    catnow("  Calculating strand groups\n")
    candLCRdf$strandGroup = 0
    refStrands = ifelse(candLCRdf[,refStrandCol] == "+", 0, 1)
    for (genome in otherGenomes)
        candLCRdf$strandGroup = 2*candLCRdf$strandGroup +
            xor(ifelse(candLCRdf[,strandCol[genome]] == "+", 0, 1), refStrands)

    ############################################################################
    # We want to identify LCRs, which are regions having these properties:
    #   1. All the N k-mers in the LCR are on the same contig of the same
    #       sequence ID within each genome.
    #   2. All the N k-mers of an LCR have consistent strands across genomes (e.g.
    #       if the reference genome's kmer's strands are +++-++----+, then all
    #       other genomes' k-mers are either the SAME strands or the opposite,
    #       ---+--++++-.).  This is the same as saying that all the N k-mers are
    #       in the same strand group.
    #   3. The difference in position of two adjacent k-mers in any genome of
    #       any LCR is no less than Dmin and no more than Dmax.
    #   4. All the N k-mers of each genome are at positions that are either
    #       monotonically increasing or decreasing (may be different for
    #       different genomes).
    #   5. The difference in minimum and maximum position of the N k-mers in any
    #       genome of any LCR is at least Lmin.
    #   6. The LCR consists of a set of N rows of candLCRdf, and N >= kMin.
    ############################################################################

    # For the reference genome, we can apply the above rules easily on the k-mers,
    # which have been sorted by reference genome position.  For the other genomes,
    # it isn't nearly so easy.  Resorting on a different genome's position does
    # not help, because once the reference genome has been done, the k-mers are
    # now part of candidate LCRs, and sorting by another genome scrambles those
    # k-mers so that each candidate LCR's k-mers are all over the place.
    # So, start by applying the above rules to the reference genome, creating sets
    # of k-mers that form candidate LCRs.  After that we'll process the candidate
    # LCRs using the other genomes.
    
    # Items 1, 2, and 3 for reference genome: Create an LCR name column in candLCRdf
    # that will provide a unique name for each LCR.  The name will incorporate
    # the contig number, sequence ID number (except reference genome, for which
    # this is unnecessary since we are looping on reference genome sequence ID),
    # strand group number, and a number that increments each time adjacent
    # k-mers exceed Dmax distance apart.  After forming this name and splitting
    # the k-mers on it, the resulting groups are candidate LCRs in the reference
    # genome.  Within each candidate LCR, test adjacent k-mers for distance
    # closer than Dmin, and eliminate k-mers so all are separated by more.

    # Item 4: this is not a consideration for the reference genome, since it is
    # already sorted by position.  Rather, the reference genome sorting creates
    # a reference order for the candidate LCRs, which the other genomes then
    # must conform to.  If the other genome has a k-mer strand the same as the
    # reference genome, its monotonicity direction must match the reference
    # genome direction, and if opposite, the direction must be opposite.

    # Items 5 and 6: these are easily handled after we create the list of
    # candidate LCRs.

    ############################################################################
    # Sort by reference genome position.
    ############################################################################

    candLCRdf = candLCRdf[order(candLCRdf[,refIdCol], candLCRdf[,refPosCol]),]

    ############################################################################
    # Item 1. All the N k-mers in the LCR are on the same contig and sequence ID
    #       within each genome.
    # Item 2. All the N k-mers of an LCR are in the same strand group.
    # Item 3: identify k-mers that are more than Dmax distance from the preceding
    #       k-mer and split LCRs at these k-mers.
    # Assign each LCR a unique LCR number, then remove k-mers so that no two
    # adjacent k-mers are separated by less than Dmin.
    ############################################################################

    catnow("  Enforcing same contig, same strand group, and Dmax on reference genome\n")

    # Make column 'exceedsDmax' containing a number that increments at each
    # k-mer that exceeds Dmax from the preceding k-mer on the same sequence ID.
    # When the sequence ID changes, we don't care what exceedsDmax is, that
    # boundary will be the start of a new candidate LCR regardless.
    candLCRdf$exceedsDmax = cumsum(c(0, candLCRdf[-1,refPosCol] - candLCRdf[-nrow(candLCRdf),refPosCol]) > Dmax)

    # Now create candidate LCR k-mer groups, which are sets of k-mers that satisfy
    # the above 6 requirements of an LCR in the reference genome, but have not yet
    # been tested in the other genomes.
    # For each row of candLCRdf, assign a name that will be different for each
    # different id, contig, strand group, and exceedsDmax column value in ANY genome.
    candLCRdf$LCRname = ""
    for (genome in genomeLtrs)
        {
        if (genome == refGenome)
            candLCRdf$LCRname = paste("R", candLCRdf[,refIdCol], candLCRdf[,refContigCol], sep=".")
        else
            candLCRdf$LCRname = paste(candLCRdf$LCRname,
                paste(candLCRdf[,idCol[genome]], candLCRdf[,contigCol[genome]], sep="."), sep="_")
        }
    candLCRdf$LCRname = paste(candLCRdf$LCRname, "_SG", candLCRdf$strandGroup, "_D", candLCRdf$exceedsDmax, sep="")
    inv(head(candLCRdf), "LCRname column to divvy by contig/strand group/Dmax separation")
    inv(dim(candLCRdf), "dim of full data")

    # Split candLCRdf into candidate LCR groups.  Column "LCR" is added containing
    # an LCR number, starting from 1, and columns "strandGroup", "exceedsDmax",
    # and "LCRname" are removed.
    LCRvals = as.integer(factor(candLCRdf$LCRname))
    candLCRdf$LCR = nextLCRnum + LCRvals - 1
    nextLCRnum = nextLCRnum + max(LCRvals)
    candLCRdf = candLCRdf[, !colnames(candLCRdf) %in% c("strandGroup", "exceedsDmax", "LCRname")]

    # Remove k-mers so as to ensure that the difference in position of two
    # adjacent k-mers of the same candidate LCR in ANY genome is no less than
    # Dmin.  For each genome, sort the k-mers by the candidate LCR number, and
    # within that by genome sequence ID, contig number, and genome position.
    # After sorting, search for adjacent k-mers with the same candidate LCR and
    # same ID and contig number (since those must be identical for all k-mers
    # within the same LCR) that have a separation distance less than Dmin, and
    # remove the second k-mer of all such pairs, repeating until there are no
    # more such pairs.
    catnow("  Enforcing Dmin\n")
    inv(dim(candLCRdf), "dim of full data before removing < Dmin")
    for (genome in genomeLtrs)
        {
        id.Col = idCol[genome]
        contig.Col = contigCol[genome]
        pos.Col = posCol[genome]
        # Note: in the reference genome there is no need to sort by ID and contig
        # because each candidate LCR number in the LCR column is already composed
        # of only k-mers with the same reference ID and contig.  In the other
        # genomes, however, this isn't true, so we must include those columns in
        # the sort.
        candLCRdf = candLCRdf[order(candLCRdf$LCR, candLCRdf[,id.Col], candLCRdf[,contig.Col], candLCRdf[,pos.Col]),]
        while (TRUE)
            {
            N = nrow(candLCRdf)
            sameLCR = (candLCRdf$LCR[-1] == candLCRdf$LCR[-N])
            sameId = (candLCRdf[-1, id.Col] == candLCRdf[-N, id.Col])
            sameContig = (candLCRdf[-1, contig.Col] == candLCRdf[-N, contig.Col])
            sepLessDmin = (abs(candLCRdf[-1, pos.Col] - candLCRdf[-N, pos.Col]) < Dmin)
            discard = (sameLCR & sameId & sameContig & sepLessDmin)
            if (!any(discard))
                break
            candLCRdf = candLCRdf[c(TRUE, !discard),] # Remove SECOND k-mer of the pair.
            }
        }
    inv(dim(candLCRdf), "dim of full data after removing < Dmin")
    catnow("   ", nrow(candLCRdf), " k-mers remaining\n")
    # If none left, move on to read more k-mers.
    if (nrow(candLCRdf) == 0)
        next

    ############################################################################
    # Items 5 and 6: Discard candidate LCR k-mers that don't satisfy kMin and Lmin.
    ############################################################################

    # Sort again by reference genome position.
    candLCRdf = candLCRdf[order(candLCRdf[,refIdCol], candLCRdf[,refPosCol]),]

    # Remove k-mers not satisfying Lmin.  Save all discarded k-mers in data frame discardDf.
    catnow("  Enforcing LMin on reference genome\n")
    discardDf = NULL
    tooSmall = tapply(candLCRdf[,refPosCol], candLCRdf$LCR, function(refPos) return(diff(range(refPos)) < Lmin))
    tooSmall = names(tooSmall)[tooSmall]
    inv(length(tooSmall), "number of candidate LCRs with bp span < Lmin in ref genome")
    if (length(tooSmall) > 0)
        {
        discardDf = rbind(discardDf, candLCRdf[candLCRdf$LCR %in% tooSmall,])
        candLCRdf = candLCRdf[!candLCRdf$LCR %in% tooSmall,]
        }
    inv(dim(candLCRdf), "dim of full data after removing < Lmin")
    inv(dim(discardDf), "dim of discarded k-mers")
    catnow("   ", nrow(candLCRdf), " k-mers remaining\n")
    # If none left, move on to read more k-mers.
    if (nrow(candLCRdf) == 0)
        next

    # Remove k-mers not satisfying kMin.  Save all discarded k-mers in data frame discardDf.
    catnow("  Enforcing kMin on reference genome\n")
    tooFewKmers = tapply(1:nrow(candLCRdf), candLCRdf$LCR, function(ii) length(ii) < kMin)
    tooFewKmers = names(tooFewKmers)[tooFewKmers]
    inv(length(tooFewKmers), "number of candidate LCRs with < kMin k-mers")
    if (length(tooFewKmers) > 0)
        {
        discardDf = rbind(discardDf, candLCRdf[candLCRdf$LCR %in% tooFewKmers,])
        candLCRdf = candLCRdf[!candLCRdf$LCR %in% tooFewKmers,]
        }
    inv(dim(candLCRdf), "dim of full data after removing < kMin")
    inv(dim(discardDf), "dim of discarded k-mers")
    catnow("   ", nrow(candLCRdf), " k-mers remaining\n")
    # If none left, move on to read more k-mers.
    if (nrow(candLCRdf) == 0)
        next

    ############################################################################
    # Reference genome finished, we have candidate LCRs in the reference genome.
    ############################################################################

    ############################################################################
    # Now move on to testing the remaining candidate LCRs in the other genomes.
    # Things that can disturb the required LCR conditions along a contig include
    # inversions, translocations, indels, too much spacing between k-mers, and
    # too-small LCRs.  As a result, a single contig can contain multiple LCRs,
    # separated at places where the above conditions are violated.  But there
    # can also be single LCRs that contain, within the middle of them, k-mers
    # and even good INNER LCRs that, once removed, leave a good outer LCR.  Note
    # that disturbances don't cause k-mers to be eliminated because they don't fit
    # the pattern, but rather they cause them to be excluded from that LCR.  As long
    # as all the k-mers included in the LCR meet the conditions, it is a valid LCR.
    # The k-mers excluded from it may be part of a different LCR, or they may end
    # up not being part of any LCR.  A given k-mer, though, is only a member of at
    # most one LCR.  The excluded k-mers are moved to a deferred k-mers list and
    # and later back to the LCR pool being tested, to retest them to see if they
    # satisfy the requirements of an LCR.
    # It may in fact possible for two LCRs to cover the same region: suppose
    # an indel is present in a region; the indel itself might be one LCR while
    # the region on both sides of it is another.  This could happen only if the
    # minimum LCR length does not exceed Dmax.  If it does, such an indel would
    # necessarily break up the surround LCR into two LCRs.
    # The challenge, then, is to find subsets of k-mers that satisfy these
    # conditions to form an LCR.  Once an LCR is formed, those k-mers are
    # eliminated from the pool of candidates.
    # Presumably, the first pass above using the reference genome took care of the
    # majority of problems, but there will still be problem areas left where the
    # genomes have diverged.  The following properties must still be tested for
    # each genome in "otherGenomes", for all candidate LCRs:
    #
    #   1. If the number of k-mers in a candidate LCR group is less than kMin,
    #       discard that candidate LCR.
    #
    #   2. If the difference in minimum and maximum position of the candidate LCR
    #       k-mers in any genome of any LCR is less than Lmin, discard it.
    #
    #   3. The k-mers are in monotonically increasing order by reference genome
    #       position because we sorted them that way.  In the other genomes, the
    #       k-mers had also better be monotonically increasing OR decreasing, if
    #       the k-mers are to be a valid LCR.  If the k-mer strand of the genome
    #       in question matches the k-mer strand of the reference genome, then
    #       the genome in question should also be monotonically increasing.  If
    #       they mismatch, it should be monotonically decreasing.  It doesn't
    #       matter which k-mer's strand is tested for this, since all k-mers of
    #       the LCR are part of the same strand group and so have the same
    #       relationship to the reference k-mer strand.
    #
    #       Suppose the monotonicity is incorrect for two adjacent k-mers of the
    #       candidate LCR in any genome.  Note that when an inversion occurs,
    #       the strand of the k-mers inverts, so we would already have separated
    #       those k-mers based on strand group.  But monotonicity can also be
    #       wrong when a non-inverting translocation occurs.  Let's look more
    #       closely at this to understand it.  For simplificity, assume two
    #       genomes:
    #
    #         non-inv xloc: genome 2:  a b c d e f g h i j A B C D E k l
    #                     ref genome:  a b A B C D E c d e f g h i j k l
    #           sign:     ref genome:   + + + + + + + + + + + + + + + +
    #   genome 2 in ref genome order:   + + + + + + - + + + + + + + + +
    #           bigger change at:         ^                         ^
    #
    #       Here, the + and - for genome 2 are taking consecutive k-mer pairs
    #       as the k-mers are shown in the ref genome line, and asking the
    #       question, what is their position difference IN GENOME 2, which
    #       can be determined by looking at where they lie in the genome 2 row.
    #       We see only a single place where two adjacent k-mers have the wrong
    #       monotonicity.  There are two other places where the size of the
    #       change may be large, and may be larger than Dmax.  That case is
    #       handled in the next step, 4, which comes after this because our
    #       resolution of monotonicity may fix the large position difference.
    #
    #       Now consider the case where the direction of movement was the
    #       opposite:
    #
    #         non-inv xloc: genome 2:  a b A B C D E c d e f g h i j k l
    #                     ref genome:  a b c d e f g h i j A B C D E k l
    #           sign:     ref genome:   + + + + + + + + + + + + + + + +
    #   genome 2 in ref genome order:   + + + + + + + + + - + + + + + +
    #           bigger change at:         ^                         ^
    #
    #       We see a similar situation, but the position of the monotonicity
    #       error is at the start rather than end of the translocated sequence.
    #
    #       Now, what should we do?  Monotonicity is funny in that it depends on
    #       two k-mers, not just one.  So, a k-mer involved in wrong monotonicity
    #       in one context may not be that way in another.  We need to be careful
    #       though, about not setting such k-mers aside and trying them later, ad
    #       infinitum.  This can be avoided by ensuring that we never recycle ALL
    #       of the k-mers of the LCR.  We will choose to keep the k-mer on the
    #       left side of the monotonicity problem, so at least one k-mer will be
    #       processed through to the end, reducing overall k-mer count by at least
    #       one.
    #
    #       Again, what to do?  In the above cases it is clear that setting a
    #       single k-mer aside is insufficient.  We need a simple solution that
    #       may not be best in all cases.  If we analyze from the first right-side
    #       k-mer of a monotonicity error rightwards, looking for the first k-mer
    #       that will not cause a monotonicity error, we would do this:
    #
    #         non-inv xloc: genome 2:  a b c d e f g h i j A B C D E k l
    #                     ref genome:  a b A B C D E c d e f g h i j k l
    #           sign:     ref genome:   + + + + + + + + + + + + + + + +
    #   genome 2 in ref genome order:   + + + + + + - + + + + + + + + +
    #           bigger change at:         ^                         ^
    #         set aside these k-mers:                c d e f g h i j
    #
    #         non-inv xloc: genome 2:  a b A B C D E c d e f g h i j k l
    #                     ref genome:  a b c d e f g h i j A B C D E k l
    #           sign:     ref genome:   + + + + + + + + + + + + + + + +
    #   genome 2 in ref genome order:   + + + + + + + + + - + + + + + +
    #           bigger change at:         ^                         ^
    #         set aside these k-mers:                      A B C D E
    #
    #       In both cases that is a satisfactory solution.  There are some cases
    #       though where we can do better.  It will often be the case that a
    #       single wacky k-mer causes a monotonicity problem.  In that case,
    #       removing that k-mer will fix the problem.  With the above solution,
    #       we will do exactly that IF the right-side k-mer is the wacky one,
    #       but not if the left side is.  So, suppose we ALSO check to see if
    #       removal of the left-side k-mer would remove the monotonicity problem,
    #       and if so, defer that k-mer and RECHECK the remaining k-mers.  We'll
    #       do that first, and if removal of the left-side k-mer doesn't fix the
    #       problem, we'll use the above solution.
    #
    #       In the second part of the solution where we look to the right, set
    #       those k-mers aside in the a deferred k-mers list, then RECHECK the
    #       remaining k-mers for a good LCR.
    #
    #       Suppose there are a bunch of k-mers that were just random common unique
    #       k-mers and their monotonicity is willy-nilly?  It doesn't really matter,
    #       we will still be tearing out some troublesome k-mers each time, retrying
    #       until either successful or the LCR becomes empty, then retrying the
    #       troublesome k-mers in the same manner.  They may cause computation to
    #       take longer, but they will eventually get kicked out.
    #
    #   4. If the difference in position of two adjacent k-mers of the candidate
    #       LCR in any genome is more than Dmax, move the right k-mer, plus any
    #       immediately following k-mers that are also more than Dmax away from
    #       the left k-mer, to a deferred k-mers list to be processed subsequently.
    #       This was made step 4, because in step 3 some large changes in position
    #       will be eliminated when monotonicity problems are corrected.
    #
    #   5. Repeat steps 1-4 over and over with the k-mers of each candidate LCR
    #       (recycling the deferred k-mers created in the process into new candidate
    #       LCRs).  Each time an LCR is modified by moving a k-mer to the auxiliary
    #       LCR, restart at step 1 for that LCR.  When all steps 1-4 are executed
    #       with no change, the LCR is good and is accepted and a sequential LCR
    #       number is assigned to it.
    ############################################################################

    ############################################################################
    # Define a function to apply the discontinuity testing.  LL is a list with
    # these members:
    #   lcr: data frame containing k-mers being tested, all assigned to same LCR,
    #       and sorted in reference genome order.
    #   deferred: NULL or data frame of k-mers that have been removed (via
    #       previous calls to this with the same LCR number) from data frame lcr
    #       because they don't satisfy requirements, but that may form a separate
    #       LCR (or more than one).
    #   discarded: NULL or data frame of k-mers that have been rejected as members
    #       of an LCR and removed from data frame lcr.
    # On return, LL$lcr has been modified, with non-compatible k-mers removed and
    # placed into either LL$deferred or LL$discarded.
    # Return value is a list with the same members and additional members:
    #   retry: TRUE if function is to be called again with same value of lcr to
    #       retest modified candidate LCR, FALSE if finished with this LCR.
    #   reason: short word giving exit status:
    #       kMin: there were fewer than kMin k-mers, k-mers discarded
    #       Lmin: maximum base pair span was less than Lmin, k-mers discarded
    #       monotonicity1: single k-mer monotonicity wrong, moved to deferred
    #       monotonicity2: k-mer monotonicity error, one or more moved to deferred
    #       Dmax: adjacent k-mer distance was greater than Dmax, moved to deferred
    #       okay: no problems detected
    # Note that even though we already tested the reference genome above, we must
    # retest it here because we keep reforming the candidate LCRs with subsets of
    # the previously tested ones.
    ############################################################################
    testCandidateLCR = function(LL)
        {
        LL$retry = FALSE
        Nkmers = nrow(LL$lcr)
        #catnow("Nkmers = ", Nkmers, "\n")

        #   1. If the number of k-mers in a candidate LCR group is less than kMin,
        #       discard that candidate LCR.
        if (Nkmers < kMin)
            {
            LL$discarded = rbind(LL$discarded, LL$lcr)
            LL$lcr = NULL
            LL$reason = "kMin"
            return(LL)
            }

        #   2. If the difference in minimum and maximum position of the candidate LCR
        #       k-mers in any genome of any LCR is less than Lmin, discard it.
        for (genome in genomeLtrs)
            if (diff(range(LL$lcr[,contigPosCol[genome]])) < Lmin)
                {
                LL$discarded = rbind(LL$discarded, LL$lcr)
                LL$lcr = NULL
                LL$reason = "Lmin"
                return(LL)
                }

        #   3. Identify places where monotonicity is incorrect in a genome, and
        #       check to see if removal of the left-side k-mer at that position
        #       fixes the problem, and if so, discard that k-mer and retry.  If
        #       not, move the first k-mer on the right side of the monotonicity
        #       error, and all immediately following ones which would also lead
        #       to a monotonicity error if the k-mers to its left were removed,
        #       to the deferred k-mers list to be processed subsequently.
        #       The expected monotonicity direction is increasing if the genome
        #       strand matches the reference genome strand, and decreasing if not.
        #       Since monotonicity is a property of two adjacent k-mers, there is
        #       one fewer monotonicity error flags than there are k-mers.
        for (genome in otherGenomes)
            {
            expectedMonotonicity = (LL$lcr[1, strandCol[genome]] == LL$lcr[1, refStrandCol])
            actualKmerMonotonicity =
                (LL$lcr[-1, posCol[genome]] > LL$lcr[-Nkmers, posCol[genome]])
            kmerMonotonicityError = (actualKmerMonotonicity != expectedMonotonicity)
            if (any(kmerMonotonicityError))
                {
                leftKmer = which(kmerMonotonicityError)[1]
                rightKmer = leftKmer + 1
                # Does removal of leftKmer fix the problem?
                leftFix = FALSE
                if (Nkmers < 3)
                    leftFix = TRUE
                else if (leftKmer == 1)
                    {
                    akm = (LL$lcr[3, posCol[genome]] > LL$lcr[2, posCol[genome]])
                    if (akm == expectedMonotonicity)
                        leftFix = TRUE
                    }
                else
                    {
                    akm = (LL$lcr[3, posCol[genome]] > LL$lcr[1, posCol[genome]])
                    if (akm == expectedMonotonicity)
                        leftFix = TRUE
                    }
                if (leftFix)
                    {
                    LL$deferred = rbind(LL$deferred, LL$lcr[leftKmer,])
                    LL$lcr = LL$lcr[-leftKmer,]
                    LL$retry = TRUE
                    LL$reason = "monotonicity1"
                    #cat("Monotonicity error 1: ", genome, " ", leftKmer, "\n", sep="")
                    return(LL)
                    }

                # Calculate monotonicity error of every k-mer to the right of leftKmer
                # as if that k-mer is the next k-mer after leftKmer.
                akm = (LL$lcr[(rightKmer:Nkmers), posCol[genome]] > LL$lcr[leftKmer, posCol[genome]])
                kme = (akm != expectedMonotonicity)

                # What is the first one that is not in error, if any?
                firstRight = leftKmer + match(FALSE, kme)
                if (is.na(firstRight))
                    lastWrong = Nkmers
                else
                    lastWrong = firstRight - 1
                #if (!expectedMonotonicity && lastWrong-rightKmer > 5 && leftKmer > 1) stop("check monotonicity error processing")
                LL$deferred = rbind(LL$deferred, LL$lcr[(rightKmer:lastWrong),])
                LL$lcr = LL$lcr[-(rightKmer:lastWrong),]
                LL$retry = TRUE
                LL$reason = "monotonicity2"
                #cat("Monotonicity error 2: ", genome, " ", rightKmer, ":", lastWrong, "\n", sep="")
                return(LL)
                }
            }

        #   4. If the difference in position of two adjacent k-mers of the candidate
        #       LCR in any genome is more than Dmax, move the right k-mer, and any
        #       following k-mers that are also more than Dmax from the left k-mer,
        #       to the deferred k-mers list to be processed subsequently.
        for (genome in otherGenomes)
            {
            tooBig = which(abs(LL$lcr[, posCol[genome]][-1] - LL$lcr[, posCol[genome]][-Nkmers]) > Dmax)
            if (length(tooBig) > 0)
                {
                leftKmer = tooBig[1]
                rightKmer = leftKmer + 1
                # Calculate difference in position of every k-mer to the right of leftKmer
                # as if that k-mer is the next k-mer after leftKmer.
                isTooBig = (abs(LL$lcr[rightKmer:Nkmers, posCol[genome]] -
                        LL$lcr[leftKmer, posCol[genome]]) > Dmax)
                firstNotTooBig = leftKmer + match(FALSE, isTooBig)
                if (is.na(firstNotTooBig))
                    lastTooBig = Nkmers
                else
                    lastTooBig = firstNotTooBig - 1
                LL$deferred = rbind(LL$deferred, LL$lcr[rightKmer:lastTooBig,])
                LL$lcr = LL$lcr[-(rightKmer:lastTooBig),]
                LL$retry = TRUE
                LL$reason = "Dmax"
                return(LL)
                }
            }

        LL$reason = "okay"
        return(LL)
        }

    ############################################################################
    # Now, for each candidate LCR in candLCRdf, test it repeatedly, removing
    # incompatible k-mers from it over to a deferred k-mers data frame,
    # repeating until no more k-mers are removed and the LCR is known to satisfy
    # the constraints (kMin, Lmin, Dmax, monotonicity), then incorporate the
    # results into new result data frame lcrDf, and if some k-mers were found
    # incompatible with the LCR and moved to the deferred k-mers data frame,
    # rbind that deferred k-mers data frame back to candLCRdf.  Thus, candLCRdf
    # will shrink as candidate LCRs are either verified and moved to lcrDf or
    # discarded and moved to discardDf.
    # This loop is slow, and in order to try to speed it up, I'm using some
    # smaller temporary data frames to store data, so we aren't passing large
    # data frames to functions.  These are:
    #   dft: small buffer of k-mers of candidate LCRs, reloaded with N.reloaod
    #       k-mers from candLCRdf whenever it becomes empty.  Process these
    #       k-mers, not the ones in candLCRdf.
    #   lcrDft: k-mers of newly confirmed LCRs, appended to lcrDf every reload.
    #   discardDft: k-mers incompatible with LCRs, appended to discardDf every reload.
    #   deferredDft: k-mers that have been deferred and have a new LCR number,
    #       and will be appended to candLCRdf every reload.
    # Keep track of stats in stats vector.
    ############################################################################

    catnow("  Testing candidate LCRs: discarding LCRs when kMin or Lmin exceeded in any genome,\n")
    catnow("  making new candidate LCRs using k-mers following one where Dmax is exceeded or all\n")
    catnow("  k-mers that violate monotonicity, and accepting LCRs when all k-Mers satisfy the\n")
    catnow("  kMin, Lmin, Dmax, and monotonicity constraints.\n")

    # Note: nextLCRnum has the next available LCR number to use when assigning new
    # LCR numbers.
    inv(nextLCRnum, "nextLCRnum")

    N.reload = 100 # Number of LCR numbers to load into dft each time.
    logEveryN = 500
    lastLog = 0
    lcrDf = NULL # # K-mers of candidate LCRs that pass the test and are declared actual LCRs.
    dft = candLCRdf[integer(0),] # This will hold k-mers for N.reload candidate LCRs being tested.
    lcrDft = NULL # K-mers of candidate LCRs that pass the test and are declared actual LCRs, transferred to lcrDf at each reload.
    discardDft = NULL # K-mers of discarded candidate LCRs, transferred to discardDf at each reload.
    deferredDft = NULL # K-mers assigned to a new LCR ID number and awaiting reprocessing, transferred back to candLCRdf at each reload.
    stats = integer() # We'll count number of times each reason code is returned by testCandidateLCR().
    LCRs = unique(candLCRdf$LCR) # These are the LCR ID numbers in candLCRdf that we will test.  New IDs assigned to deferred LCRs are added to this.
    idx = 1 # Index of LCR ID in LCRs[] that is currently being tested.
    catnow("  Number of candidate LCRs: ", length(LCRs), " Number of k-mers: ", nrow(candLCRdf), "\n")

    # Loop processing one LCR at a time, the LCR whose ID (in column "LCR") is LCRs[idx].
    while (idx <= length(LCRs))
        {
        #catnow("idx = ", idx, "length(LCRs) = ", length(LCRs), "\n")
        # Reload temporary data frames if dft is empty.
        if (nrow(dft) == 0)
            {
            #stop()
            # Transfer deferred k-mers back to candLCRdf.
            candLCRdf = rbind(candLCRdf, deferredDft)
            deferredDft = NULL
                            # candLCRdf must be sorted again!!!
                            # Sort by reference genome position.
                            #candLCRdf = candLCRdf[order(candLCRdf[,refPosCol]),]
            # Transfer discarded k-mers to discardDf.
            discardDf = rbind(discardDf, discardDft)
            discardDft = NULL
            # Transfer k-mers of validated LCRs to lcrDf.
            lcrDf = rbind(lcrDf, lcrDft)
            lcrDft = NULL
            # Compute index into LCRs[] of last LCR whose k-mers are to be loaded into dft during this round.
            lastIdx = idx+N.reload-1
            if (lastIdx > length(LCRs))
                lastIdx = length(LCRs)
            # Determine which k-mers in candLCRdf are to be transferred, then transfer them.
            xferKmers = (candLCRdf$LCR %in% LCRs[idx:lastIdx])
            dft = candLCRdf[xferKmers,]
            candLCRdf = candLCRdf[!xferKmers,]
                            # if (nrow(dft) > 0 && any(order(dft[,refPosCol]) != 1:nrow(dft))) stop("Order")
            # Logging.
            lenLCRs = length(LCRs)
            if (idx - lastLog >= logEveryN)
                {
                catnow("    At candidate LCR #", idx, " of ", lenLCRs, " leaving ",
                    lenLCRs-idx, " with ", nrow(candLCRdf)+nrow(dft), " k-mers, have ",
                    length(unique(lcrDf$LCR)), " good LCRs with ", nrow(lcrDf), " k-mers.\n", sep="")
                lastLog = idx
                }
            }
        # Get candidate LCR k-mers for the LCR with LCR ID = LCRs[idx] into list LL and remove them from dft.
        #catnow("Get LCR\n")
        lcrKmers = (dft$LCR == LCRs[idx])
        LL = list(lcr=dft[lcrKmers,], deferred=NULL, discarded=NULL)
        dft = dft[!lcrKmers,]
        #catnow("Loop calling testCandidateLCR\n")
        repeat
            {
            #catnow("Call\n")
            LL = testCandidateLCR(LL)
            # Update stats.
            #catnow(" end call, stat: ", LL$reason, "\n")
            if (is.na(stats[LL$reason]))
                stats[LL$reason] = 1
            else
                stats[LL$reason] = stats[LL$reason] + 1
            # Exit if retry not flagged.
            if (!LL$retry)
                break
            }
        # If the candidate k-mer list is not NULL, append it to lcrDft.
        if (!is.null(LL$lcr))
            {
            #catnow(" lcrDft\n")
            lcrDft = rbind(lcrDft, LL$lcr)
            }
        # If some k-mers appear in deferred, append them to deferredDft using a new LCR number.
        if (!is.null(LL$deferred))
            {
            #catnow(" deferredDft\n")
            LCRs = c(LCRs, nextLCRnum)
            LL$deferred$LCR = nextLCRnum
            deferredDft = rbind(deferredDft, LL$deferred)
            nextLCRnum = nextLCRnum + 1
            }
        # If some k-mers were discarded, append them to discardDft.
        if (!is.null(LL$discarded))
            {
            #catnow(" discardDft\n")
            discardDft = rbind(discardDft, LL$discarded)
            }
        # Next LCR index.
        idx = idx + 1
        }
    inv(stats, "LCR verification statistics")

    # Transfer final results from temporary data frames to permanent ones.
    if (nrow(dft) != 0) stop("Expected empty dft")
    if (nrow(candLCRdf) != 0) stop("Expected empty candLCRdf")
    if (!is.null(deferredDft)) stop("Expected null deferredDft")
    discardDf = rbind(discardDf, discardDft)
    lcrDf = rbind(lcrDf, lcrDft)
    rm(dft, lcrDft, discardDft)

    catnow("\n  Found", length(unique(lcrDf$LCR)), "LCRs containing", nrow(lcrDf), "k-mers\n")

    # Accumulate results.
    cumLcbDf = rbind(cumLcbDf, lcrDf)
    cumDiscardDf = rbind(cumDiscardDf, discardDf)
    rm(candLCRdf, lcrDf, discardDf)
    catnow("\n  Cumulative total of", length(unique(cumLcbDf$LCR)), "LCRs containing", nrow(cumLcbDf), "k-mers\n")
    }

# Close the input file.
close(inFile)

# Check LCR count.
if (nrow(cumLcbDf) == 0)
    stop("No LCRs found.")

# Sort by LCR number and reference genome position.
cumLcbDf = cumLcbDf[order(cumLcbDf$LCR, cumLcbDf[, refPosCol]),]

################################################################################
# Write results to files.
################################################################################

write.table(cumLcbDf, outLcbFile, quote=FALSE, sep="\t")
# cumLcbDf = read.table(outLcbFile, header=TRUE, row.names=1, sep="\t", stringsAsFactors=FALSE)
write.table(cumDiscardDf, outBadKmers, quote=FALSE, sep="\t")
# cumDiscardDf = read.table(outBadKmers, header=TRUE, row.names=1, sep="\t", stringsAsFactors=FALSE)

catnow("\nFound a total of", length(unique(cumLcbDf$LCR)), "LCRs containing",
    nrow(cumLcbDf), "k-mers common to the ", Ngenomes, "genomes\n")
}
