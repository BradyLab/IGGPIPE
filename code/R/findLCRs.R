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
testing = FALSE
#testing = TRUE # For testing only.
{
if (!testing)
    args = commandArgs(TRUE)
else
    {
    args = c("~/Documents/UCDavis/BradyLab/Genomes/kmers/IGGPIPE",
        "outTestHP11/Kmers/Split/Genome_", "HP", 2, 100, 10, 2000,
        "outTestHP11/LCRs_K11Km2Lm100Dm10Dx2000.tsv",
        "outTestHP11/BadKmers_K11Km2Lm100Dm10Dx2000.tsv", TRUE)
    }
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
        "Usage: Rscript findLCRs.R <wd> <inPfx> <ltrs> <kmin> <Lmin> <Dmin> <Dmax> \\",
        "       <outLcbFile> <outBadKmers> <investigate>",
        "",
        "Arguments:",
        "   <wd>    : Path of R working directory, specify other file paths relative to this.",
        "   <inPfx> : Prefix, including directory, of input files of k-mers and their positions.",
        "   <ltrs>  : String of genome designator letters, one letter each, first letter is",
        "             for genome 1, the reference genome, second letter for genome 2, etc.",
        "   <kmin>  : Minimum number of sequential k-mers to create an LCR.",
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

inPfx = args[2]
catnow("  inPfx: ", inPfx, "\n")

genomeLtrs = args[3]
catnow("  genomeLtrs: ", genomeLtrs, "\n")
if (is.na(genomeLtrs))
    stop("genomeLtrs must be specified")

kmin = as.integer(args[4])
catnow("  kmin: ", kmin, "\n")
if (is.na(kmin) || kmin < 2)
    stop("kmin must be > 1")

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

# Split the inPfx argument into a directory name and a filename prefix.
inDir = sub("/[^/]*$", "", inPfx)
inPfx = sub("^.*/([^/]*$)", "\\1", inPfx)
inv(inDir, "inDir")
inv(inPfx, "inPfx")

# Get a list of all filenames in the input directory whose names begin with the
# input prefix followed by genome number 1 (=reference genome), "_", a sequence
# ID name, and finally end in ".isect.split".
allFiles = list.files(inDir)
RE = paste("^", inPfx, "1_(.*)\\.isect.split$", sep="")
refFiles = allFiles[grepl(RE, allFiles)]
if (length(refFiles) == 0)
    stop("No reference genome files found in ", inDir, ", with names starting with ", inPfx, "_1_")
# Extract the reference genome IDs from the reference genome file names.
refIDs = sub(RE, "\\1", refFiles)
refIDs = sort(refIDs)

# For each reference ID, make a vector of input filenames, one per genome, in the
# same order as the 'genomeLtrs' vector, with the reference genome first.
inFiles = list()
for (ID in refIDs)
    {
    inFiles[[ID]] = paste(inDir, PATHSEP, inPfx, 1:Ngenomes, "_", ID, ".isect.split", sep="")
    names(inFiles[[ID]]) = genomeLtrs
    }

# Create vectors of column names in df, indexed by genome letter, for the .seqID, .pos,
# .strand, .contig, and .contigPos columns.
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

# Loop for each reference sequence ID, read the data files, and find LCRs.
# Accumulate results in data frames cumLcbDf and cumDiscardDf.
cumLcbDf = NULL
cumDiscardDf = NULL
for (ID in refIDs)
    {
    catnow("\nDoing reference ID", ID, "\n")

    # Read k-mer data for each genome for this ID into data frames and bind them
    # together column-wise to create data frame df.  Make sure each genome has
    # exactly the same list of k-mers in it.  The column names have the genome
    # name prepended, with a "." separator.
    catnow("  Reading data\n")
    df = NULL
    for (genome in genomeLtrs)
        {
        inv(genome, "Genome")
        filename = inFiles[[ID]][[genome]]
        inv(filename, "Input file")
        dfTmp = read.table(filename, header=TRUE, row.names=1, sep="\t", stringsAsFactors=FALSE)
        # Sort the rows by k-mer sequence, which are the row names.
        dfTmp = dfTmp[order(rownames(dfTmp)),]
        inv(dim(dfTmp), "input data dim")
        inv(colnames(dfTmp), "input data columns")
        inv(head(dfTmp), "input data head")
        colnames(dfTmp) = paste(genome, colnames(dfTmp), sep=".")
        if (genome == refGenome)
            {
            saveDim = dim(dfTmp)
            df = dfTmp
            }
        else
            {
            if (any(dim(dfTmp) != saveDim))
                stop("Genome ", genome, " data not same size as reference genome as expected")
            if (any(rownames(df) != rownames(dfTmp)))
                stop("Genome ", genome, " kmer row names not same as reference genome as expected")
            df = data.frame(df, dfTmp, stringsAsFactors=FALSE)
            }
        rm(dfTmp)
        # Note: row names are k-mers, column names are:
        # "seqID"     "pos"       "strand"    "contig"    "contigPos"
        }
    inv(dim(df), "All genome data dim")
    inv(colnames(df), "All genome data columns")
    inv(head(df), "All genome data head")

    # Find out just what k is, i.e. how big are the k-mers?
    k = nchar(rownames(df)[1])
    inv(k, "k-mer size k")

    # findMers.cpp reports a k-mer position that is the position of the 5' end of
    # the k-mer on the strand where the k-mer is found, and it reports which strand
    # that is.  To simplify processing, we will convert all positions to be the
    # position of the 5' end of the + strand side of the k-mer, even if the k-mer
    # itself was on the - strand.
    catnow("  Canonicalizing positions\n")
    for (genome in genomeLtrs)
        {
        minusStrand = (df[,strandCol[genome]] == "-")
        df[minusStrand, posCol[genome]] = df[minusStrand, posCol[genome]] - (k-1)
        }

    # Calculate the "strand group" of each k-mer and add column "strandGroup".
    # A strand group is a particular strand direction for each of the genomes.
    # If there are N genomes, there are 2^(N-1) possible strand groups, since
    # N +/-'s have 2^N possible arrangements (e.g. +++, ++-, +-+, +--, ... ---)
    # and two arrangements that are mutual complements are strand-consistent with
    # one another (because it does not matter if some reference genome k-mers are
    # on one strand and some on the other, since each k-mer is on both strands.
    # Its canonical sequence (lower in alphabetical order) is the one whose
    # strand is given and used to form the strand groups.)  To deal with the
    # complementary pairs, complement all strands if the reference genome is "-"
    # strand, so that reference is always "+" strand.  We use the function xor()
    # (exclusive-or function) to invert the strand when ref strand is '-' = 1.
    # We form a number from the strands of each kmer using 0 for '-' strand
    # and 1 for '+' strand, with the number being represented as a binary number.
    # It is that binary number that is stored in column "strandGroup".
    catnow("  Calculating strand groups\n")
    df$strandGroup = 0
    refStrands = ifelse(df[,refStrandCol] == "+", 0, 1)
    for (genome in otherGenomes)
        df$strandGroup = 2*df$strandGroup + xor(ifelse(df[,strandCol[genome]] == "+", 0, 1), refStrands)

    # Remove k-mers so as to ensure that the difference in position of two
    # adjacent k-mers of the same LCR in ANY genome is no less than Dmin.  For
    # each genome, sort the k-mers by the strand group number, and within that
    # by the genome contig number, and within that by the genome.  After sorting,
    # search for adjacent k-mers with the same strand group number and same
    # contig number (since those must be identical for all k-mers within the
    # same LCR) that have a separation distance less than Dmin, and remove the
    # second k-mer of all such pairs, repeating until there are no more such
    # pairs.
    catnow("  Enforcing Dmin\n")
    inv(dim(df), "dim of full data before removing < Dmin")
    for (genome in genomeLtrs)
        {
        contig.Col = contigCol[genome]
        pos.Col = posCol[genome]
        df = df[order(df$strandGroup, df[,contig.Col], df[,pos.Col]),]
        while (TRUE)
            {
            N = nrow(df)
            sameGroup = (df$strandGroup[-1] == df$strandGroup[-N])
            sameContig = (df[-1, contig.Col] == df[-N, contig.Col])
            sepLessDmin = (abs(df[-1, pos.Col] - df[-N, pos.Col]) < Dmin)
            discard = (sameGroup & sameContig & sepLessDmin)
            if (!any(discard))
                break
            df = df[c(TRUE, !discard),] # Remove SECOND k-mer of the pair.
            }
        }
    inv(dim(df), "dim of full data after removing < Dmin")

    # Sort by reference genome position.
    df = df[order(df[,refPosCol]),]

    ############################################################################
    # We want to identify LCRs, which are regions having these properties:
    #   1. All k-mers in the LCR are on the same contig within each genome.
    #   2. All the N k-mers of an LCR have consistent strands across genomes (e.g.
    #       if the reference genome's kmer's strands are +++-++----+, then all
    #       other genome's kmer's are either the SAME strands or the opposite,
    #       ---+--++++-.).
    #   3. The difference in position of two adjacent k-mers in any genome of
    #       any LCR is no less than Dmin and no more than Dmax.
    #   4. All the k-mers of each genome are at positions that are either
    #       monotonically increasing or decreasing (may be different for
    #       different genomes).
    #   5. The difference in minimum and maximum position of the N k-mers in any
    #       genome of any LCR is at least Lmin.
    #   6. The LCR consists of a set of N rows of df such that, within each genome,
    #       there are N >= kmin k-mers all on the same contig.
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
    
    # Items 1 and 3 for reference genome: Create an LCR name column in df that
    # will provide a unique name for each LCR.  The name will incorporate the
    # contig number, strand group number, and a number that increments each time
    # adjacent k-mers exceed Dmax distance apart.  After forming this name and
    # splitting the k-mers on it, the resulting groups are candidate LCRs in the
    # reference genome.  The Dmin criteria has already been handled earlier.

    # Item 4: this is not a consideration for the reference genome, since it is
    # already sorted by position.  Rather, the reference genome sorting creates
    # a reference order for the candidate LCRs, which the other genomes then must
    # conform to.

    # Items 5 and 6: these are easily handled after we create the list of
    # candidate LCRs.

    ############################################################################
    # Item 3: identify k-mers that are more than Dmax distance from the preceding
    # k-mer and split LCRs at these k-mers.  Assign each LCR a unique LCR number.
    ############################################################################

    catnow("  Enforcing Dmax on reference genome\n")

    # Make column 'exceedsDmax' containing a number that increments at each
    # k-mer that exceeds Dmax from the preceding k-mer.
    df$exceedsDmax = cumsum(c(0, df[-1,refPosCol] - df[-nrow(df),refPosCol]) > Dmax)

    # Now create candidate LCR k-mer groups, which are sets of k-mers that satisfy
    # the above 6 requirements of an LCR in the reference genome, but have not yet
    # been tested in the other genomes.
    # For each row of df, assign a name that will be different for each different
    # contig, strand group, and exceedsDmax column value.
    df$LCRname = ""
    for (genome in genomeLtrs)
        {
        if (genome == refGenome)
            df$LCRname = paste("R", df[,refContigCol], sep=".")
        else
            df$LCRname = paste(df$LCRname, paste(df[,idCol[genome]], df[,contigCol[genome]], sep="."), sep="_")
        }
    df$LCRname = paste(df$LCRname, "_SG", df$strandGroup, "_D", df$exceedsDmax, sep="")
    inv(head(df), "LCRname column to divvy by contig/strand group/Dmax separation")
    inv(dim(df), "dim of full data")

    # Split df into candidate LCR groups.  Column "LCR" is added containing an
    # LCR number, starting from 1, and columns "strandGroup", "exceedsDmax", and
    # "LCRname" are removed.
    df$LCR = as.integer(factor(df$LCRname))
    df = df[, !colnames(df) %in% c("strandGroup", "exceedsDmax", "LCRname")]

    ############################################################################
    # Items 5 and 6: Discard candidate LCR k-mers that don't satisfy kMin and Lmin.
    ############################################################################

    catnow("  Enforcing kMin and Lmin on reference genome\n")

    # Save all discarded k-mers in data frame discardDf.
    discardDf = NULL
    tooSmall = tapply(df[,refPosCol], df$LCR, function(refPos) return(diff(range(refPos)) < Lmin))
    tooSmall = names(tooSmall)[tooSmall]
    inv(length(tooSmall), "number of candidate LCRs with bp span < Lmin in ref genome")
    if (length(tooSmall) > 0)
        {
        discardDf = rbind(discardDf, df[df$LCR %in% tooSmall,])
        df = df[!df$LCR %in% tooSmall,]
        }
    inv(dim(df), "dim of full data after removing < Lmin")
    inv(dim(discardDf), "dim of discarded k-mers")
    tooFewKmers = tapply(1:nrow(df), df$LCR, function(ii) length(ii) < kmin)
    tooFewKmers = names(tooFewKmers)[tooFewKmers]
    inv(length(tooFewKmers), "number of candidate LCRs with < kmin k-mers")
    if (length(tooFewKmers) > 0)
        {
        discardDf = rbind(discardDf, df[df$LCR %in% tooFewKmers,])
        df = df[!df$LCR %in% tooFewKmers,]
        }
    inv(dim(df), "dim of full data after removing < kmin")
    inv(dim(discardDf), "dim of discarded k-mers")

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
    # most one LCR.  The excluded k-mers are moved to an auxiliary LCR list and are
    # retested by themselves to see if they satisfy the requirements of an LCR.
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
    #   1. If the number of k-mers in a candidate LCR group is less than kmin,
    #       discard that candidate LCR.
    #   2. If the difference in minimum and maximum position of the candidate LCR
    #       k-mers in any genome of any LCR is less than Lmin, discard it.
    #   3. If the difference in position of two adjacent k-mers of the candidate
    #       LCR in any genome is more than Dmax, move the right k-mer, plus any
    #       immediately following k-mers that are also more than Dmax away from
    #       the left k-mer, to an auxiliary LCR to be processed subsequently
    #       (see discussion in next item).
    #   4. Identify places where monotonicity changes in a genome, and in the same
    #       way as with strand consistency, separate the k-mers into monotonicity
    #       pattern groups.  For example, say there are three genomes and one of
    #       the three contains an inversion, so its monotonicity flips.  Those
    #       k-mers in the region where it flipped are moved to a separate group
    #       and analyzed independently.  A challenge here is that monotonicity
    #       flips can be created not only by inversions but also by non-inverting
    #       translocations, which can create a single-k-mer monotonicity flip.
    #       In such a case, we want to break the candidate LCR into two at that
    #       flip.  Let's look more closely at this to understand it.
    #       For simplificity, assume two genomes.  Consider 2 cases: an inversion
    #       and a short-distance non-inverting translocation:
    #
    #         inversion:    genome 2:  . . . . F E D C B A . . . .
    #                     ref genome:  . . . . A B C D E F . . . .
    #           sign:     ref genome:   + + + + + + + + + + + + + +
    #                       genome 2:   + + + + - - - - - + + + +
    #           bigger change at:             ^           ^
    #           what we'd like at ref: LCR1--- LCR2------- LCR1---
    #
    #         non-inv xloc: genome 2:  . . . . . . . . . . A B C D E . .
    #                     ref genome:  . . A B C D E . . . . . . . . . .
    #           sign:     ref genome:   + + + + + + + + + + + + + + + +
    #                       genome 2:   + + + + + + - + + + + + + + + +
    #           bigger change at:         ^                         ^
    #           what we'd like at ref: LCR1LCR2------LCR1---------------
    #       But, we don't want to get into the ugliness of trying to analyze and
    #       recognize all the possible cases.  Instead, we want a simple algorithm
    #       for identifying k-mers that clearly don't below in the current group
    #       and eliminating them.  A different monotonicity group (here, the - signs)
    #       clearly should be separated out.  But in the non-inverting translocation,
    #       removing the one '-' will cause a new '-' to appear at the following k-mer.
    #       Also, if the first bigger-change + marked with ^ has a change bigger
    #       than Dmin and we move all following k-mers to a different LCR, the '-'
    #       will bring the following k-mers back so that a big change no longer
    #       occurs from the first LCR, and so the k-mers following the '-' should
    #       really remain part of the first LCR.  This is what makes this tough.
    #       This method works: when a discontinuity such as change of monotonicity
    #       or gap > Dmin occurs, don't split the LCR at that point, but instead,
    #       move that single disrupting k-mer out of the LCR to be handled in a
    #       new auxiliary LCR, and RECHECK for discontinuities.  By this means
    #       successive k-mers will be moved to the auxiliary LCR until either all
    #       remaining k-mers of the main LCR have been moved, or until a k-mer is
    #       reached at which there is no longer a discontinuity.
    #       Let's use an improved version of that technique: if montonicity
    #       sign changes, move the k-mer at the sign change, and all immediately
    #       following k-mers that would also introduce a sign change, to the
    #       auxiliary LCR and then RECHECK remaining k-mers for discontinuities.
    #   5. Repeat steps 1-4 over and over with each candidate LCR (including the
    #       auxiliary LCRs created in the process).  Each time an LCR is modified
    #       by moving a k-mer to the auxiliary LCR, restart at step 1 for that LCR.
    #       When all steps 1-4 are executed with no change, the LCR is good and is
    #       accepted and a sequential LCR number is assigned to it.
    ############################################################################

    ############################################################################
    # Since we need to apply the discontinuity testing over and over, define a
    # function for it.  LL is a list with these members:
    #   lcb: data frame containing k-mers being tested, all assigned to same LCR.
    #   aux: NULL or data frame of k-mers that have been removed (via previous
    #       calls to this with the same LCR number) from data frame lcb because
    #       they don't satisfy requirements, but that may form a separate LCR (or
    #       more than one).
    #   discarded: NULL or data frame of k-mers that have been rejected as members
    #       of an LCR and removed from data frame lcb.
    # On return, LL$lcb has been modified, with non-compatible k-mers removed and
    # placed into either LL$aux or LL$discarded.
    # Return value is a list with the same members and additional members:
    #   retry: TRUE if function is to be called again with same value of lcb to
    #       retest modified candidate LCR.
    #   reason: short word giving exit status:
    #       kmin: there were fewer than kmin k-mers, k-mers discarded
    #       Lmin: maximum base pair span was less than Lmin, k-mers discarded
    #       Dmax: adjacent k-mer distance was greater than Dmax, moved to aux
    #       monotonicity: k-mer monotonicity change, moved to aux
    #       okay: no problems detected
    # Note that even though we already tested the reference genome above, we must
    # retest it here because we keep reforming the candidate LCRs with subsets of
    # the previously tested ones.
    ############################################################################
    testCandidateLCR = function(LL)
        {
        LL$retry = FALSE
        Nkmers = nrow(LL$lcb)
        #catnow("Nkmers = ", Nkmers, "\n")

        #   1. If the number of k-mers in a candidate LCR group is less than kmin,
        #       discard that candidate LCR.
        if (Nkmers < kmin)
            {
            LL$discarded = rbind(LL$discarded, LL$lcb)
            LL$lcb = NULL
            LL$reason = "kmin"
            return(LL)
            }

        #   2. If the difference in minimum and maximum position of the candidate LCR
        #       k-mers in any genome of any LCR is less than Lmin, discard it.
        for (genome in genomeLtrs)
            if (diff(range(LL$lcb[,contigPosCol[genome]])) < Lmin)
                {
                LL$discarded = rbind(LL$discarded, LL$lcb)
                LL$lcb = NULL
                LL$reason = "Lmin"
                return(LL)
                }

        #   3. If the difference in position of two adjacent k-mers of the candidate
        #       LCR in any genome is more than Dmax, move the right k-mer, and any
        #       following k-mers that are also more than Dmax from the left k-mer,
        #       to the auxiliary LCR to be processed subsequently.
        for (genome in otherGenomes)
            {
            tooBig = which(abs(LL$lcb[, posCol[genome]][-1] - LL$lcb[, posCol[genome]][-Nkmers]) > Dmax)
            if (length(tooBig) > 0)
                {
                leftKmer = tooBig[1]
                rightKmer = leftKmer + 1
                # Calculate difference in position of every k-mer to the right of leftKmer
                # as if that k-mer is the next k-mer after leftKmer.
                isTooBig = (abs(LL$lcb[rightKmer:Nkmers, posCol[genome]] -
                        LL$lcb[leftKmer, posCol[genome]]) > Dmax)
                firstNotTooBig = leftKmer + match(FALSE, isTooBig)
                if (is.na(firstNotTooBig))
                    lastTooBig = Nkmers
                else
                    lastTooBig = firstNotTooBig - 1
                LL$aux = rbind(LL$aux, LL$lcb[rightKmer:lastTooBig,])
                LL$lcb = LL$lcb[-(rightKmer:lastTooBig),]
                LL$retry = TRUE
                LL$reason = "Dmax"
                return(LL)
                }
            }

        #   4. Identify places where monotonicity changes in a genome, and move
        #       the first inconsistent k-mer, and all immediately following
        #       inconsistent ones, to the auxiliary LCR to be processed subsequently.
        #       We define a monotonicity group much like the strand group above.
        #       Each genome is a binary digit position, and each binary digit is 0
        #       if genome position is increasing, 1 if decreasing (and the ref genome
        #       is always increasing since we started with it sorted, and have not
        #       changed the sort order since, so it acts as a reference).  Since
        #       monotonicity is a property of two adjacent k-mers, there are one
        #       fewer group numbers than there are k-mers.
        monotonicityGroup = rep(0, Nkmers-1)
        for (genome in otherGenomes)
            monotonicityGroup = 2*monotonicityGroup +
                ((LL$lcb[-1,posCol[genome]] - LL$lcb[-Nkmers,posCol[genome]]) < 0)
        montonicityTbl = table(monotonicityGroup)
        if (length(montonicityTbl) > 1)
            {
            mostCommonGroup = names(which.max(montonicityTbl))[1]
            leftKmer = which(monotonicityGroup != mostCommonGroup)[1]
            rightKmer = leftKmer + 1
            # Calculate monotonicity group of every k-mer to the right of leftKmer
            # as if that k-mer is the next k-mer after leftKmer.
            mg = rep(0, Nkmers-leftKmer)
            for (genome in otherGenomes)
                mg = 2*mg +
                    ((LL$lcb[(rightKmer:Nkmers),posCol[genome]] - LL$lcb[leftKmer,posCol[genome]]) < 0)
            # What is the first one that is equal to mostCommonGroup, if any?
            firstRight = leftKmer + match(mostCommonGroup, mg)
            if (is.na(firstRight))
                lastWrong = Nkmers
            else
                lastWrong = firstRight - 1
            LL$aux = rbind(LL$aux, LL$lcb[(rightKmer:lastWrong),])
            LL$lcb = LL$lcb[-(rightKmer:lastWrong),]
            LL$retry = TRUE
            LL$reason = "monotonicity"
            return(LL)
            }
        LL$reason = "okay"
        return(LL)
        }

    ############################################################################
    # Now, for each candidate LCR in df, test it repeatedly, removing incompatible
    # k-mers from it over to a holding area for creation of an auxiliary LCR,
    # repeating until no more k-mers are removed and the LCR is known to satisfy
    # the constraints (kMin, Lmin, Dmax, monotonicity), then incorporate the
    # results into new result data frame lcbDf, and if some k-mers were found
    # incompatible with the LCR and moved to an aux data frame, rbind that aux
    # data frame back to df.  Thus, df will shrink as candidate LCRs are either
    # verified and moved to lcbDf or discarded and moved to discardDf.
    # This loop is slow, and in order to try to speed it up, I'm using some smaller
    # temporary data frames to store data.  These are:
    #   dft: k-mers of candidate LCRs from df, reloaded with k-mers for the next
    #           N.reload candidate LCRs whenever it becomes empty.
    #   lcbDft: k-mers of newly confirmed LCRs, appended to lcbDf every reload.
    #   discardDft: k-mers incompatible with LCRs, appended to discardDf every reload.
    #   auxDft: k-mers newly added to the candidate k-mers with a new LCR number,
    #       appended to df every reload.
    # Keep track of stats in stats vector.
    ############################################################################

    catnow("  Testing candidate LCRs: discarding LCRs when kMin or Lmin exceeded in any genome,\n")
    catnow("  making new candidate LCRs using k-mers following one where Amax is exceeded or all\n")
    catnow("  k-mers that violate monotonicity, and accepting LCRs when all k-Mers satisfy the\n")
    catnow("  kMin, Lmin, Dmax, and monotonicity constraints.\n")

    N.reload = 100 # Number of LCR numbers to load into dft each time.
    logEveryN = 500
    lastLog = 0
    dft = df[integer(0),]
    lcbDft = NULL
    discardDft = NULL
    auxDft = NULL
    #
    idx = 1
    stats = integer()
    LCRs = unique(df$LCR)
    maxLCR = max(LCRs)
    lcbDf = NULL
    while (idx <= length(LCRs))
        {
        #catnow("idx = ", idx, "length(LCRs) = ", length(LCRs), "\n")
        # Reload temporary data frames if dft is empty.
        if (nrow(dft) == 0)
            {
            lenLCRs = length(LCRs)
            if (lenLCRs - lastLog >= logEveryN)
                {
                catnow("    At candidate LCR #", idx, " Total # of candidate LCRs to be tested: ", lenLCRs, "\n")
                lastLog = lenLCRs
                }
            df = rbind(df, auxDft)
            auxDft = NULL
            discardDf = rbind(discardDf, discardDft)
            discardDft = NULL
            lcbDf = rbind(lcbDf, lcbDft)
            lcbDft = NULL
            lastIdx = idx+N.reload-1
            if (lastIdx > length(LCRs))
                lastIdx = length(LCRs)
            xferKmers = (df$LCR %in% LCRs[idx:lastIdx])
            dft = df[xferKmers,]
            df = df[!xferKmers,]
            }
        # Get candidate LCR k-mers into list LL and remove them from df.
        #catnow("Get LCR\n")
        lcbKmers = (dft$LCR == LCRs[idx])
        LL = list(lcb=dft[lcbKmers,], aux=NULL, discarded=NULL)
        dft = dft[!lcbKmers,]
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
        # If the candidate k-mer list is not NULL, append it to lcbDft.
        if (!is.null(LL$lcb))
            {
            #catnow(" lcbDft\n")
            lcbDft = rbind(lcbDft, LL$lcb)
            }
        # If some k-mers appear in aux, append them to auxDft using a new LCR number.
        if (!is.null(LL$aux))
            {
            #catnow(" auxDft\n")
            maxLCR = maxLCR+1
            LCRs = c(LCRs, maxLCR)
            LL$aux$LCR = maxLCR
            auxDft = rbind(auxDft, LL$aux)
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

    # Append final results from temporary data frames to permanent ones.
    if (nrow(dft) != 0) stop("Expected empty dft")
    if (nrow(df) != 0) stop("Expected empty df")
    if (!is.null(auxDft)) stop("Expected null auxDft")
    discardDf = rbind(discardDf, discardDft)
    lcbDf = rbind(lcbDf, lcbDft)
    rm(dft, lcbDft, discardDft)

    catnow("\nFor reference ID", ID, "found", length(unique(lcbDf$LCR)), "LCRs containing", nrow(lcbDf), "k-mers\n")

    # Accumulate results.
    cumLcbDf = rbind(cumLcbDf, lcbDf)
    cumDiscardDf = rbind(cumDiscardDf, discardDf)
    rm(df, lcbDf, discardDf)
    }

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
