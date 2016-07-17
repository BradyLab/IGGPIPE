################################################################################
# See usage below for description.
# Author: Ted Toal
# Date: 2015
# Brady Lab, UC Davis
################################################################################

# Enclose everything in braces so stop statements will work correctly.
{

# Pathname separator.
#PATHSEP = ifelse(grepl("/", Sys.getenv("HOME")), "/", "\\")
PATHSEP = "/"

# Get directory where this file resides.
XSEP = ifelse(PATHSEP == "\\", "\\\\", PATHSEP)
RE = paste("^.*--file=(([^", XSEP, "]*", XSEP, ")*)[^", XSEP, "]+$", sep="")
args = commandArgs(FALSE)
thisDir = sub(RE, "\\1", args[grepl("--file=", args)])
#thisDir = "~/Documents/UCDavis/BradyLab/IGGPIPE/IGGPIPE/code/R/" # For testing only.

# Source the necessary include files from the same directory containing this file.
source(paste(thisDir, "Include_Common.R", sep=""))

# Get arguments.
testing = 0
#testing = 1 # For testing only, outTestHP11 LCRs
#testing = 2 # For testing only, outTaCW15 LCRs
{
if (testing == 0)
    {
    args = commandArgs(TRUE)
    investigate = FALSE
    }
else if (testing == 1)
    {
    args = c("~/Documents/UCDavis/BradyLab/IGGPIPE/IGGPIPE",
    "outTestHP11/LCRs_K11k2L100D10_2000.tsv", "outTestHP11/LCBs_K11k2L100D10_2000.tsv")
    investigate = TRUE
    }
else if (testing == 2)
    {
    args = c("~/Documents/UCDavis/BradyLab/IGGPIPE/IGGPIPE",
    "outTaCW15/LCRs_K15k2L200D5_1000.tsv", "outTaCW15/LCBs_K15k2L200D5_1000.tsv")
    investigate = TRUE
    }
else stop("Unknown value for 'testing'")
}

Nexpected = 3
if (length(args) != Nexpected)
    {
    usage = c(
        "Read a data frame of LCRs and remove overlapping regions by breaking the LCRs into",
        "multiple pieces which we will call 'locally conserved blocks'.  When the common",
        "unique k-mers are sorted in order by any genome, the consecutive k-mers that are",
        "part of the same LCR become an LCB, and when a different LCR interrupts the first",
        "LCR (i.e. is located in the middle of it; overlaps it), a new LCB is started.  Write",
        "the LCBs to the output file, containing the positions of the LCB k-mers in each",
        "genome, along with their relative phases (+ or - strand relative to reference genome)",
        "and an LCB ID number.  The contig information that is part of the input file is",
        "dropped.",
        "",
        "Usage: Rscript LCRsToNonoverlappingLCBs.R <wd> <inputFile> <outputFile>",
        "",
        "Arguments:",
        "   <wd> : Path of R working directory, specify other file paths relative to this.",
        "   <inputFile>   : Input file containing LCRs.",
        "   <outputFile>  : Output file name to which to write LCBs.",
        )
    for (S in usage)
        catnow(S, "\n", sep="")
    stop("Try again with correct number of arguments")
    }

catnow("LCRsToNonoverlappingLCBs.R arguments:\n")
workingDirectory = args[1]
catnow("  workingDirectory: ", workingDirectory, "\n")
if (!dir.exists(workingDirectory))
    stop("Directory doesn't exist: ", workingDirectory)
setwd(workingDirectory)

inputFile = args[2]
catnow("  inputFile: ", inputFile, "\n")
if (!file.exists(inputFile))
    stop("File doesn't exist: ", inputFile)

outputFile = args[3]
catnow("  outputFile: ", outputFile, "\n")

########################################
# Initialize.
########################################

catnow("Removing overlaps from input LCRs file data\n")

# Read input data.
catnow("Reading input data.\n")
df = read.table(inputFile, header=TRUE, row.names=1, sep="\t", stringsAsFactors=FALSE)
if (nrow(df) == 0)
    stop("There is no input LCR data.")
inv(dim(df), "input data dim")
inv(colnames(df), "input data columns")
inv(head(df), "input data head")

# Convert df into a data frame with columns 'LCB', 'X.seqID', 'X.pos', and
# 'X.strand', where X=genome letter.  When two LCRs overlap in any genome after
# sorting in that genome's order, the k-mer at the lower address at the overlap
# junction is marked, and then the data is sorted by LCRs and within them by
# that genome position, and then the LCRs are split at each marked k-mer by
# adding a number to the end of the LCR number.  The split-up LCRs are called
# LCBs.  At the end, any LCB with just one k-mer is discarded.

# Input:
#                   H.seqID	H.pos	H.strand	H.contig	H.contigPos	P.seqID	P.pos	P.strand	P.contig	P.contigPos	LCR
# GATTGCGTGCATCC	SL2.40ch00	21658	+	1	21671	Spenn-ch02	3548436	-	143	14555	81
# AGGCTGATAATTGC	SL2.40ch00	21729	+	1	21742	Spenn-ch02	3548286	-	143	14405	81
# CGACTGTATTGAAC	SL2.40ch00	22014	+	1	22027	Spenn-ch02	3548021	-	143	14140	81

# Output:
#                   LCB	phases  Hid	Hpos	Pid	Ppos
# GATTGCGTGCATCC	81_1_2    "+-"    SL2.40ch00	21658	Spenn-ch02	3548436

# The k-mer length we are working with.
catnow("Computing k-mer length\n")
kmerLen = nchar(rownames(df)[1])
inv(kmerLen, "k")
if (is.na(kmerLen) || kmerLen < 1 || kmerLen > 50)
    stop("Bad LCR file, can't find k-mer length from row names.\n")

# Get genome letters.
catnow("Getting genome letters\n")
genomeLtrs = sub("\\.seqID$", "", colnames(df)[grepl("\\.seqID$", colnames(df))])
if (length(genomeLtrs) < 2)
    stop("Expected at least two id columns in <inputFile>")
inv(genomeLtrs, "genome letters")
otherGenomeLtrs = genomeLtrs[-1]
catnow("Genome letters: ", paste(genomeLtrs, collapse=","), "\n")

# Drop "contig" columns.
catnow("Dropping contig columns\n")
df = df[, !grepl("contig", colnames(df))]

catnow("Renaming columns\n")
# Remove "." from column names.
colnames(df) = sub("\\.", "", colnames(df))
# Change XseqID to Xid.
colnames(df) = sub("seqID", "id", colnames(df))
# Change LCR to LCB.
colnames(df) = sub("LCR", "LCB", colnames(df))

# Get vectors of names of X columns.
catnow("Creating column name vectors\n")
idCols = colnames(df)[grepl("id$", colnames(df))]
posCols = colnames(df)[grepl("pos$", colnames(df))]
strandCols = colnames(df)[grepl("strand$", colnames(df))]
names(idCols) = genomeLtrs
names(posCols) = genomeLtrs
names(strandCols) = genomeLtrs

# Add column "phases".
catnow("Adding column 'phases'\n")
df$phases = ""
for (genome in genomeLtrs)
    df$phases = paste(df$phases, c("+", "-")[1+(df[,strandCols[genome]] != df[,strandCols[1]])], sep="")

# Remove strand columns, don't need them any more.
catnow("Removing strand columns\n")
df = df[, !colnames(df) %in% strandCols]

# For each genome, sort in that genome's order, then mark the k-mer where the
# LCB number changes.  Following that, re-sort in LCB order and within the LCB
# in that genome's order, then append a sequence number to the LCB, starting at
# 1 and incrementing at each position where a k-mer is marked to require an LCB
# number change.
for (ltr in genomeLtrs)
    {
    catnow("Removing overlaps from genome ", ltr, "\n")

    # Sort in genome order.
    catnow("  Sorting in genome order\n")
    df = df[order(df[, idCols[ltr]], df[, posCols[ltr]]),]
    N = nrow(df)

    # Get flags for next k-mer having a different LCB number.
    catnow("  Getting flags for LCB change\n")
    LCBpairDifferent = (df$LCB[-1] != df$LCB[-N])

    # Mark k-mers where LCB number changes, marking the SECOND k-mer of the pair
    # at which the number changes.  Do this without regard for phase of the LCB.
    # The meaning of a mark is that the marked k-mer must be the first k-mer of
    # an LCB, when the k-mers are sorted in LCB order and within that in genome
    # order.
    catnow("  Marking LCB changes\n")
    df$breakLCB = c(TRUE, LCBpairDifferent) # Mark second one of pair.  Very first k-mer is considered different LCB.

    # Re-sort in LCB order, and within that in genome order.
    catnow("  Resorting in LCB/genome order\n")
    df = df[order(df$LCB, df[, posCols[ltr]]),]

    # Get flags for next k-mer having a different LCB sequence number.
    catnow("  Getting flags for LCB change\n")
    LCBpairDifferent = (df$LCB[-1] != df$LCB[-N])

    # Form the cumulative sum of the marked k-mers.
    catnow("  Cumulative summing LCB change flags\n")
    cumMarks = cumsum(df$breakLCB)

    # Get number of k-mers for each LCB.
    catnow("  Getting number of k-mers per LCB\n")
    firstKmerOfLCB = c(TRUE, LCBpairDifferent)
    firstKmerIdx = which(c(firstKmerOfLCB, TRUE))
    numKmersPerLCB = firstKmerIdx[-1] - firstKmerIdx[-length(firstKmerIdx)]
    numLCBs = length(numKmersPerLCB)

    # Get vector of offsets to be subtracted from cumMarks to get new LCB suffix.
    catnow("  Getting LCB suffix offsets\n")
    subtractFromOffset = cumMarks[firstKmerOfLCB]-1
    if (length(subtractFromOffset) != numLCBs) stop("Number of LCBs computed wrong")
    subtractFromOffset = unlist(sapply(1:numLCBs, function(i)
        rep(subtractFromOffset[i], numKmersPerLCB[i])), use.names=FALSE)
    if (length(subtractFromOffset) != nrow(df)) stop("Incorrect offset vector")

    # Determine what to append to each LCB number.  For the first LCB, the number
    # to append is cumMarks.  For the second one, it is cumMarks+1-cumMarks[first
    # k-mer of second one], etc.  Then append a suffix to the LCB number with the
    # offset as the suffix value.
    catnow("  Suffixing LCB number\n")
    LCBoffset = cumMarks - subtractFromOffset
    df$LCB = paste(df$LCB, LCBoffset, sep="_")
    }

# Remove df rows of LCBs containing just one k-mer.
catnow("  Removing stand-alone k-mer rows\n")
NkmersPerLCB = table(df$LCB)
keepLCBs = names(NkmersPerLCB)[NkmersPerLCB > 1]
df = df[df$LCB %in% keepLCBs,]
N = nrow(df)
catnow("  ", sum(NkmersPerLCB <= 1), " stand-alone k-mer rows removed\n")

# Do a final double-check.  There should not be any overlapping LCBs.
for (ltr in genomeLtrs)
    {
    catnow("Double-checking genome ", ltr, "\n")

    # Sort in genome order.
    catnow("  Sorting in genome ", ltr, " order\n")
    df = df[order(df[, idCols[ltr]], df[, posCols[ltr]]),]
    N = nrow(df)
    catnow("  Computing breakLCB flags\n")
    df$breakLCB = c(df$LCB[-1] != df$LCB[-N], TRUE)

    # Re-sort in LCB order, and within that in genome order.
    catnow("  Sorting in LCB/genome order\n")
    df = df[order(df$LCB, df[, posCols[ltr]]),]

    # Get flags for next k-mer having a different LCB number.
    catnow("  Computing next LCB different flags\n")
    LCBpairDifferent = c(df$LCB[-1] != df$LCB[-N], TRUE)

    # Every case where breakLCB is TRUE should also have LCBpairDifferent.
    catnow("  Comparing actual to expected flags.\n")
    if (any(df$breakLCB & !LCBpairDifferent))
        stop("Genome ", ltr, " error")
    }

# Put the LCB and phases column first and remove LCB and breakLCB columns.
catnow("Selecting and re-ordering columns.\n")
cols = c("LCB", "phases")
for (ltr in genomeLtrs)
    cols = c(cols, idCols[ltr], posCols[ltr])
df = df[, cols]

catnow("  Sorting in reference genome order\n")
df = df[order(df[, idCols[1]], df[, posCols[1]]),]

########################################
# Save the output file.
########################################

write.table.winSafe(df, outputFile, row.names=FALSE, quote=FALSE, sep="\t")
# df = read.table(outputFile, header=TRUE, row.names=NULL, sep="\t", stringsAsFactors=FALSE)
catnow("Number of non-overlapping LCBs:", nrow(df), "\n")
catnow("Finished converting LCRs to LCBs, output file:\n", outputFile, "\n")
}

################################################################################
# End of file.
################################################################################
