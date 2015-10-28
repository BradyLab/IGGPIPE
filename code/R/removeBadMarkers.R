################################################################################
# See usage below for description.
# Author: Ted Toal
# Date: 2013-2015
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
#thisDir = "~/Documents/UCDavis/BradyLab/Genomes/IGGPIPE/code/R/" # For testing only.

# Source the necessary include files from the same directory containing this file.
source(paste(thisDir, "Include_Common.R", sep=""))
source(paste(thisDir, "Include_RemoveOverlapping.R", sep=""))

# Get arguments.
testing = 0
#testing = 1 # For testing only. outTestHP11
#testing = 2 # For testing only. outHPT14_400_1500_50_300
{
if (testing == 0)
    args = commandArgs(TRUE)
else if (testing == 1)
    {
    args = c("~/Documents/UCDavis/BradyLab/Genomes/IGGPIPE", "MIN",
        "outTestHP11/NonvalidatedMarkers_K11k2L100D10_2000A100_2000d10_100N2F0X20.tsv",
        "outTestHP11/MarkerErrors_K11k2L100D10_2000A100_2000d10_100N2F0X20V3000W8M3G1",
        "outTestHP11/MarkersOverlapping_K11k2L100D10_2000A100_2000d10_100N2F0X20V3000W8M3G1.tsv",
        "outTestHP11/MarkersNonoverlapping_K11k2L100D10_2000A100_2000d10_100N2F0X20V3000W8M3G1.tsv",
        TRUE)
    }
else if (testing == 2)
    {
    args = c("~/Documents/UCDavis/BradyLab/Genomes/IGGPIPE", "MIN",
        "outHPT14_400_1500_50_300/NonvalidatedMarkers_K14k4L400D1_1500A400_1500d50_300N2F0X10.tsv",
        "outHPT14_400_1500_50_300/MarkerErrors_K14k4L400D1_1500A400_1500d50_300N2F0X10V2500W8M0G0",
        "outHPT14_400_1500_50_300/MarkersOverlapping_K14k4L400D1_1500A400_1500d50_300N2F0X10V2500W8M0G0.tsv",
        "outHPT14_400_1500_50_300/MarkersNonoverlapping_K14k4L400D1_1500A400_1500d50_300N2F0X10V2500W8M0G0.tsv",
        TRUE)
    }
else stop("Unknown value for 'testing'")
}

NexpectedMin = 7
if (length(args) < NexpectedMin)
    {
    usage = c(
        "Read a data frame of candidate IGG marker primer pairs and, for each genome, another data",
        "frame of bad primer pairs, then remove the bad pairs from the candidate IGG marker file",
        "and write the cleaned-up IGG marker data to a new file.",
        "",
        "Usage: Rscript removeBadMarkers.R <wd> <minMax> <tsvMarkerFile> <badMarkerPfx> \\",
        "       <overlappingFile> <nonoverlappingFile> <investigate>",
        "",
        "Arguments:",
        "   <wd>                 : Path of R working directory, specify other file paths relative to this.",
        "   <minMax>             : Either MIN or MAX to indicate the method to be used to remove overlapping",
        "                          Indel Groups.  MIN means that the smallest Indel Groups are retained over",
        "                          larger ones, while MAX means the larger are retained.",
        "   <tsvMarkerFile>      : Input file containing the candidate markers with primer pairs.",
        "   <badMarkerPfx>       : Prefix, including directory, of input files containing bad markers to",
        "                           be removed.  The filename suffix is '_<genome letter>.bad.tsv'.",
        "   <overlappingFile>    : Output file to be written containing overlapping markers with",
        "                          e-PCR-checked primers.",
        "   <nonoverlappingFile> : Output file to be written containing non-overlapping markers with",
        "                          e-PCR-checked primers.",
        "   <investigate>        : FALSE for normal operation, TRUE for more verbose debugging output."
        )
    for (S in usage)
        catnow(S, "\n", sep="")
    stop("Try again with correct number of arguments")
    }

catnow("removeBadMarkers.R arguments:\n")
workingDirectory = args[1]
catnow("  workingDirectory: ", workingDirectory, "\n")
if (!dir.exists(workingDirectory))
    stop("Directory doesn't exist: ", workingDirectory)
setwd(workingDirectory)

minMax = args[2]
catnow("  minMax: ", minMax, "\n")
if (!minMax %in% c("MIN", "MAX"))
    stop("minMax must be MIN or MAX")

tsvMarkerFile = args[3]
catnow("  tsvMarkerFile: ", tsvMarkerFile, "\n")
if (!file.exists(tsvMarkerFile))
    stop("File doesn't exist: ", tsvMarkerFile)

badMarkerPfx = args[4]
catnow("  badMarkerPfx: ", badMarkerPfx, "\n")

overlappingFile = args[5]
catnow("  overlappingFile: ", overlappingFile, "\n")

nonoverlappingFile = args[6]
catnow("  nonoverlappingFile: ", nonoverlappingFile, "\n")

investigate = as.logical(args[7])
catnow("  investigate: ", investigate, "\n")
if (is.na(investigate))
    stop("investigate must be TRUE or FALSE")

########################################
# Initialize.
########################################

# Read candidate marker data including primers.
dfMarkers = read.table(tsvMarkerFile, header=TRUE, row.names=NULL, sep="\t", stringsAsFactors=FALSE)
catnow("Number of candidate markers read from input file:", nrow(dfMarkers), "\n")
if (nrow(dfMarkers) == 0)
    stop("There are no candidate markers.")
inv(dim(dfMarkers), "input data dim")
inv(colnames(dfMarkers), "input data columns")
inv(head(dfMarkers), "input data head")
NoriginalMarkers = nrow(dfMarkers)
rownames(dfMarkers) = paste(dfMarkers$kmer1, dfMarkers$kmer2, sep="_")

# The k-mer length we are working with.
kmerLen = nchar(dfMarkers$kmer1[1])
inv(kmerLen, "k")

# Get the genome letters.
genomeLtrs = substr(colnames(dfMarkers)[grepl("^.id$", colnames(dfMarkers))], 1, 1)
if (length(genomeLtrs) < 2)
    stop("Expected at least two id columns in <tsvMarkerFile>")
inv(genomeLtrs, "genome letter")
Ngenomes = length(genomeLtrs)
inv(Ngenomes, "Ngenomes")

# Create vectors of column names in dfMarkers, indexed by genome, for the id,
# ampPos1, ampPos2, and phase columns.
makeColVec = function(S)
    {
    V = paste(genomeLtrs, S, sep="")
    names(V) = genomeLtrs
    return(V)
    }
idCol = makeColVec("id")
ampPos1Col = makeColVec("ampPos1")
ampPos2Col = makeColVec("ampPos2")

########################################
# Read bad marker files, remove the markers from dfMarkers, and accumulate stats.
########################################

reasonCounts = list() # Indexed by reasonDiscarded column value.
for (i in 1:Ngenomes)
    {
    badMarkerFile = paste(badMarkerPfx, "_", i, ".bad.tsv", sep="")
    df = read.table(badMarkerFile, header=TRUE, row.names=NULL, sep="\t", stringsAsFactors=FALSE)
    reasons = unique(df$reasonDiscarded)
    reasonCounts[setdiff(reasons, names(reasonCounts))] = 0
    counts = table(df$reasonDiscarded)
    for (reason in names(counts))
        reasonCounts[[reason]] = reasonCounts[[reason]] + counts[reason]
    dfMarkers = dfMarkers[!rownames(dfMarkers) %in% paste(df$kmer1, df$kmer2, sep="_"),]
    }
NgoodMarkers = nrow(dfMarkers)
catnow("\n")
catnow("Totals of reasons for removal (some markers have multiple reasons):\n")
catnow("    ", format("Reason", width=20), "   Count\n")
catnow("    ", format("------", width=20), "   -----\n")
for (reason in names(reasonCounts))
    catnow("    ", format(reason, width=20), "   ", reasonCounts[[reason]], "\n")
catnow("\n")
catnow("Number of markers to start with:", NoriginalMarkers, "\n")
catnow("Good markers after removing bad:", NgoodMarkers, "\n")
catnow("\n")
if (nrow(dfMarkers) == 0)
    stop("There are no markers left after removing bad markers!!")
catnow("Finished removing e-PCR-identified bad markers from file\n")

# Put the data frame in order by reference genome position.
catnow("Sorting by reference genome position...")
dfMarkers = dfMarkers[order(dfMarkers[, idCol[1]], dfMarkers[, ampPos1Col[1]]),]
rownames(dfMarkers) = NULL
catnow("\n")

########################################
# Write the overlapping markers to a file.
########################################

rownames(dfMarkers) = NULL
write.table.winSafe(dfMarkers, overlappingFile, col.names=TRUE, row.names=FALSE, quote=FALSE, sep="\t")
# dfMarkers = read.table(overlappingFile, header=TRUE, sep="\t", stringsAsFactors=FALSE)
catnow(nrow(dfMarkers), "overlapping markers output to file:\n", overlappingFile, "\n")

########################################
# Remove overlapping markers, guided by the value of minMax, which is either
# MIN or MAX.
########################################

inv(nrow(dfMarkers), "Number of markers including overlapping markers")

# We will use the same data frame that holds the overlapping markers, dfMarkers,
# to hold the non-overlapping markers, since there might be a LOT of them (so
# copying the data frame would be costly) and we need to start out this algorithm
# with the data frame containing the overlapping markers, which dfMarkers already
# does.  The data frame will be transformed into one containing only non-overlapping
# markers.  We no longer need the overlapping markers (except for this).

dfMarkers = removeOverlappingRows(dfMarkers, "Marker", "genome", genomeLtrs, idCol, ampPos1Col, ampPos2Col, verbose=TRUE)

########################################
# Write the non-overlapping markers to a file.
########################################

rownames(dfMarkers) = NULL
write.table.winSafe(dfMarkers, nonoverlappingFile, col.names=TRUE, row.names=FALSE, quote=FALSE, sep="\t")
# dfMarkers = read.table(nonoverlappingFile, header=TRUE, sep="\t", stringsAsFactors=FALSE)
catnow(nrow(dfMarkers), "non-overlapping markers output to file:\n", nonoverlappingFile, "\n")

}

################################################################################
# End of file.
################################################################################
