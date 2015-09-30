################################################################################
# See usage below for description.
# Author: Ted Toal
# Date: 2013-2015
# Brady Lab, UC Davis
################################################################################

# Enclose everything in braces so stop statements will work correctly.
{

# Pathname separator.
PATHSEP = ifelse(grepl("/", Sys.getenv("HOME")), "/", "\\")

# Get directory where this file resides.
XSEP = ifelse(PATHSEP == "\\", "\\\\", PATHSEP)
RE = paste("^.*--file=(([^", XSEP, "]*", XSEP, ")*)[^", XSEP, "]+$", sep="")
args = commandArgs(FALSE)
thisDir = sub(RE, "\\1", args[grepl("--file=", args)])
#thisDir = "~/Documents/UCDavis/BradyLab/Genomes/kmers/IGGPIPE/code/R/" # For testing only.

# Source the necessary include files from the same directory containing this file.
source(paste(thisDir, "Include_Common.R", sep=""))

# Get arguments.
testing = 0
#testing = 1 # For testing only.
{
if (testing == 0)
    args = commandArgs(TRUE)
else if (testing == 1)
    {
    args = c("~/Documents/UCDavis/BradyLab/Genomes/kmers/IGGPIPE", "MIN",
        "outTestHP11/NonvalidatedMarkers_K11k2L100D10_2000A100_2000d10_100N2F0X20.tsv",
        "outTestHP11/MarkerErrors_K11k2L100D10_2000A100_2000d10_100N2F0X20V3000W8M3G1_",
        "outTestHP11/MarkersOverlapping_K11k2L100D10_2000A100_2000d10_100N2F0X20V3000W8M3G1.tsv",
        "outTestHP11/MarkersNonoverlapping_K11k2L100D10_2000A100_2000d10_100N2F0X20V3000W8M3G1.tsv",
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
dfMarkers = dfMarkers[dfMarkers(dfMarkers[, idCol[refGenome]], dfMarkers[, pos1Col[refGenome]]),]
rownames(dfMarkers) = NULL
catnow("\n")

########################################
# Write the overlapping markers to a file.
########################################

rownames(dfMarkers) = NULL
write.table(dfMarkers, overlappingFile, col.names=TRUE, row.names=FALSE, quote=FALSE, sep="\t")
# dfMarkers = read.table(overlappingFile, header=TRUE, sep="\t", stringsAsFactors=FALSE)
catnow(nrow(dfMarkers), "overlapping markers output to file:\n", overlappingFile, "\n")

########################################
# Remove overlapping markers, guided by the value of minMax, which is either
# MIN or MAX.  This is a copy of the same code from findIndelGroups.R.
########################################

inv(nrow(dfMarkers), "Number of markers including overlapping markers")

# We will use the same data frame that holds the overlapping markers, dfMarkers,
# to hold the non-overlapping markers, since there might be a LOT of them (so
# copying the data frame would be costly) and we need to start out this algorithm
# with the data frame containing the overlapping markers, which dfMarkers already
# does.  The data frame will be transformed into one containing only non-overlapping
# markers.  We no longer need the overlapping markers (except for this).

# Test each genome, one by one, for marker overlaps, and remove markers to get
# rid of them.
for (genome in genomeLtrs)
    {
    inv(genome, "Remove overlaps")
    id.Col = idCol[genome]
    pos1.Col = ampPos1Col[genome]
    pos2.Col = ampPos2Col[genome]

    # Refer to more comments in findIndelGroups.R for more information.

    # Loop until no more markers are found to overlap in this genome.
    while (TRUE)
        {
        inv(nrow(dfMarkers), "Loop with # markers remaining")

        # Set new column "start" to the position of the 5' end of the 5' k-mer
        # of the amplicon, and new column "end" to the position of the 5' end of
        # the 3' k-mer of the amplicon.  There are three issues to be dealt with:
        #   1. The amplicon end position for "start" may be either ampPos1 or
        #       ampPos2 depending on the strand polarity; "start" is always the
        #       smaller of the two, "end" the larger.
        #   2. We want "start" and "end" to be k-mer ends, not amplicon ends.
        #       Columns kmer1offset and kmer2offset give us the offset to add,
        #       but the sign is complicated.  These offset columns give offset
        #       positive towards center of amplicon, negative away from center.
        #       We need to negate the "end" offset depending on strand polarity.
        #   3. kmer1offset and kmer2offset give offset to outside edge of the
        #       k-mer, but we want offset to 5' side of k-mer so we may need an
        #       additional "end" offset of kmerLen-1 depending on strand polarity.
        dfMarkers$start = dfMarkers[, pos1.Col]
        dfMarkers$end = dfMarkers[, pos2.Col]
        pos2IsSmaller = (dfMarkers[,pos1.Col] > dfMarkers[,pos2.Col])
        dfMarkers$start[pos2IsSmaller] = dfMarkers[pos2IsSmaller, pos2.Col]
        dfMarkers$end[pos2IsSmaller] = dfMarkers[pos2IsSmaller, pos1.Col]
        startOffset = dfMarkers$kmer1offset
        endOffset = -dfMarkers$kmer2offset - kmerLen + 1
        startOffset[pos2IsSmaller] = dfMarkers$kmer2offset[pos2IsSmaller]
        endOffset[pos2IsSmaller] = -dfMarkers$kmer1offset[pos2IsSmaller] - kmerLen + 1
        dfMarkers$start = dfMarkers$start + startOffset
        dfMarkers$end = dfMarkers$end + endOffset

        # Add new column "len" equal to length of the marker segments.
        dfMarkers$len = dfMarkers$end - dfMarkers$start + 1

        # Sort by ID and start position.
        dfMarkers = dfMarkers[order(dfMarkers[, id.Col], dfMarkers$start),]
        N = nrow(dfMarkers)

        # Find index of marker that each marker overlaps through and put it in column thruIdx.
        dfMarkers$thruIdx = NA
        for (id in unique(dfMarkers[, id.Col]))
            {
            thisId = (dfMarkers[, id.Col] == id)
            dfMarkers$thruIdx[thisId] = match(TRUE, thisId) - 1 + findInterval(dfMarkers$end[thisId], dfMarkers$start[thisId]+1)
            }

        # Set thruIdx of markers that do not overlap even the next marker to 0.
        dfMarkers$thruIdx[dfMarkers$thruIdx == 1:nrow(dfMarkers)] = 0

        # Get the set of indexes of all markers which overlap at least one other marker.
        overlapIdxs = sapply(1:nrow(dfMarkers), function(i)
            {
            if (dfMarkers$thruIdx[i] == 0)
                return(0)
            return(i:dfMarkers$thruIdx[i])
            })
        overlapIdxs = unique(unlist(overlapIdxs))
        # Remove index 0, which comes from non-overlapping markers.
        overlapIdxs = overlapIdxs[overlapIdxs != 0]

        # If no markers overlap, break out of loop.
        inv(length(overlapIdxs), "Number of overlapping markers")
        if (length(overlapIdxs) == 0)
            break

        # Set column "overlap" TRUE for each of those overlapping markers.
        dfMarkers$overlap = FALSE
        dfMarkers$overlap[overlapIdxs] = TRUE

        # Get the start and end index of each group of mutually overlapping markers.
        startOverlap = which(!c(FALSE, dfMarkers$overlap[-N]) & dfMarkers$overlap)
        endOverlap = which(dfMarkers$overlap & !c(dfMarkers$overlap[-1], FALSE))
        if (length(startOverlap) != length(endOverlap)) stop("Expected equal start/end overlap vectors")

        # Find the index within each group of the marker with the smallest or largest length.
        # Then get the indexes within the group of the markers that overlap that marker with
        # the smallest or largest length.
        idxsToRemove = sapply(1:length(startOverlap), function(i)
            {
            idxs = startOverlap[i]:endOverlap[i]
            len.SL = ifelse(minMax == "MIN", min(dfMarkers$len[idxs]), max(dfMarkers$len[idxs]))
            idx.SL = idxs[dfMarkers$len[idxs] == len.SL][1] # If more than one, pick the first.
            return(idxs[!(dfMarkers$end[idxs] <= dfMarkers$start[idx.SL]) & !(dfMarkers$start[idxs] >= dfMarkers$end[idx.SL]) & idxs != idx.SL])
            })
        idxsToRemove = unlist(idxsToRemove)

        # Remove those markers.
        dfMarkers = dfMarkers[-idxsToRemove,]
        inv(length(idxsToRemove), "Number of markers deleted")
        inv(nrow(dfMarkers), "Number of markers remaining")
        }
    }
dfNoOverlaps = dfNoOverlaps[, !colnames(dfNoOverlaps) %in% c("start", "end", "len", "thruIdx", "overlap")]
inv(nrow(dfNoOverlaps), "Number of markers with overlapping markers removed")

# Put the data frame in order by reference genome position.
catnow("Sorting by reference genome position...")
dfNoOverlaps = dfNoOverlaps[dfNoOverlaps(dfNoOverlaps[, idCol[refGenome]], dfNoOverlaps[, pos1Col[refGenome]]),]
rownames(dfNoOverlaps) = NULL
 if (nrow(dfNoOverlaps) == 0)
    stop("There are no markers left after removing overlapping markers!!")
catnow("\n")

########################################
# Write the non-overlapping markers to a file.
########################################

rownames(dfNoOverlaps) = NULL
write.table(dfNoOverlaps, nonoverlappingFile, col.names=TRUE, row.names=FALSE, quote=FALSE, sep="\t")
# dfNoOverlaps = read.table(nonoverlappingFile, header=TRUE, sep="\t", stringsAsFactors=FALSE)
catnow(nrow(dfNoOverlaps), "non-overlapping markers output to file:\n", nonoverlappingFile, "\n")

}

################################################################################
# End of file.
################################################################################
