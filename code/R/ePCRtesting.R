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
#testing = 1 # For testing only, outTestHP11, genome 1.
#testing = 2 # For testing only, outTestHP11, genome 2.
{
if (testing == 0)
    args = commandArgs(TRUE)
else if (testing == 1)
    {
    args = c("~/Documents/UCDavis/BradyLab/Genomes/kmers/IGGPIPE",
        "outTestHP11/NonvalidatedMarkers_K11k2L100D10_2000A100_2000d10_100N2F0X20.tsv",
        1, "outTestHP11/MarkerErrors_K11k2L100D10_2000A100_2000d10_100N2F0X20V3000W8M3G1_H.bad.tsv",
        "outTestHP11/GenomeData", "/Users/tedtoal/bin/e-PCR", 3000, 8, 3, 1,
        "testFASTA/ITAG2.4_test.fasta", TRUE)
    }
else if (testing == 2)
    {
    args = c("~/Documents/UCDavis/BradyLab/Genomes/kmers/IGGPIPE",
        "outTestHP11/NonvalidatedMarkers_K11k2L100D10_2000A100_2000d10_100N2F0X20.tsv",
        2, "outTestHP11/MarkerErrors_K11k2L100D10_2000A100_2000d10_100N2F0X20V3000W8M3G1_P.bad.tsv",
        "outTestHP11/GenomeData", "/Users/tedtoal/bin/e-PCR", 3000, 8, 3, 1,
        "testFASTA/Spenn2.0_test.fasta", TRUE)
    }
else stop("Unknown value for 'testing'")
}

Nexpected = 12
if (length(args) != Nexpected)
    {
    usage = c(
        "Read a data frame of primer pairs for candidate IGG markers and test each pair using e-PCR.",
        "",
        "Usage: Rscript ePCRtesting.R <wd> <tsvMarkerFile> <genomeNum> <badMarkerFile> <pcrInfoDir> \\",
        "       <ePCRpath> <maxDeviation> <wordSize> <maxMismatches> <maxGaps> <fastaFile> <investigate>",
        "",
        "Arguments:",
        "   <wd> : Path of R working directory, specify other file paths relative to this.",
        "   <tsvMarkerFile> : Input file containing the candidate markers with primer pairs.",
        "   <genomeNum>     : The genome number whose primers are to be tested, 1=reference genome, etc.",
        "   <badMarkerFile> : Output file to be written containing info on any markers failing the ePCR.",
        "   <pcrInfoDir>    : Directory to hold e-PCR input and output files.",
        "   <ePCRpath>      : Full path of the e-PCR program.",
        "   <maxDeviation>  : Maximum deviation from target amplicon size and still report it.",
        "   <wordSize>      : Word size of e-PCR hash table in base-pairs, also sets number of 3' primer",
        "                     end base-pairs that must match exactly.",
        "   <maxMismatches> : Maximum mismatches in an off-target primer hit and still report it.",
        "   <maxGaps>       : Maximum gaps in an off-target primer hit and still report it.",
        "   <fastaFile>     : FASTA file containing genomic sequences for genome 1.",
        "   <investigate>   : FALSE for normal operation, TRUE for more verbose debugging output."
        )
    for (S in usage)
        catnow(S, "\n", sep="")
    stop("Try again with correct number of arguments")
    }

catnow("ePCRtesting.R arguments:\n")
workingDirectory = args[1]
catnow("  workingDirectory: ", workingDirectory, "\n")
if (!dir.exists(workingDirectory))
    stop("Directory doesn't exist: ", workingDirectory)
setwd(workingDirectory)

tsvMarkerFile = args[2]
catnow("  tsvMarkerFile: ", tsvMarkerFile, "\n")
if (!file.exists(tsvMarkerFile))
    stop("File doesn't exist: ", tsvMarkerFile)

genomeNum = as.integer(args[3])
catnow("  genomeNum: ", genomeNum, "\n")
if (genomeNum < 1)
    stop("genomeNum must be >= 1")

badMarkerFile = args[4]
catnow("  badMarkerFile: ", badMarkerFile, "\n")

pcrInfoDir = args[5]
catnow("  pcrInfoDir: ", pcrInfoDir, "\n")
if (!dir.exists(pcrInfoDir))
    stop("Directory doesn't exist: ", pcrInfoDir)

ePCRpath = args[6]
catnow("  ePCRpath: ", ePCRpath, "\n")

maxDeviation = as.integer(args[7])
catnow("  maxDeviation: ", maxDeviation, "\n")
if (is.na(maxDeviation) || maxDeviation < 10 || maxDeviation > 100000)
    stop("maxDeviation must be between 10 and 100000")

wordSize = as.integer(args[8])
catnow("  wordSize: ", wordSize, "\n")
if (is.na(wordSize) || wordSize < 5 || wordSize > 25)
    stop("wordSize must be between 5 and 25")

maxMismatches = as.integer(args[9])
catnow("  maxMismatches: ", maxMismatches, "\n")
if (is.na(maxMismatches) || maxMismatches < 0 || maxMismatches > 10)
    stop("maxMismatches must be between 0 and 10")

maxGaps = as.integer(args[10])
catnow("  maxGaps: ", maxGaps, "\n")
if (is.na(maxGaps) || maxGaps < 0 || maxGaps > 10)
    stop("maxGaps must be between 0 and 10")

fastaFile = args[11]
catnow("  fastaFile: ", fastaFile, "\n")
if (!file.exists(fastaFile))
    stop("File doesn't exist: ", fastaFile)

investigate = as.logical(args[12])
catnow("  investigate: ", investigate, "\n")
if (is.na(investigate))
    stop("investigate must be TRUE or FALSE")

catnow("Preparing to do electronic PCR using primers of genome", genomeNum, "on its FASTA file\n")

########################################
# Initialize.
########################################

# Read candidate marker data including primers.
df = read.table(tsvMarkerFile, header=TRUE, row.names=NULL, sep="\t", stringsAsFactors=FALSE)
if (nrow(df) == 0)
    stop("There are no candidate markers.")
inv(dim(df), "input data dim")
inv(colnames(df), "input data columns")
inv(head(df), "input data head")

# The k-mer length we are working with.
kmerLen = nchar(df$kmer1[1])
inv(kmerLen, "k")

# Get the genome letter for this genome.
genomeLtr = substr(colnames(df)[grepl("^.id$", colnames(df))], 1, 1)
Ngenomes = length(genomeLtr)
if (Ngenomes < 2)
    stop("Expected at least two id columns in <tsvMarkerFile>")
refGenomeLtr = genomeLtr[1]
if (genomeNum > Ngenomes)
    stop("genome number is larger than the number of genomes in <tsvMarkerFile>")
genomeLtr = genomeLtr[genomeNum]
inv(genomeLtr, "genome letter")

############################################################################
# Create a ".epcr.in" tab-separated text file to use as input to e-PCR,
# giving the primer pairs to be tested.  Example:
#
# 1 TGAAGATTACGACTTGAAGCTCC TGCTCACGAGTTCAGGAGAT    1506    H=1506,P=691
#
# The last column is info only and can contain anything.  We place the df
# row number in the first column.
############################################################################

# Amplicon length columns.
allAmpCols = colnames(df)[grepl("ampLen$", colnames(df))]
ampCol = paste(genomeLtr, "ampLen", sep="")

dfInfo = sapply(allAmpCols, function(ampCol) paste(substr(ampCol, 1, 1), df[, ampCol], sep="="))

dfEPCR = data.frame(rowNum=1:nrow(df), prmSeqL=df$prmSeqL, prmSeqR=df$prmSeqR,
    ampLen=df[, ampCol], info=apply(dfInfo, 1, paste, collapse=","),
    stringsAsFactors=FALSE)

ePCRinputFile = paste("Genome_", genomeNum, ".epcr.in", sep="")
ePCRinputFile = paste(pcrInfoDir, ePCRinputFile, sep=PATHSEP)
cat("ePCRinputFile='", ePCRinputFile, "'\n", sep="")
write.table(dfEPCR, ePCRinputFile, row.names=FALSE, col.names=FALSE, quote=FALSE, sep="\t")

########################################
# Now run the e-PCR command.  This can take a long time, because the entire FASTA
# file must be read and searched for each primer.
########################################

ePCRoutputFile = paste("Genome_", genomeNum, ".epcr.out", sep="")
ePCRoutputFile = paste(pcrInfoDir, ePCRoutputFile, sep=PATHSEP)
cmdLine = paste(ePCRpath, "-v+", "-p+", "-t3",
    paste("-m", maxDeviation, sep=""),
    paste("-w", wordSize, sep=""),
    paste("-n", maxMismatches, sep=""),
    paste("-g", maxGaps, sep=""),
    "-o", ePCRoutputFile, ePCRinputFile, fastaFile)
inv(nrow(dfEPCR), "nrow(dfEPCR)")
inv(cmdLine, "e-PCR command line")
catnow("Running e-PCR command to search for marker primer sequences from FASTA file\n")
catnow("   ", cmdLine, "\n")
system(cmdLine)

########################################
# Read the output file and check to see that each pair had one and only one
# amplicon, at the correct position.
########################################

if (!file.exists(ePCRoutputFile)) stop("e-PCR output file not found: ", ePCRoutputFile)
dfAmps = read.table(ePCRoutputFile, header=FALSE, row.names=NULL, sep="\t", stringsAsFactors=FALSE)
colnames(dfAmps) = c("id", "rowNum", "strand", "posL", "posR", "lenExpLens", "mismatches", "gaps", "info")

# Group the results by df row number.
L = split(1:nrow(dfAmps), dfAmps$rowNum)

# Convert character row number names(L) to integers and re-order L so it can be
# indexed by a df row number.
n = as.integer(names(L))
L = L[match(1:nrow(df), n)]

# Find those found zero times, those found more than once, and those found once
# with wrong sequence id, posL or posR.  Note: ampPos1 and ampPos2 may be in
# reverse order so that ampPos1 > ampPos2.  However, it is always the case that
# posL < posR.
len = sapply(L, length)
foundZero = (len == 0)
foundOnce = (len == 1)
foundMultiple = (len > 1)
idx1 = sapply(L, '[', 1)
idCol = paste(genomeLtr, "id", sep="")
ampPos1Col = paste(genomeLtr, "ampPos1", sep="")
ampPos2Col = paste(genomeLtr, "ampPos2", sep="")
foundOnceIdWrong = (foundOnce & df[, idCol] != dfAmps$id[idx1])
posL = df[, ampPos1Col]
posR = df[, ampPos2Col]
swap = (posL > posR)
tmp = posL[swap]
posL[swap] = posR[swap]
posR[swap] = tmp
foundOncePosLwrong = (foundOnce & posL != dfAmps$posL[idx1])
foundOncePosRwrong = (foundOnce & posR!= dfAmps$posR[idx1])
rmvMarker = (foundZero | foundMultiple | foundOnceIdWrong | foundOncePosLwrong | foundOncePosRwrong)
{
cat("Number candidate markers:        ", nrow(df), "\n")
cat("No amplicon found:               ", sum(foundZero), "\n")
cat("Multiple amplicons found:        ", sum(foundMultiple), "\n")
cat("One amplicon but wrong id:       ", sum(foundOnceIdWrong), "\n")
cat("One amplicon but wrong left pos: ", sum(foundOncePosLwrong), "\n")
cat("One amplicon but wrong right pos:", sum(foundOncePosRwrong), "\n")
cat("----------------------------\n")
cat("Total markers to be removed:", sum(rmvMarker), "\n")
}

########################################
# Write output file containing markers to be removed, with "reason" column.
########################################

dfRemove = df[rmvMarker,]
foundZero = foundZero[rmvMarker]
foundMultiple = foundMultiple[rmvMarker]
foundOnceIdWrong = foundOnceIdWrong[rmvMarker]
foundOncePosLwrong = foundOncePosLwrong[rmvMarker]
foundOncePosRwrong = foundOncePosRwrong[rmvMarker]
if (nrow(dfRemove) > 0)
    {
    dfRemove = data.frame(reasonDiscarded="unknown", dfRemove, stringsAsFactors=FALSE)
    dfRemove$reasonDiscarded[foundZero] = "not found"
    dfRemove$reasonDiscarded[foundMultiple] = "found multiple"
    dfRemove$reasonDiscarded[foundOnceIdWrong] = "wrong seq id"
    dfRemove$reasonDiscarded[!foundOnceIdWrong & foundOncePosLwrong & foundOncePosRwrong] = "wrong pos"
    dfRemove$reasonDiscarded[!foundOnceIdWrong & foundOncePosLwrong & !foundOncePosRwrong] = "wrong posL"
    dfRemove$reasonDiscarded[!foundOnceIdWrong & !foundOncePosLwrong & foundOncePosRwrong] = "wrong posR"
    }
else
    {
    catnow("No bad markers to be removed.\n")
    dfRemove = data.frame(reasonDiscarded=character(), df[c(),], stringsAsFactors=FALSE)
    }
# Write the data frame to the output file.
write.table(dfRemove, badMarkerFile, row.names=FALSE, quote=FALSE, sep="\t")
# dfRemove = read.table(badMarkerFile, header=TRUE, row.names=NULL, sep="\t", stringsAsFactors=FALSE)

catnow("Finished testing primers of candidate markers for genome", genomeNum, "\n")
catnow("Markers to remove are in file:\n", badMarkerFile, "\n")
}
