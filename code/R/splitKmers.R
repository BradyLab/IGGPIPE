#######################################################################################
# See usage below for description.
# Author: Ted Toal
# Date: 2013-2014
# Brady Lab, UC Davis
#######################################################################################

# Enclose everything in braces so stop statements will work correctly.
{

# Pathname separator.
PATHSEP = ifelse(grepl("/", Sys.getenv("HOME")), "/", "\\")

# Get arguments.
testing = FALSE ; testGenome1 = FALSE
#testing = TRUE ; testGenome1 = TRUE # For testing only.
{
if (!testing)
    args = commandArgs(TRUE)
# For testing only:
else
    {
    if (testGenome1)
        args = c("~/Documents/UCDavis/BradyLab/Genomes/kmers/SCARF",
           "outTestHP11/Kmers/Split/Genome_1_", "outTestHP11/Kmers/Genome_1.isect.kmers")
    else
        args = c("~/Documents/UCDavis/BradyLab/Genomes/kmers/SCARF",
           "outTestHP11/Kmers/Split/Genome_2_", "outTestHP11/Kmers/Genome_2.isect.kmers",
           "outTestHP11/Kmers/Genome_1.isect.kmers")
    }
}

NexpectedMin = 3
if (length(args) < NexpectedMin)
    {
    usage = c(
        "Read a file of unique-k-mer positions for a genome and split it up into multiple",
        "files based on the sequence ID where the k-mer was found or based upon the sequence",
        "ID where it was found within another genome.  A different output file is assigned",
        "to each sequence ID.  The filename includes both the genome name and the sequence",
        "ID.  Generally the IDs will be chromosome IDs in ONE genome.  For example, the first",
        "run might be for genome 'Heinz', output files include the Heinz chromosome IDs, and",
        "no <guideKmersFile> argument is specified.  The second run might be for genome 'Penn',",
        "the <guideKmersFile> argument is specified as the Heinz k-mer positions file, and",
        "that file is read and used to sort k-mer position lines from <kmersFile> into output",
        "files, so the output files again include the Heinz, not the Penn, chromosome IDs.",
        "This allows the total set of k-mers to be split into several (=number of sequence",
        "IDs in first genome) smaller pieces, with files named with a given sequence ID all",
        "containing the SAME k-mers, so that LCRs on that sequence ID in the first genome may",
        "later be constructed by using the common unique k-mers in each genome.",
        "The reason for this splitting of the k-mer position data is that the input k-mer",
        "file size is too large to work with, especially when multiple genomes, each with",
        "its own k-mer file, are being processed together.  By splitting, they can be",
        "processed one sequence ID at a time (using the sequence IDs from just one of the",
        "genomes, say genome 1).  The input and output files are tab-separated files.",
        "",
        "Usage: Rscript splitKmers.R <wd> <outPrefix> <kmersFile> [<guideKmersFile>]",
        "",
        "Arguments:",
        "   <wd> : Path of R working directory, specify other file paths relative to this.",
        "   <outPrefix> : Prefix for output files, including a directory path if desired.",
        "                 This must be different for each genome.  The suffix added on",
        "                 to this is the sequence ID followed by '.isect.split'.",
        "   <kmersFile> : Pathname of input k-mer file to be split, listing k-mers with position.",
        "   <guideKmersFile> : (optional) Pathname of another input k-mer file, to be used",
        "                      to guide the sorting of <kmersFile> lines to output files.",
        )
    for (S in usage)
        cat(S, "\n", sep="")
    stop("Try again with correct number of arguments")
    }

cat("splitKmers.R arguments:\n")
workingDirectory = args[1]
cat("  workingDirectory: ", workingDirectory, "\n")
if (!dir.exists(workingDirectory))
    stop("Directory doesn't exist: ", workingDirectory)
setwd(workingDirectory)

outPrefix = args[2]
cat("  outPrefix: ", outPrefix, "\n")

kmersFile = args[3]
cat("  kmersFile: ", kmersFile, "\n")
if (!file.exists(kmersFile))
    stop("File doesn't exist: ", kmersFile)
guideKmersFile = NULL

if (length(args) >= 4)
    {
    guideKmersFile = args[4]
    cat("  guideKmersFile: ", guideKmersFile, "\n")
    if (!file.exists(guideKmersFile))
        stop("File doesn't exist: ", guideKmersFile)
    }

########################################
# Initialize.
########################################

# Expected names of columns in output file.
expectedColNames = c("kmer", "seqID", "pos", "strand", "contig", "contigPos")

# Number of lines to read at one shot from large files.
nLinesAtOnce = 100000

########################################
# If a guide file was specified, read it and save the mapping of k-mers to IDs
# in vector kmerSeqIDs.
########################################
if (!is.null(guideKmersFile))
    {
    cat("Read file ", guideKmersFile, "\n")

    # Open guide file, read a chunk of nLinesAtOnce lines, including the header,
    # and confirm header length and column names.  We will read the file in
    # chunks of nLinesAtOnce, because the file is LARGE.
    inFile = file(guideKmersFile, open="r")
    df = read.table(inFile, header=TRUE, sep="\t", nrows=nLinesAtOnce, stringsAsFactors=FALSE)
    if (ncol(df) != length(expectedColNames))
        stop("Expected ", expectedNumCols, " columns in input file ", guideKmersFile, " but there were ", ncol(df))
    if (any(colnames(df) != expectedColNames))
        stop("Actual column names in input file ", guideKmersFile, " are not as expected.")

    # Loop reading data from the input file and saving data in kmerSeqIDs.
    N = 0 # Count lines.
    dfKmerSeqIDs = NULL
    while (nrow(df) != 0)
        {
        dfKmerSeqIDs = rbind(dfKmerSeqIDs, df[, c("kmer", "seqID")])
        N = N + nrow(df)
        cat(round(N/1e6, 2), "M lines\n")
        df = read.table(inFile, header=FALSE, sep="\t", nrows=nLinesAtOnce, col.names=expectedColNames,
            stringsAsFactors=FALSE)
        }
    close(inFile)

    # Make a vector of seq IDs, whose names are the k-mers.
    kmerSeqIDs = dfKmerSeqIDs$seqID
    names(kmerSeqIDs) = dfKmerSeqIDs$kmer
    rm(dfKmerSeqIDs)
    }

########################################
# Read genome k-mer file and split k-mers into different output files.
########################################

# Initialize for reading the genome k-mer file and writing the output to an output
# file.  As needed, create a new output file, and accumulate output filenames in
# vector outFilenames, whose values are the output filenames and whose names are
# the sequence IDs whose k-mers are to be written to that output file.
outFilenames = c()
getFileNames = function(IDs) return(paste(outPrefix, IDs, ".isect.split", sep=""))

# Open input file, read a chunk of nLinesAtOnce lines, including the header,
# and confirm header length and column names.  We will read the file in
# chunks of nLinesAtOnce, because the file is LARGE.
cat("Read and split file ", kmersFile, "\n")
inFile = file(kmersFile, open="r")
df = read.table(inFile, header=TRUE, sep="\t", nrows=nLinesAtOnce, stringsAsFactors=FALSE)
if (ncol(df) != length(expectedColNames))
    stop("Expected ", expectedNumCols, " columns in input file ", kmersFile, " but there were ", ncol(df))
if (any(colnames(df) != expectedColNames))
    stop("Actual column names in input file ", kmersFile, " are not as expected.")

# Loop reading data from the input file and writing them to the appropriate
# output files.
N = 0 # Count lines.
while (nrow(df) != 0)
    {
    # Get the set of IDs for each k-mer.  If guideKmersFile was not specified, we
    # the k-mer's own ID.  If guideKmersFile was specified, we use kmerSeqIDs[]
    # to provide the ID for each k-mer.
    if (is.null(guideKmersFile))
        {
        seqIDs = df$seqID
        names(seqIDs) = df$kmer
        }
    else
        {
        if (any(!df$kmer %in% names(kmerSeqIDs)))
            stop("Expected every k-mer in input file to be listed in <guideKmersFile>, but not all are")
        seqIDs = kmerSeqIDs[df$kmer]
        }

    # Get the set of unique IDs and make sure we have an output file and k-mer
    # vector for each one.  For new output files, create the file by writing only
    # the header line.
    newIDs = setdiff(unique(seqIDs), names(outFilenames))
    if (length(newIDs) > 0)
        {
        newFilenames = getFileNames(newIDs)
        names(newFilenames) = newIDs
        outFilenames = c(outFilenames, newFilenames)
        for (ID in newIDs)
            write.table(df[c(),], outFilenames[[ID]], row.names=FALSE, quote=FALSE, sep="\t")
        }

    # For each unique output ID, write the data to the corresponding output file.
    invisible(tapply(1:length(seqIDs), seqIDs, function(ii)
        {
        write.table(df[ii,], outFilenames[[seqIDs[ii[1]]]], row.names=FALSE,
            col.names=FALSE, append=TRUE, quote=FALSE, sep="\t")
        }))
    N = N + nrow(df)
    cat(round(N/1e6, 2), "M lines\n")
    df = read.table(inFile, header=FALSE, sep="\t", nrows=nLinesAtOnce, col.names=expectedColNames,
        stringsAsFactors=FALSE)
    }
# Close the input file.
close(inFile)
cat("Split completed, done.\n")
}
