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
#testing = 1 # For testing only, outTestHP11 LCBs
#testing = 2 # For testing only, outTestHP11 IndelGroupsNonoverlapping
#testing = 3 # For testing only, outTestHP11 MarkersNonoverlapping
#testing = 4 # For testing only, outTaCW15 LCBs
{
if (testing == 0)
    args = commandArgs(TRUE)
else if (testing == 1)
    {
    args = c("~/Documents/UCDavis/BradyLab/IGGPIPE/IGGPIPE", 11,
        "outTestHP11/LCBs_K11k2L100D10_2000.tsv",
        "outTestHP11/LCBs_K11k2L100D10_2000.withseqs.tsv",
        "outTestHP11/GenomeData",
        "/Users/tedtoal/perl5/perlbrew/perls/perl-5.14.2/bin/perl",
        "code/perl/getSeqsFromFasta.pl",
        TRUE,
        "testFASTA/ITAG2.4_test.fasta", "testFASTA/SpennV2.0_test.fasta")
    }
else if (testing == 2)
    {
    args = c("~/Documents/UCDavis/BradyLab/IGGPIPE/IGGPIPE", 11,
        "outTestHP11/IndelGroupsNonoverlapping_K11k2L100D10_2000A100_2000d10_100N2F0.tsv",
        "outTestHP11/IndelGroupsNonoverlapping_K11k2L100D10_2000A100_2000d10_100N2F0.withseqs.tsv",
        "outTestHP11/GenomeData",
        "/Users/tedtoal/perl5/perlbrew/perls/perl-5.14.2/bin/perl",
        "code/perl/getSeqsFromFasta.pl",
        TRUE,
        "testFASTA/ITAG2.4_test.fasta", "testFASTA/SpennV2.0_test.fasta")
    }
else if (testing == 3)
    {
    args = c("~/Documents/UCDavis/BradyLab/IGGPIPE/IGGPIPE", 11,
        "outTestHP11/MarkersNonoverlapping_K11k2L100D10_2000A100_2000d10_100N2F0X20V3000W8M3G1.tsv",
        "outTestHP11/MarkersNonoverlapping_K11k2L100D10_2000A100_2000d10_100N2F0X20V3000W8M3G1.withseqs.tsv",
        "outTestHP11/GenomeData",
        "/Users/tedtoal/perl5/perlbrew/perls/perl-5.14.2/bin/perl",
        "code/perl/getSeqsFromFasta.pl",
        TRUE,
        "testFASTA/ITAG2.4_test.fasta", "testFASTA/SpennV2.0_test.fasta")
    }
else if (testing == 4)
    {
    args = c("~/Documents/UCDavis/BradyLab/IGGPIPE/IGGPIPE", 15,
        "outTaCW15/LCBs_K15k2L200D5_1000.tsv",
        "outTaCW15/LCBs_K15k2L200D5_1000.withseqs.tsv",
        "outTaCW15/GenomeData",
        "/Users/tedtoal/perl5/perlbrew/perls/perl-5.14.2/bin/perl",
        "code/perl/getSeqsFromFasta.pl",
        TRUE,
        "/Users/tedtoal/Documents/UCDavis/BradyLab/Genomes/Ta/Genome_IWGSC1.0/Triticum_aestivum.IWGSC1.0+popseq.28.dna.genome.fa",
        "/Users/tedtoal/Documents/UCDavis/BradyLab/Genomes/Ta/Genome_W7984_scaffoldsMar28/w7984.meraculous.scaffolds.Mar28_contamination_removed.fa")
    }
else stop("Unknown value for 'testing'")
}

NexpectedMin = 10
if (length(args) < NexpectedMin)
    {
    usage = c(
        "Read a data frame of LCBs, Indel Groups, or markers, extract DNA sequences from all",
        "genomes between the start and end position of each one and add the sequences to the",
        "data frame, then write the data frame to the output file.",
        "",
        "Usage: Rscript getDNAseqsForIndelsSNPs.R <wd> <inputFile> <outputFile> \\",
        "       <extractDir> <perlPath> <getSeqsFromFasta> <investigate> \\",
        "       <fastaFile1> <fastaFile2> ...",
        "",
        "Arguments:",
        "   <wd> : Path of R working directory, specify other file paths relative to this.",
        "   <k>  : Value of k, i.e. k-mer size.",
        "   <inputFile>   : Input file containing LCBs, Indel groups, or markers.",
        "   <outputFile>  : Output file to which to write prepped input data with DNA sequences.",
        "   <extractDir>  : Directory to hold DNA sequence extraction files.",
        "   <perlPath>    : Full path of the Perl language interpreter program",
        "   <getSeqsFromFasta> : Full path of the Perl script getSeqsFromFasta.pl",
        "   <investigate> : FALSE for normal operation, TRUE for more verbose debugging output.",
        "   <fastaFile1>  : Genome 1 FASTA file.",
        "   <fastaFile2>  : Genome 2 FASTA file.",
        "   ...           : FASTA files of other genomes, for all the genomes."
        )
    for (S in usage)
        catnow(S, "\n", sep="")
    stop("Try again with correct number of arguments")
    }

catnow("getDNAseqsForIndelsSNPs.R arguments:\n")
workingDirectory = args[1]
catnow("  workingDirectory: ", workingDirectory, "\n")
if (!dir.exists(workingDirectory))
    stop("Directory doesn't exist: ", workingDirectory)
setwd(workingDirectory)

kmerLen = as.integer(args[2])
catnow("  kmerLen: ", kmerLen, "\n")
if (is.na(kmerLen) || kmerLen < 1)
    stop("kmerLen must be > 0")

inputFile = args[3]
catnow("  inputFile: ", inputFile, "\n")
if (!file.exists(inputFile))
    stop("File doesn't exist: ", inputFile)

outputFile = args[4]
catnow("  outputFile: ", outputFile, "\n")

extractDir = args[5]
catnow("  extractDir: ", extractDir, "\n")
if (!dir.exists(extractDir))
    stop("Directory doesn't exist: ", extractDir)

perlPath = args[6]
catnow("  perlPath: ", perlPath, "\n")

getSeqsFromFasta = args[7]
catnow("  getSeqsFromFasta: ", getSeqsFromFasta, "\n")
if (!file.exists(getSeqsFromFasta))
    stop("File doesn't exist: ", getSeqsFromFasta)

investigate = as.logical(args[8])
catnow("  investigate: ", investigate, "\n")
if (is.na(investigate))
    stop("investigate must be TRUE or FALSE")

fastaFiles = args[-(1:8)]
if (length(fastaFiles) < 2)
    stop("Must have at least two FASTA files specified")
catnow("  fastaFiles:\n")
for (f in fastaFiles)
    {
    catnow("   ", f, "\n")
    if (!file.exists(f))
        stop("File doesn't exist: ", f)
    }

########################################
# Initialize.
########################################

catnow("Preparing to retrieve DNA sequences from FASTA files\n")

# Read input data.
df = read.table(inputFile, header=TRUE, row.names=NULL, sep="\t", stringsAsFactors=FALSE)
if (nrow(df) == 0)
    stop("There is no input data.")
inv(dim(df), "input data dim")
inv(colnames(df), "input data columns")
inv(head(df), "input data head")

# Convert df into a data frame with columns 'ID', 'phases', 'Xid', 'Xpos1',
# and 'Xpos2', where X=genome letter.  This conversion involves doing different
# things depending on whether the input file was LCBs, Indel Groups, or markers.
# Set the Xpos1 and Xpos2 values to be the exact start and end positions of the
# outside k-mers, so a DNA extraction will start and end with the k-mers.  For
# LCBs, each LCB may become one or more output data frame rows.  (The entire LCB
# could be one row, but the size can be too large for alignment, so it is broken
# up into smaller segments at k-mer boundaries).

# If this is an LCB file (has an "LCB" column), convert each LCB into one or more
#       output rows, trying to keep length no more than maxOutLen, and ensuring
#       that pos1 and pos2 are set correctly and offset correctly so they are
#       PRECISELY the positions of the first bases on the outside side of each
#       k-mer that defines that output row.
{
if (any(colnames(df) == "LCB"))
    {
    # Get vectors of names of X columns.
    idCols = colnames(df)[grepl("id$", colnames(df))]
    posCols = colnames(df)[grepl("pos$", colnames(df))]
    genomeLtrs = substring(idCols, 1, 1)
    otherGenomeLtrs = genomeLtrs[-1]
    if (length(genomeLtrs) < 2)
        stop("Expected at least two id columns in <inputFile>")
    inv(genomeLtrs, "genome letters")
    names(idCols) = genomeLtrs
    names(posCols) = genomeLtrs
    refIdCol = idCols[1]
    refPosCol = posCols[1]

    # Break the LCB up into one or more output rows and build the rows.  A single
    # LCB can be so long that we can't align it.  Keep the maximum length of the
    # output rows short enough that we will be able to align them.
    maxOutLen = 500 # Try to keep output length less than this.
    df = df[order(df$LCB, df[, refIdCol], df[, refPosCol]),] # Make sure df is sorted by reference genome within each LCB.
    # Get indexes indicating starting and ending row of each LCB.  Note that it
    # is always true that length(startLCB) = length(endLCB)
    N = nrow(df)
    startLCB = which(c(TRUE, df$LCB[-1] != df$LCB[-N]))
    endLCB = which(c(df$LCB[-1] != df$LCB[-N], TRUE))
    # Accumulate df row indexes indicating start and end k-mers of each LCB span
    # that is to become an output row.  Note that it is always true that
    # length(startRows) = length(endRows)
    startRows = c()
    endRows = c()
    # curStarts are df row indexes of rows that will be the start k-mers of the
    # next output data frame rows, initially the first df row of each LCB, and
    # subsequently the same row that was the end k-mer of the previous output row.
    # (The output rows of one LCB will overlap by one k-mer).
    curStarts = startLCB
    # curEnds are df row indexes of the rows currently being tested to see if they
    # should be end-of-span k-mers, initially the SECOND df row of each LCB.
    # Note that it is always true that length(curStarts) = length(curEnds)
    # It will also always be maintained true that all(curEnds-curStarts >= 1)
    curEnds = curStarts+1
    # Loop moving curEnds along the df rows and watching for end of LCB or
    # maxOutLen being exceeded.  Record a span and start a new span when
    # maxOutLen is exceeded.  Record the final span of an LCB when the end of
    # the LCB is reached.
    while (length(curEnds) > 0)
        {
        # If curEnds is the last row of the LCB, move those curStarts/curEnds to
        # startRows/endRows.
        isEndLCB = (curEnds %in% endLCB)
        startRows = c(startRows, curStarts[isEndLCB])
        endRows = c(endRows, curEnds[isEndLCB])

        # Remove them from curStarts/curEnds.
        curStarts = curStarts[!isEndLCB]
        curEnds = curEnds[!isEndLCB]

        # If curEnds+1 exceeds maxOutLen, move those curStarts/curEnds to
        # startRows/endRows.
        # Check to see if the k-mer following the curEnds k-mer would cause the next
        # output row length to exceed the maximum.  If so, make a new output row
        # ending at the current k-mer.  
        nextCurEnds = curEnds+1
        nextExceedsMax = (df[nextCurEnds, refPosCol] - df[curStarts, refPosCol] > maxOutLen)
        startRows = c(startRows, curStarts[nextExceedsMax])
        endRows = c(endRows, curEnds[nextExceedsMax])

        # And, remove the start k-mer from curStarts and replace it with the end k-mer
        # in curEnds, so that the next row will start with the same k-mer that this row
        # ends with.  Move the maxed-out curEnds to the end of the curEnds vector.
        curStarts = c(curStarts[!nextExceedsMax], curEnds[nextExceedsMax])
        curEnds = c(curEnds[!nextExceedsMax], curEnds[nextExceedsMax])

        # Finally, advance curEnds by one position.
        curEnds = curEnds+1
        }
    # head(sort(startRows))
    # head(sort(endRows))
    # tail(sort(startRows))
    # tail(sort(endRows))

    # Collect the new rows in dfx.
    pos1Cols = paste(genomeLtrs, "pos1", sep="")
    pos2Cols = paste(genomeLtrs, "pos2", sep="")
    names(pos1Cols) = genomeLtrs
    names(pos2Cols) = genomeLtrs
    refPos1Col = pos1Cols[1]
    refPos2Col = pos2Cols[1]
    dfx = df[startRows, colnames(df) != "LCB"]
    colnames(dfx) = sub("pos$", "pos1", colnames(dfx))
    for (genome in genomeLtrs)
        dfx[, pos2Cols[genome]] = df[endRows, posCols[genome]]
    # Make the ID.
    dfx$ID = paste(dfx[, refIdCol], dfx[, refPos1Col], dfx[, refPos2Col], sep="_")

    # Now change pos1 and pos2 so that instead of being the 5' end of the k-mer, they
    # are the positions of the first bases on the outside side of each k-mer.
    for (i in 1:length(genomeLtrs))
        {
        genome = genomeLtrs[i]
        pos1Col = pos1Cols[genome]
        pos2Col = pos2Cols[genome]
        minusPhase = (substr(dfx$phases, i, i) == "-")

        # Adjustment for a "+" phase genome to go from 5' end of k-mer to outside base of k-mer:
        #   pos1 = pos1
        #   pos2 = pos2 + (k-1)
        # Adjustment for a "-" phase genome to go from 5' end of k-mer to outside base of k-mer:
        #   pos1 = pos1 + (k-1)
        #   pos2 = pos2

        dfx[!minusPhase, pos2Col] = dfx[!minusPhase, pos2Col] + rep(kmerLen-1, sum(!minusPhase))
        dfx[minusPhase, pos1Col] = dfx[minusPhase, pos1Col] + rep(kmerLen-1, sum(minusPhase))

        # Also, swap the positions of "-" phase positions so that pos1 < pos2 always.
        t = dfx[minusPhase, pos1Col]
        dfx[minusPhase, pos1Col] = dfx[minusPhase, pos2Col]
        dfx[minusPhase, pos2Col] = t

        rm(minusPhase)
        }
    df = dfx
    rm(dfx)
    }

# If this is an IndelGroups file (has an "Xctg1" column), convert it by:
#   (1) remove row names, if any
#   (2) retain only the needed columns: 'Xid', 'Xpos1', 'Xpos2'
#   (3) add an "ID" column by concatenating reference Xid with reference Xpos1
#       and Xpos2 numbers
#   (4) create a "phase" column which has "-" if Xpos1 > Xpos2, swapping Xpos1/Xpos2
#       columns when that is the case so Xpos1 < Xpos2 always, and adjusting Xpos1/Xpos2
#       to change the positions from the 5' end of the "+" strand k-mer to being
#       PRECISELY the positions of the first bases on the outside side of each k-mer.
else if (any(grepl("ctg1$", colnames(df))))
    {
    # Make sure kmerLen is consistent.
    if (kmerLen != nchar(df$kmer1[1]))
        stop("kmerLen of ", kmerLen, " did not match input file k-mer size")
    # Remove row names.
    rownames(df) = NULL
    # Get vectors of names of X columns.
    idCols = colnames(df)[grepl("id$", colnames(df))]
    pos1Cols = colnames(df)[grepl("pos1$", colnames(df))]
    pos2Cols = colnames(df)[grepl("pos2$", colnames(df))]
    genomeLtrs = substring(idCols, 1, 1)
    otherGenomeLtrs = genomeLtrs[-1]
    if (length(genomeLtrs) < 2)
        stop("Expected at least two id columns for two genomes in <inputFile>")
    inv(genomeLtrs, "genome letters")
    names(idCols) = genomeLtrs
    names(pos1Cols) = genomeLtrs
    names(pos2Cols) = genomeLtrs
    # Retain only needed columns.
    df = df[, c(idCols, pos1Cols, pos2Cols)]
    # Create ID column.
    df$ID = paste(df[,idCols[1]], df[,pos1Cols[1]], df[,pos2Cols[1]], sep="_")
    # Create phases column, swap Xpos1/Xpos2 when backwards (means k-mer "-" strand), add k-1 to Xpos2.
    df$phases = ""

    # Adjustment for a "+" phase genome to go from 5' end of k-mer to outside base of k-mer:
    #   pos1 = pos1
    #   pos2 = pos2 + (k-1)
    # Adjustment for a "-" phase genome to go from 5' end of k-mer to outside base of k-mer:
    #   pos1 = pos1 + (k-1)
    #   pos2 = pos2

    # In addition to these adjustments, we want to swap pos1 and pos2 if pos1 > pos2 ("-" phase).
    for (genome in genomeLtrs)
        {
        pos1Col = pos1Cols[genome]
        pos2Col = pos2Cols[genome]

        minusPhase = (df[, pos1Col] > df[, pos2Col])
        df[!minusPhase, pos2Col] = df[!minusPhase, pos2Col] + rep(kmerLen-1, sum(!minusPhase))
        df[minusPhase, pos1Col] = df[minusPhase, pos1Col] + rep(kmerLen-1, sum(minusPhase))

        t = df[minusPhase, pos1Col]
        df[minusPhase, pos1Col] = df[minusPhase, pos2Col]
        df[minusPhase, pos2Col] = t
        phases = rep("+", nrow(df))
        phases[minusPhase] = "-"
        df$phases = paste(df$phases, phases, sep="")

        rm(minusPhase, t, phases)
        }
    }

# Else this must be a Markers file, convert it by:
#   (1) remove row names if any
#   (2) convert "XampPos1/2" columns to "Xpos1/2"
#   (3) add an "ID" column by concatenating reference Xid with reference Xpos1
#       and Xpos2 numbers
#   (4) swap Xpos1/2 columns if necessary so that Xpos1 < Xpos2
#   (5) if no swap, add kmer1offset to Xpos1 and subtract kmer2offset from Xpos2,
#       if swap, add kmer2offset to Xpos1 and subtract kmer1offset from Xpos2,
#       placing Xpos1/Xpos2 at the outside ends of the k-mers.
#   (6) adjust Xpos1/Xpos2 to change the positions from the 5' end of the "+"
#       strand k-mer to being PRECISELY the positions of the first bases on the
#       outside side of each k-mer.
#   (7) retain only the needed columns: 'Xid', 'XYphase', 'Xpos1', 'Xpos2'
#   (8) concatenat columns "XYphase" to form column "phases" and delet the
#       XYphase columns.
else
    {
    # Make sure kmerLen is consistent.
    if (kmerLen != nchar(df$kmer1[1]))
        stop("kmerLen of ", kmerLen, " did not match input file k-mer size")
    # Remove row names.
    rownames(df) = NULL
    # Change XampPos to pos.
    colnames(df) = sub("ampPos", "pos", colnames(df))
    # Get vectors of names of X columns.
    idCols = colnames(df)[grepl("id$", colnames(df))]
    pos1Cols = colnames(df)[grepl("pos1$", colnames(df))]
    pos2Cols = colnames(df)[grepl("pos2$", colnames(df))]
    phaseCols = colnames(df)[grepl("phase$", colnames(df))]
    genomeLtrs = substring(idCols, 1, 1)
    otherGenomeLtrs = genomeLtrs[-1]
    if (length(genomeLtrs) < 2)
        stop("Expected at least two id columns for two genomes in <inputFile>")
    inv(genomeLtrs, "genome letters")
    names(idCols) = genomeLtrs
    names(pos1Cols) = genomeLtrs
    names(pos2Cols) = genomeLtrs
    names(phaseCols) = otherGenomeLtrs
    # Create ID column.
    df$ID = paste(df[,idCols[1]], df[,pos1Cols[1]], df[,pos2Cols[1]], sep="_")

    # Here we need to reverse the changes made by findPrimers.R to convert positions
    # from 5' end of each k-mer on + strand, to ampPos columns.  More than that, we
    # want the positions to be the OUTSIDE BASE of each k-mer, not the 5' end.

    # Adjustment for a "+" phase genome to return position to 5' end of k-mer:
    #   pos1 = ampPos1 + kmer1offset
    #   pos2 = ampPos2 - kmer2offset - (k-1)
    # Adjustment for a "-" phase genome to return position to 5' end of k-mer:
    #   pos1 = ampPos1 - kmer1offset - (k-1)
    #   pos2 = ampPos2 + kmer2offset

    # Adjustment for a "+" phase genome to go from 5' end of k-mer to outside base of k-mer:
    #   pos1 = pos1
    #   pos2 = pos2 + (k-1)
    # Adjustment for a "-" phase genome to go from 5' end of k-mer to outside base of k-mer:
    #   pos1 = pos1 + (k-1)
    #   pos2 = pos2

    # Overall adjustment for a "+" phase genome:
    #   pos1 = ampPos1 + kmer1offset
    #   pos2 = ampPos2 - kmer2offset
    # Overall adjustment for a "-" phase genome:
    #   pos1 = ampPos1 - kmer1offset
    #   pos2 = ampPos2 + kmer2offset

    # In addition to these adjustments, we want to swap pos1 and pos2 if pos1 > pos2 ("-" phase).
    for (genome in genomeLtrs)
        {
        pos1Col = pos1Cols[genome]
        pos2Col = pos2Cols[genome]

        minusPhase = (df[, pos1Col] > df[, pos2Col])
        df[!minusPhase, pos1Col] = df[!minusPhase, pos1Col] + df$kmer1offset[!minusPhase]
        df[!minusPhase, pos2Col] = df[!minusPhase, pos2Col] - df$kmer2offset[!minusPhase]
        df[minusPhase, pos1Col] = df[minusPhase, pos1Col] - df$kmer1offset[minusPhase]
        df[minusPhase, pos2Col] = df[minusPhase, pos2Col] + df$kmer2offset[minusPhase]
        t = df[minusPhase, pos1Col]
        df[minusPhase, pos1Col] = df[minusPhase, pos2Col]
        df[minusPhase, pos2Col] = t

        rm(minusPhase, t)
        }

    # Retain only needed columns.
    df = df[, c("ID", idCols, pos1Cols, pos2Cols, phaseCols)]
    # Create phases column and remove individual phase columns.
    df$phases = "+"
    for (genome in otherGenomeLtrs)
        df$phases = paste(df$phases, df[, phaseCols[genome]], sep="")
    df = df[, !colnames(df) %in% phaseCols]
    rm(phaseCols)
    }
}

# Make genome column name vectors for df.
idCols = paste(genomeLtrs, "id", sep="")
names(idCols) = genomeLtrs
pos1Cols = paste(genomeLtrs, "pos1", sep="")
names(pos1Cols) = genomeLtrs
pos2Cols = paste(genomeLtrs, "pos2", sep="")
names(pos2Cols) = genomeLtrs

# Put columns in desired order.
df = df[, c("ID", "phases", c(rbind(idCols, pos1Cols, pos2Cols)))]
catnow("Number of LCBs/Indel Groups/Markers for which to retrieve DNA sequences:", nrow(df), "\n")
if (nrow(df) == 0)
    stop("No regions to be processed.")

# Check to see if the correct number of genome FASTA files were specified.
Ngenomes = length(genomeLtrs)
if (Ngenomes != length(fastaFiles))
    stop("Expected one FASTA file for each of the ", Ngenomes, " genomes but got " , length(fastaFiles), " of them")
names(fastaFiles) = genomeLtrs

################################################################################
# Create a text file used as a sequence specification file for getSeqsFromFasta.pl
# to extract the DNA sequence between pos1 and pos2 of each genome for each row
# of df.  The extraction direction must be properly chosen: use the strand given
# by the "phases" column.  The extraction direction is specified by including
# (for the - strand) or not (for the + strand) a "!" character in front of the
# sequence position specification, which determines whether or not it reverse
# complements the extracted + strand sequence.
################################################################################

# Create extraction directory if it doesn't exist.
if (!file.exists(extractDir))
    dir.create(extractDir)

# Create filenames for extraction position files.
extractPosFiles = paste(extractDir, paste("extract_IndelGroupRgns_", 1:Ngenomes, ".txt", sep=""), sep=PATHSEP)
names(extractPosFiles) = genomeLtrs

# Put sequence extraction strings in list seqExtStrs, with one vector of strings
# member per genome and keep them for later when we read back the sequences, to
# check them.  Write the strings to text files, one per genome.
seqExtStrs = list()
Ngenome = 0
for (genome in genomeLtrs)
    {
    Ngenome = Ngenome + 1
    revComp = (substr(df$phases, Ngenome, Ngenome) == "-")
    revComp = 1+as.integer(revComp)
    revComp = c("", "!")[revComp]
    seqExtStrs[[genome]] = paste(revComp, df[, idCols[genome]], ":",
        as.integer(df[, pos1Cols[genome]]), "..", as.integer(df[, pos2Cols[genome]]), sep="")
    inv(length(seqExtStrs[[genome]]), "length(seqExtStrs[[genome]])")
    writeLines.winSafe(seqExtStrs[[genome]], extractPosFiles[genome])
    rm(revComp)
    }

################################################################################
# A command line to run getSeqsFromFasta.pl for hypothetical genome G and put
# the entire extracted sequence on one line looks like this:
#   perl blah/getSeqsFromFasta.pl -l 0 blah/genome.fasta -i extractDir/extract1.txt -o extractDir/seqs1.txt
#
# Create the command line for each genome, then run them one genome at a time.
# This can take a long time, because the genome sequence file is very large,
# and if the sequence IDs are entire chromosomes, getSeqsFromFasta.pl reads an
# entire sequence into memory before processing it.
################################################################################

seqFiles = paste(extractDir, paste("seqs_IndelGroupRgns_", 1:Ngenomes, ".txt", sep=""), sep=PATHSEP)
names(seqFiles) = genomeLtrs
getSeqs_args = sapply(genomeLtrs, function(genome)
    c(getSeqsFromFasta, fastaFiles[genome], "-l", "0", "-i", extractPosFiles[genome], "-o", seqFiles[genome]), simplify=FALSE)

catnow("Retrieving sequences at each LCB or marker for each genome.\n")
for (genome in genomeLtrs)
    {
    inv(extractPosFiles[genome], "extraction positions file")
    inv(seqFiles[genome], "DNA sequences file")
    catnow("  Genome ", genome, "command line:\n")
    catnow("   ", perlPath, " ", paste(getSeqs_args[[genome]], collapse=" "), "\n")
    stat = system2(perlPath, getSeqs_args[[genome]])
    if (stat != 0)
        stop("Perl program ", getSeqsFromFasta, " exited with error status ", stat)
    }

################################################################################
# Read the output files, extract the sequence extraction strings from the
# sequence ID lines and use those strings as the NAMES of the extracted
# sequences.  Then, make sure that the extracted sequences include ALL the
# expected sequences.  Save the sequences in the df data frame, using column
# names "Xseq".
################################################################################

# Read sequences and make sure all are present.
catnow("Processing sequence data\n")
for (genome in genomeLtrs)
    {
    catnow("  Genome ", genome, "\n")
    seqs = readLines(seqFiles[genome])
    dft = as.data.frame(matrix(seqs, ncol=2, byrow=TRUE), stringsAsFactors=FALSE)
    colnames(dft) = c("seqExtStr", "seq")
    dft$seqExtStr = sub("^>([^ ]+) revcomp:no sub_region(:[0-9]+)-([0-9]+).*$", "\\1\\2..\\3", dft$seqExtStr)
    dft$seqExtStr = sub("^>([^ ]+) revcomp:yes sub_region(:[0-9]+)-([0-9]+).*$", "!\\1\\2..\\3", dft$seqExtStr)
    seqs = toupper(dft$seq)
    names(seqs) = dft$seqExtStr
    numMissing = sum(!seqExtStrs[[genome]] %in% dft$seqExtStr)
    if (numMissing > 0)
        stop(numMissing, " sequences are missing from the FASTA extraction of genome ", genome,
            ", check regular expression patterns in sub() calls to assign to dft$seqExtStr")
    seqs = seqs[seqExtStrs[[genome]]]

    # Append the sequences to the df data frame.
    df = data.frame(df, seq=seqs, stringsAsFactors=FALSE)
    colnames(df)[ncol(df)] = paste(genome, "seq", sep="")
    rownames(df) = NULL

    rm(seqs, dft)
    }
rm(seqExtStrs)
seqCols = paste(genomeLtrs, "seq", sep="")
names(seqCols) = genomeLtrs
refSeqCol = seqCols[1]

################################################################################
# Confirm that all sequences start and end with the same kmerLen characters,
# which are the common unique k-mer.
#
# This is where we find out if the Xpos1 and Xpos2 adjustments we made after
# reading the input file were correct.
################################################################################

ref.startKmers = substr(df[, refSeqCol], 1, kmerLen)
lens = nchar(df[, refSeqCol])
ref.endKmers = substr(df[, refSeqCol], lens-kmerLen+1, lens)
startKmerMatches = rep(TRUE, nrow(df))
endKmerMatches = startKmerMatches
for (genome in otherGenomeLtrs)
    {
    seqCol = seqCols[genome]
    g.startKmers = substr(df[, seqCol], 1, kmerLen)
    lens = nchar(df[, seqCol])
    g.endKmers = substr(df[, seqCol], lens-kmerLen+1, lens)
    startKmerMatches = startKmerMatches & (ref.startKmers == g.startKmers)
    endKmerMatches = endKmerMatches & (ref.endKmers == g.endKmers)
    }
if (any(!startKmerMatches) || any(!endKmerMatches))
    stop("Error, expected all sequences to start and end with the same k-mer but ",
        sum(!startKmerMatches), " start seqs and ", sum(!endKmerMatches),
        " end seqs do not match out of ", nrow(df))

################################################################################
# Remove all df rows for which the DNA sequences are identical for all genomes,
# since there will be no Indels or SNPs in such sequences.
################################################################################

allIdentical = rep(TRUE, nrow(df))
for (genome in otherGenomeLtrs)
    allIdentical = allIdentical & (df[, refSeqCol] == df[, seqCols[genome]])
df = df[!allIdentical,]
rownames(df) = NULL

########################################
# Save the output file.
########################################

write.table.winSafe(df, outputFile, row.names=FALSE, quote=FALSE, sep="\t")
# df = read.table(outputFile, header=TRUE, row.names=NULL, sep="\t", stringsAsFactors=FALSE)
catnow("Number of sequence sets for alignment written to output file:", nrow(df), "\n")
catnow("Finished retrieving DNA sequences, output file:\n", outputFile, "\n")
}

################################################################################
# End of file.
################################################################################
