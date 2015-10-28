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
#thisDir = "~/Documents/UCDavis/BradyLab/Genomes/IGGPIPE/code/R/" # For testing only.

# Source the necessary include files from the same directory containing this file.
source(paste(thisDir, "Include_Common.R", sep=""))

# Get arguments.
testing = 0
#testing = 1 # For testing only, outTestHP11 LCRs
#testing = 2 # For testing only, outTestHP11 IndelGroupsNonoverlapping
#testing = 3 # For testing only, outTestHP11 MarkersNonoverlapping
#testing = 4 # For testing only, outTaCW15 LCRs
{
# Which aligner?  I had problems with ClustalW2, see comments later.
useClustal = FALSE
aligner = ifelse(useClustal,
    "/Users/tedtoal/bin/clustalw-2.1-macosx/clustalw2",
    "/Users/tedtoal/bin/muscle3.8.31_i86darwin64")

if (testing == 0)
    args = commandArgs(TRUE)
else if (testing == 1)
    {
    args = c("~/Documents/UCDavis/BradyLab/Genomes/IGGPIPE",
        "outTestHP11/LCRs_K11k2L100D10_2000.tsv",
        "outTestHP11/LCRs_K11k2L100D10_2000.indels.tsv",
        "outTestHP11/LCRs_K11k2L100D10_2000.snps.tsv",
        "outTestHP11/GenomeData",
        "/Users/tedtoal/perl5/perlbrew/perls/perl-5.14.2/bin/perl",
        "code/perl/getSeqsFromFasta.pl",
        aligner, 0.1, TRUE,
        "testFASTA/ITAG2.4_test.fasta", "testFASTA/SpennV2.0_test.fasta")
    }
else if (testing == 2)
    {
    args = c("~/Documents/UCDavis/BradyLab/Genomes/IGGPIPE",
        "outTestHP11/IndelGroupsNonoverlapping_K11k2L100D10_2000A100_2000d10_100N2F0.tsv",
        "outTestHP11/IndelGroupsNonoverlapping_K11k2L100D10_2000A100_2000d10_100N2F0.indels.tsv",
        "outTestHP11/IndelGroupsNonoverlapping_K11k2L100D10_2000A100_2000d10_100N2F0.snps.tsv",
        "outTestHP11/GenomeData",
        "/Users/tedtoal/perl5/perlbrew/perls/perl-5.14.2/bin/perl",
        "code/perl/getSeqsFromFasta.pl",
        aligner, 0.1, TRUE,
        "testFASTA/ITAG2.4_test.fasta", "testFASTA/SpennV2.0_test.fasta")
    }
else if (testing == 3)
    {
    args = c("~/Documents/UCDavis/BradyLab/Genomes/IGGPIPE",
        "outTestHP11/MarkersNonoverlapping_K11k2L100D10_2000A100_2000d10_100N2F0X20V3000W8M3G1.tsv",
        "outTestHP11/MarkersNonoverlapping_K11k2L100D10_2000A100_2000d10_100N2F0X20V3000W8M3G1.indels.tsv",
        "outTestHP11/MarkersNonoverlapping_K11k2L100D10_2000A100_2000d10_100N2F0X20V3000W8M3G1.snps.tsv",
        "outTestHP11/GenomeData",
        "/Users/tedtoal/perl5/perlbrew/perls/perl-5.14.2/bin/perl",
        "code/perl/getSeqsFromFasta.pl",
        aligner, 0.1, TRUE,
        "testFASTA/ITAG2.4_test.fasta", "testFASTA/SpennV2.0_test.fasta")
    }
else if (testing == 4)
    {
    args = c("~/Documents/UCDavis/BradyLab/Genomes/IGGPIPE",
        "outTaCW15/LCRs_K15k2L200D5_1000.tsv",
        "outTaCW15/LCRs_K15k2L200D5_1000.indels.tsv",
        "outTaCW15/LCRs_K15k2L200D5_1000.snps.tsv",
        "outTaCW15/GenomeData",
        "/Users/tedtoal/perl5/perlbrew/perls/perl-5.14.2/bin/perl",
        "code/perl/getSeqsFromFasta.pl",
        aligner, 0.1, TRUE,
        "/Users/tedtoal/Documents/UCDavis/BradyLab/Genomes/Ta/Genome_IWGSC1.0/Triticum_aestivum.IWGSC1.0+popseq.28.dna.genome.fa",
        "/Users/tedtoal/Documents/UCDavis/BradyLab/Genomes/Ta/Genome_W7984_scaffoldsMar28/w7984.meraculous.scaffolds.Mar28_contamination_removed.fa")
    }
else stop("Unknown value for 'testing'")
}

NexpectedMin = 11
if (length(args) < NexpectedMin)
    {
    usage = c(
        "Read a data frame of LCRs, Indel Groups, or markers, extract DNA sequences from all",
        "genomes between the start and end position of each one, align the sequences, and",
        "extract from the alignments the size of each Indel and SNP, writing them to output",
        "files.  The Indel output file contains the Indel positions, with columns 'ID',",
        "'phases', 'idx', 'Xdel', 'Xid', 'Xstart', and 'Xend', and with one row per Indel.",
        "The 'ID' column contains a unique ID composed of the reference genome ID followed by",
        "'_', the reference genome 'pos1' number, an '_', and the reference genome 'pos2' number",
        "from the input data frame (or in the case of LCRs, of the two k-mer 'pos' values).",
        "The 'phases' column gives the 'phase' of each genome including the reference genome,",
        "relative to the reference genome, as a string of '+' and '-' characters, where '+'",
        "means the sequences in the two genomes run in the same direction, and '-' means they",
        "run in the opposite direction.",
        "The 'idx' column contains integers that indicate which Indel this is among all Indels",
        "with the same 'ID' value. The 'idx' value will start at 1 and count each Indel within",
        "an 'ID'. In cases with more than two genomes where the alignment shows a complex Indel",
        "with insertions and deletions in various genomes overlapping one another, the entire",
        "region where one or more genomes shows an indel at some base position is counted as",
        "one Indel.  The number of Indels found in the region for a given input row is equal to",
        "the maximum value found in the 'idx' column of all rows with that 'ID'.",
        "There is one column of 'Xdel', 'Xid', 'Xstart', and 'Xend' for each genome, with X",
        "replaced by the genome letter.  These columns give the total number of deleted bps",
        "within the Indel in each genome (Xdel), the sequence ID of the Indel in each genome (Xid)",
        "and the overall Indel starting and ending position in each genome (Xstart and Xend).",
        "Xstart and Xend are the positions of the bps just before the Indel (before the first gap",
        "in the alignment) and just after the Indel (after the last gap in the alignment), so that",
        "Xstart/Xend refer to the same two bps in all genomes.  However, it is always true that",
        "Xstart < Xend, but if the 'phases' value for that genome is '-', then Xstart is the bp",
        "just AFTER the Indel, not just before it, and likewise for Xend, i.e. the Indel region was",
        "reverse-complemented in order to align it to the reference genome sequence (which always",
        "has phase '+' and is never reverse-complemented).  The length of the Indel region in each",
        "genome is therefore Xend-Xstart-1.  Xdel is always the total number of deleted bp within the",
        "Indel in the genome.  With 2 genomes, when Xdel is 0 that is the genome with the insertion",
        "(no gaps), and the length of it is Xend-Xstart-1, which will be equal to Xdel of the genome",
        "with the deletion (whose Xend-Xstart-1 will be 0).  With >2 genomes, Xdel can be non-zero",
        "for all genomes.  A genome has only insertions in the Indel if Xdel is 0, and it has only",
        "deletions if Xend-Xstart-1 = 0, and otherwise it has a mixture of at least one insertion and",
        "one deletion within the Indel interval.",
        "",
        "The SNP output file contains the SNP positions and values, with columns 'ID', 'phases',",
        "'idx', 'Xid', 'Xpos', and 'Xval', and with one row per SNP.  The first three of those are",
        "defined exactly the same as with the Indel output file as described above.  The 'Xid' and",
        "'Xpos' columns (one set per genome) give the sequence ID and position of the SNP, and the",
        "'Xval' column gives the SNP value as a single base letter ATCG.",
        "",
        "Usage: Rscript alignAndGetIndelsSNPs.R <wd> <inputFile> <outIndelFile> <outSNPfile> \\",
        "       <extractDir> <perlPath> <getSeqsFromFasta> <alignerPath> <maxSNPsFrac> \\",
        "       <investigate> <fastaFile1> <fastaFile2> ...",
        "",
        "Arguments:",
        "   <wd> : Path of R working directory, specify other file paths relative to this.",
        "   <inputFile>    : Input file containing LCRs, Indel groups, or markers.",
        "   <outIndelFile> : Output file to which to write data frame of Indel information.",
        "   <outSNPfile>   : Output file to which to write data frame of SNP information.",
        "   <extractDir>   : Directory to hold DNA sequence extraction files.",
        "   <perlPath>     : Full path of the Perl language interpreter program",
        "   <getSeqsFromFasta> : Full path of the Perl script getSeqsFromFasta.pl",
        ifelse(useClustal,
            "   <alignerPath>  : Full path of the sequence alignment program 'clustalW2'",
            "   <alignerPath>  : Full path of the sequence alignment program 'muscle'"),
        "   <maxSNPsFrac>  : Maximum allowed fraction of total non-insert sequence that can be SNPs.",
        "                    More than this many SNPs causes that LCR/IndelGroup/Marker to be",
        "                    ignored.  Use 1 to allow any number of SNPs.",
        "   <investigate>  : FALSE for normal operation, TRUE for more verbose debugging output.",
        "   <fastaFile1>   : Genome 1 FASTA file.",
        "   <fastaFile2>   : Genome 2 FASTA file.",
        "   ...            : FASTA files of other genomes, for all the genomes."
        )
    for (S in usage)
        catnow(S, "\n", sep="")
    stop("Try again with correct number of arguments")
    }

catnow("alignAndGetIndelsSNPs.R arguments:\n")
workingDirectory = args[1]
catnow("  workingDirectory: ", workingDirectory, "\n")
if (!dir.exists(workingDirectory))
    stop("Directory doesn't exist: ", workingDirectory)
setwd(workingDirectory)

inputFile = args[2]
catnow("  inputFile: ", inputFile, "\n")
if (!file.exists(inputFile))
    stop("File doesn't exist: ", inputFile)

outIndelFile = args[3]
catnow("  outIndelFile: ", outIndelFile, "\n")

outSNPfile = args[4]
catnow("  outSNPfile: ", outSNPfile, "\n")

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

alignerPath = args[8]
catnow("  alignerPath: ", alignerPath, "\n")

maxSNPsFrac = as.numeric(args[9])
catnow("  maxSNPsFrac: ", maxSNPsFrac, "\n")
if (is.na(maxSNPsFrac) || maxSNPsFrac < 0 || maxSNPsFrac > 1)
    stop("maxSNPsFrac must be between 0 and 1")

investigate = as.logical(args[10])
catnow("  investigate: ", investigate, "\n")
if (is.na(investigate))
    stop("investigate must be TRUE or FALSE")

fastaFiles = args[-(1:10)]
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
# things depending on whether the input file was LCRs, Indel Groups, or markers.
# Set the Xpos1 and Xpos2 values to be the exact start and end positions of the
# outside k-mers, so a DNA extraction will start and end with the k-mers.  For
# LCRs, each LCR may become one or more output data frame rows.  (The entire LCR
# could be one row, but the size can be too large for alignment, so it is broken
# up into smaller segments at k-mer boundaries).

# If this is an LCR file (has an "LCR" column), convert it by:
#   (1) drop "contig" columns
#   (2) remove row names, which should actually be in column "row.names"
#       because we read the file w/o row names
#   (3) convert "X." column name prefix to "X"
#   (4) change "XseqID" column names to "Xid"
#   (5) add column "phases" that has "-" when the direction is reversed relative
#       to the reference strand.
#   (6) convert each LCR into one or more output rows, trying to keep length
#       no more than maxOutLen, and ensuring that pos1 and pos2 are set correctly
#       and offset correctly so they are PRECISELY the positions of the first
#       bases on the outside side of each k-mer that defines that output row.
{
if (any(colnames(df) == "LCR"))
    {
    # The k-mer length we are working with.
    kmerLen = nchar(df[1,"row.names"])
    inv(kmerLen, "k")
    # Drop "contig" columns.
    df = df[, !grepl("contig", colnames(df))]
    # Get rid of row names.
    rownames(df) = NULL
    df = df[, colnames(df) != "row.names"]
    # Remove "." from column names.
    colnames(df) = sub("\\.", "", colnames(df))
    # Change XseqID to Xid.
    colnames(df) = sub("seqID", "id", colnames(df))
    # Get vectors of names of X columns.
    idCols = colnames(df)[grepl("id$", colnames(df))]
    posCols = colnames(df)[grepl("pos$", colnames(df))]
    strandCols = colnames(df)[grepl("strand$", colnames(df))]
    genomeLtrs = substring(idCols, 1, 1)
    otherGenomeLtrs = genomeLtrs[-1]
    if (length(genomeLtrs) < 2)
        stop("Expected at least two id columns in <inputFile>")
    inv(genomeLtrs, "genome letters")
    names(idCols) = genomeLtrs
    names(posCols) = genomeLtrs
    names(strandCols) = genomeLtrs
    refIdCol = idCols[1]
    refPosCol = posCols[1]
    refStrandCol = strandCols[1]
    # Add column "phases".
    df$phases = ""
    for (genome in genomeLtrs)
        df$phases = paste(df$phases, c("+", "-")[1+(df[,strandCols[genome]] != df[,refStrandCol])], sep="")
    # Remove strand columns, don't need them any more.
    df = df[, !colnames(df) %in% strandCols]

    # Break the LCR up into one or more output rows and build the rows.  A single
    # LCR can be so long that we can't align it.  Keep the maximum length of the
    # output rows short enough that we will be able to align them.
    maxOutLen = 500 # Try to keep output length less than this.
    df = df[order(df$LCR, df[, refIdCol], df[, refPosCol]),] # Make sure df is sorted by reference genome within each LCR.
    # Get indexes indicating starting and ending row of each LCR.  Note that it
    # is always true that length(startLCR) = length(endLCR)
    N = nrow(df)
    startLCR = which(c(TRUE, df$LCR[-1] != df$LCR[-N]))
    endLCR = which(c(df$LCR[-1] != df$LCR[-N], TRUE))
    # Accumulate df row indexes indicating start and end k-mers of each LCR span
    # that is to become an output row.  Note that it is always true that
    # length(startRows) = length(endRows)
    startRows = c()
    endRows = c()
    # curStarts are df row indexes of rows that will be the start k-mers of the
    # next output data frame rows, initially the first df row of each LCR, and
    # subsequently the same row that was the end k-mer of the previous output row.
    # (The output rows of one LCR will overlap by one k-mer).
    curStarts = startLCR
    # curEnds are df row indexes of the rows currently being tested to see if they
    # should be end-of-span k-mers, initially the SECOND df row of each LCR.
    # Note that it is always true that length(curStarts) = length(curEnds)
    # It will also always be maintained true that all(curEnds-curStarts >= 1)
    curEnds = curStarts+1
    # Loop moving curEnds along the df rows and watching for end of LCR or
    # maxOutLen being exceeded.  Record a span and start a new span when
    # maxOutLen is exceeded.  Record the final span of an LCR when the end of
    # the LCR is reached.
    while (length(curEnds) > 0)
        {
        # If curEnds is the last row of the LCR, move those curStarts/curEnds to
        # startRows/endRows.
        isEndLCR = (curEnds %in% endLCR)
        startRows = c(startRows, curStarts[isEndLCR])
        endRows = c(endRows, curEnds[isEndLCR])
        # Remove them from curStarts/curEnds.
        curStarts = curStarts[!isEndLCR]
        curEnds = curEnds[!isEndLCR]
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
        # ends with.
        curStarts = c(curStarts[!nextExceedsMax], curEnds[nextExceedsMax])
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
    dfx = df[startRows, colnames(df) != "LCR"]
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
    # The k-mer length we are working with.
    kmerLen = nchar(df$kmer1[1])
    inv(kmerLen, "k")
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
    # The k-mer length we are working with.
    kmerLen = nchar(df$kmer1[1])
    inv(kmerLen, "k")
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
catnow("Number of LCRs/Indel Groups/Markers to be searched for Indels/SNPs:", nrow(df), "\n")
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

catnow("Retrieving sequences at each LCR or marker for each genome.\n")
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
# since there will be no Indels or SNPs in such sequences.  There will be many
# of these if the input file was an LCR file, none if it was an Indel Group or
# marker file.
################################################################################

allIdentical = rep(TRUE, nrow(df))
for (genome in otherGenomeLtrs)
    allIdentical = allIdentical & (df[, refSeqCol] == df[, seqCols[genome]])
df = df[!allIdentical,]

################################################################################
# Now perform and analyze alignments.  For each row of df, write a FASTA file
# containing the sequences for each genome, then invoke the sequence alignment
# program in "alignerPath" to perform a multiple alignment.  Read the alignment
# back, parse it to extract the positions of all Indels and SNPs, and create new
# data frames dfIndels and dfSNPs containing them.
################################################################################

# Doing a for loop here may take a LONG time.
outputType = "FASTA" # ClustalW2.  Muscle default is FASTA.
tempFastaFileName = paste(extractDir, "align.fa", sep=PATHSEP)
tempAlignFileName = paste(extractDir, "align.fasta", sep=PATHSEP)
tempStdoutFileName = ifelse(useClustal, paste(extractDir, "align.stdout", sep=PATHSEP), "") # ClustalW2.  Muscle is quiet.

# Output data frames, initially empty.
dfIndels = NULL
dfSNPs = NULL

# Make genome column name vectors for dfIndels and dfSNPs for those X columns we
# don't already have vectors for.
#    ID, phases, idx, Xdel, Xid, Xstart, Xend
#    ID, phases, idx, Xid, Xpos, Xval
delCols = paste(genomeLtrs, "del", sep="")
names(delCols) = genomeLtrs
startCols = paste(genomeLtrs, "start", sep="")
names(startCols) = genomeLtrs
endCols = paste(genomeLtrs, "end", sep="")
names(endCols) = genomeLtrs
posCols = paste(genomeLtrs, "pos", sep="")
names(posCols) = genomeLtrs
valCols = paste(genomeLtrs, "val", sep="")
names(valCols) = genomeLtrs

# Count number of df rows in which there were too many SNPs in the alignment (> maxSNPsFrac).
tooManySNPs = 0

# Prepare for assembling sequence lines into FASTA file.
idLines = paste(">", genomeLtrs, sep="")
names(idLines) = genomeLtrs
maxLineChars = 80

# Open the output files and write results to them as we go along.  Use "wb" so
# it will write Unix line ends even under Windows.
indelsFile = file(outIndelFile, "wb")
SNPsFile = file(outSNPfile, "wb")
writeColNamesIndels = TRUE
writeColNamesSNPs = TRUE

# Write data to output file after every outEveryN lines of df have been aligned.
outEveryN = 100
catnow("Performing", nrow(df), "alignments and extracting Indels and SNPs\n")

# Do garbage collection every gcEveryN, which must be multiple of outEveryN.
gcEveryN = 10000

# Timing note: with Muscle, and with 459 alignments, total time without anything
# below EXCEPT the Muscle alignment was 1:05 minutes, and with everything below
# was 1:11 minutes.  So most of the time is the alignment time.  ClustalW2 was
# probably even slower.

# Start the loop.  This may take forever!
Nindels = 0
Nsnps = 0
NalignerFails = 0
NalignerGarbage = 0
for (alignCount in 1:nrow(df))
    {
    #catnow("alignCount =", alignCount, "of", nrow(df), "\n")
    d = df[alignCount,]

    # Write sequences to FASTA file.  Split lines so each is no longer than
    # maxLineChars bases, both for convenience examining and because aligner
    # might have trouble with very long lines.
    outLines = c()
    for (genome in genomeLtrs)
        {
        sequence = d[, seqCols[genome]]
        Len = nchar(sequence)
        startChar = seq(1, Len, by=maxLineChars)
        endChar = startChar + maxLineChars - 1
        endChar[length(endChar)] = Len
        sequence = substring(sequence, startChar, endChar)
        outLines = c(outLines, idLines[genome], sequence)
        }
    writeLines.winSafe(outLines, tempFastaFileName)
    rm(sequence, startChar, endChar, outLines)

    # Create and execute a command to align the sequences.
    {
    if (useClustal)
        {
        # Originally the extracted sequence included the two common unique k-mers.
        # Originally I tried to use ClustalW2.  A gap extension penalty of 0 should
        # ensure that the alignment will align the two common unique k-mers at either end
        # with one another, EXCEPT if the gap open penalty is more than the penalty of a
        # possible single bp match in the k-mer.  So, set the gap open penalty to something
        # less than the penalty for any base mismatch.  What are those?  The default DNA
        # weight matrix for ClustalW2 is IUB:
        #   X's and N's are treated as matches to any IUB ambiguity symbol. All matches
        #   score 1.9; all mismatches for IUB symbols score 0.
        # Therefore, we should set the gap open penalty to no less than -1.9, say -1.
        # This did not work.  Even with a POSITIVE value for the gap extension penalty,
        # ClustalW2 did not align the ends of a particular sequence when it found a 1bp
        # mismatch farther down.  That's why I switched to Muscle.  Most of the time it
        # DID align the two k-mers, but then I encountered a case where it didn't, I
        # don't know why.  So then I switched to NOT include the flanking k-mers in the
        # alignment, and adjusted the code below to allow for gaps at the start or end
        # of the alignment.  I continued using Muscle, however.

        alignerArgs = c("-quiet", "-quicktree", "-gapext=0", "-gapopen=-1",
            "-outorder=INPUT", "-type=DNA", paste("-output=", outputType, sep=""),
            paste("-infile=", tempFastaFileName, sep=""))
        }
    else
        {
        alignerArgs = c("-quiet", "-in", tempFastaFileName, "-out", tempAlignFileName)
        }
    }
    #inv(paste(alignerArgs, collapse=" "), "alignerArgs")
    stat = system2(alignerPath, alignerArgs, stdout=tempStdoutFileName)
    if (stat != 0)
        {
        NalignerFails = NalignerFails + 1
        next
        }

    # Read the alignment output file and split the lines by genome, paste
    # together the sequence lines into one long character vector.
    align = readLines(tempAlignFileName)
    genomeIdLines = unlist(sapply(idLines, function(idLine) which(align == idLine)))
    if (length(genomeIdLines) != length(idLines))
        stop("Expected all ", length(idLines), " ID lines in alignment output")
    names(genomeIdLines) = genomeLtrs
    genomeIdLines = sort(genomeIdLines)
    firstLines = genomeIdLines + 1
    lastLines = c(firstLines[-1]-2, length(align))
    names(lastLines) = names(firstLines)
    genomeIdLines = genomeIdLines[genomeLtrs]
    firstLines = firstLines[genomeLtrs]
    lastLines = lastLines[genomeLtrs]
    alignSeqs = sapply(1:length(firstLines), function(i) paste(align[firstLines[i]:lastLines[i]], collapse=""))
    if (any(grepl("[^-A-Z]", alignSeqs)))
        {
        NalignerGarbage = NalignerGarbage + 1
        next
        }
    N = nchar(alignSeqs)[1]
    if (!all(nchar(alignSeqs) == N))
        stop("Expected all sequences to be the same size")
    alignMtx = matrix(unlist(strsplit(alignSeqs, "", fixed=TRUE)), ncol=Ngenomes, dimnames=list(NULL, genomeLtrs))

    # Check to see if the number of SNPs exceeds maxSNPsFrac in any genome.  If
    # so, skip this one.  We calculate seqLens as the total number of bases in
    # that genome that have non-"-" for both it and the reference genome.
    refV = alignMtx[, 1, drop=TRUE]
    seqLens = apply(alignMtx[(refV != "-"), -1, drop=FALSE], 2, function(V) sum(V != "-"))
    maxSNPsAllowed = max(as.integer(seqLens * maxSNPsFrac))
    isSNP = apply(alignMtx[, -1, drop=FALSE], 2, function(V) (V != "-" & refV != "-" & V != refV))
    if (ncol(isSNP) == 1)
        isSNP = isSNP[,1] # Converts to vector.
    else
        isSNP = apply(isSNP, 1, any)
    numSNPs = sum(isSNP)
    if (numSNPs > maxSNPsAllowed)
        {
        tooManySNPs = tooManySNPs+1
        next
        }

    # Create a data frame of SNPs and append it to dfSNPs.
    if (numSNPs > 0)
        {
        phases = d$phases
        dfs = data.frame(ID=d$ID, phases=phases, idx=1:numSNPs, stringsAsFactors=FALSE)

        # We will need the position offset of each SNP in each genome.  This offset ignores
        # the gaps in the alignment.
        offsetMtx = apply(alignMtx, 2, function(V) (cumsum(V != "-") - 1)) # Subtract 1 so first offset is 0.

        # Get SNP columns for each genome.
        for (j in 1:Ngenomes)
            {
            genome = genomeLtrs[j]
            phase = substring(phases, j, j)
            dfs[, idCols[genome]] = d[, idCols[genome]]

            # For positions, we must add the offset in that genome to the starting position in it.
            SNPoffset.genome = offsetMtx[isSNP, j]

            # Genomes with negative phase require that the position account for
            # the fact that the aligned sequence was reverse complemented after
            # being extracted from the genome position given by pos1Cols/pos2Cols.
            if (phase == "+")
                dfs[, posCols[genome]] = d[, pos1Cols[genome]] + SNPoffset.genome
            else
                dfs[, posCols[genome]] = d[, pos2Cols[genome]] - SNPoffset.genome

            # And finally, the SNP value.
            dfs[, valCols[genome]] = alignMtx[isSNP, j]
            }
        dfSNPs = rbind.fast(dfSNPs, dfs)
        }

    # Here is how we will find Indels.  If a base position has no "-" gap characters
    # in any genome, that base position is not part of an Indel, else it is.  The
    # Indel spans all consecutive bases where one or more genomes has a "-".
    gaps = apply(alignMtx, 1, function(V) any(V == "-"))
    indelStarts = which(c(TRUE, !gaps[-N]) & gaps) # Previous bp not a gap and this bp is a gap
    indelEnds = which(c(!gaps[-1], TRUE) & gaps) # Next bp not a gap and this bp is a gap
    numIndels = length(indelStarts)
    if (length(indelEnds) != numIndels)
        stop("Programming error, expected equal number of starts/ends")
    if (numIndels > 0)
        {
        # Each Indel in each genome has a number of gaps ("-" characters) between the
        # Indel start and end position, and we need to count these.
        gapCount = list()
        for (genome in genomeLtrs)
            gapCount[[genome]] = sapply(1:numIndels, function(i) sum(alignMtx[indelStarts[i]:indelEnds[i], genome] == "-"))

        # We use the base BEFORE the Indel as its starting position and the base AFTER the
        # Indel as its ending position.  If the alignment includes gaps at the start or end,
        # the base before or after may be at an offset of 0 or (seq length + 1) from the
        # extraction position start/end (Xpos1/Xpos2).
        indelBaseBeforeStart = indelStarts - 1
        indelBaseAfterEnd = indelEnds + 1

        # Create a data frame of Indels and append it to dfIndels.
        phases = d$phases
        dfi = data.frame(ID=d$ID, phases=phases, idx=1:numIndels, stringsAsFactors=FALSE)

        # Get Indel columns for each genome.
        for (j in 1:Ngenomes)
            {
            genome = genomeLtrs[j]
            phase = substring(phases, j, j)
            dfi[, delCols[genome]] = gapCount[[genome]]
            dfi[, idCols[genome]] = d[, idCols[genome]]

            # For start and end positions, we must SUBTRACT THE NUMBER OF GAPS IN EACH GENOME.
            # This number is held in "gapCount", but we must do it properly.  The first count
            # is not subtracted from the first start position but is subtracted from all
            # subsequent positions.  The second count is not subtracted from the first or second
            # start position of first end position, but is subtracted from the rest.  Etc.
            gapSubtract = cumsum(gapCount[[genome]])
            indelBaseBeforeStart.genome = indelBaseBeforeStart - c(0, gapSubtract[-numIndels])
            indelBaseAfterEnd.genome = indelBaseAfterEnd - gapSubtract

            # Genomes with negative phase require that the start/end positions be
            # adjusted for the fact that the aligned sequence was reverse complemented
            # after being extracted from the genome position given by pos1Cols/pos2Cols.
            if (phase == "+")
                {
                dfi[, startCols[genome]] = d[, pos1Cols[genome]] + indelBaseBeforeStart.genome - 1
                dfi[, endCols[genome]] = d[, pos1Cols[genome]] + indelBaseAfterEnd.genome - 1
                }
            else
                {
                dfi[, startCols[genome]] = d[, pos2Cols[genome]] - indelBaseBeforeStart.genome + 1
                dfi[, endCols[genome]] = d[, pos2Cols[genome]] - indelBaseAfterEnd.genome + 1
                }
            }
        dfIndels = rbind.fast(dfIndels, dfi)
        }

    # Output the data.
    if (alignCount %% outEveryN == 0)
        {
        cat(" N =", alignCount, "of", nrow(df))
        dfIndels = rbind.fast.finish(dfIndels)
        if (!is.null(dfIndels) && nrow(dfIndels) > 0)
            {
            Nindels = Nindels + nrow(dfIndels)
            write.table(dfIndels, indelsFile, row.names=FALSE, col.names=writeColNamesIndels, quote=FALSE, sep="\t")
            writeColNamesIndels = FALSE
            dfIndels = NULL
            }
        dfSNPs = rbind.fast.finish(dfSNPs)
        if (!is.null(dfSNPs) && nrow(dfSNPs) > 0)
            {
            Nsnps = Nsnps + nrow(dfSNPs)
            write.table(dfSNPs, SNPsFile, row.names=FALSE, col.names=writeColNamesSNPs, quote=FALSE, sep="\t")
            writeColNamesSNPs = FALSE
            dfSNPs = NULL
            }
        if (alignCount %% gcEveryN == 0)
            {
            cat(" GC...")
            gc()
            }
        cat("\n")
        }
    }

# Output the remaining data if any.
dfIndels = rbind.fast.finish(dfIndels)
if (!is.null(dfIndels) && nrow(dfIndels) > 0)
    {
    Nindels = Nindels + nrow(dfIndels)
    write.table(dfIndels, indelsFile, row.names=FALSE, col.names=writeColNamesIndels, quote=FALSE, sep="\t")
    }
dfSNPs = rbind.fast.finish(dfSNPs)
if (!is.null(dfSNPs) && nrow(dfSNPs) > 0)
    {
    Nsnps = Nsnps + nrow(dfSNPs)
    write.table(dfSNPs, SNPsFile, row.names=FALSE, col.names=writeColNamesSNPs, quote=FALSE, sep="\t")
    }
close(indelsFile)
close(SNPsFile)

# Display statistics.
catnow("Finished aligning sequences for Indel Groups and locating Indels and SNPs\n")
catnow("Total number of Indels is ", Nindels, "\n")
catnow("Total number of SNPs is ", Nsnps, "\n")
catnow("Number of LCR/IndelGroup/Marker alignments that had too many SNPs and were ignored:", tooManySNPs, "\n")
catnow("Aligner program failed ", NalignerFails, " times\n")
catnow("Aligner program output not recognizable ", NalignerGarbage, " times\n")
catnow("Indel output file:\n", outIndelFile, "\n")
catnow("SNP output file:\n", outSNPfile, "\n")
}

################################################################################
# End of file.
################################################################################
