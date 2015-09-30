################################################################################
# See usage below for description.
# Author: Ted Toal
# Date: 2015
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
#testing = 1 # For testing only, outTestHP11 LCRs
#testing = 2 # For testing only, outTestHP11 IndelGroupsNonoverlapping
#testing = 3 # For testing only, outTestHP11 MarkersNonoverlapping
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
    args = c("~/Documents/UCDavis/BradyLab/Genomes/kmers/IGGPIPE",
        "outTestHP11/LCRs_K11k2L100D10_2000.tsv",
        "outTestHP11/LCRs_K11k2L100D10_2000.indels.tsv",
        "outTestHP11/GenomeData",
        "/Users/tedtoal/perl5/perlbrew/perls/perl-5.14.2/bin/perl",
        "code/perl/getSeqsFromFasta.pl",
        aligner, TRUE,
        "testFASTA/ITAG2.4_test.fasta", "testFASTA/SpennV2.0_test.fasta")
    }
else if (testing == 2)
    {
    args = c("~/Documents/UCDavis/BradyLab/Genomes/kmers/IGGPIPE",
        "outTestHP11/IndelGroupsNonoverlapping_K11k2L100D10_2000A100_2000d10_100N2F0.tsv",
        "outTestHP11/IndelGroupsNonoverlapping_K11k2L100D10_2000A100_2000d10_100N2F0.indels.tsv",
        "outTestHP11/GenomeData",
        "/Users/tedtoal/perl5/perlbrew/perls/perl-5.14.2/bin/perl",
        "code/perl/getSeqsFromFasta.pl",
        aligner, TRUE,
        "testFASTA/ITAG2.4_test.fasta", "testFASTA/SpennV2.0_test.fasta")
    }
else if (testing == 3)
    {
    args = c("~/Documents/UCDavis/BradyLab/Genomes/kmers/IGGPIPE",
        "outTestHP11/MarkersNonoverlapping_K11k2L100D10_2000A100_2000d10_100N2F0X20V3000W8M3G1.tsv",
        "outTestHP11/MarkersNonoverlapping_K11k2L100D10_2000A100_2000d10_100N2F0X20V3000W8M3G1.indels.tsv",
        "outTestHP11/GenomeData",
        "/Users/tedtoal/perl5/perlbrew/perls/perl-5.14.2/bin/perl",
        "code/perl/getSeqsFromFasta.pl",
        aligner, TRUE,
        "testFASTA/ITAG2.4_test.fasta", "testFASTA/SpennV2.0_test.fasta")
    }
else stop("Unknown value for 'testing'")
}

NexpectedMin = 10
if (length(args) < NexpectedMin)
    {
    usage = c(
        "Read a data frame of LCRs, Indel Groups, or markers, extract DNA sequences from all",
        "genomes between the start and end position of each one, align the sequences, and",
        "extract from the alignments the size of each InDel that larger than the specified",
        "minimum size.  Write a new file containing the InDel positions, with columns 'ID',",
        "'phases', 'idx', 'Xdel', 'Xid', 'Xstart', and 'Xend', and with one row per InDel.",
        "The 'ID' column contains a unique ID that ties the row back to the input file row(s)",
        "from which it originated.  For LCR input files, the ID is the LCR number.  For InDel",
        "group and marker input files, the ID is the reference genome ID followed by '_', the",
        "reference genome 'pos1' number, an '_', and the reference genome 'pos2' number.",
        "The 'phases' column gives the 'phase' of each genome including the reference genome,",
        "relative to the reference genome, as a string of '+' and '-' characters, where '+'",
        "means the sequences in the two genomes run in the same direction, and '-' means they",
        "run in the opposite direction.",
        "The 'idx' column contains integers that indicate which InDel this is among all InDels",
        "with the same 'ID' value. The 'idx' value will start at 1 and count each InDel within",
        "an 'ID'. In cases with more than two genomes where the alignment shows a complex InDel",
        "with insertions and deletions in various genomes overlapping one another, the entire",
        "region where one or more genomes shows an indel at some base position is counted as",
        "one InDel.  The number of InDels found in the region for a given input row is equal to",
        "the maximum value found in the 'idx' column of all rows with that 'ID'.",
        "The 'Xdel', 'Xid', 'Xstart', and 'Xend' columns give the total number of deleted bps",
        "within the InDel in each genome (Xdel), the sequence ID of the InDel in each genome (Xid)",
        "and the overall InDel starting and ending position in each genome (Xstart and Xend).",
        "Xstart and Xend are the positions of the bps just before the InDel (before the first gap",
        "in the alignment) and just after the InDel (after the last gap in the alignment), so that",
        "Xstart/Xend refer to the same two bps in all genomes.  However, it is always true that",
        "Xstart < Xend, but if the 'phases' value for that genome is '-', then Xstart is the bp",
        "just AFTER the InDel, not just before it, and likewise for Xend, i.e. the InDel region was",
        "reverse-complemented in order to align it to the reference genome sequence (which always",
        "has phase '+' and is never reverse-complemented).  The length of the InDel region in each",
        "genome is therefore Xend-Xstart-1.  Xdel is always the total number of deleted bp within the",
        "InDel in the genome.  With 2 genomes, when Xdel is 0 that is the genome with the insertion",
        "(no gaps), and the length of it is Xend-Xstart-1, which will be equal to Xdel of the genome",
        "with the deletion (whose Xend-Xstart-1 will be 0).  With >2 genomes, Xdel can be non-zero",
        "for all genomes.  A genome has only insertions in the InDel if Xdel is 0, and it has only",
        "deletions if Xend-Xstart-1 = 0, and otherwise it has a mixture of at least one insertion and",
        "one deletion within the InDel interval.",
        "",
        "Usage: Rscript alignAndGetIndels.R <wd> <inputFile> <outputFile> <extractDir> \\",
        "       <perlPath> <getSeqsFromFasta> <alignerPath> <investigate> \\",
        "       <fastaFile1> <fastaFile2> ...",
        "",
        "Arguments:",
        "   <wd> : Path of R working directory, specify other file paths relative to this.",
        "   <inputFile>    : Input file containing LCRs, InDel groups, or markers.",
        "   <outputFile>   : Output file to which to write data frame of InDel information.",
        "   <extractDir>   : Directory to hold DNA sequence extraction files.",
        "   <perlPath>     : Full path of the Perl language interpreter program",
        "   <getSeqsFromFasta> : Full path of the Perl script getSeqsFromFasta.pl",
        ifelse(useClustal,
            "   <alignerPath>  : Full path of the sequence alignment program 'clustalW2'",
            "   <alignerPath>  : Full path of the sequence alignment program 'muscle'"),
        "   <investigate>  : FALSE for normal operation, TRUE for more verbose debugging output.",
        "   <fastaFile1>   : Genome 1 FASTA file.",
        "   <fastaFile2>   : Genome 2 FASTA file.",
        "   ...            : FASTA files of other genomes, for all the genomes."
        )
    for (S in usage)
        catnow(S, "\n", sep="")
    stop("Try again with correct number of arguments")
    }

catnow("alignAndGetIndels.R arguments:\n")
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

extractDir = args[4]
catnow("  extractDir: ", extractDir, "\n")
if (!dir.exists(extractDir))
    stop("Directory doesn't exist: ", extractDir)

perlPath = args[5]
catnow("  perlPath: ", perlPath, "\n")

getSeqsFromFasta = args[6]
catnow("  getSeqsFromFasta: ", getSeqsFromFasta, "\n")
if (!file.exists(getSeqsFromFasta))
    stop("File doesn't exist: ", getSeqsFromFasta)

alignerPath = args[7]
catnow("  alignerPath: ", alignerPath, "\n")

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

# If this is an LCR file (has an "LCR" column), convert it by:
#   (1) dropping "contig" columns
#   (2) removing row names, which should actually be in column "row.names"
#       because we read the file w/o row names
#   (3) converting "X." column name prefix to "X"
#   (4) changing "XseqID" column names to "Xid"
#   (5) adding an "ID" column by concatenating reference ID with LCR number and
#       deleting the LCR column
#   (6) merging rows with the same ID, unsuring that pos1 and pos2 are set
#       correctly and offset correctly so they are PRECISELY the positions of
#       the first bases INSIDE of the two flanking k-mers, and adding column
#       "phases" that has "-" when the direction is reversed relative to the
#       reference strand.
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
    # Rename LCR column to ID column.
    colnames(df) = sub("^LCR$", "ID", colnames(df))

    # Merge: collapse all the identical ID rows together, compute correct Xpos1,
    # Xpos2, and phase.  Add phases column.
    df = df[order(df$ID, df[, idCols[1]], df[, posCols[1]]),] # Make sure df is sorted by reference genome within each LCR.
    dfx = data.frame(ID=unique(df$ID), phases="", row.names=unique(df$ID), stringsAsFactors=FALSE)
    pos1Cols = paste(genomeLtrs, "pos1", sep="")
    names(pos1Cols) = genomeLtrs
    pos2Cols = paste(genomeLtrs, "pos2", sep="")
    names(pos2Cols) = genomeLtrs
    for (genome in genomeLtrs)
        {
        idCol = idCols[genome]
        posCol = posCols[genome]
        strandCol = strandCols[genome]
        refStrandCol = strandCols[1]
        pos1Col = pos1Cols[genome]
        pos2Col = pos2Cols[genome]
        i1 = tapply(1:nrow(df), df$ID, function(ii) return(ii[1]))
        i2 = tapply(1:nrow(df), df$ID, function(ii) return(ii[length(ii)]))
        IDs = df$ID[i1]
        revPhase = (df[i1, strandCol] != df[i1, refStrandCol])
        dfx[IDs, idCol] = df[i1, idCol]
        dfx[IDs, pos1Col] = df[i1, posCol]
        dfx[IDs, pos2Col] = df[i2, posCol]
        dfx[IDs[revPhase], pos1Col] = df[i2[revPhase], posCol]
        dfx[IDs[revPhase], pos2Col] = df[i1[revPhase], posCol]
        dfx[IDs, "phases"] = paste(dfx[IDs, "phases"], c("+", "-")[1+revPhase], sep="")

        # Adjust positions so that instead of being the 5' end of the k-mer, they
        # are the positions of the first bases TO THE INSIDE of the k-mer.
        dfx[, pos1Col] = dfx[, pos1Col] + kmerLen
        dfx[, pos2Col] = dfx[, pos2Col] - 1
        }
    df = dfx
    rm(dfx)
    }

# If this is an IndelGroups file (has an "Xctg1" column), convert it by:
#   (1) removing row names if any
#   (2) retaining only the needed columns: 'Xid', 'Xpos1', 'Xpos2'
#   (3) adding an "ID" column by concatenating reference Xid with reference Xpos1
#       and Xpos2 numbers
#   (4) creating a "phase" column which has "-" if Xpos1 > Xpos2, swapping Xpos1/Xpos2
#       columns when that is the case so Xpos1 < Xpos2 always, and adjusting Xpos1/Xpos2
#       to change the positions from the 5' end of the "+" strand k-mer to being
#       PRECISELY the positions of the first bases to the inside of the k-mer.
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
    for (genome in genomeLtrs)
        {
        isReversed = (df[, pos1Cols[genome]] > df[, pos2Cols[genome]])
        t = df[isReversed, pos1Cols[genome]]
        df[isReversed, pos1Cols[genome]] = df[isReversed, pos2Cols[genome]]
        df[isReversed, pos2Cols[genome]] = t
        phases = rep("+", nrow(df))
        phases[isReversed] = "-"
        df$phases = paste(df$phases, phases, sep="")

        # Adjust positions so that instead of being the 5' end of the k-mer, they
        # are the positions of the first bases TO THE INSIDE of the k-mer.
        df[, pos1Cols[genome]] = df[, pos1Cols[genome]] + kmerLen
        df[, pos2Cols[genome]] = df[, pos2Cols[genome]] - 1
        }
    }

# Else this must be a Markers file, convert it by:
#   (1) removing row names if any
#   (2) converting "XampPos1/2" columns to "Xpos1/2"
#   (3) adding an "ID" column by concatenating reference Xid with reference Xpos1
#       and Xpos2 numbers
#   (4) swapping Xpos1/2 columns if necessary so that Xpos1 < Xpos2
#   (5) if no swap, adding kmer1offset to Xpos1 and subtracting kmer2offset from Xpos2,
#       if swap, adding kmer2offset to Xpos1 and subtracting kmer1offset from Xpos2,
#       placing Xpos1/Xpos2 at the outside ends of the k-mers.
#   (6) adjusting Xpos1/Xpos2 to change the positions from the 5' end of the "+"
#       strand k-mer to being PRECISELY the positions of the first bases to the
#       inside of the k-mer.
#   (7) retaining only the needed columns: 'Xid', 'XYphase', 'Xpos1', 'Xpos2'
#   (8) concatenating columns "XYphase" to form column "phases" and deleting the
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

    # Swap pos1 and pos2 if they are backwards ("-" strand, but we don't care).
    # If no swap, add kmer1offset to Xpos1 and subtracting kmer2offset from Xpos2,
    # If swap, add kmer2offset to Xpos1 and subtracting kmer1offset from Xpos2,
    # placing Xpos1/Xpos2 at the outside ends of the k-mers.  Then adjust Xpos1/Xpos2
    # to change the positions from the 5' end of the "+" strand k-mer to the positions
    # of the first base to the inside of the k-mer.
    for (genome in genomeLtrs)
        {
        isReversed = (df[, pos1Cols[genome]] > df[, pos2Cols[genome]])
        t = df[isReversed, pos1Cols[genome]]
        df[isReversed, pos1Cols[genome]] = df[isReversed, pos2Cols[genome]]
        df[isReversed, pos2Cols[genome]] = t
        df[!isReversed,pos1Cols[genome]] = df[!isReversed,pos1Cols[genome]] + df$kmer1offset[!isReversed]
        df[!isReversed,pos2Cols[genome]] = df[!isReversed,pos2Cols[genome]] - df$kmer2offset[!isReversed]
        df[isReversed,pos1Cols[genome]] = df[isReversed,pos1Cols[genome]] + df$kmer2offset[isReversed]
        df[isReversed,pos2Cols[genome]] = df[isReversed,pos2Cols[genome]] - df$kmer1offset[isReversed]

        # Adjust positions so that instead of being the outside end of the k-mer,
        # they are the positions of the first bases TO THE INSIDE of the k-mer.
        df[, pos1Cols[genome]] = df[, pos1Cols[genome]] + kmerLen
        df[, pos2Cols[genome]] = df[, pos2Cols[genome]] - kmerLen
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

# Put columns in desired order.
df = df[, c("ID", "phases", c(rbind(idCols, pos1Cols, pos2Cols)))]
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
    writeLines(seqExtStrs[[genome]], extractPosFiles[genome])
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
cmdLines = paste(perlPath, getSeqsFromFasta, fastaFiles, "-l 0", "-i", extractPosFiles, "-o", seqFiles)
names(cmdLines) = genomeLtrs

catnow("Retrieving sequences at each LCR or marker for each genome.\n")
for (genome in genomeLtrs)
    {
    inv(cmdLines[genome], "Perl command line")
    inv(extractPosFiles[genome], "extraction positions file")
    inv(seqFiles[genome], "DNA sequences file")
    catnow("  Genome ", genome, "command line:\n")
    catnow("   ", cmdLines[genome], "\n")
    system(cmdLines[genome])
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
    }
seqCols = paste(genomeLtrs, "seq", sep="")

################################################################################
# Now perform and analyze alignments.  For each row of df, write a FASTA file
# containing the sequences for each genome, then invoke the sequence alignment
# program in "alignerPath" to perform a multiple alignment.  Read the alignment
# back, parse it to extract the positions of all InDels, and create new data
# frame dfIndels containing them.
################################################################################

# Doing a for loop here may take a LONG time.
outputType = "FASTA" # ClustalW2.  Muscle default is FASTA.
tempFastaFileName = paste(extractDir, "align.fa", sep=PATHSEP)
tempAlignFileName = paste(extractDir, "align.fasta", sep=PATHSEP)
tempStdoutFileName = paste(extractDir, "align.stdout", sep=PATHSEP) # ClustalW2.  Muscle is quiet.

# Make missing genome column name vectors for dfIndels.
#    ID, phases, idx, Xdel, Xid, Xstart, Xend
delCols = paste(genomeLtrs, "del", sep="")
names(delCols) = genomeLtrs
startCols = paste(genomeLtrs, "start", sep="")
names(startCols) = genomeLtrs
endCols = paste(genomeLtrs, "end", sep="")
names(endCols) = genomeLtrs

# Start the loop.  This may take forever!
dfIndels = NULL
logEveryN = 100
logCount = 0
catnow("Performing", nrow(df), "alignments and extracting InDel positions\n")

# Timing note: with Muscle, and with 459 alignments, total time without anything
# below EXCEPT the Muscle alignment was 1:05 minutes, and with everything below
# was 1:11 minutes.  So most of the time is the alignment time.  ClustalW2 was
# probably even slower.

for (i in 1:nrow(df))
    {
    #catnow("i =", i, "of", nrow(df), "\n")

    # Write sequences to FASTA file.
    idLines = paste(">", genomeLtrs, sep="")
    writeLines(unlist(c(rbind(idLines, df[i, seqCols]))), tempFastaFileName)

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

        cmdLine = paste(alignerPath, " -quiet", " -quicktree", " -gapext=0", " -gapopen=-1",
            " -outorder=INPUT", " -type=DNA", " -output=", outputType, " -infile=",
            tempFastaFileName, " >", tempStdoutFileName, sep="")
        }
    else
        {
        cmdLine = paste(alignerPath, " -quiet", " -in ", tempFastaFileName, " -out ", tempAlignFileName, sep="")
        }
    }
    #inv(cmdLine)
    system(cmdLine)

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
    N = nchar(alignSeqs)[1]
    if (!all(nchar(alignSeqs) == N))
        stop("Expected all sequences to be the same size")
    alignMtx = matrix(unlist(strsplit(alignSeqs, "", fixed=TRUE)), ncol=Ngenomes, dimnames=list(NULL, genomeLtrs))

    # Here is how we will find InDels.  If a base position has no "-" gap characters
    # in any genome, that base position is not part of an InDel, else it is.  The
    # InDel spans all consecutive bases where one or more genomes has a "-".
    gaps = apply(alignMtx, 1, function(V) any(V == "-"))
    indelStarts = which(c(TRUE, !gaps[-N]) & gaps) # Previous bp not a gap and this bp is a gap
    indelEnds = which(c(!gaps[-1], TRUE) & gaps) # Next bp not a gap and this bp is a gap
    Nindels = length(indelStarts)
    if (length(indelEnds) != Nindels)
        stop("Programming error, expected equal number of starts/ends")
    if (Nindels > 0)
        {
        # Each InDel in each genome has a number of gaps ("-" characters) between the
        # InDel start and end position, and we need to count these.
        gapCount = list()
        for (genome in genomeLtrs)
            gapCount[[genome]] = sapply(1:Nindels, function(i) sum(alignMtx[indelStarts[i]:indelEnds[i], genome] == "-"))

        # We use the base BEFORE the InDel as its starting position and the base AFTER the
        # InDel as its ending position.  If the alignment includes gaps at the start or end,
        # the base before or after may be at an offset of 0 or (seq length + 1) from the
        # extraction position start/end (Xpos1/Xpos2).
        indelBaseBeforeStart = indelStarts - 1
        indelBaseAfterEnd = indelEnds + 1

        # Create a data frame of InDels and append it to dfIndels.
        phases = df$phases[i]
        dfi = data.frame(ID=df$ID[i], phases=phases, idx=1:Nindels, stringsAsFactors=FALSE)
        for (i in 1:Ngenomes)
            {
            genome = genomeLtrs[i]
            phase = substring(phases, i, i)
            dfi[, delCols[genome]] = gapCount[[genome]]
            dfi[, idCols[genome]] = df[i, idCols[genome]]
            # For start and end positions, we must SUBTRACT THE NUMBER OF GAPS IN EACH GENOME.
            # This number is held in "gapCount", but we must do it properly.  The first count
            # is not subtracted from the first start position but is subtracted from all
            # subsequent positions.  The second count is not subtracted from the first or second
            # start position of first end position, but is subtracted from the rest.  Etc.
            gapSubtract = cumsum(gapCount[[genome]])
            indelBaseBeforeStart.genome = indelBaseBeforeStart - c(0, gapSubtract[-Nindels])
            indelBaseAfterEnd.genome = indelBaseAfterEnd - gapSubtract
            # Genomes with negative phase require that the start/end positions be
            # adjusted for the fact that the aligned sequence was reverse complemented
            # after being extracted from the genome position given by pos1Cols/pos2Cols.
            if (phase == "+")
                {
                dfi[, startCols[genome]] = df[i, pos1Cols[genome]] + indelBaseBeforeStart.genome - 1
                dfi[, endCols[genome]] = df[i, pos1Cols[genome]] + indelBaseAfterEnd.genome - 1
                }
            else
                {
                dfi[, startCols[genome]] = df[i, pos2Cols[genome]] - indelBaseBeforeStart.genome + 1
                dfi[, endCols[genome]] = df[i, pos2Cols[genome]] - indelBaseAfterEnd.genome + 1
                }
            }
        dfIndels = rbind(dfIndels, dfi)
        }

    # Logging.
    logCount = logCount + 1
    if (logCount %% logEveryN == 0)
        cat(" N =", logCount, "of", nrow(df), "\n")
    }
catnow("Finished alignments.\n")
inv(dim(dfIndels), "dim(dfIndels)")

########################################
# Finish up.
########################################

# Write the data frame to the output file.
write.table(dfIndels, outputFile, row.names=FALSE, quote=FALSE, sep="\t")
# dfIndels = read.table(outputFile, header=TRUE, row.names=NULL, sep="\t", stringsAsFactors=FALSE)

catnow("Finished aligning sequences for Indel Groups and locating InDels, output file:\n", outputFile, "\n")
}

################################################################################
# End of file.
################################################################################
