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
#testing = 5 # For testing only, outHP14 IndelGroupsNonoverlapping
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
    args = c("~/Documents/UCDavis/BradyLab/IGGPIPE/IGGPIPE",
        "outTestHP11/MarkersNonoverlapping_K11k2L100D10_2000A100_2000d10_100N2F0X20V3000W8M3G1.withseqs.tsv",
        "outTestHP11/MarkersNonoverlapping_K11k2L100D10_2000A100_2000d10_100N2F0X20V3000W8M3G1.indels.tsv",
        "outTestHP11/MarkersNonoverlapping_K11k2L100D10_2000A100_2000d10_100N2F0X20V3000W8M3G1.snps.tsv",
        aligner, 10, 30, FALSE, TRUE)
    }
else if (testing == 2)
    {
    args = c("~/Documents/UCDavis/BradyLab/IGGPIPE/IGGPIPE",
        "outTestHP11/IndelGroupsNonoverlapping_K11k2L100D10_2000A100_2000d10_100N2F0.withseqs.tsv",
        "outTestHP11/IndelGroupsNonoverlapping_K11k2L100D10_2000A100_2000d10_100N2F0.indels.tsv",
        "outTestHP11/IndelGroupsNonoverlapping_K11k2L100D10_2000A100_2000d10_100N2F0.snps.tsv",
        aligner, 10, 30, FALSE, TRUE)
    }
else if (testing == 3)
    {
    args = c("~/Documents/UCDavis/BradyLab/IGGPIPE/IGGPIPE",
        "outTestHP11/MarkersNonoverlapping_K11k2L100D10_2000A100_2000d10_100N2F0X20V3000W8M3G1.withseqs.tsv",
        "outTestHP11/MarkersNonoverlapping_K11k2L100D10_2000A100_2000d10_100N2F0X20V3000W8M3G1.indels.tsv",
        "outTestHP11/MarkersNonoverlapping_K11k2L100D10_2000A100_2000d10_100N2F0X20V3000W8M3G1.snps.tsv",
        aligner, 10, 30, FALSE, TRUE)
    }
else if (testing == 4)
    {
    args = c("~/Documents/UCDavis/BradyLab/IGGPIPE/IGGPIPE",
        "outTaCW15/LCBs_K15k2L200D5_1000.withseqs.tsv",
        "outTaCW15/LCBs_K15k2L200D5_1000.indels.tsv",
        "outTaCW15/LCBs_K15k2L200D5_1000.snps.tsv",
        aligner, 10, 30, FALSE, TRUE)
    }
else if (testing == 5)
    {
    args = c("~/Documents/UCDavis/BradyLab/IGGPIPE/IGGPIPE",
        "outHP14/IndelGroupsNonoverlapping_K14k4L100D1_3000A100_3000d100_100N2F0.withseqs.tsv",
        "outHP14/IndelGroupsNonoverlapping_K14k4L100D1_3000A100_3000d100_100N2F0.indels.tsv",
        "outHP14/IndelGroupsNonoverlapping_K14k4L100D1_3000A100_3000d100_100N2F0.snps.tsv",
        aligner, 20, 120, FALSE, TRUE)
    }
else stop("Unknown value for 'testing'")
}

NexpectedMin = 9
if (length(args) < NexpectedMin)
    {
    usage = c(
        "Read a data frame of LCBs, Indel Groups, or markers with DNA sequences for all",
        "genomes between the start and end position of each one, align the sequences, and",
        "extract from the alignments the size of each Indel and SNP, writing them to output",
        "files.  The Indel output file contains the Indel positions, with columns 'ID',",
        "'phases', 'idx', 'Xdel', 'Xid', 'Xstart', and 'Xend', and with one row per Indel.",
        "The 'ID' column contains a unique ID composed of the reference genome ID followed by",
        "'_', the reference genome 'pos1' number, an '_', and the reference genome 'pos2' number",
        "from the input data frame (or in the case of LCBs, of the two k-mer 'pos' values).",
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
        "       <alignerPath> <maxSNPsPerKbp> <maxIndelsPerKbp> <scramble> <investigate>",
        "",
        "Arguments:",
        "   <wd>              : Path of R working directory, specify other file paths relative to this.",
        "   <inputFile>       : Input file containing LCBs, Indel groups, or markers, with DNA sequences.",
        "   <outIndelFile>    : Output file to which to write data frame of Indel information.",
        "   <outSNPfile>      : Output file to which to write data frame of SNP information.",
        ifelse(useClustal,
            "   <alignerPath>     : Full path of the sequence alignment program 'clustalW2'",
            "   <alignerPath>     : Full path of the sequence alignment program 'muscle'"),
        "   <maxIndelsPerKbp> : Maximum allowed Indels per Kbp of reference genome sequence.",
        "                       More than this many Indels causes that LCB/IndelGroup/Marker",
        "                       to be ignored.  Typical minimum of random sequence is 10.",
        "   <maxSNPsPerKbp>   : Maximum allowed SNPs per Kbp of reference genome sequence.",
        "                       More than this many SNPs causes that LCB/IndelGroup/Marker",
        "                       to be ignored.  Typical minimum of random sequence is 30.",
        "   <scramble>        : FALSE for normal operation, TRUE to scramble non-ref seqs before alignment.",
        "   <investigate>     : FALSE for normal operation, TRUE for more verbose debugging output."
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

alignerPath = args[5]
catnow("  alignerPath: ", alignerPath, "\n")

maxIndelsPerKbp = as.numeric(args[6])
catnow("  maxIndelsPerKbp: ", maxIndelsPerKbp, "\n")
if (is.na(maxIndelsPerKbp) || maxIndelsPerKbp < 1 || maxIndelsPerKbp > 500)
    stop("maxIndelsPerKbp must be between 1 and 500")

maxSNPsPerKbp = as.numeric(args[7])
catnow("  maxSNPsPerKbp: ", maxSNPsPerKbp, "\n")
if (is.na(maxSNPsPerKbp) || maxSNPsPerKbp < 1 || maxSNPsPerKbp > 1000)
    stop("maxSNPsPerKbp must be between 1 and 1000")

scramble = as.logical(args[8])
catnow("  scramble: ", scramble, "\n")
if (is.na(scramble))
    stop("scramble must be TRUE or FALSE")

investigate = as.logical(args[9])
catnow("  investigate: ", investigate, "\n")
if (is.na(investigate))
    stop("investigate must be TRUE or FALSE")

########################################
# Initialize.
########################################

catnow("Reading data to be aligned, including DNA sequences, from input file\n")

# Read input data.
df = read.table(inputFile, header=TRUE, row.names=NULL, sep="\t", stringsAsFactors=FALSE)
if (nrow(df) == 0)
    stop("There is no input data.")
inv(dim(df), "input data dim")
inv(colnames(df), "input data columns")
inv(head(df), "input data head")

# Get genome letters and number of genomes.
genomeLtrs = sub("id$", "", colnames(df)[grepl("id$", colnames(df))])
Ngenomes = length(genomeLtrs)
if (Ngenomes < 2)
    stop("Expected at least two genomes")

# Make genome column name vectors for df.
idCols = paste(genomeLtrs, "id", sep="")
names(idCols) = genomeLtrs
pos1Cols = paste(genomeLtrs, "pos1", sep="")
names(pos1Cols) = genomeLtrs
pos2Cols = paste(genomeLtrs, "pos2", sep="")
names(pos2Cols) = genomeLtrs
seqCols = paste(genomeLtrs, "seq", sep="")
names(seqCols) = genomeLtrs

catnow("Number of LCBs/Indel Groups/Markers to be searched for Indels/SNPs:", nrow(df), "\n")
if (nrow(df) == 0)
    stop("No regions to be processed.")

################################################################################
# Now perform and analyze alignments.  For each row of df, write a FASTA file
# containing the sequences for each genome, then invoke the sequence alignment
# program in "alignerPath" to perform a multiple alignment.  Read the alignment
# back, parse it to extract the positions of all Indels and SNPs, and create new
# data frames dfIndels and dfSNPs containing them.
################################################################################

# Minimum reference sequence length in bp for gathering stats.
minRefLenStats = 100

# Doing a for loop here may take a LONG time.
outputType = "FASTA" # ClustalW2.  Muscle default is FASTA.
tempFastaFileName = "align.fa"
tempAlignFileName = "align.fasta"
tempStdoutFileName = ifelse(useClustal, "align.stdout", "") # ClustalW2.  Muscle is quiet.

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

# Count number of df rows in which there were too many Indels in the alignment (> maxIndelsPerKbp).
tooManyIndels = 0

# Count number of df rows in which there were too many SNPs in the alignment (> maxSNPsPerKbp).
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

# Keep track of min, max, mean, and standard deviation of number of indels and
# number of SNPs per Kbp of reference sequence, using these vectors.
indelsPerKbp = numeric()
snpsPerKbp = numeric()

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
        if (scramble)
            {
            x = unlist(strsplit(sequence, ""), use.names=FALSE)
            sequence = paste(x[sample.int(Len)], collapse="")
            }
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

    # Check for unexpected characters or wrong-size sequences.
    if (any(grepl("[^-A-Z]", alignSeqs)))
        {
        NalignerGarbage = NalignerGarbage + 1
        next
        }
    N = nchar(alignSeqs)[1]
    if (!all(nchar(alignSeqs) == N))
        stop("Expected all sequences to be the same size")

    # Put all characters of alignment into a matrix.
    alignMtx = matrix(unlist(strsplit(alignSeqs, "", fixed=TRUE)), ncol=Ngenomes, dimnames=list(NULL, genomeLtrs))

    # Here is how we will find Indels.  If a base position has no "-" gap characters
    # in any genome, that base position is not part of an Indel, else it is.  The
    # Indel spans all consecutive bases where one or more genomes has a "-".
    gaps = apply(alignMtx, 1, function(V) any(V == "-"))
    indelStarts = which(c(TRUE, !gaps[-N]) & gaps) # Previous bp not a gap and this bp is a gap
    indelEnds = which(c(!gaps[-1], TRUE) & gaps) # Next bp not a gap and this bp is a gap
    numIndels = length(indelStarts)
    if (length(indelEnds) != numIndels)
        stop("Programming error, expected equal number of starts/ends")

    # Here is how we will find SNPs.  Compare each alignment position in each genome
    # to the reference genome value at that position, and if they are not equal and
    # neither is "-", that is a SNP.
    refV = alignMtx[, 1, drop=TRUE]
    isSNP = apply(alignMtx[, -1, drop=FALSE], 2, function(V) (V != "-" & refV != "-" & V != refV))
    if (ncol(isSNP) == 1)
        isSNP = isSNP[,1] # Converts to vector.
    else
        isSNP = apply(isSNP, 1, any)
    numSNPs = sum(isSNP)

    # Check to see if the number of Indels exceeds maxIndelsPerKbp or the number
    # of SNPs exceeds maxSNPsPerKbp.  If so, skip this one.  We calculate seqLen as the total
    # number of bases in reference genome that have non-"-".
    seqLen = sum(refV != "-")
    cur.IndelsPerKbp = 1000*numIndels/seqLen
    cur.SNPsPerKbp = 1000*numSNPs/seqLen
    if (cur.IndelsPerKbp > maxIndelsPerKbp || cur.SNPsPerKbp > maxSNPsPerKbp)
        {
        if (cur.IndelsPerKbp > maxIndelsPerKbp)
            tooManyIndels = tooManyIndels+1
        if (cur.SNPsPerKbp > maxSNPsPerKbp)
            tooManySNPs = tooManySNPs+1
        next
        }

    # Create a data frame of Indels and append it to dfIndels.
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
                dfi[, startCols[genome]] = d[, pos2Cols[genome]] - indelBaseAfterEnd.genome + 1
                dfi[, endCols[genome]] = d[, pos2Cols[genome]] - indelBaseBeforeStart.genome + 1
                }
            if (any(dfi[, startCols[genome]] >= dfi[, endCols[genome]]))
                stop("Software error: expected start < end always")
            }
        dfIndels = rbind.fast(dfIndels, dfi)

        # Update Indel stats.
        if (seqLen >= minRefLenStats)
            indelsPerKbp = c(indelsPerKbp, 1000*numIndels/seqLen)
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

        # Update SNP stats.
        if (seqLen >= minRefLenStats)
            snpsPerKbp = c(snpsPerKbp, 1000*numSNPs/seqLen)
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

# Compute min, max, mean, and standard deviation of SNP and Indel frequency.
if (length(snpsPerKbp) > 0)
    {
    min.SNPsperKbp = signif(min(snpsPerKbp), 2)
    max.SNPsperKbp = signif(max(snpsPerKbp), 2)
    mean.SNPsperKbp = signif(mean(snpsPerKbp), 2)
    sd.SNPsperKbp = signif(sd(snpsPerKbp), 2)
    }

if (length(indelsPerKbp) > 0)
    {
    min.IndelsperKbp = signif(min(indelsPerKbp), 2)
    max.IndelsperKbp = signif(max(indelsPerKbp), 2)
    mean.IndelsperKbp = signif(mean(indelsPerKbp), 2)
    sd.IndelsperKbp = signif(sd(indelsPerKbp), 2)
    }

# Display statistics.
catnow("Finished aligning sequences for Indel Groups and locating Indels and SNPs\n")
catnow("Total number of Indels is ", Nindels, "\n")
catnow("Total number of SNPs is ", Nsnps, "\n")
catnow("Number of LCB/IndelGroup/Marker alignments that had too many Indels and were ignored:", tooManyIndels, "\n")
catnow("Number of LCB/IndelGroup/Marker alignments that had too many SNPs and were ignored:", tooManySNPs, "\n")
catnow("Aligner program failed ", NalignerFails, " times\n")
catnow("Aligner program output not recognizable ", NalignerGarbage, " times\n")
if (length(indelsPerKbp) > 0)
    catnow("Indels per Kbp: min:", min.IndelsperKbp, " max:", max.IndelsperKbp,
        " mean:", mean.IndelsperKbp, " stddev:", sd.IndelsperKbp, "\n", sep="")
if (length(snpsPerKbp) > 0)
    catnow("SNPs per Kbp: min:", min.SNPsperKbp, " max:", max.SNPsperKbp,
        " mean:", mean.SNPsperKbp, " stddev:", sd.SNPsperKbp, "\n", sep="")
catnow("Indel output file:\n", outIndelFile, "\n")
catnow("SNP output file:\n", outSNPfile, "\n")
}

################################################################################
# End of file.
################################################################################
