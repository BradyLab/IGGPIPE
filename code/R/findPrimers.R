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
#thisDir = "~/Documents/UCDavis/BradyLab/Genomes/kmers/IGGPIPE/code/R/" # For testing only.

# Source the necessary include files from the same directory containing this file.
source(paste(thisDir, "Include_Common.R", sep=""))

# Get arguments.
testing = 0
#testing = 1 # For testing only.
#testing = 2 # For testing only.
#testing = 3 # For testing only.
{
if (testing == 0)
    args = commandArgs(TRUE)
else if (testing == 1)
    {
    args = c("~/Documents/UCDavis/BradyLab/Genomes/kmers/IGGPIPE",
        "outTestHP11/IndelGroupsOverlapping_K11k2L100D10_2000A100_2000d10_100N2F0.tsv",
        "outTestHP11/GenomeData/DNAseqs_K11k2L100D10_2000A100_2000d10_100N2F0X20",
        "outTestHP11/NonvalidatedMarkers_K11k2L100D10_2000A100_2000d10_100N2F0X20.tsv",
        "~/bin/primer3_core", "primer3settings.txt",
        "/Users/tedtoal/src/primer3-2.3.6/primer3_config",
        "outTestHP11/Primers/Primer3In.txt", "outTestHP11/Primers/Primer3Out.txt", TRUE)
    }
else if (testing == 2)
    {
    args = c("~/Documents/UCDavis/BradyLab/Genomes/kmers/IGGPIPE",
        "outHP14/IndelGroupsOverlapping_K14k2L400D10_1500A400_1500d50_300N2F0.tsv",
        "outHP14/GenomeData/DNAseqs_K14k2L400D10_1500A400_1500d50_300N2F0X20",
        "outHP14/NonvalidatedMarkers_K14k2L400D10_1500A400_1500d50_300N2F0X20.tsv",
        "~/bin/primer3_core", "primer3settings.txt",
        "/Users/tedtoal/src/primer3-2.3.6/primer3_config",
        "outHP14/Primers/Primer3In.txt", "outHP14/Primers/Primer3Out.txt", TRUE)
    }
else if (testing == 3)
    {
    args = c("~/Documents/UCDavis/BradyLab/Genomes/kmers/IGGPIPE",
        "outHPT14/IndelGroupsOverlapping_K14k2L300D5_1500A300_1500d50_300N2F0.tsv",
        "outHPT14/GenomeData/DNAseqs_K14k2L300D5_1500A300_1500d50_300N2F0X15",
        "outHPT14/NonvalidatedMarkers_K14k2L300D5_1500A300_1500d50_300N2F0X15.tsv",
        "~/bin/primer3_core", "primer3settings.txt",
        "/Users/tedtoal/src/primer3-2.3.6/primer3_config",
        "outHPT14/Primers/Primer3In.txt", "outHPT14/Primers/Primer3Out.txt", TRUE)
    }
else stop("Unknown value for 'testing'")
}

Nexpected = 10
if (length(args) != Nexpected)
    {
    usage = c(
        "Read a data frame of Indel Groups bounded by k-mer pairs and, for each genome, another data frame",
        "of genome DNA sequences around the k-mer pairs, then run Primer3 to create primers from those",
        "sequences.  Write the new candidate IGG marker data including primers to a new file.",
        "",
        "Usage: Rscript findPrimers.R <wd> <indelFile> <seqInPfx> <tsvMarkerFile> \\",
        "       <primer3core> <primer3settings> <primer3config> <primer3DataFile> <primer3OutFile> <investigate>",
        "",
        "Arguments:",
        "   <wd>              : Path of R working directory, specify other file paths relative to this.",
        "   <indelFile>       : Input file containing the indel k-mer pairs.",
        "   <seqInPfx>        : Prefix, including directory, of input files containing the DNA sequences.  The filename",
        "                       suffix is '_<genome number>.dnaseqs'",
        "   <tsvMarkerFile>   : Output file to be written containing the candidate marker k-mer pairs with DNA sequences.",
        "   <primer3core>     : Full path of the primer3_core program file",
        "   <primer3settings> : File containing primer3 settings edited for desired IGG marker primer characteristics.",
        "   <primer3config>   : Path of primer3 config directory containing .ds and .dh thermodynamic parameter files.",
        "   <primer3DataFile> : File in which to place sequence data for primer3 to use to design primers.",
        "   <primer3OutFile>  : File into which primer3 should place its primer design results.",
        "   <investigate>     : FALSE for normal operation, TRUE for more verbose debugging output."
        )
    for (S in usage)
        catnow(S, "\n", sep="")
    stop("Try again with correct number of arguments")
    }

catnow("findPrimers.R arguments:\n")
workingDirectory = args[1]
catnow("  workingDirectory: ", workingDirectory, "\n")
if (!dir.exists(workingDirectory))
    stop("Directory doesn't exist: ", workingDirectory)
setwd(workingDirectory)

indelFile = args[2]
catnow("  indelFile: ", indelFile, "\n")
if (!file.exists(indelFile))
    stop("File doesn't exist: ", indelFile)

seqInPfx = args[3]
catnow("  seqInPfx: ", seqInPfx, "\n")

tsvMarkerFile = args[4]
catnow("  tsvMarkerFile: ", tsvMarkerFile, "\n")

primer3core = args[5]
catnow("  primer3core: ", primer3core, "\n")

primer3settings = args[6]
catnow("  primer3settings: ", primer3settings, "\n")

primer3config = args[7]
catnow("  primer3config: ", primer3config, "\n")

primer3DataFile = args[8]
catnow("  primer3DataFile: ", primer3DataFile, "\n")

primer3OutFile = args[9]
catnow("  primer3OutFile: ", primer3OutFile, "\n")

investigate = as.logical(args[10])
catnow("  investigate: ", investigate, "\n")
if (is.na(investigate))
    stop("investigate must be TRUE or FALSE")

########################################
# Initialize.
########################################

# Read indel group k-mer pairs data, which become the initial marker candidates.
catnow("Reading InDel Group file...")
dfMarkers = read.table(indelFile, header=TRUE, row.names=NULL, sep="\t", stringsAsFactors=FALSE)
catnow("\n")
if (nrow(dfMarkers) == 0)
    stop("There are no Indel Groups.")
inv(dim(dfMarkers), "input markers dim")
inv(colnames(dfMarkers), "input markers columns")
inv(head(dfMarkers), "input markers head")

# Remove contig columns.
dfMarkers = dfMarkers[, !grepl("^ctg", colnames(dfMarkers))]

# The k-mer length we are working with.
kmerLen = nchar(dfMarkers$kmer1[1])
inv(kmerLen, "k")

# Get the genome identifying letters from the dfMarkers kkLen columns names.
genomeLtrs = sub("kkLen", "", colnames(dfMarkers)[grepl("kkLen", colnames(dfMarkers))])
Ngenomes = length(genomeLtrs)
if (Ngenomes < 2)
    stop("Expected to recognize at least two genomes in the data column names")
refGenomeLtr = genomeLtrs[1]
otherGenomeLtrs = genomeLtrs[-1]
inv(genomeLtrs, "genomeLtrs")

########################################
# Read indel k-mer DNA sequence data files and cbind the data into data frame df.
########################################

df = NULL
for (i in 1:Ngenomes)
    {
    catnow("Reading DNA sequence data for genome #", i, "...")
    seqInFile = paste(seqInPfx, "_", i, ".dnaseqs", sep="")
    dft = read.table(seqInFile, header=TRUE, row.names=NULL, sep="\t", stringsAsFactors=FALSE)
    dft = dft[order(dft$kmer),]
    # Add genome letter prefix.
    pfxNames = (colnames(dft) %in% c("id", "pos", "s", "seq"))
    colnames(dft)[pfxNames] = paste(genomeLtrs[i], colnames(dft)[pfxNames], sep="")
    if (is.null(df))
        {
        expectedCols = colnames(dft)
        df = dft
        }
    else
        {
        if (ncol(dft) != length(expectedCols)) stop("Expected all DNA sequence files to have same number of columns.")
        if (nrow(dft) != nrow(df)) stop("Expected all DNA sequence files to have same number of rows.")
        sameCols = colnames(dft)[colnames(dft) == expectedCols]
        for (col in sameCols)
            {
            noMatch = (dft[,col] != df[,col])
            if (any(noMatch))
                stop("Expected matching columns in DNA sequence files, ",
                    sum(noMatch), " mismatches out of ", length(noMatch),
                    " in column '", col, "' at i=", i)
            }
        dft = dft[, !colnames(dft) %in% sameCols]
        df = cbind(df, dft)
        }
    rm(dft)
    catnow("\n")
    }
if (nrow(df) == 0)
    stop("There is no indel DNA sequence data.")
inv(dim(df), "input DNA dim")
inv(colnames(df), "input DNA columns")
inv(head(df), "input DNA head")

########################################
# Create vectors of column names in df and dfMarkers, indexed by genome.
########################################

makeColVec = function(S)
    {
    V = paste(genomeLtrs, S, sep="")
    names(V) = genomeLtrs
    return(V)
    }

idCol = makeColVec("id")
refIdCol = idCol[refGenomeLtr]
posCol = makeColVec("pos")
refPosCol = posCol[refGenomeLtr]
pos1Col = makeColVec("pos1")
refPos1Col = pos1Col[refGenomeLtr]
pos2Col = makeColVec("pos2")
refPos2Col = pos2Col[refGenomeLtr]
strandCol = makeColVec("s")
refStrandCol = strandCol[refGenomeLtr]
strand1Col = makeColVec("s1")
refStrand1Col = strand1Col[refGenomeLtr]
strand2Col = makeColVec("s2")
refStrand2Col = strand2Col[refGenomeLtr]
kkLenCol = makeColVec("kkLen")
refKkLenCol = kkLenCol[refGenomeLtr]
pctCol = makeColVec("pct")
refPctCol = pctCol[refGenomeLtr]
dnaSeqCol = makeColVec("seq")
refDnaSeqCol = dnaSeqCol[refGenomeLtr]
dnaSeq1Col = makeColVec("seq1")
refDnaSeq1Col = dnaSeq1Col[refGenomeLtr]
dnaSeq2Col = makeColVec("seq2")
refDnaSeq2Col = dnaSeq2Col[refGenomeLtr]

########################################
# Replace the sequence bases in all but the reference genome with "." if the
# base matches the corresponding reference genome base.
########################################

catnow("Replacing bases with '.'...")

# Length of each DNA sequence.
DNAseqLen = nchar(df[1, refDnaSeqCol])

# The extensionLen value used to extract the DNA sequences.
extensionLen = (DNAseqLen - kmerLen)/2

refBases = strsplit(df[, refDnaSeqCol], "", fixed=TRUE)
for (genome in otherGenomeLtrs)
    {
    genomeBases = strsplit(df[, dnaSeqCol[genome]], "", fixed=TRUE)
    newBases = sapply(1:length(genomeBases), function(i)
        {
        v = genomeBases[[i]]
        v[v == refBases[[i]]] = "."
        return(list(v))
        })
    rm(genomeBases)
    df[, dnaSeqCol[genome]] = sapply(newBases, paste, collapse="")
    rm(newBases)
    }
rm(refBases)
catnow("\n")

# Put the data frame rows in order by reference genome position.
catnow("Sorting by reference genome position...")
df = df[order(df[, refIdCol], df[, refPosCol]),]
rownames(df) = NULL
catnow("\n")

# Put the data frame columns in order by genome.
catnow("Reordering columns...")
df = df[, c(setdiff(colnames(df), "kmer"), "kmer")]
for (genome in genomeLtrs)
    {
    genomeCols = colnames(df)[grepl(paste("^", genome, sep=""), colnames(df))]
    df = df[, c(setdiff(colnames(df), genomeCols), genomeCols)]
    }
catnow("\n")

# The extensionLen value used to extract the DNA sequences.
extensionLen = (nchar(df[1, refDnaSeqCol]) - kmerLen)/2

########################################
# Now merge the DNA sequence data in columns dnaSeqCol of df into dfMarkers
# columns dnaSeq1Col and dnaSeq2Col.  Use the kmer1 and kmer2 columns in
# dfMarkers as keys to match up with the df$kmer column.
########################################

catnow("Matching k-mers 1...")
kmer1idxs = match(dfMarkers$kmer1, df$kmer)
catnow("\n")
catnow("Matching k-mers 1...")
kmer2idxs = match(dfMarkers$kmer2, df$kmer)
catnow("\n")

if (any(is.na(kmer1idxs))) stop("Missing kmer1 in primer output")
if (any(is.na(kmer2idxs))) stop("Missing kmer2 in primer output")

catnow("Merging sequence data...")
dfMarkers[, dnaSeq1Col] = df[kmer1idxs, dnaSeqCol]
dfMarkers[, dnaSeq2Col] = df[kmer2idxs, dnaSeqCol]
rm(df, kmer1idxs, kmer2idxs)
inv(dim(dfMarkers), "Marker data frame with DNA sequence data dim")
catnow("\n")

########################################
# You may be wondering how primers are made if the extracted sequences in the
# extension areas differ between genomes (the k-mer regions themselves are
# obviously identical).  The answer is that bases that differ are set to N and
# the primer program will not use them for primers, so only indel k-mers with
# sufficient identity on one or both sides of the k-mers will result in primers
# being generated.
########################################

# Concatenate the two DNA sequences at kmer1 and kmer2 for each marker.  The
# primers will be designed using the concatenated sequence.  Note that when the
# amplicon is on the opposite strand of the reference genome, it was extracted
# in reverse complement order, so concatenating the k-mer #1 and k-mer #2
# sequences is the correct thing to do; k-mer #1 position in that case is >
# k-mer #2 position.

# The following can be extremely intensive in memory usage.  We can easily
# break this up and do it in pieces to avoid that.  Let's do a million rows of
# dfMarkers at a time.
NatOnce = 1000000
Npasses = as.integer((nrow(dfMarkers)-1)/NatOnce) + 1
curPass = 0
curStart = 1
seqs = c() # Result sequences.
catnow("Converting '.' to 'N' in DNA sequences...\n")
while (curStart <= nrow(dfMarkers))
    {
    curEnd = curStart + NatOnce - 1
    if (curEnd > nrow(dfMarkers))
        curEnd = nrow(dfMarkers)
    rng = curStart:curEnd

    #catnow("Concatenating DNA sequences...")
    tmpSeqs = paste(dfMarkers[rng, refDnaSeq1Col], dfMarkers[rng, refDnaSeq2Col], sep="")
    #catnow("\n")

    # Split the sequences apart into vectors of characters.
    #catnow("Splitting sequences...")
    tmpSeqs = strsplit(tmpSeqs, "", fixed=TRUE)
    #catnow("\n")

    # Convert reference genome sequence bases to N if corresponding base in other
    # genome does not match (is not ".").
    #catnow("Changing bases to N...\n")
    for (genome in otherGenomeLtrs)
        {
        #catnow("  Concatenating DNA sequences for genome", genome, "...")
        tmpSeqsO = paste(dfMarkers[rng, dnaSeq1Col[genome]], dfMarkers[rng, dnaSeq2Col[genome]], sep="")
        #catnow("\n")
        #catnow("  Splitting sequences for genome", genome, "...")
        tmpSeqsO = strsplit(tmpSeqsO, "", fixed=TRUE)
        #catnow("\n")
        #catnow("  Changing '.' to 'N'...")
        tmpSeqs = sapply(1:length(tmpSeqs), function(i)
            {
            v = tmpSeqs[[i]]
            v[tmpSeqsO[[i]] != "."] = "N"
            return(list(v))
            })
        rm(tmpSeqsO)
        #catnow("\n")
        }
    #catnow("Done\n")

    # Rejoin the split-up sequence characters.
    catnow("Rejoining sequence characters...")
    tmpSeqs = sapply(tmpSeqs, paste, collapse="")
    seqs = c(seqs, tmpSeqs)
    rm(tmpSeqs)
    catnow("\n")

    # Next NatOnce rows of dfMarkers.
    curStart = curEnd + 1
    curPass = curPass+1
    catnow(" ", curPass, " of ", Npasses, "\n", sep="")
    }

########################################
# Create an input file for primer3 to make a pair of primers for each marker.
# For SEQUENCE_ID, we will paste together kmer1 and kmer2 columns with "_".
# All SEQUENCE_PRIMER_PAIR_OK_REGION_LIST records are identical.
# Sample file record:
#   SEQUENCE_ID=CATCCGGACCC_CGTCACGAGTC
#   SEQUENCE_TEMPLATE=TGGTACAANGGGGTTCATCCGGACCCCNTNNNNNNAAAANTAAAAGAAAAAGACATGACTCGTGACGACAATTTAAATAATT
#   SEQUENCE_PRIMER_PAIR_OK_REGION_LIST=1,41,42,41
#   =
#
# We write a first line to the input file that gives a path to the primer
# thermodynamic settings files with PRIMER_THERMODYNAMIC_PARAMETERS_PATH.
########################################

catnow("Creating primer3 input file...")
thermo = paste("PRIMER_THERMODYNAMIC_PARAMETERS_PATH=", primer3config, PATHSEP, sep="")
seqIDs = paste("SEQUENCE_ID=", dfMarkers$kmer1, "_", dfMarkers$kmer2, sep="")
primerTemplates = paste("SEQUENCE_TEMPLATE=", seqs, sep="")
rm(seqs)
seqLen = 2*extensionLen+kmerLen
primerRegions = paste("SEQUENCE_PRIMER_PAIR_OK_REGION_LIST=1,", seqLen, ",", seqLen+1,",", seqLen, sep="")
primerRegions = rep(primerRegions, nrow(dfMarkers))
recordSeps = rep("=", nrow(dfMarkers))
records = c(rbind(thermo, seqIDs, primerTemplates, primerRegions, recordSeps))
numSeqIDs = length(seqIDs)
rm(seqIDs, primerTemplates, primerRegions, recordSeps)
writeLines(records, primer3DataFile)
rm(records)
catnow("\n")

########################################
# Run primer3 on the sequence data file we just created.
########################################

if (!file.exists(primer3settings)) stop("Primer3 settings file ", primer3settings, " not found")
primer3_args = c("-p3_settings_file", primer3settings)
catnow("Running primer3 to design primers\n")
estimateTime = FALSE
if (estimateTime)
    {
    maxSecondsPerPrimer = 0.025
    totalSeconds = round(maxSecondsPerPrimer * numSeqIDs)
    maxMinutes = ceiling(totalSeconds/60)
    catnow("Expect this to take up to", maxMinutes, "minutes on a slower computer.\n")
    }
stat = system2(primer3core, primer3_args, stdout=primer3OutFile, stdin=primer3DataFile)
if (stat != 0)
    stop("Program ", primer3core, " exited with error status ", stat)

########################################
# Read and parse the primer3 output file to extract the primers, and include
# them as columns in dfDNA.
########################################

# Sample output.
#   PRIMER_THERMODYNAMIC_PARAMETERS_PATH=/Users/tedtoal/src/primer3-2.3.6/primer3_config/
#   SEQUENCE_ID=GGGGCACAATTGGA_GTCGATCAAAGCAC
#   SEQUENCE_TEMPLATE=CAATTTAAAATCCAATTGTGCCCCAAATTCTCTTGTCTGAAATTGTGCTTTGATCGACAAAAGTNTCG
#   SEQUENCE_PRIMER_PAIR_OK_REGION_LIST=1,34,35,34
#   PRIMER_LEFT_NUM_RETURNED=1
#   PRIMER_RIGHT_NUM_RETURNED=1
#   PRIMER_INTERNAL_NUM_RETURNED=0
#   PRIMER_PAIR_NUM_RETURNED=1
#   PRIMER_PAIR_0_PENALTY=7.768749
#   PRIMER_LEFT_0_PENALTY=3.982241
#   PRIMER_RIGHT_0_PENALTY=3.786507
#   PRIMER_LEFT_0_SEQUENCE=TTTAAAATCCAATTGTGCCCCA
#   PRIMER_RIGHT_0_SEQUENCE=TGTCGATCAAAGCACAATTTCAG
#   PRIMER_LEFT_0=4,22
#   PRIMER_RIGHT_0=59,23
#   PRIMER_LEFT_0_TM=57.018
#   PRIMER_RIGHT_0_TM=58.213
#   PRIMER_LEFT_0_GC_PERCENT=36.364
#   PRIMER_RIGHT_0_GC_PERCENT=39.130
#   PRIMER_LEFT_0_SELF_ANY_TH=0.00
#   PRIMER_RIGHT_0_SELF_ANY_TH=0.00
#   PRIMER_LEFT_0_SELF_END_TH=0.00
#   PRIMER_RIGHT_0_SELF_END_TH=0.00
#   PRIMER_LEFT_0_HAIRPIN_TH=0.00
#   PRIMER_RIGHT_0_HAIRPIN_TH=0.00
#   PRIMER_LEFT_0_END_STABILITY=4.9600
#   PRIMER_RIGHT_0_END_STABILITY=3.0200
#   PRIMER_LEFT_0_TEMPLATE_MISPRIMING=4.0000
#   PRIMER_RIGHT_0_TEMPLATE_MISPRIMING=6.0000
#   PRIMER_PAIR_0_COMPL_ANY_TH=10.00
#   PRIMER_PAIR_0_COMPL_END_TH=0.57
#   PRIMER_PAIR_0_PRODUCT_SIZE=56
#   PRIMER_PAIR_0_TEMPLATE_MISPRIMING=10.00

# Discard unneeded lines.  We wish to retain lines:
#   SEQUENCE_ID
#   PRIMER_LEFT_0_SEQUENCE
#   PRIMER_RIGHT_0_SEQUENCE
#   PRIMER_LEFT_0
#   PRIMER_RIGHT_0
#   PRIMER_LEFT_0_TM
#   PRIMER_RIGHT_0_TM
#   =

# It is possible that the primer3 output file is huge.  We reduce the size of the
# data vastly by removing unwanted parameters.  So, let's read the file in chunks
# of NlinesAtOnce until end of file.
NlinesAtOnce = 1000000

primers = c()   # Collect filtered and simplified data lines here.
catnow("Processing primer data\n")
isEOF = FALSE
inFile = file(primer3OutFile, open="r")
while (!isEOF)
    {
    lines = readLines(inFile, n=NlinesAtOnce)
    if (length(lines) < NlinesAtOnce)
        isEOF = TRUE
    if (length(lines) == 0)
        break

    # Remove lines we don't care about.
    lines = lines[!grepl("PRIMER_THERMODYNAMIC_PARAMETERS_PATH", lines, fixed=TRUE)]
    lines = lines[!grepl("SEQUENCE_PRIMER_PAIR_OK_REGION_LIST", lines, fixed=TRUE)]
    lines = lines[!grepl("PRIMER_INTERNAL_NUM_RETURNED", lines, fixed=TRUE)]
    lines = lines[!grepl("_0_GC_PERCENT", lines, fixed=TRUE)]
    lines = lines[!grepl("_0_SELF_END_TH", lines, fixed=TRUE)]
    lines = lines[!grepl("_0_END_STABILITY", lines, fixed=TRUE)]
    lines = lines[!grepl("SEQUENCE_TEMPLATE", lines, fixed=TRUE)]
    lines = lines[!grepl("_NUM_RETURNED", lines, fixed=TRUE)]
    lines = lines[!grepl("_0_SELF_ANY_TH", lines, fixed=TRUE)]
    lines = lines[!grepl("_0_HAIRPIN_TH", lines, fixed=TRUE)]
    lines = lines[!grepl("_0_TEMPLATE_MISPRIMING", lines, fixed=TRUE)]
    lines = lines[!grepl("_0_PENALTY", lines, fixed=TRUE)]
    lines = lines[!grepl("PRIMER_PAIR_0_", lines, fixed=TRUE)]
    lines = lines[!grepl("SEQUENCE_ID", lines, fixed=TRUE)] # We previously were retaining this but not using it, then discarding it.

    # Reduce the size of remaining lines by trimming away some "fat".
    lines = sub("SEQUENCE", "SEQ", lines)
    lines = sub("PRIMER_", "", lines)
    lines = sub("LEFT_0", "L", lines)
    lines = sub("RIGHT_0", "R", lines)
    lines = sub("^L=", "L_POS=", lines)
    lines = sub("^R=", "R_POS=", lines)

    primers = c(primers, lines)
    rm(lines)
    }
close(inFile)

# Find starting and ending lines of each primer group.
endGroup = which(primers == "=")
startGroup = c(1, endGroup[-length(endGroup)]+1)

# Make sure each group's size is either 0 (no primer found) or 6 (left and right
# sequence, position, and Tm) lines.
Nlines = endGroup-startGroup
if (!all(Nlines == 0 | Nlines == 6))
    stop("Expected all primers to have 6 data lines but something else was found")
rm(Nlines)

# Split them up into separate groups, one per primer.
primers = sapply(1:length(startGroup), function(i)
    {
    if (startGroup[i] == endGroup[i])
        v = c("", "", "", "", "", "")
    else
        {
        v = primers[startGroup[i]:(endGroup[i]-1)]
        v = unlist(strsplit(v, "=", fixed=TRUE))
        v = v[seq(2, length(v), 2)]
        }
    return(list(v))
    })
rm(startGroup, endGroup)

catnow("Binding primer data...")
primers = matrix(unlist(primers, use.names=FALSE), ncol=length(primers[[1]]), byrow=TRUE)
colnames(primers) = c("seqL", "seqR", "poslenL", "poslenR", "tmL", "tmR")
primers = as.data.frame(primers, stringsAsFactors=FALSE)
catnow("\n")

catnow("Editing primer data...")
primers$tmL = round(as.numeric(primers$tmL), 2)
primers$tmR = round(as.numeric(primers$tmR), 2)
primers$posL = as.integer(sub(",[0-9]+$", "", primers$poslenL))
primers$lenL = as.integer(sub("^[0-9]+,", "", primers$poslenL))
primers = primers[, colnames(primers) != "poslenL"]
primers$posR = as.integer(sub(",[0-9]+$", "", primers$poslenR))
primers$lenR = as.integer(sub("^[0-9]+,", "", primers$poslenR))
primers = primers[, colnames(primers) != "poslenR"]
primers = primers[, c("seqL", "tmL", "posL", "lenL", "seqR", "tmR", "posR", "lenR")]
colnames(primers) = c("prmSeqL",  "prmTmL",  "prmPosL",  "prmLenL", "prmSeqR",  "prmTmR",  "prmPosR",  "prmLenR")
if (nrow(primers) != nrow(dfMarkers))
    stop("Wrong number of primers returned by primer3")
dfMarkers = cbind(dfMarkers, primers, stringsAsFactors=FALSE)
rm(primers)
catnow("\n")

########################################
# There will be many markers for which primers could not be found.  Get rid of these.
########################################

{
catnow("Number of markers initially: ", nrow(dfMarkers), "\n")
catnow("No left primer found:        ", sum(dfMarkers$prmSeqL == "" & dfMarkers$prmSeqR != ""), "\n")
catnow("No right primer found:       ", sum(dfMarkers$prmSeqL != "" & dfMarkers$prmSeqR == ""), "\n")
catnow("Neither left nor right found:", sum(dfMarkers$prmSeqL == "" & dfMarkers$prmSeqR == ""), "\n")
dfMarkers = dfMarkers[dfMarkers$prmSeqL != "" & dfMarkers$prmSeqR != "",]
catnow("Number of markers remaining: ", nrow(dfMarkers), "\n")
}

########################################
# If no markers remain, exit with error.
########################################

if (nrow(dfMarkers) == 0)
    stop("No markers remaining.")

########################################
# There may be markers with identical primer-pairs.  This happens when two k-mers
# are very near each other and the same primer is generated for both, on both the
# start and end sides of the amplicon.  Get rid of duplicates.  We don't care which
# ones we discard (losing certain k-mers is not an issue).
########################################

catnow("Removing duplicates...")
dupes = duplicated(paste(dfMarkers$prmSeqL, dfMarkers$prmSeqR, sep="_"))
dfMarkers = dfMarkers[!dupes,]
catnow(sum(dupes), "duplicate primer-pairs removed\n")
rm(dupes)

########################################
# Add "phase" columns giving the genome strand which matches the reference genome
# "+" strand.  This column is "-" if a marker's region is running in the opposite
# direction as the reference genome.  For the reference genome it is always "+",
# so we will remove it from the data frame we write to a file.
########################################

catnow("Adding phase columns...")
phaseCol = makeColVec(paste(refGenomeLtr, "phase", sep=""))
for (genome in genomeLtrs)
    dfMarkers[, phaseCol[genome]] = ifelse(dfMarkers[, refStrand1Col] == dfMarkers[, strand1Col[genome]], "+", "-")
catnow("\n")

########################################
# Change the "prmPosL" and "prmPosR" columns, which give the offset+1 of the 5'
# end of the primer from the 5' end of the DNA sequence, to columns "kmer1offset"
# and "kmer2offset", which give the offset+0 of the outside end of the k-mer from
# the end of the amplicon (the amplicon end that is closest to that k-mer).  A
# value of 0 means the amplicon and k-mer ends correspond.  Values may be positive
# or negative: positive if the k-mer lies inside the amplicon and negative if it
# lies outside.  The same values apply to all genomes.
# Note: The DNA seqs from which the primers were made include extensionLen extra
# base pairs in each direction beyond the k-mer, and the kmer2 DNA sequence was
# concatenated to the kmer1 DNA sequence.
########################################

catnow("Adding kmer offset columns...")
dfMarkers$kmer1offset = extensionLen - (dfMarkers$prmPosL-1)
seqLen = 2*(2*extensionLen+kmerLen)
kmer2End = seqLen - extensionLen
dfMarkers$kmer2offset = dfMarkers$prmPosR - kmer2End
catnow("\n")

########################################
# Change the "*pos1" and "*pos2" columns, which give the position of the 5' end
# of each k-mer on the + strand, to columns "*ampPos1" and "*ampPos2", giving
# the positions of the ends of the amplicon.
#
# The position change is different for each genome depending on whether the k-mer
# "phase" is + or -.  The reference genome is always "+" phase for all k-mers,
# while the other genomes have both phases.  For a "+" phase genome, the pos1
# values are adjusted in one manner, and the pos2 values are adjusted in a
# different manner, to give the ampPos1 and ampPos2 values.  For a "-" phase
# genome, the adjustment manners are switched.
# Adjustment for a "+" phase genome:
#   - ampPos1 = pos1 - kmer1offset
#   - ampPos2 = pos2 + kmer2offset + k-1
# Adjustment for a "-" phase genome:
#   - ampPos1 = pos1 + kmer1offset + k-1
#   - ampPos2 = pos2 - kmer2offset
########################################

catnow("Adding amplicon position columns...")
ampPos1Col = makeColVec("ampPos1")
ampPos2Col = makeColVec("ampPos2")
refAmpPos1Col = ampPos1Col[refGenomeLtr]
refAmpPos2Col = ampPos2Col[refGenomeLtr]

for (genome in genomeLtrs)
    {
    posPhase = (dfMarkers[, phaseCol[genome]] == "+")
    negPhase = !posPhase
    dfMarkers[, ampPos1Col[genome]] = dfMarkers[, pos1Col[genome]] - dfMarkers$kmer1offset
    dfMarkers[, ampPos2Col[genome]] = dfMarkers[, pos2Col[genome]] - dfMarkers$kmer2offset
    dfMarkers[negPhase, ampPos1Col[genome]] =
        dfMarkers[negPhase, pos1Col[genome]] + dfMarkers[negPhase, "kmer1offset"] + (kmerLen-1)
    dfMarkers[posPhase, ampPos2Col[genome]] =
        dfMarkers[posPhase, pos2Col[genome]] + dfMarkers[posPhase, "kmer2offset"] + (kmerLen-1)
    rm(posPhase, negPhase)
    }
catnow("\n")

########################################
# Remove the *kkLen columns, which contain k-mer-to-k-mer lengths, and replace
# them with *ampLen columns, which contain amplicon lengths that are simply
# abs(ampPos1-ampPos2+1).
########################################

catnow("Adding amplicon length columns...")
ampLenCol = makeColVec("ampLen")
refAmpLenCol = ampLenCol[refGenomeLtr]
for (genome in genomeLtrs)
    {
    dfMarkers[, ampLenCol[genome]] = abs(dfMarkers[, ampPos1Col[genome]] - dfMarkers[, ampPos2Col[genome]]) + 1
    dfMarkers = dfMarkers[, colnames(dfMarkers) != kkLenCol[genome]]
    }
catnow("\n")

########################################
# Add "dif" columns giving difference of amplicon length of each other genome
# to the reference genome.  There is no such column for the reference genome
# since it would always be 0, but columns exist for each other genome.
########################################

catnow("Adding dif columns...")
difToRefCol = makeColVec(paste(refGenomeLtr, "dif", sep=""))
for (genome in otherGenomeLtrs)
    dfMarkers[, difToRefCol[genome]] = dfMarkers[, ampLenCol[genome]] - dfMarkers[, refAmpLenCol]
catnow("\n")

########################################
# Concatenate the k-mer strand signs of all genomes into a single column.
########################################

catnow("Concatenating strand signs...")
dfMarkers$kmer1strands = apply(dfMarkers[, strand1Col], 1, paste, collapse="")
dfMarkers$kmer2strands = apply(dfMarkers[, strand2Col], 1, paste, collapse="")
catnow("\n")

########################################
# Put the data in a nice order and write it to the output file.
########################################

catnow("Setting column order...")

# Put columns in a nice order.  Put NDA column first.
colOrder = c("NDA")

# Important marker position columns first.
for (genome in genomeLtrs)
    {
    colOrder = c(colOrder, idCol[genome], pctCol[genome], ampLenCol[genome])
    if (genome != refGenomeLtr)
        colOrder = c(colOrder, difToRefCol[genome], phaseCol[genome])
    }

# Primer columns next, in left-right pairs, most important primer columns first.
colOrder = c(colOrder, "prmSeqL", "prmSeqR", "prmTmL", "prmTmR", "prmLenL", "prmLenR")

# Amplicon position.
for (genome in genomeLtrs)
    colOrder = c(colOrder, ampPos1Col[genome], ampPos2Col[genome])

# K-mer columns.
colOrder = c(colOrder, "kmer1", "kmer1strands", "kmer1offset", "kmer2", "kmer2strands", "kmer2offset")

# DNA sequence columns.
for (genome in genomeLtrs)
    colOrder = c(colOrder, dnaSeq1Col[genome], dnaSeq2Col[genome])

# Put the columns in order and discard unwanted columns.
dfMarkers = dfMarkers[, colOrder[]]
rownames(dfMarkers) = NULL

catnow("\n")

########################################
# Write the dfMarkers data frame to the output file.
########################################

catnow("Writing output file...")
write.table(dfMarkers, tsvMarkerFile, col.names=TRUE, row.names=FALSE, quote=FALSE, sep="\t")
# dfMarkers = read.table(tsvMarkerFile, header=TRUE, sep="\t", stringsAsFactors=FALSE)
catnow("\n\n")
catnow("Finished adding primer sequences to Indel Groups, candidate marker output file:\n", tsvMarkerFile, "\n")
}

################################################################################
# End of file.
################################################################################
