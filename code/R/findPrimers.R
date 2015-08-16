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

# Cat that immediately flushes to console.
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
testing = FALSE
#testing = TRUE # For testing only.
{
if (!testing)
    args = commandArgs(TRUE)
else
    {
    args = c("~/Documents/UCDavis/BradyLab/Genomes/kmers/SCARF",
        "outTestHP11/IndelsOverlapping_K11Km2Lm100Dm10Dx2000Am100Ax2000ADm10ADx100ND2mF0.tsv",
        "outTestHP11/GenomeData/Genome_",
        "outTestHP11/CandidateMarkers_K11Km2Lm100Dm10Dx2000Am100Ax2000ADm10ADx100ND2mF0XL20.tsv",
        "~/bin/primer3_core", "Primer3Settings.txt",
        "outTestHP11/Primers/Primer3Data.txt", "outTestHP11/Primers/Primer3Out.txt", TRUE)
    }
}

Nexpected = 9
if (length(args) != Nexpected)
    {
    usage = c(
        "Read a data frame of indels bounded by k-mer pairs and, for each genome, another data frame",
        "of genome DNA sequences around the k-mer pairs, then run Primer3 to create primers from those",
        "sequences.  Write the new candidate SCAR marker data including primers to a new file.",
        "",
        "Usage: Rscript findPrimers.R <wd> <indelFile> <seqInPfx> <tsvMarkerFile> \\",
        "       <primer3core> <primer3Settings> <primer3DataFile> <primer3OutFile> <investigate>",
        "",
        "Arguments:",
        "   <wd>              : Path of R working directory, specify other file paths relative to this.",
        "   <indelFile>       : Input file containing the indel k-mer pairs.",
        "   <seqInPfx>        : Prefix, including directory, of input files containing the DNA sequences.  The filename",
        "                       suffix is the genome number followed by '.dnaseqs'",
        "   <tsvMarkerFile>   : Output file to be written containing the candidate marker k-mer pairs with DNA sequences.",
        "   <primer3core>     : Full path of the primer3_core program file",
        "   <primer3Settings> : File containing primer3 settings edited for desired SCAR marker primer characteristics.",
        "   <primer3DataFile> : File in which to place sequence data for primer3 to use to design primers.",
        "   <primer3OutFile>  : File into which primer3 should place its primer design results.",
        "   <investigate>     : FALSE for normal operation, TRUE for more verbose debugging output."
        )
    for (S in usage)
        catnow(S, "\n", sep="")
    stop("Try again with correct number of arguments")
    }

catnow("findPrimers arguments:\n")
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
if (!file.exists(primer3core))
    stop("File doesn't exist: ", primer3core)

primer3Settings = args[6]
catnow("  primer3Settings: ", primer3Settings, "\n")

primer3DataFile = args[7]
catnow("  primer3DataFile: ", primer3DataFile, "\n")

primer3OutFile = args[8]
catnow("  primer3OutFile: ", primer3OutFile, "\n")

investigate = as.logical(args[9])
catnow("  investigate: ", investigate, "\n")
if (is.na(investigate))
    stop("investigate must be TRUE or FALSE")

########################################
# Initialize.
########################################

# Read indel k-mer pairs data.
dfInDels = read.table(indelFile, header=TRUE, row.names=NULL, sep="\t", stringsAsFactors=FALSE)
if (nrow(dfInDels) == 0)
    stop("There are no INDELs.")
inv(dim(dfInDels), "input markers dim")
inv(colnames(dfInDels), "input markers columns")
inv(head(dfInDels), "input markers head")

# Remove contig columns.
dfInDels = dfInDels[, !grepl("^ctg", colnames(dfInDels))]

# The k-mer length we are working with.
kmerLen = nchar(dfInDels$kmer1[1])
inv(kmerLen, "k")

# Get the genome identifying letters from the dfInDels kkLen columns names.
genomeLtrs = sub("kkLen", "", colnames(dfInDels)[grepl("kkLen", colnames(dfInDels))])
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
    seqInFile = paste(seqInPfx, i, ".dnaseqs", sep="")
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
    }
if (nrow(df) == 0)
    stop("There is no indel DNA sequence data.")
inv(dim(df), "input DNA dim")
inv(colnames(df), "input DNA columns")
inv(head(df), "input DNA head")

########################################
# Create vectors of column names in df and dfInDels, indexed by genome.
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
    df[, dnaSeqCol[genome]] = sapply(newBases, paste, collapse="")
    }

# Put the data frame rows in order by reference genome position.
df = df[order(df[, refIdCol], df[, refPosCol]),]
rownames(df) = NULL

# Put the data frame columns in order by genome.
df = df[, c(setdiff(colnames(df), "kmer"), "kmer")]
for (genome in genomeLtrs)
    {
    genomeCols = colnames(df)[grepl(paste("^", genome, sep=""), colnames(df))]
    df = df[, c(setdiff(colnames(df), genomeCols), genomeCols)]
    }

# The extensionLen value used to extract the DNA sequences.
extensionLen = (nchar(df[1, refDnaSeqCol]) - kmerLen)/2

########################################
# Now merge the DNA sequence data in columns dnaSeqCol of df into dfInDels
# columns dnaSeq1Col and dnaSeq2Col.  Use the kmer1 and kmer2 columns in
# dfInDels as keys to match up with the df$kmer column.
########################################

kmer1idxs = match(dfInDels$kmer1, df$kmer)
kmer2idxs = match(dfInDels$kmer2, df$kmer)

if (any(is.na(kmer1idxs))) stop("Missing kmer1 in primer output")
if (any(is.na(kmer2idxs))) stop("Missing kmer2 in primer output")

dfInDels[, dnaSeq1Col] = df[kmer1idxs, dnaSeqCol]
dfInDels[, dnaSeq2Col] = df[kmer2idxs, dnaSeqCol]
inv(dim(dfInDels), "Marker data frame with DNA sequence data dim")

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
seqs = paste(dfInDels[, refDnaSeq1Col], dfInDels[, refDnaSeq2Col], sep="")

# Split the sequences apart into vectors of characters.
seqs = strsplit(seqs, "", fixed=TRUE)

# Convert reference genome sequence bases to N if corresponding base in other
# genome does not match (is not ".").
for (genome in otherGenomeLtrs)
    {
    seqsO = paste(dfInDels[, dnaSeq1Col[genome]], dfInDels[, dnaSeq2Col[genome]], sep="")
    seqsO = strsplit(seqsO, "", fixed=TRUE)
    seqs = sapply(1:length(seqs), function(i)
        {
        v = seqs[[i]]
        v[seqsO[[i]] != "."] = "N"
        return(list(v))
        })
    }

# Rejoin the split-up sequence characters.
seqs = sapply(seqs, paste, collapse="")

########################################
# Create an input file for primer3 to make a pair of primers for each marker.
# For SEQUENCE_ID, we will paste together kmer1 and kmer2 columns with "_".
# All SEQUENCE_PRIMER_PAIR_OK_REGION_LIST records are identical.
# Sample file record:
#   SEQUENCE_ID=CATCCGGACCC_CGTCACGAGTC
#   SEQUENCE_TEMPLATE=TGGTACAANGGGGTTCATCCGGACCCCNTNNNNNNAAAANTAAAAGAAAAAGACATGACTCGTGACGACAATTTAAATAATT
#   SEQUENCE_PRIMER_PAIR_OK_REGION_LIST=1,41,42,41
#   =
########################################

seqIDs = paste("SEQUENCE_ID=", dfInDels$kmer1, "_", dfInDels$kmer2, sep="")
primerTemplates = paste("SEQUENCE_TEMPLATE=", seqs, sep="")
seqLen = 2*extensionLen+kmerLen
primerRegions = paste("SEQUENCE_PRIMER_PAIR_OK_REGION_LIST=1,", seqLen, ",", seqLen+1,",", seqLen, sep="")
primerRegions = rep(primerRegions, nrow(dfInDels))
recordSeps = rep("=", nrow(dfInDels))
records = c(rbind(seqIDs, primerTemplates, primerRegions, recordSeps))
writeLines(records, primer3DataFile)

########################################
# Run primer3 on the sequence data file we just created.
########################################

if (!file.exists(primer3Settings)) stop("Primer3 settings file ", primer3Settings, " not found")
primer3_cmd = paste(primer3core, "-p3_settings_file", primer3Settings, "<", primer3DataFile, ">", primer3OutFile)
catnow("Running primer3 to design primers:\n   ", primer3_cmd, "\n")
estimateTime = FALSE
if (estimateTime)
    {
    maxSecondsPerPrimer = 0.025
    totalSeconds = round(maxSecondsPerPrimer * length(seqIDs))
    maxMinutes = ceiling(totalSeconds/60)
    catnow("Expect this to take up to", maxMinutes, "minutes on a slower computer.\n")
    }
system(primer3_cmd)

########################################
# Read and parse the primer3 output file to extract the primers, and include
# them as columns in dfDNA.
########################################

# Sample output.
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

catnow("Processing primer data\n")
primers = readLines(primer3OutFile)

primers = primers[!grepl("SEQUENCE_PRIMER_PAIR_OK_REGION_LIST", primers, fixed=TRUE)]
primers = primers[!grepl("PRIMER_INTERNAL_NUM_RETURNED", primers, fixed=TRUE)]
primers = primers[!grepl("_0_GC_PERCENT", primers, fixed=TRUE)]
primers = primers[!grepl("_0_SELF_END_TH", primers, fixed=TRUE)]
primers = primers[!grepl("_0_END_STABILITY", primers, fixed=TRUE)]
primers = primers[!grepl("SEQUENCE_TEMPLATE", primers, fixed=TRUE)]
primers = primers[!grepl("_NUM_RETURNED", primers, fixed=TRUE)]
primers = primers[!grepl("_0_SELF_ANY_TH", primers, fixed=TRUE)]
primers = primers[!grepl("_0_HAIRPIN_TH", primers, fixed=TRUE)]
primers = primers[!grepl("_0_TEMPLATE_MISPRIMING", primers, fixed=TRUE)]
primers = primers[!grepl("_0_PENALTY", primers, fixed=TRUE)]
primers = primers[!grepl("PRIMER_PAIR_0_", primers, fixed=TRUE)]

# Split them up into separate groups, one per primer.
endGroup = which(primers == "=")
startGroup = c(1, endGroup[-length(endGroup)]+1)
groups = sapply(1:length(startGroup), function(i)
    {
    v = primers[startGroup[i]:(endGroup[i]-1)]
    v = unlist(strsplit(v, "=", fixed=TRUE))
    vnames = v[seq(1, length(v), 2)]
    vnames = sub("PRIMER_(LEFT|RIGHT)_0", "\\1", vnames)
    vnames[vnames == "LEFT"] = "LEFT_POSLEN"
    vnames[vnames == "RIGHT"] = "RIGHT_POSLEN"
    v = v[seq(2, length(v), 2)]
    names(v) = vnames
    if (length(v) != 7)
        {
        if (length(v) != 1) stop("Unexpected primer3 output data, expected these edited names:\n",
            "SEQUENCE_ID, and on good primers: LEFT_SEQUENCE, RIGHT_SEQUENCE, LEFT_POS, RIGHT_POS,\n",
            "LEFT_TM, and RIGHT_TM.  Instead we see these:\n",
            paste(vnames, collapse=" "))
        # Following must be in same order as cases where a primer WAS found.
        v["LEFT_SEQUENCE"] = ""
        v["RIGHT_SEQUENCE"] = ""
        v["LEFT_POSLEN"] = ""
        v["RIGHT_POSLEN"] = ""
        v["LEFT_TM"] = ""
        v["RIGHT_TM"] = ""
        }
    return(list(v))
    })
dfPrimers = do.call(rbind, groups)
rm(groups)
colnames(dfPrimers) = c("name", "seqL", "seqR", "poslenL", "poslenR", "tmL", "tmR")
dfPrimers = as.data.frame(dfPrimers, stringsAsFactors=FALSE)
dfPrimers$tmL = round(as.numeric(dfPrimers$tmL), 2)
dfPrimers$tmR = round(as.numeric(dfPrimers$tmR), 2)
dfPrimers$posL = as.integer(sub(",[0-9]+$", "", dfPrimers$poslenL))
dfPrimers$posR = as.integer(sub(",[0-9]+$", "", dfPrimers$poslenR))
dfPrimers$lenL = as.integer(sub("^[0-9]+,", "", dfPrimers$poslenL))
dfPrimers$lenR = as.integer(sub("^[0-9]+,", "", dfPrimers$poslenR))
dfPrimers = dfPrimers[, !colnames(dfPrimers) %in% c("name", "poslenL", "poslenR")]
dfPrimers = dfPrimers[, c("seqL", "tmL", "posL", "lenL", "seqR", "tmR", "posR", "lenR")]
colnames(dfPrimers) = c("prmSeqL",  "prmTmL",  "prmPosL",  "prmLenL",
    "prmSeqR",  "prmTmR",  "prmPosR",  "prmLenR")
if (nrow(dfPrimers) != nrow(dfInDels))
    stop("Wrong number of primers returned by primer3")
dfInDels = cbind(dfInDels, dfPrimers, stringsAsFactors=FALSE)
rm(dfPrimers)

########################################
# There will be many markers for which primers could not be found.  Get rid of these.
########################################

{
catnow("Number of markers initially: ", nrow(dfInDels), "\n")
catnow("No left primer found:        ", sum(dfInDels$prmSeqL == "" & dfInDels$prmSeqR != ""), "\n")
catnow("No right primer found:       ", sum(dfInDels$prmSeqL != "" & dfInDels$prmSeqR == ""), "\n")
catnow("Neither left nor right found:", sum(dfInDels$prmSeqL == "" & dfInDels$prmSeqR == ""), "\n")
dfInDels = dfInDels[dfInDels$prmSeqL != "" & dfInDels$prmSeqR != "",]
catnow("Number of markers remaining: ", nrow(dfInDels), "\n")
}

########################################
# If no markers remain, exit with error.
########################################

if (nrow(dfInDels) == 0)
    stop("No markers remaining.")

########################################
# Add "phase" columns giving the genome strand which matches the reference genome
# "+" strand.  This column is "-" if a marker's region is running in the opposite
# direction as the reference genome.  For the reference genome it is always "+",
# so we will remove it from the data frame we write to a file.
########################################

phaseCol = makeColVec(paste(refGenomeLtr, "phase", sep=""))
for (genome in genomeLtrs)
    dfInDels[, phaseCol[genome]] = ifelse(dfInDels[, refStrand1Col] == dfInDels[, strand1Col[genome]], "+", "-")

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

dfInDels$kmer1offset = extensionLen - (dfInDels$prmPosL-1)
seqLen = 2*(2*extensionLen+kmerLen)
kmer2End = seqLen - extensionLen
dfInDels$kmer2offset = dfInDels$prmPosR - kmer2End

########################################
# Change the "*pos1" and "*pos2" columns, which give the position of the 5' end
# of each k-mer on the + strand, to columns "*ampPos1" and "*ampPos2", giving
# the positionn of the endn of the amplicon.
#
# The position change is different for each genome depending on whether the k-mer
# "phase" is + or -.  The reference genome is always "+" phase for all k-mers,
# while the other genomes have k-mers of both phases.  For a "+" phase genome,
# the pos1 values are adjusted in one manner, and the pos2 values are adjusted
# in a different manner, to give the ampPos1 and ampPos2 values.  For a "-" phase
# genome, the adjustment manners are switched.
# Adjustment for a "+" phase genome:
#   - ampPos1 = pos1 - kmer1offset
#   - ampPos2 = pos2 + kmer2offset + k-1
# Adjustment for a "-" phase genome:
#   - ampPos1 = pos1 + kmer1offset + k-1
#   - ampPos2 = pos2 - kmer2offset
########################################

ampPos1Col = makeColVec("ampPos1")
ampPos2Col = makeColVec("ampPos2")
refAmpPos1Col = ampPos1Col[refGenomeLtr]
refAmpPos2Col = ampPos2Col[refGenomeLtr]

for (genome in genomeLtrs)
    {
    posPhase = (dfInDels[, phaseCol[genome]] == "+")
    negPhase = !posPhase
    dfInDels[, ampPos1Col[genome]] = dfInDels[, pos1Col[genome]] - dfInDels$kmer1offset
    dfInDels[, ampPos2Col[genome]] = dfInDels[, pos2Col[genome]] - dfInDels$kmer2offset
    dfInDels[negPhase, ampPos1Col[genome]] =
        dfInDels[negPhase, pos1Col[genome]] + dfInDels[negPhase, "kmer1offset"] + (kmerLen-1)
    dfInDels[posPhase, ampPos2Col[genome]] =
        dfInDels[posPhase, pos2Col[genome]] + dfInDels[posPhase, "kmer2offset"] + (kmerLen-1)
    }

########################################
# Remove the *kkLen columns, which contain k-mer-to-k-mer lengths, and replace
# them with *ampLen columns, which contain amplicon lengths that are simply
# abs(ampPos1-ampPos2+1).
########################################

ampLenCol = makeColVec("ampLen")
refAmpLenCol = ampLenCol[refGenomeLtr]
for (genome in genomeLtrs)
    {
    dfInDels[, ampLenCol[genome]] = abs(dfInDels[, ampPos1Col[genome]] - dfInDels[, ampPos2Col[genome]]) + 1
    dfInDels = dfInDels[, colnames(dfInDels) != kkLenCol[genome]]
    }

########################################
# Add "dif" columns giving difference of amplicon length of each other genome
# to the reference genome.  There is no such column for the reference genome
# since it would always be 0, but columns exist for each other genome.
########################################

difToRefCol = makeColVec(paste(refGenomeLtr, "dif", sep=""))
for (genome in otherGenomeLtrs)
    dfInDels[, difToRefCol[genome]] = dfInDels[, ampLenCol[genome]] - dfInDels[, refAmpLenCol]

########################################
# Concatenate the k-mer strand signs of all genomes into a single column.
########################################

dfInDels$kmer1strands = apply(dfInDels[, strand1Col], 1, paste, collapse="")
dfInDels$kmer2strands = apply(dfInDels[, strand2Col], 1, paste, collapse="")

########################################
# Put the data in a nice order and write it to the output file.
########################################

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
dfInDels = dfInDels[, colOrder[]]
rownames(dfInDels) = NULL

########################################
# Write the dfInDels data frame to the output file.
########################################

write.table(dfInDels, tsvMarkerFile, col.names=TRUE, row.names=FALSE, quote=FALSE, sep="\t")
# dfInDels = read.table(tsvMarkerFile, header=TRUE, sep="\t", stringsAsFactors=FALSE)

catnow("Finished adding primer sequences to indels, candidate marker output file:\n", tsvMarkerFile, "\n")
}
