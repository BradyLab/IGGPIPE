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

# Get arguments.
testing = 0
#testing = 1 # For testing only, outTestHP11, genome 1.
#testing = 2 # For testing only, outTestHP11, genome 2.
{
if (testing == 0)
    args = commandArgs(TRUE)
else if (testing == 1)
    {
    args = c("~/Documents/UCDavis/BradyLab/Genomes/IGGPIPE",
        "outTestHP11/IndelGroupsOverlapping_K11k2L100D10_2000A100_2000d10_100N2F0.tsv",
        1, "outTestHP11/GenomeData/IndelGroupsOverlapping_K11k2L100D10_2000A100_2000d10_100N2F0X20_1.dnaseqs",
        "outTestHP11/GenomeData", 20,
        "/Users/tedtoal/perl5/perlbrew/perls/perl-5.14.2/bin/perl",
        "code/perl/getSeqsFromFasta.pl", "testFASTA/ITAG2.4_test.fasta",
        "outTestHP11/GenomeData/Genome_1.contigs", TRUE)
    }
else if (testing == 2)
    {
    args = c("~/Documents/UCDavis/BradyLab/Genomes/IGGPIPE",
        "outTestHP11/IndelGroupsOverlapping_K11Km2Lm100Dm10Dx2000Am100Ax2000ADm10ADx100ND2mF0.tsv",
        2, "outTestHP11/GenomeData/IndelGroupsOverlapping_K11Km2Lm100Dm10Dx2000Am100Ax2000ADm10ADx100ND2mF0X20_2.dnaseqs",
        "outTestHP11/GenomeData", 20,
        "/Users/tedtoal/perl5/perlbrew/perls/perl-5.14.2/bin/perl",
        "code/perl/getSeqsFromFasta.pl", "testFASTA/SpennV2.0_test.fasta",
        "outTestHP11/GenomeData/Genome_2.contigs", TRUE)
    }
else stop("Unknown value for 'testing'")
}

Nexpected = 11
if (length(args) != Nexpected)
    {
    usage = c(
        "Read a data frame of k-mer pairs that have been identified as bounding Indel Groups,",
        "then read a genomic FASTA file for one genome and extract DNA sequences at each",
        "k-mer position, and finally, write the resulting data frame to a new file.",
        "",
        "Usage: Rscript getDNAseqsForPrimers.R <wd> <indelFile> <genomeNum> <seqOutFile> <extractDir> \\",
        "       <extensionLen> <perlPath> <getSeqsFromFasta> <fastaFile> <contigsFile> <investigate>",
        "",
        "Arguments:",
        "   <wd> : Path of R working directory, specify other file paths relative to this.",
        "   <indelFile>    : Input file containing the indel k-mer pairs.",
        "   <genomeNum>    : The genome number whose DNA sequences are to be retrieved, 1=reference genome, etc.",
        "   <seqOutFile>   : Output file to be written containing all unique indel k-mer pairs with DNA sequences.",
        "   <extractDir>   : Directory to hold DNA sequence extraction files.",
        "   <extensionLen> : Number of base pairs of additional sequence to get on each side of k-mer (e.g. 10)",
        "   <perlPath>     : Full path of the Perl language interpreter program",
        "   <getSeqsFromFasta> : Full path of the Perl script getSeqsFromFasta.pl",
        "   <fastaFile>    : FASTA file containing genomic sequences for genome 1.",
        "   <contigsFile>  : File of contig lengths/positions for genome 1.",
        "   <investigate>  : FALSE for normal operation, TRUE for more verbose debugging output."
        )
    for (S in usage)
        catnow(S, "\n", sep="")
    stop("Try again with correct number of arguments")
    }

catnow("getDNAseqsForPrimers.R arguments:\n")
workingDirectory = args[1]
catnow("  workingDirectory: ", workingDirectory, "\n")
if (!dir.exists(workingDirectory))
    stop("Directory doesn't exist: ", workingDirectory)
setwd(workingDirectory)

indelFile = args[2]
catnow("  indelFile: ", indelFile, "\n")
if (!file.exists(indelFile))
    stop("File doesn't exist: ", indelFile)

genomeNum = as.integer(args[3])
catnow("  genomeNum: ", genomeNum, "\n")
if (genomeNum < 1)
    stop("genomeNum must be >= 1")

seqOutFile = args[4]
catnow("  seqOutFile: ", seqOutFile, "\n")

extractDir = args[5]
catnow("  extractDir: ", extractDir, "\n")
if (!dir.exists(extractDir))
    stop("Directory doesn't exist: ", extractDir)

extensionLen = as.integer(args[6])
catnow("  extensionLen: ", extensionLen, "\n")
if (is.na(extensionLen) || extensionLen < 0 || extensionLen > 50)
    stop("extensionLen must be between 0 and 50")

perlPath = args[7]
catnow("  perlPath: ", perlPath, "\n")

getSeqsFromFasta = args[8]
catnow("  getSeqsFromFasta: ", getSeqsFromFasta, "\n")
if (!file.exists(getSeqsFromFasta))
    stop("File doesn't exist: ", getSeqsFromFasta)

fastaFile = args[9]
catnow("  fastaFile: ", fastaFile, "\n")
if (!file.exists(fastaFile))
    stop("File doesn't exist: ", fastaFile)

contigsFile = args[10]
catnow("  contigsFile: ", contigsFile, "\n")
if (!file.exists(contigsFile))
    stop("File doesn't exist: ", contigsFile)

investigate = as.logical(args[11])
catnow("  investigate: ", investigate, "\n")
if (is.na(investigate))
    stop("investigate must be TRUE or FALSE")

########################################
# Initialize.
########################################

catnow("Preparing to retrieve k-mer-vicinity sequences from FASTA file\n")

# Read indel k-mer pairs data.
df = read.table(indelFile, header=TRUE, row.names=NULL, sep="\t", stringsAsFactors=FALSE)
catnow("Number of input file Indel Groups read:", nrow(df), "\n")
if (nrow(df) == 0)
    stop("There are no Indel Groups.")
inv(dim(df), "input data dim")
inv(colnames(df), "input data columns")
inv(head(df), "input data head")

# The k-mer length we are working with.
kmerLen = nchar(df$kmer1[1])
inv(kmerLen, "k")

# Get the genome letter for this genome.
genomeLtr = substr(colnames(df)[grepl("^.id$", colnames(df))], 1, 1)
if (length(genomeLtr) < 2)
    stop("Expected at least two id columns in <indelFile>")
refGenomeLtr = genomeLtr[1]
if (genomeNum > length(genomeLtr))
    stop("genome number is larger than the number of genomes in <indelFile>")
genomeLtr = genomeLtr[genomeNum]
inv(genomeLtr, "genome letter")

############################################################################
# Create a text file used as a sequence specification file for getSeqsFromFasta.pl
# to extract an extended DNA sequence around each k-mer.  You might think we
# would extract the entire sequence between k-mer 1 and 2, but we don't because
# this sequence is by definition DIFFERENT between each genome, so we can't make
# primers using a single sequence like that.  Instead, we only extract sequence
# around EACH of the two k-mers, and make primers using the concatenated sequence
# with no intervening sequence.
#
# We will extract the DNA sequence of the k-mer itself plus "extensionLen" base
# pairs on each side of it, UNLESS that extends beyond the contig containing the
# k-mer (in which case the extraction is reduced to extend only to the contig
# boundary).
#
# The extraction direction must be properly chosen to be the one that is
# appropriate, i.e. so that the sequences of the other genomes match the one
# extracted for the reference genome (which is always the + strand sequence).
# Each k-mer position was previously adjusted to be the + strand position of the
# left-most (most 5') k-mer base pair.  Also, the reference strand k-mer #1 was
# chosen as the one at the lower position, and k-mer #2 at higher position.
# Given a k-mer pair, the sequence to extract for the reference genome is the
# + strand sequence for both k-mers (regardless of which strand the k-mers
# themselves are on).  For the other genomes, the extraction strand depends on
# whether their k-mers are in phase or out of phase with the reference genome
# k-mers.  If the k-mer 1 strand of the other genome equals the k-mer 1 strand
# of the reference genome, then the other genome should also be extracted on the
# + strand.  If not, it should be on the minus strand.  Note that since the
# k-mer positions are the positions relative to the 5' end of the + strand, and
# since extensionLen base pairs are extracted on BOTH sides of the k-mer, even
# if the extraction is on the minus strand, the positions for the extraction are
# kmerPos-extensionLen thru kmerPos+k+extensionLen.
#
# The extraction direction is specified to getSeqsFromFasta.pl by including (for
# the - strand) or not (for the + strand) a "!" character in front of the
# sequence position specification, which determines whether or not it reverse-
# complements the extracted + strand sequence.
#
# Many Indel Groups will share the same left-side or same right-side k-mer.
# To reduce time required, we will not extract the same area multiple times.  We
# will find all the unique left-side and right-side k-mers and extract sequences
# only for those.
############################################################################

############################################################################
# Keep only the data for genome genomeNum plus the k-mer and the reference
# genome strand columns.  Strip the genome letter from the start of the
# column name, and copy the reference genome strand columns to column names
# "refS1" and "refS2".
############################################################################

df$refS1 = df[, paste(refGenomeLtr, "s1", sep="")]
df$refS2 = df[, paste(refGenomeLtr, "s2", sep="")]
RE = paste("^(", genomeLtr, "|kmer|refS)", sep="")
df = df[, grepl(RE, colnames(df))]
colnames(df) = sub(paste("^", genomeLtr, sep=""), "", colnames(df))

############################################################################
# Split df into two parts: the left-side k-mer portion (suffix 1) and the
# right-side k-mer portion (suffix 2).  Leave out the "kkLen" and "pct"
# columns, we don't need them.  Make the column names of df1 and df2 be
# the same.  Then, row-bind the two data frames.
############################################################################

df1 = df[, grepl("(1|id)$", colnames(df))]
df2 = df[, grepl("(2|id)$", colnames(df))]
rm(df)
colnames(df1) = sub("1$", "", colnames(df1))
colnames(df2) = sub("2$", "", colnames(df2))
df = rbind(df1, df2)
rm(df1, df2)

############################################################################
# Remove duplicate k-mer rows.  The entire row is identical in each case.
# We don't care whether the duplicates occurred in kmer1 or kmer2 columns
# of the original data frame, the extracted sequence is the same in either
# case.
############################################################################

inv(nrow(df), "before removing duplicate k-mer rows") 
df = df[!duplicated(df$kmer),]
inv(nrow(df), "after removing duplicate k-mer rows") 

############################################################################
# Create an extraction positions data frame dfXP.
# Columns:
#   pos.left, pos.right: left and right genome positions for k-mer extraction
#   pos.recomp: whether to extract the - strand, i.e. reverse-complement the
#                + strand sequence
############################################################################

refKmerStrand = (df$refS == "+")

dfXP = data.frame(pos.left=df$pos, stringsAsFactors=FALSE)
dfXP$pos.right = dfXP$pos.left + kmerLen-1
dfXP$pos.left = dfXP$pos.left - extensionLen
dfXP$pos.right = dfXP$pos.right + extensionLen
# The following also works correctly when our genome is the reference genome.
# In that case, it sets pos.revcomp FALSE for the entire vector.
dfXP$pos.revcomp = (df$s == "-")
dfXP$pos.revcomp[!refKmerStrand] = !dfXP$pos.revcomp[!refKmerStrand]

########################################
# Read the contigs file, we need to know the starting position and how long each contig is.
########################################

dfContigs = read.table(contigsFile, header=TRUE, sep="\t", stringsAsFactors=FALSE)
rownames(dfContigs) = paste(dfContigs$seqID, dfContigs$contig, sep="_")
dfContigs = dfContigs[, c("pos", "len")]

########################################
# Make adjustments to the pos.left/pos.right positions in the case where they lie
# off the ends of the contig.
########################################

# Get the contig start and end positions for the contigs in df, which are the
# same ones as in dfXP.
ctgID = paste(df$id, df$ctg, sep="_")
ctgStart = dfContigs[ctgID, "pos"]
ctgEnd = ctgStart + dfContigs[ctgID, "len"] - 1

# FIX LEFT
# It is possible (if the k-mer lies near the start of the contig) for .left
# values computed above to be left of the start of the contig, and if so, we
# need to move these .left values to the position of the start of the contig.
# We will save the amount of adjustment to each position in variables Hpos1.Lcut
# etc., for later adjustment of retrieved sequences (N's will be added to those
# whose .left was adjusted so they will all be the same length).

dfXP$pos.Lcut = (ctgStart - dfXP$pos.left) # Distance that .left lies left of start of contig
dfXP$pos.Lcut[dfXP$pos.Lcut < 0] = 0
dfXP$pos.left = dfXP$pos.left + dfXP$pos.Lcut

# FIX RIGHT
# Likewise it is possible (if the k-mer lies near the end of the contig)
# for .right values to be beyond the length of contig (leading into N's or
# beyond the actual sequence length), and if so, we need to move these .right
# values to the last valid base character position in the contig.

dfXP$pos.Rcut = (dfXP$pos.right - ctgEnd) # Distance that .right lies right of end of contig
dfXP$pos.Rcut[dfXP$pos.Rcut < 0] = 0
dfXP$pos.right = dfXP$pos.right - dfXP$pos.Rcut

########################################
# Create a sequence extraction positions data file for use by getSeqsFromFasta.pl.
#
# The format of a sequence extraction string (an argument to getSeqsFromFasta.pl)
# is id:left..right and this can be preceded by a "!" to reverse complement the
# extracted sequence.
#
# Put sequence extraction strings in vector seqExtStrs and write it to a text file.
#
# A command line to run getSeqsFromFasta.pl for hypothetical genome G looks like this:
#   perl blah/getSeqsFromFasta.pl blah/genome.fasta -i extractDir/extractG.txt -o extractDir/seqsG.txt
########################################

if (!file.exists(extractDir))
    dir.create(extractDir)

seqExtStrs = paste(c("", "!")[1+dfXP$pos.revcomp], df$id,
    ":", as.integer(dfXP$pos.left),
    "..", as.integer(dfXP$pos.right), sep="")

extractPosFile = paste(extractDir, paste("extract_KmerRgns_", genomeNum, ".txt", sep=""), sep=PATHSEP)
writeLines.winSafe(seqExtStrs, extractPosFile)

seqFile = paste(extractDir, paste("seqs_KmerRgns_", genomeNum, ".txt", sep=""), sep=PATHSEP)
getSeqs_args = c(getSeqsFromFasta, fastaFile, "-i", extractPosFile, "-o", seqFile)

inv(length(seqExtStrs), "length(seqExtStrs)")
inv(paste(getSeqs_args, collapse=" "), "Perl getSeqsFromFasta arguments")
inv(seqFile, "DNA sequence file")

########################################
# Now run the command.  This can take a long time, because the genome sequence
# file is very large, and if the sequence IDs are entire chromosomes,
# getSeqsFromFasta.pl reads an entire sequence into memory before processing it.
########################################

catnow("Running command to retrieve k-mer-vicinity sequences from FASTA file:\n")
catnow("   ", perlPath, " ", paste(getSeqs_args, collapse=" "), "\n")
stat = system2(perlPath, getSeqs_args)
if (stat != 0)
    stop("Perl program ", getSeqsFromFasta, " exited with error status ", stat)

########################################
# Read the output file, extract the sequence extraction strings from the
# sequence ID lines, and use the sequence extraction command strings as the
# names of the extracted sequences themselves.  Then, make sure that the
# sequence file has the expected sequences in it.  Save the sequences in the
# df data frame, using column name "seq".
########################################

# Extract sequences and make sure all are present.
catnow("Processing sequence data\n")
seqs = readLines(seqFile)
dft = as.data.frame(matrix(seqs, ncol=2, byrow=TRUE), stringsAsFactors=FALSE)
colnames(dft) = c("seqExtStr", "seq")
dft$seqExtStr = sub("^>([^ ]+) revcomp:no sub_region(:[0-9]+)-([0-9]+).*$", "\\1\\2..\\3", dft$seqExtStr)
dft$seqExtStr = sub("^>([^ ]+) revcomp:yes sub_region(:[0-9]+)-([0-9]+).*$", "!\\1\\2..\\3", dft$seqExtStr)
seqs = toupper(dft$seq)
names(seqs) = dft$seqExtStr
numMissing = sum(!seqExtStrs %in% dft$seqExtStr)
if (numMissing > 0)
    stop(numMissing, " sequences are missing from the FASTA extraction of genome ", genomeLtr,
        ", check regular expression patterns in sub() calls to assign to dft$seqExtStr")
seqs = seqs[seqExtStrs]

# Append the sequences to the df data frame.
df = data.frame(df, seq=seqs, stringsAsFactors=FALSE)
rownames(df) = NULL

########################################
# We must deal with cases where the sequence lengths are different because
# one or the other or both sides of an extracted sequence ran up against the
# left or right edge of the contig.  We will do this by adding "N" characters
# to sequences where they were truncated, to make them all the same length.
# This is where we make use of the earlier-saved pos.Lcut/pos.Rcut values.
########################################

Ns = sapply(1:extensionLen, function(i) paste(rep("N", i), collapse=""))
adjSeqs = function(seqs, amts, LR)
    {
    adjWhich = which(amts > 0) # It won't be many
    amts = amts[adjWhich]
    if (LR == "left")
        seqs[adjWhich] = paste(Ns[amts], seqs[adjWhich], sep="")
    else
        seqs[adjWhich] = paste(seqs[adjWhich], Ns[amts], sep="")
    return(seqs)
    }

df$seq = adjSeqs(df$seq, dfXP$pos.Lcut, "left")
df$seq = adjSeqs(df$seq, dfXP$pos.Rcut, "right")

########################################
# Save the output file.
########################################

# Remove refS and ctg columns.
df = df[, colnames(df) != "refS"]
df = df[, colnames(df) != "ctg"]

# Write the data frame to the output file.
write.table.winSafe(df, seqOutFile, row.names=FALSE, quote=FALSE, sep="\t")
# df = read.table(seqOutFile, header=TRUE, row.names=NULL, sep="\t", stringsAsFactors=FALSE)
catnow("Number of DNA sequences written to output file:", nrow(df), "\n")
catnow("Finished retrieving DNA sequences for Indel Groups, output file:\n", seqOutFile, "\n")
}

################################################################################
# End of file.
################################################################################
