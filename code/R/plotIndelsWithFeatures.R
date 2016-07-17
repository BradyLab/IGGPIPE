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
#testing = 1 # For testing only.  .test
#testing = 2 # For testing only.  .HP14longITAG2.4
#testing = 3 # For testing only.  .CL13long
{
if (testing == 0)
    args = commandArgs(TRUE)
else if (testing == 1)
    {
    args = c("~/Documents/UCDavis/BradyLab/IGGPIPE/IGGPIPE",
        "outTestHP11/MarkersNonoverlappingWithInNearFeatures_K11k2L100D10_2000A100_2000d10_100N2F0X20V3000W8M3G1.indels.tsv",
        "outTestHP11/MarkersNonoverlappingWithInNearFeatures_K11k2L100D10_2000A100_2000d10_100N2F0X20V3000W8M3G1.indels.pdf",
        1000, 2000, '(^[^:]+):', "S.lycopersicum,S.pennellii")
    }
else if (testing == 2)
    {
    args = c("~/Documents/UCDavis/BradyLab/IGGPIPE/IGGPIPE",
        "outHP14/IndelGroupsNonoverlappingWithInNearFeatures_K14k4L100D1_3000A100_3000d100_100N2F0.indels.tsv",
        "outHP14/IndelGroupsNonoverlappingWithInNearFeatures_K14k4L100D1_3000A100_3000d100_100N2F0.indels.pdf",
        1000, 3000, '(^[^:]+):', "S.lycopersicum,S.pennellii")
    }
else if (testing == 3)
    {
    args = c("~/Documents/UCDavis/BradyLab/IGGPIPE/IGGPIPE",
        "outCL13/IndelGroupsNonoverlappingWithInNearFeatures_K13k4L100D1_3000A100_3000d100_100N2F0.indels.tsv",
        "outCL13/IndelGroupsNonoverlappingWithInNearFeatures_K13k4L100D1_3000A100_3000d100_100N2F0.indels.pdf",
        1000, 3000, '^[^-]+-([^-(]+)[-(]', "Col-0,Ler-0")
    }
else stop("Unknown value for 'testing'")
}

Nexpected = 7
if (length(args) != Nexpected)
    {
    usage = c(
        "Read a data frame of Indel information that has had gene feature information added to it, and",
        "make plots of various aspects of the data.",
        "See annotate/HPlongITAG2.4_addFeatureIDs.nonoverlappingIndelGroups for parameters used to create",
        "the input file for this.",
        "",
        "Usage: Rscript plotIndelsWithFeatures.R <wd> <IndelsFile> <pdfPlotFile> <nearInBp> <maxIndelLen> <featureRE> <genomeNames>",
        "",
        "Arguments:",
        "   <wd>                 : Path of R working directory, specify other file paths relative to this.",
        "   <IndelsFile>         : Name of file containing Indel position information with a column named 'isInNear'.",
        "                          See annotate/test_add_isInNearColumn.markers",
        "   <pdfPlotFile>        : Name of pdf file to be created.",
        "   <nearInBp>           : Number of base pairs used to define 'near' when running annotateFile.R",
        "   <maxIndelLen>        : Maximum possible size of any indel that is detected with the method used to detect them.",
        "   <featureRE>          : Perl-compatible regular expression that matches a feature name within the isInNear",
        "                          column of <IndelsFile>, but does not match any part of the file name, nor does it",
        "                          match the portion of the feature name that is typically a numeric suffix that is used",
        "                          to discriminate for example between intron1 and intron2.  The expression must have one",
        "                          parenthesized () subexpression that matches the feature name, and which is used to",
        "                          extract the feature names from the isInNear column.  The isInNear column value is",
        "                          split on comma and <featureRE> is applied to the resulting strings.",
        "   <genomeNames>        : Comma-separated names of genomes in the input file, same order as input file."
        )
    for (S in usage)
        catnow(S, "\n", sep="")
    stop("Try again with correct number of arguments")
    }

catnow("plotIndelsWithFeatures.R arguments:\n")
workingDirectory = args[1]
catnow("  workingDirectory: ", workingDirectory, "\n")
if (!dir.exists(workingDirectory))
    stop("Directory doesn't exist: ", workingDirectory)
setwd(workingDirectory)

IndelsFile = args[2]
catnow("  IndelsFile: ", IndelsFile, "\n")
if (!file.exists(IndelsFile))
    stop("File doesn't exist: ", IndelsFile)

pdfPlotFile = args[3]
catnow("  pdfPlotFile: ", pdfPlotFile, "\n")

nearInBp = as.integer(args[4])
catnow("  nearInBp: ", nearInBp, "\n")

maxIndelLen = as.integer(args[5])
catnow("  maxIndelLen: ", maxIndelLen, "\n")

featureRE = args[6]
catnow("  featureRE: ", featureRE, "\n")

genomeNames = unlist(strsplit(args[7], ",", fixed=TRUE))
catnow("  genomeNames: ", paste(genomeNames, collapse=","), "\n")

########################################
# Initialize.
########################################

# Let's not plot grid lines.
plotGridLines = FALSE

# Read Indels data.
dfIndels = read.table(IndelsFile, header=TRUE, row.names=NULL, sep="\t", stringsAsFactors=FALSE)
if (nrow(dfIndels) == 0)
    stop("There are no Indels.")

# Get the genome identifying letters from the dfIndels id columns names.
genomeLtrs = sub("id", "", colnames(dfIndels)[grepl("id$", colnames(dfIndels))])
N.genomes = length(genomeLtrs)
if (N.genomes < 2)
    stop("Expected to recognize at least two genomes in the data column names")
if (N.genomes != length(genomeNames))
    stop("Expected ", N.genomes, " genome names but only got ", length(genomeNames), " in the argument list")
names(genomeNames) = genomeLtrs
refGenomeLtr = genomeLtrs[1]
otherGenomeLtrs = genomeLtrs[-1]

# Get all pairs of genome letters.
genomePairs = t(combn(genomeLtrs, 2))
Npairs = nrow(genomePairs)
cat("Npairs = ", Npairs, "\n")

makeColVec = function(S)
    {
    V = paste(genomeLtrs, S, sep="")
    names(V) = genomeLtrs
    return(V)
    }

delCol = makeColVec("del")
refDelCol = delCol[refGenomeLtr]
idCol = makeColVec("id")
refIdCol = idCol[refGenomeLtr]
startCol = makeColVec("start")
refStartCol = startCol[refGenomeLtr]
endCol = makeColVec("end")
refEndCol = endCol[refGenomeLtr]

catnow("Number of Indels:", nrow(dfIndels), "\n")

# Add a "len" column that gives length of the insertion when delCol is 0.  The
# length will be 0.
lenCol = makeColVec("len")
refLenCol = lenCol[refGenomeLtr]
# Note that it is end-start-1, not end-start+1, see alignAndGetIndelsSNPs.R comments.
for (genome in genomeLtrs)
    dfIndels[, lenCol[genome]] = dfIndels[, endCol[genome]] - dfIndels[, startCol[genome]] - 1

########################################
# Simplify things.  For each row of dfIndels that overlaps more than one feature,
# split and duplicate each such row.
########################################
whichMultiple = grepl(",", dfIndels$isInNear)
dfIndels.Multiple = dfIndels[whichMultiple,,drop=FALSE]
dfIndels.Single = dfIndels[!whichMultiple,,drop=FALSE]
isInNear = strsplit(dfIndels.Multiple$isInNear, ",", fixed=TRUE)
lens = sapply(isInNear, length)
idxs = sapply(1:length(lens), function(i) rep(i, lens[i]))
dfIndels.Duped = dfIndels.Multiple[unlist(idxs, use.names=FALSE),]
dfIndels.Duped$isInNear = unlist(isInNear, use.names=FALSE)
dfIndels = rbind(dfIndels.Single, dfIndels.Duped)
rm(whichMultiple, dfIndels.Multiple, dfIndels.Single, isInNear, lens, idxs, dfIndels.Duped)

########################################
# Parse the isInNear column.  We don't care about the gene ID, discard it, but
# part of the ID field includes the type of feature (intron, etc.) and we do
# want that.  It is not predictable how to find the feature name.  Examples:
#   Soly: "intron:Solyc00g006830.2.1.2(@-1344:-)"
#   Arabidopsis: "AT1G12938.1-CDS-1(-70:+)"
# That's why we have the <featureRE> argument.
# Get the strand in featStrand, and the flag indicating overlap (@), upstream (-),
# or downstream(+).  We don't care about the amount of overlap, but get the amount
# upstream or downstream.
########################################

# Extract the feature name into new column "feature" using featureRE.
m = regexec(featureRE, dfIndels$isInNear)
hasMatch = sapply(m, function(V) V[1] != -1)
startPos = sapply(m[hasMatch], function(V) V[2])
matchLen = sapply(m[hasMatch], function(V) attr(V, "match.length")[2])
stopPos = startPos + matchLen - 1
features = substr(dfIndels$isInNear[hasMatch], startPos, stopPos)
u.features = unique(features)
N.features = length(u.features)
if (N.features < 2)
    stop("There were ", N.features, " unique feature types found using <featureRE>, and this is too few")
if (N.features > 20)
    stop("There were ", N.features, " unique feature types found using <featureRE>,",
        "and this is too many, first few were:\n", paste(u.features[1:20], collapse=" "), "\n")
dfIndels$feature = ""
dfIndels[hasMatch, "feature"] = features
rm(m, hasMatch, startPos, matchLen, stopPos, features)

# Extract the strand, + or -, into new column "featStrand".
dfIndels$featStrand = sub("^[^(]*\\([^:]*:([^)]*).*$", "\\1", dfIndels$isInNear)

# Extract the flag character after "(", either @, +, or -, into new column "flag".
dfIndels$flag = sub("^[^(]*\\((.).*$", "\\1", dfIndels$isInNear)

# Remove columns we don't need.
dfIndels = dfIndels[, c("feature", "flag", "featStrand", delCol, lenCol)]

# Make some clean-up edits to the feature names.
dfIndels$feature = sub("five_prime_UTR", "5'UTR", dfIndels$feature)
dfIndels$feature = sub("5UTR", "5'UTR", dfIndels$feature)
dfIndels$feature = sub("three_prime_UTR", "3'UTR", dfIndels$feature)
dfIndels$feature = sub("3UTR", "3'UTR", dfIndels$feature)

# Make sure some feature names are as we expect.
u.features = unique(dfIndels$feature)
if (!any(u.features == "CDS"))
    stop("There were no CDS features and we expected some")
if (!any(u.features == "5'UTR"))
    stop("There were no 5'UTR features and we expected some")
if (!any(u.features == "3'UTR"))
    stop("There were no 3'UTR features and we expected some")

# If there was no feature that was near the indel, it is in an intergenic region.
dfIndels$feature[dfIndels$feature == ""] = "intergenic"

# We need to add "upstream" and "downstream" features.  Any indel whose feature is CDS,
# 5'UTR, or 3'UTR and has a flag that is "+" or "-" rather than "@" indicates an indel
# lying upstream or downstream of the feature.  If the feature lies on the "+" strand,
# then a "+" flag means downstream and a "-" flag means upstream.  If the feature lies
# on the "-" strand, then a "+" flag means upstream and a "-" flag means downstream.
# (Double-check that!)
isNear_CDS = (dfIndels$feature == "CDS") & (dfIndels$flag != "@")
isNear_5UTR = (dfIndels$feature == "5'UTR") & (dfIndels$flag != "@")
isNear_3UTR = (dfIndels$feature == "3'UTR") & (dfIndels$flag != "@")
isNear_CDS_UTR = isNear_CDS | isNear_5UTR | isNear_3UTR
strand = dfIndels$featStrand[isNear_CDS_UTR]
flag = dfIndels$flag[isNear_CDS_UTR]
isUpstream = (strand == flag)
dfIndels$feature[isNear_CDS_UTR][isUpstream] = "upstream"
dfIndels$feature[isNear_CDS_UTR][!isUpstream] = "downstream"
#table(dfIndels$feature[isNear_CDS]) # Yes!
#table(dfIndels$feature[isNear_5UTR]) # Yes!
#table(dfIndels$feature[isNear_3UTR]) # Yes!

u.features = unique(dfIndels$feature)
N.features = length(u.features)
ourFeatures = c("upstream", "5'UTR", "CDS", "intron", "3'UTR", "downstream", "intergenic")
skipped = (!u.features %in% ourFeatures)
if (any(skipped))
    stop("Features not recognized and skipped: ", paste(u.features[skipped], collapse=" "))
missing = (!ourFeatures %in% u.features)
if (any(missing))
    stop("Missing features that we expected to see: ", paste(ourFeatures[missing], collapse=" "))

########################################
# Plot colors for each genome and each feature.
########################################

genomeCols = colorBlind.8[1:N.genomes]
if (N.genomes > 8)
    genomeCols = rainbow(N.genomes)
names(genomeCols) = genomeLtrs

featureCols = colorBlind.8[1:N.features]
if (N.features > 8)
    featureCols = rainbow(N.features)
names(featureCols) = ourFeatures

########################################
# Grid line color.
########################################

gridCol = "gray"

########################################
# Line types for each genome and each feature.
########################################

ltys = c("solid", "dashed", "dotted", "dotdash", "longdash", "twodash")
genomeLtys = rep(ltys, length.out=N.genomes)
names(genomeLtys) = genomeLtrs
featureLtys = rep(ltys, length.out=N.features)
names(featureLtys) = ourFeatures

########################################
# Create pdf file.
########################################

pdf(pdfPlotFile, height=8, width=6)
par(mar=c(9,8,6,3))

########################################
# Make a bar plot (one bar per genome) showing the number of indels occurring in
# each feature type.
########################################

# Count number of indels for each feature.
countMtx = matrix(0, ncol=N.features, nrow=N.genomes, dimnames=list(genomeLtrs, ourFeatures))

for (genome in genomeLtrs)
    {
    # Get counts of deletions in genome 'genome' within each feature.
    deletions = (dfIndels[,delCol[genome]] > 0)
    feats = dfIndels[deletions, "feature"]
    x = table(feats)
    missing = is.na(x[ourFeatures])
    if (any(missing))
        x[ourFeatures[missing]] = 0
    countMtx[genome,] = x[ourFeatures]

    # Double-check that #insertions in OTHER genome equals number of deletions
    # in THIS genome.  This test showed a problem with alignAndGetIndelsSNPs initially.
    otherGenome = ifelse(genome == genomeLtrs[1], genomeLtrs[2], genomeLtrs[1])
    insertions = (dfIndels[,delCol[otherGenome]] == 0)
    feats = dfIndels[insertions, "feature"]
    x = table(feats)
    missing = is.na(x[ourFeatures])
    if (any(missing))
        x[ourFeatures[missing]] = 0
    if (any(countMtx[genome,] != x[ourFeatures]))
        stop("Expected insertions in one genome to equal deletions in other (1)")
    }

# Compute axis labels.
ypretty = pretty.log(c(countMtx))
ylim = range(ypretty)
yat = ypretty
ylab = "Number of Deletions or Insertions"
if (ypretty[1] >= 1000)
    {
    ypretty = ypretty / 1000
    ylab = paste(ylab, "  (x1000)", sep="")
    }
ypretty = format(ypretty, scientific=FALSE)

# Make the bar plot.
barplot(countMtx, beside=TRUE, col=genomeCols, axes=FALSE, log="y", las=2, cex.names=1.5,
    main="Distribution of indel regions",
    ylab="")
title(paste("upstream/downstream defined as within ", nearInBp, " bp of 5'UTR/CDS/3'UTR", sep=""), line=1.5, cex.main=0.8)
axis(2, at=yat, labels=ypretty, las=2, cex.axis=1.5)
mtext(ylab, 2, line=6, cex=1.5)

# Legend.
txt = paste(genomeNames, " deletion, ", genomeNames[2:1], " insertion", sep="")
legend("topleft", txt, fill=genomeCols, cex=1)

########################################
# Make a line plot (one line per genome) where x-axis positions where points are
# plotted are bins containing a range of indel deletion lengths, and the y-axis
# is the number of indels with a deletion length in that range.
########################################

# Indel length bins.
minIndelLen = c(1, 5, 25, 100, 250, 1000)
Nbins = length(minIndelLen)
binRanges = paste(minIndelLen, "..", c((minIndelLen-1)[-1], "+"), sep="")
binRanges = sub("0..0", "0", binRanges)
binRanges = sub("1000..+", "1000+", binRanges)

# Count number of indels in each indel size range.
countMtx = matrix(0, ncol=Nbins, nrow=N.genomes, dimnames=list(genomeLtrs, minIndelLen))
for (genome in genomeLtrs)
    {
    # Deletions.
    deletions = (dfIndels[,delCol[genome]] > 0)
    x = findInterval(dfIndels[deletions, delCol[genome]], minIndelLen)
    x = table(x)
    missing = setdiff(1:Nbins, as.integer(names(x)))
    if (length(missing) > 0)
        x[as.character(missing)] = 0
    names(x) = as.character(minIndelLen[as.integer(names(x))])
    countMtx[genome,] = x[as.character(minIndelLen)]

    # Double-check that #insertions in OTHER genome equals number of deletions
    # in THIS genome.
    otherGenome = ifelse(genome == genomeLtrs[1], genomeLtrs[2], genomeLtrs[1])
    insertions = (dfIndels[,delCol[otherGenome]] == 0)
    x = findInterval(dfIndels[insertions, lenCol[otherGenome]], minIndelLen)
    x = table(x)
    missing = setdiff(1:Nbins, as.integer(names(x)))
    if (length(missing) > 0)
        x[as.character(missing)] = 0
    names(x) = as.character(minIndelLen[as.integer(names(x))])
    if (any(countMtx[genome,] != x[as.character(minIndelLen)]))
        stop("Expected insertions in one genome to equal deletions in other (2)")
    }

# Compute axis labels.
xlim = c(1, Nbins)
xlab = paste("Indel size range (bp)\nmax detectable ~ ", maxIndelLen, " bp", sep="")
ypretty = pretty.log(c(countMtx))
ylim = range(ypretty)
yat = ypretty
ylab = "Number of Deletions or Insertions"
if (ypretty[1] >= 1000)
    {
    ypretty = ypretty / 1000
    ylab = paste(ylab, "  (x1000)", sep="")
    }
ypretty = format(ypretty, scientific=FALSE)

# Start the plot.
plot(NA, type="n", xlim=xlim, ylim=ylim, axes=FALSE, log="y", xlab="", ylab="",
    main="Distribution of indel sizes")
mtext(binRanges, 1, at=1:Nbins, line=0, las=2, cex=1.5)
mtext(xlab, 1, line=7.5, cex=1.5)
axis(2, at=yat, labels=ypretty, las=2, cex.axis=1.5)
mtext(ylab, 2, line=5, cex=1.5)

# Loop and plot line/points for each genome.
for (genome in genomeLtrs)
    {
    lines(1:Nbins, countMtx[genome,], col=genomeCols[genome], lwd=3, lty=genomeLtys[genome])
    points(1:Nbins, countMtx[genome,], col=genomeCols[genome], pch=20, cex=2)
    }

# Make some grid lines.
if (plotGridLines)
    {
    segments(1:Nbins, ylim[1], 1:Nbins, ylim[2], col=gridCol)
    segments(0, yat, xlim[2], yat, col=gridCol)
    }

# Legend.
txt = paste(genomeNames, " deletion, ", genomeNames[2:1], " insertion", sep="")
legend("topright", txt, pch=20, col=genomeCols, pt.cex=2, lwd=3, lty=genomeLtys,
    seg.len=5, bty="n", cex=1)

########################################
# Similar to the previous, except show a line for each feature type, and make a
# separate plot for each genome.
########################################

# Count number of indels for each feature within each indel size range.
counts = array(0, c(N.genomes, nrow=N.features, Nbins), dimnames=list(genomeLtrs, ourFeatures, minIndelLen))
for (genome in genomeLtrs)
    {
    deletions = (dfIndels[,delCol[genome]] > 0)

    for (feat in ourFeatures)
        {
        x = findInterval(dfIndels[dfIndels$feature == feat & deletions, delCol[genome]], minIndelLen)
        x = table(x)
        missing = setdiff(1:Nbins, as.integer(names(x)))
        if (length(missing) > 0)
            x[as.character(missing)] = 0
        names(x) = as.character(minIndelLen[as.integer(names(x))])
        counts[genome, feat,] = x[as.character(minIndelLen)]
        }
    }

# Compute axis labels.
xlim = c(1, Nbins)
xlab = paste("Indel size range (bp)\nmax detectable ~ ", maxIndelLen, " bp", sep="")
ypretty = pretty.log(1+c(counts))
ylim = range(ypretty)
yat = ypretty
ylab = "Number of Deletions or Insertions"
if (ypretty[1] >= 1000)
    {
    ypretty = ypretty / 1000
    ylab = paste(ylab, "  (x1000)", sep="")
    }
ypretty = format(ypretty, scientific=FALSE)

# Loop and plot feature lines/points for each genome.
for (genome in genomeLtrs)
    {
    otherGenome = ifelse(genome == genomeLtrs[1], genomeLtrs[2], genomeLtrs[1])

    # Start the plot.
    plot(NA, type="n", xlim=xlim, ylim=ylim, axes=FALSE, log="y", xlab="", ylab="",
        main=paste("Distribution of indel sizes\nDeletions in ",
            genomeNames[genome], ", insertions in ", genomeNames[otherGenome], sep=""))
    title(paste("upstream/downstream defined as within ", nearInBp, " bp of 5'UTR/CDS/3'UTR", sep=""), line=0.5, cex.main=0.8)
    mtext(binRanges, 1, at=1:Nbins, line=0, las=2, cex=1.5)
    mtext(xlab, 1, line=7.5, cex=1.5)
    axis(2, at=yat, labels=ypretty, las=2)
    mtext(ylab, 2, line=4, cex=1.5)

    # Loop for each feature type.
    for (feat in ourFeatures)
        {
        lines(1:Nbins, counts[genome, feat,], col=featureCols[feat], lwd=2, lty=featureLtys[feat])
        points(1:Nbins, counts[genome, feat,], col=featureCols[feat], pch=20, cex=2)
        }

    # Label points left of x=1.  When distance is close in y-direction, spread them.
    if (FALSE)
        {
        y = counts[genome, ourFeatures, 1]
        txt = ourFeatures
        ord = order(y)
        y = y[ord]
        txt = txt[ord]
        near = (diff(y)/y[-1] < 0.2)
        while (any(near))
            {
            y[c(FALSE,near)] = y[c(FALSE,near)]*1.05
            near = (diff(y)/y[-1] < 0.2)
            }
        text(1, y, txt, pos=2, cex=0.8)
        }

    # Make some grid lines.
    if (plotGridLines)
        {
        segments(1:Nbins, ylim[1], 1:Nbins, ylim[2], col=gridCol)
        segments(0, yat, xlim[2], yat, col=gridCol)
        }

    # Legend.
    legend("bottomright", ourFeatures, pch=20, col=featureCols, pt.cex=2, lwd=2, lty=featureLtys, seg.len=5, bty="n")
    }

########################################
# Close pdf file.
########################################

dev.off()
catnow("Finished making plots of Indel information, output file:\n", pdfPlotFile, "\n")
}

################################################################################
# End of file.
################################################################################
