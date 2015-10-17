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
#testing = 1 # For testing only.  .test
#testing = 2 # For testing only.  .HPlongITAG2.4
{
if (testing == 0)
    args = commandArgs(TRUE)
else if (testing == 1)
    {
    args = c("~/Documents/UCDavis/BradyLab/Genomes/IGGPIPE",
        "outTestHP11/MarkersNonoverlappingWithInNearFeatures_K11k2L100D10_2000A100_2000d10_100N2F0X20V3000W8M3G1.indels.tsv",
        "outTestHP11/MarkersNonoverlappingWithInNearFeatures_K11k2L100D10_2000A100_2000d10_100N2F0X20V3000W8M3G1.indels.pdf",
        1000, 2000, "S.lycopersicum,S.pennellii")
    }
else if (testing == 2)
    {
    args = c("~/Documents/UCDavis/BradyLab/Genomes/IGGPIPE",
        "outHP14/IndelGroupsNonoverlappingWithInNearFeatures_K14k2L100D1_3000A100_3000d100_100N2F0.indels.tsv",
        "outHP14/MarkersNonoverlappingWithInNearFeatures_K14k2L400D10_1500A400_1500d50_300N2F0X20V3000W8M3G1.indels.pdf",
        1000, 3000, "S.lycopersicum,S.pennellii")
    }
else stop("Unknown value for 'testing'")
}

Nexpected = 6
if (length(args) != Nexpected)
    {
    usage = c(
        "Read a data frame of Indel information that has had gene feature information added to it, and",
        "make plots of various aspects of the data.  NOTE: This is currently specific to S. lycopersicum.",
        "See annotate/HPlongITAG2.4_addFeatureIDs.nonoverlappingIndelGroups for parameters used to create the",
        "input file for this.",
        "",
        "Usage: Rscript plotIndelsWithFeatures.R <wd> <IndelsFile> <pdfPlotFile> <nearInBp> <maxIndelLen> <genomeNames>",
        "",
        "Arguments:",
        "   <wd>                 : Path of R working directory, specify other file paths relative to this.",
        "   <IndelsFile>         : Name of file containing Indel position information with a column named 'isInNear'.",
        "                          See annotate/test_add_isInNearColumn.markers",
        "   <pdfPlotFile>        : Name of pdf file to be created.",
        "   <nearInBp>           : Number of base pairs used to define 'near' when running annotateFile.R",
        "   <maxIndelLen>        : Maximum possible size of any indel that is detected with the method used to detect them.",
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

genomeNames = unlist(strsplit(args[6], ",", fixed=TRUE))
catnow("  genomeNames: ", paste(genomeNames, collapse=","), "\n")

########################################
# Initialize.
########################################

# Read Indels data.
dfIndels = read.table(IndelsFile, header=TRUE, row.names=NULL, sep="\t", stringsAsFactors=FALSE)
if (nrow(dfIndels) == 0)
    stop("There are no Indels.")

# Get the genome identifying letters from the dfIndels id columns names.
genomeLtrs = sub("id", "", colnames(dfIndels)[grepl("id$", colnames(dfIndels))])
Ngenomes = length(genomeLtrs)
if (Ngenomes < 2)
    stop("Expected to recognize at least two genomes in the data column names")
if (Ngenomes != length(genomeNames))
    stop("Expected ", Ngenomes, " genome names but only got ", length(genomeNames), " in the argument list")
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

########################################
# What fraction overlap more than one?
########################################
sum(grepl(",", dfIndels$isInNear))/nrow(dfIndels)

# With such a low frequency, simplify things.  Let's keep the extras, but split
# and duplicate each such row.
whichMultiple = which(grepl(",", dfIndels$isInNear))
isInNear = strsplit(dfIndels[whichMultiple, "isInNear"], ",", fixed=TRUE)
dupes = sapply(1:length(whichMultiple), function(i)
    {
    N = length(isInNear[[i]])
    df = dfIndels[rep(whichMultiple[i], N),]
    df$isInNear = isInNear[[i]]
    return(df)
    }, simplify=FALSE)
dupes = do.call.rbind.fast(dupes)
dfIndels = rbind(dfIndels[-whichMultiple,], dupes)
rm(whichMultiple, isInNear, dupes)

########################################
# Modify the isInNear column to simplify it.  We don't care about the gene, only
# about intron, etc.  We don't care about strand, but we do care whether it is
# overlap (@), upstream (-), or downstream(+).  We don't care about the amount
# of overlap, but we do care about the amount upstream or downstream.
########################################

dfIndels$isInNear = gsub("Solyc[0-9]+g[0-9]+\\.[0-9+]\\.[0-9]+\\.", "", dfIndels$isInNear)
dfIndels$featStrand = gsub("(^|,)[^(]*\\([^:]*:([^)]*)\\)", "\\1\\2", dfIndels$isInNear)
dfIndels$flags = gsub("(^|,)[^(]*\\(([^:]*):[^)]*\\)", "\\1\\2", dfIndels$isInNear)
dfIndels$featDist = gsub("@", "", dfIndels$flags)
dfIndels$flags = gsub("(^|,)(.)[^,]+", "\\1\\2", dfIndels$flags)
dfIndels$isInNear = gsub("\\([^)]*\\)", "", dfIndels$isInNear)
dfIndels$whichFeature = gsub("(^|,)[^:]*:", "\\1", dfIndels$isInNear)
dfIndels$isInNear = gsub(":[0-9]+", "", dfIndels$isInNear)
dfIndels$isInNear = gsub("five_prime_UTR", "5'UTR", dfIndels$isInNear)
dfIndels$isInNear = gsub("three_prime_UTR", "3'UTR", dfIndels$isInNear)
dfIndels$isInNear[dfIndels$isInNear == ""] = "intergenic"

# We need to add "upstream" and "downstream" features.  Any indel whose feature is CDS,
# 5'UTR, or 3'UTR and has a flag that is "+" or "-" rather than "@" indicates an indel
# lying upstream or downstream of the feature.  If the feature lies on the "+" strand,
# then a "+" flag means downstream and a "-" flag means upstream.  If the feature lies
# on the "-" strand, then a "+" flag means upstream and a "-" flag means downstream.
# (Double-check that!)  Note: all features are assumed to be S. lycopersicum features,
# and it is assumed that is the reference, or first, genome.
isNear_CDS = (dfIndels$isInNear == "CDS") & (dfIndels$flags != "@")
isNear_5UTR = (dfIndels$isInNear == "5'UTR") & (dfIndels$flags != "@")
isNear_3UTR = (dfIndels$isInNear == "3'UTR") & (dfIndels$flags != "@")
isNear_CDS_UTR = isNear_CDS | isNear_5UTR | isNear_3UTR
strands = dfIndels$featStrand[isNear_CDS_UTR]
flags = dfIndels$flags[isNear_CDS_UTR]
isUpstream = (strands == flags)
dfIndels$isInNear[isNear_CDS_UTR][isUpstream] = "upstream"
dfIndels$isInNear[isNear_CDS_UTR][!isUpstream] = "downstream"
#table(dfIndels$isInNear[isNear_CDS]) # Yes!
#table(dfIndels$isInNear[isNear_5UTR]) # Yes!
#table(dfIndels$isInNear[isNear_3UTR]) # Yes!

features = c("upstream", "5'UTR", "CDS", "intron", "3'UTR", "downstream", "intergenic")
Nfeatures = length(features)
if (any(!features %in% unique(dfIndels$isInNear)))
    stop("Features skipped (this is set up for S. lycopersicum only)")
if (any(!unique(dfIndels$isInNear) %in% features))
    stop("Missing features (this is set up for S. lycopersicum only)")

########################################
# Plot colors for each genome and each feature.
########################################

genomeCols = colorBlind.8[1:Ngenomes]
if (Ngenomes > 8)
    genomeCols = rainbow(Ngenomes)
names(genomeCols) = genomeLtrs

featureCols = colorBlind.8[1:Nfeatures]
if (Nfeatures > 8)
    featureCols = rainbow(Nfeatures)
names(featureCols) = features

########################################
# Grid line color.
########################################

gridCol = "gray"

########################################
# Line types for each genome and each feature.
########################################

ltys = c("solid", "dashed", "dotted", "dotdash", "longdash", "twodash")
genomeLtys = rep(ltys, length.out=Ngenomes)
names(genomeLtys) = genomeLtrs
featureLtys = rep(ltys, length.out=Nfeatures)
names(featureLtys) = features

########################################
# Create pdf file.
########################################

pdf(pdfPlotFile, height=8, width=6)
par(mar=c(9,8,6,3))

########################################
# Make a bar plot (one bar per genome) showing the number of indels occurring
# in each feature type.
########################################

# Count number for each feature.
countMtx = matrix(0, ncol=Nfeatures, nrow=Ngenomes, dimnames=list(genomeLtrs, features))
for (genome in genomeLtrs)
    {
    # Ignore insertions.
    feats = dfIndels[dfIndels[,delCol[genome]] > 0, "isInNear"]
    x = table(feats)
    missing = is.na(x[features])
    if (any(missing))
        x[features[missing]] = 0
    countMtx[genome,] = x
    }

# Compute axis labels.
ypretty = pretty.log(countMtx)
ylim = range(ypretty)
yat = ypretty
ylab="Number of Indels"
if (ypretty[1] >= 1000)
    {
    ypretty = ypretty / 1000
    ylab="Number of Indels  (x1000)"
    }
ypretty = format(ypretty, scientific=FALSE)

# Make the bar plot.
barplot(countMtx, beside=TRUE, col=genomeCols, axes=FALSE, log="y", las=2, cex.names=1.5,
    main="Distribution of indel locations",
    ylab="")
title(paste("upstream/downstream defined as within ", nearInBp, " bp of 5'UTR/CDS/3'UTR", sep=""), line=1.5, cex.main=0.8)
axis(2, at=yat, labels=ypretty, las=2, cex.axis=1.5)
mtext(ylab, 2, line=6, cex=1.5)

# Legend.
legend("topleft", genomeNames, fill=genomeCols, cex=2)

########################################
# Make a line plot (one line per genome) where x-axis positions where points are
# plotted are bins containing a range of indel deletion lengths, and the y-axis
# is the number of indels with a deletion length in that range.
########################################

# Indel deletion length bins.
minIndelLen = c(0, 1, 5, 25, 100, 250, 1000)
Nbins = length(minIndelLen)
binRanges = paste(minIndelLen, "..", c((minIndelLen-1)[-1], "+"), sep="")
binRanges = sub("0..0", "0", binRanges)
binRanges = sub("1000..+", "1000+", binRanges)

# Count number of indels in each indel size range.
countMtx = matrix(0, ncol=Nbins, nrow=Ngenomes, dimnames=list(genomeLtrs, minIndelLen))
for (genome in genomeLtrs)
    {
    x = findInterval(dfIndels[, delCol[genome]], minIndelLen)
    x = table(x)
    missing = setdiff(1:Nbins, as.integer(names(x)))
    if (length(missing) > 0)
        x[as.character(missing)] = 0
    names(x) = as.character(minIndelLen[as.integer(names(x))])
    countMtx[genome,] = x[as.character(minIndelLen)]
    }

# Compute axis labels.
xlim = c(1, Nbins)
xlab = paste("Deletion size range (bp)\n0=insertion, max detectable ~ ", maxIndelLen, " bp", sep="")
ypretty = pretty.log(countMtx)
ylim = range(ypretty)
yat = ypretty
ylab="Number of Indels"
if (ypretty[1] >= 1000)
    {
    ypretty = ypretty / 1000
    ylab="Number of Indels  (x1000)"
    }
ypretty = format(ypretty, scientific=FALSE)

# Start the plot.
plot(NA, type="n", xlim=xlim, ylim=ylim, axes=FALSE, log="y", xlab="", ylab="",
    main="Distribution of indel deletion sizes")
mtext(binRanges, 1, at=1:Nbins, line=0, las=2, cex=1.5)
mtext(xlab, 1, line=7.5, cex=1.5)
axis(2, at=yat, labels=ypretty, las=2, cex.axis=1.5)
mtext(ylab, 2, line=4, cex=1.5)

# Loop and plot line/points for each genome.
for (genome in genomeLtrs)
    {
    lines(2:Nbins, countMtx[genome, -1], col=genomeCols[genome], lwd=3, lty=genomeLtys[genome])
    points(1:Nbins, countMtx[genome,], col=genomeCols[genome], pch=20, cex=2)
    }

# Make some grid lines.
segments(1:Nbins, ylim[1], 1:Nbins, ylim[2], col=gridCol)
segments(0, yat, xlim[2], yat, col=gridCol)

# Legend.
legend("topright", genomeNames, pch=20, col=genomeCols, pt.cex=2, lwd=3, lty=genomeLtys, seg.len=5, bty="n", cex=2)

########################################
# Similar to the previous, except show a line for each feature type, and make a
# separate plot for each genome.
########################################

# Count number of indels for each feature within each indel size range.
counts = array(0, c(Ngenomes, nrow=Nfeatures, Nbins), dimnames=list(genomeLtrs, features, minIndelLen))
for (genome in genomeLtrs)
    {
    for (feat in features)
        {
        x = findInterval(dfIndels[dfIndels$isInNear == feat, delCol[genome]], minIndelLen)
        x = table(x)
        missing = setdiff(1:Nbins, as.integer(names(x)))
        if (length(missing) > 0)
            x[as.character(missing)] = 0
        names(x) = as.character(minIndelLen[as.integer(names(x))])
        counts[genome,feat,] = x[as.character(minIndelLen)]
        }
    }

# Compute axis labels.
xlim = c(1, Nbins)
xlab = paste("Deletion size range (bp)\n0=insertion, max detectable ~ ", maxIndelLen, " bp", sep="")
ypretty = pretty.log(counts)
ylim = range(ypretty)
yat = ypretty
ylab="Number of Indels"
if (ypretty[1] >= 1000)
    {
    ypretty = ypretty / 1000
    ylab="Number of Indels  (x1000)"
    }
ypretty = format(ypretty, scientific=FALSE)

# Loop and plot feature lines/points for each genome.
for (genome in genomeLtrs)
    {
    # Start the plot.
    plot(NA, type="n", xlim=xlim, ylim=ylim, axes=FALSE, log="y", xlab="", ylab="",
        main=paste("Distribution of indel sizes in genome ",
            genomeNames[genome], " features\n(features are from S. lycopersicum gene model file)", sep=""))
    title(paste("upstream/downstream defined as within ", nearInBp, " bp of 5'UTR/CDS/3'UTR", sep=""), line=0.5, cex.main=0.8)
    mtext(binRanges, 1, at=1:Nbins, line=0, las=2, cex=1.5)
    mtext(xlab, 1, line=7.5, cex=1.5)
    axis(2, at=yat, labels=ypretty, las=2)
    mtext(ylab, 2, line=4, cex=1.5)

    # Loop for each feature type.
    for (feat in features)
        {
        lines(2:Nbins, counts[genome, feat, -1], col=featureCols[feat], lwd=2, lty=featureLtys[feat])
        points(1:Nbins, counts[genome, feat,], col=featureCols[feat], pch=20, cex=2)
        }

    # Label points left of x=2.  When distance is close in y-direction, spread them.
    y = counts[genome, features, 2]
    txt = features
    ord = order(y)
    y = y[ord]
    txt = txt[ord]
    near = (diff(y)/y[-1] < 0.2)
    while (any(near))
        {
        y[c(FALSE,near)] = y[c(FALSE,near)]*1.05
        near = (diff(y)/y[-1] < 0.2)
        }
    text(2, y, txt, pos=2, cex=0.8)

    # Make some grid lines.
    segments(1:Nbins, ylim[1], 1:Nbins, ylim[2], col=gridCol)
    segments(0, yat, xlim[2], yat, col=gridCol)

    # Legend.
    legend("topright", features, pch=20, col=featureCols, pt.cex=2, lwd=2, lty=featureLtys, seg.len=5, bty="n")
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
