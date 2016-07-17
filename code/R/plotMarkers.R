################################################################################
# See usage below for description.
# Author: Ted Toal
# Date: 2013-2016
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
#testing = 2 # For testing only.  .HP
#testing = 3 # For testing only.  .TaCW
#testing = 4 # For testing only.  .AtCL
#testing = 5 # For testing only.  .HPT
{
if (testing == 0)
    args = commandArgs(TRUE)
else if (testing == 1)
    {
    args = c("~/Documents/UCDavis/BradyLab/IGGPIPE/IGGPIPE", 2, 0.3,
        "outTestHP11/MarkersOverlapping_K11k2L100D10_2000A100_2000d10_100N2F0X20V3000W8M3G1.tsv",
        "outTestHP11/MarkersNonoverlapping_K11k2L100D10_2000A100_2000d10_100N2F0X20V3000W8M3G1.tsv",
        "outTestHP11/MarkerCounts_K11k2L100D10_2000A100_2000d10_100N2F0X20V3000W8M3G1",
        "outTestHP11/MarkerDensity_K11k2L100D10_2000A100_2000d10_100N2F0X20V3000W8M3G1",
        "outTestHP11/GenomeData/Genome_1.idlens", "outTestHP11/GenomeData/Genome_2.idlens")
    }
else if (testing == 2)
    {
    args = c("~/Documents/UCDavis/BradyLab/IGGPIPE/IGGPIPE", 2, 0.25,
        "outHP14/MarkersOverlapping_K14k2L400D10_1500A400_1500d50_300N2F0X20V3000W8M3G1.tsv",
        "outHP14/MarkersNonoverlapping_K14k2L400D10_1500A400_1500d50_300N2F0X20V3000W8M3G1.tsv",
        "outHP14/MarkerCounts_K14k2L400D10_1500A400_1500d50_300N2F0X20V3000W8M3G1",
        "outHP14/MarkerDensity_K14k2L400D10_1500A400_1500d50_300N2F0X20V3000W8M3G1",
        "outHP14/GenomeData/Genome_1.idlens", "outHP14/GenomeData/Genome_2.idlens")
    }
else if (testing == 3)
    {
    args = c("~/Documents/UCDavis/BradyLab/IGGPIPE/IGGPIPE", 2, 0.5,
        "outTaCW15/MarkersOverlapping_K15k2L200D5_1000A200_1000d50_200N2F0X20V3000W9M3G1.tsv",
        "outTaCW15/MarkersNonoverlapping_K15k2L200D5_1000A200_1000d50_200N2F0X20V3000W9M3G1.tsv",
        "outTaCW15/MarkerCounts_K15k2L200D5_1000A200_1000d50_200N2F0X20V3000W9M3G1",
        "outTaCW15/MarkerDensity_K15k2L200D5_1000A200_1000d50_200N2F0X20V3000W9M3G1",
        "outTaCW15/GenomeData/Genome_1.idlens", "outTaCW15/GenomeData/Genome_2.idlens")
    }
else if (testing == 4)
    {
    args = c("~/Documents/UCDavis/BradyLab/IGGPIPE/IGGPIPE", 2, 0.85,
        "outCL13/MarkersOverlapping_K13k4L400D15_1500A400_1500d50_300N2F2X20V5000W8M0G0.tsv",
        "outCL13/MarkersNonoverlapping_K13k4L400D15_1500A400_1500d50_300N2F2X20V5000W8M0G0.tsv",
        "outCL13/MarkerCounts_K13k4L400D15_1500A400_1500d50_300N2F2X20V5000W8M0G0",
        "outCL13/MarkerDensity_K13k4L400D15_1500A400_1500d50_300N2F2X20V5000W8M0G0",
        "outCL13/GenomeData/Genome_1.idlens", "outCL13/GenomeData/Genome_2.idlens")
    }
else if (testing == 5)
    {
    args = c("~/Documents/UCDavis/BradyLab/IGGPIPE/IGGPIPE", 2, 0.8,
        "outHPT14/MarkersOverlapping_K14k2L300D5_1500A300_1500d50_300N2F0X15V2500W8M0G0.tsv",
        "outHPT14/MarkersNonoverlapping_K14k2L300D5_1500A300_1500d50_300N2F0X15V2500W8M0G0.tsv",
        "outHPT14/MarkerCounts_K14k2L300D5_1500A300_1500d50_300N2F0X15V2500W8M0G0",
        "outHPT14/MarkerDensity_K14k2L300D5_1500A300_1500d50_300N2F0X15V2500W8M0G0",
        "outHPT14/GenomeData/Genome_1.idlens", "outHPT14/GenomeData/Genome_2.idlens",
        "outHPT14/GenomeData/Genome_3.idlens")
    }
else stop("Unknown value for 'testing'")
}

NexpectedMin = 9
if (length(args) < NexpectedMin)
    {
    usage = c(
        "Read a data frame of good candidate IGG markers that have passed all tests, and make",
        "marker bar plots and density plots showing numbers of markers on each chromosome and",
        "distribution of amplicon size differences.",
        "",
        "Usage: Rscript plotMarkers.R <wd> <plotNDAmin> <alpha> <overlappingFile> <nonoverlappingFile> \\",
        "           <pdfCountsFilePfx> <pngDensityFilePfx> <idlensFile1> ...",
        "",
        "Arguments:",
        "   <wd>                 : Path of R working directory, specify other file paths relative to this.",
        "   <plotNDAmin>         : Minimum number of distinct amplicons for marker to be plotted.",
        "   <alpha>              : Value of alpha channel of color for plotting marker density lines.",
        "   <overlappingFile>    : Input file containing overlapping markers with e-PCR-checked primers.",
        "   <nonoverlappingFile> : Input file containing non-overlapping markers with e-PCR-checked primers.",
        "   <pdfCountsFilePfx>   : Prefix, including path, of output .pdf file to write, containing plots of",
        "                          number of markers per seq ID and histogram of marker amplicon size differences.",
        "                          The suffix '.plot.pdf' is appended to this prefix to form the full file name.",
        "   <pngDensityFilePfx>  : Prefix, including path, of output .png files to write, containing plots of",
        "                          density of markers per seq ID.  A separate .png file is produced for each",
        "                          genome, appending '_<genome letter>.plot.png' to this prefix to form the",
        "                          full file name.",
        "   <idlensFile1>        : File containing sequence IDs/lengths for genome 1.",
        "   ...                  : Additional sequence IDs/lengths files for genomes 2..Ngenomes."
        )
    for (S in usage)
        catnow(S, "\n", sep="")
    stop("Try again with correct number of arguments")
    }

catnow("plotMarkers.R arguments:\n")
workingDirectory = args[1]
catnow("  workingDirectory: ", workingDirectory, "\n")
if (!dir.exists(workingDirectory))
    stop("Directory doesn't exist: ", workingDirectory)
setwd(workingDirectory)

plotNDAmin = as.integer(args[2])
if (is.na(plotNDAmin))
    stop("plotNDAmin must be an integer")
catnow("  plotNDAmin: ", plotNDAmin, "\n")

alpha = as.numeric(args[3])
if (is.na(alpha) || alpha <= 0 || alpha > 1)
    stop("alpha must be 0 < alpha <= 1")
catnow("  alpha: ", alpha, "\n")

overlappingFile = args[4]
catnow("  overlappingFile: ", overlappingFile, "\n")
if (!file.exists(overlappingFile))
    stop("File doesn't exist: ", overlappingFile)

nonoverlappingFile = args[5]
catnow("  nonoverlappingFile: ", nonoverlappingFile, "\n")
if (!file.exists(nonoverlappingFile))
    stop("File doesn't exist: ", nonoverlappingFile)

pdfCountsFilePfx = args[6]
catnow("  pdfCountsFilePfx: ", pdfCountsFilePfx, "\n")

pngDensityFilePfx = args[7]
catnow("  pngDensityFilePfx: ", pngDensityFilePfx, "\n")

idlensFiles = args[8:length(args)]
for (idlensFile in idlensFiles)
    {
    if (!file.exists(idlensFile))
        stop("File doesn't exist: ", idlensFile)
    catnow("  idlensFile: ", idlensFile, "\n")
    }

########################################
# Initialize.
########################################

# Read overlapping marker data.
dfMarkers = read.table(overlappingFile, header=TRUE, row.names=NULL, sep="\t", stringsAsFactors=FALSE)
if (nrow(dfMarkers) == 0)
    stop("There are no overlapping IGG markers.")

# The k-mer length we are working with.
kmerLen = nchar(dfMarkers$kmer1[1])

# Get the genome identifying letters from the dfMarkers id columns names.
genomeLtrs = sub("id", "", colnames(dfMarkers)[grepl("id$", colnames(dfMarkers))])
Ngenomes = length(genomeLtrs)
if (Ngenomes < 2)
    stop("Expected to recognize at least two genomes in the data column names")
refGenomeLtr = genomeLtrs[1]
otherGenomeLtrs = genomeLtrs[-1]
if (length(idlensFiles) != Ngenomes)
    stop("Number of <idlensFile>'s must be equal to the number of genomes")
names(idlensFiles) = genomeLtrs

makeColVec = function(S)
    {
    V = paste(genomeLtrs, S, sep="")
    names(V) = genomeLtrs
    return(V)
    }

idCol = makeColVec("id")
refIdCol = idCol[refGenomeLtr]
pctCol = makeColVec("pct")
refPctCol = pctCol[refGenomeLtr]
ampPos1Col = makeColVec("ampPos1")
refAmpPos1Col = ampPos1Col[refGenomeLtr]
ampPos2Col = makeColVec("ampPos2")
refAmpPos2Col = ampPos2Col[refGenomeLtr]
ampDifCol = makeColVec(paste(refGenomeLtr, "dif", sep=""))
ampDifCol = ampDifCol[names(ampDifCol) != refGenomeLtr]

catnow("Number of overlapping markers:", nrow(dfMarkers), "\n")

########################################
# Read non-overlapping marker data.
########################################

dfNoOverlaps = read.table(nonoverlappingFile, header=TRUE, row.names=NULL, sep="\t", stringsAsFactors=FALSE)
if (nrow(dfNoOverlaps) == 0)
    stop("There are no non-overlapping IGG markers.")

# The k-mer length should be the same.
if (nchar(dfNoOverlaps$kmer1[1]) != kmerLen)
    stop("Expected non-overlapping markers to use same k-mer length as overlapping markers")

# The genome identifying letters from the dfNoOverlaps id columns names should be the same.
ltrs = sub("id", "", colnames(dfNoOverlaps)[grepl("id$", colnames(dfNoOverlaps))])
if (length(ltrs) != Ngenomes)
if (Ngenomes < 2)
    stop("Expected non-overlapping markers have same number of genomes as overlapping markers")
if (ltrs[1] != refGenomeLtr)
    stop("Expected non-overlapping markers have same number reference genome as overlapping markers")

catnow("Number of non-overlapping markers:", nrow(dfNoOverlaps), "\n")

########################################
# Get sequence IDs/lengths.
########################################

# Read the sequence ID/length files for each genome to get a list of all possible
# sequence IDs and their lengths.  Remove IDs that are not in the dfMarkers.
idlens = list()
for (genome in genomeLtrs)
    {
    id.Col = idCol[genome]
    theseIDs = sort(unique(dfMarkers[,id.Col]))
    idlens[[genome]] = read.table(idlensFiles[genome], header=TRUE, row.names=1, sep="\t", stringsAsFactors=FALSE)
    idlens[[genome]] = idlens[[genome]][theseIDs,,drop=FALSE]
    idlens[[genome]]$len = idlens[[genome]]$len * 1e-6 # Convert to Mbp.
    }

########################################
# Remove letters that are common to the beginning of all IDs to form shorter ID names
# that can be used for plotting names on axes.  Make lists of values for las, cex, line,
# and xlab that depend on how much the names are shortened.
########################################

las = list()
cex = list()
line = list()
xlab = list()
for (genome in genomeLtrs)
    {
    las[[genome]] = 2
    cex[[genome]] = 2
    line[[genome]] = -2
    xlab[[genome]] = ""
    S = rownames(idlens[[genome]])
    N = max(sapply(S, nchar))
    if (N > 3)
        {
        for (i in 1:N)
            {
            x = substr(S, 1, i)
            if (!all(x == x[1]))
                break
            }
        idlens[[genome]]$shortIDs = substr(S, i, N)

        # If the longest shortID is no more than 3 characters long, use las=0,
        # use an x-axis label, and use larger cex for displaying the IDs.
        if (max(sapply(idlens[[genome]]$shortIDs, nchar)) <= 3)
            {
            las[[genome]] = 0
            cex[[genome]] = 5
            line[[genome]] = 0
            xlab[[genome]] = "Sequence ID"
            }
        }
    }

########################################
# Remove markers not satisfying plotNDAmin.
########################################

dfMarkers = dfMarkers[dfMarkers$NDA >= plotNDAmin,]
dfNoOverlaps = dfNoOverlaps[dfNoOverlaps$NDA >= plotNDAmin,]
if (nrow(dfMarkers) == 0) stop("No overlapping markers.")
if (nrow(dfNoOverlaps) == 0) stop("No non-overlapping markers.")
catnow("After NDAmin applied, number of overlapping markers:", nrow(dfMarkers), "\n")
catnow("After NDAmin applied, number of non-overlapping markers:", nrow(dfNoOverlaps), "\n")

########################################
# Create pdf file.
########################################

pdfCountsFile = paste(pdfCountsFilePfx, ".plot.pdf", sep="")
pdf(pdfCountsFile, height=8, width=11)
par(mar=par("mar")+c(4,2,1,0))

########################################
# Plot counts of markers on chromosomes as either a bar or line plot.
# Any chromosome or scaffold that has 0 markers is not plotted.
########################################

for (genome in genomeLtrs)
    {
    id.Col = idCol[genome]
    pct.Col = pctCol[genome]
    pos1.Col = ampPos1Col[genome]
    NeachID = table(dfMarkers[,id.Col])
    V = NeachID/(idlens[[genome]][names(NeachID), "len"])
    names(V) = idlens[[genome]][names(V), "shortIDs"]
    NeachID2 = table(dfNoOverlaps[,id.Col])
    V2 = NeachID2/(idlens[[genome]][names(NeachID2), "len"])
    V2 = V2[names(NeachID)]
    # Make either a bar plot or a line plot, depending on number of chromosomes/scaffolds.
    if (length(V) <= 50) # Reasonable maximum number of chromosomes.
        {
        ylim = range(pretty(c(0, max(V)*1.2)))
        bp = barplot(V, col="gray", cex.main=2.5, cex.lab=3, cex.axis=2, names.arg=NA,
            ylim=ylim, yaxp=c(ylim, 10), mgp=c(4,1,0),
            xlab=xlab[[genome]], ylab="IGG markers per 1 Mbp",
            main=paste("IGG markers density in genome '", genome, "'", sep=""))
        barplot(V2, width=0.5, space=c(1.4, 0.9), col="blue", add=TRUE, axes=FALSE, axisnames=FALSE, main="")
        mtext(names(V), 1, at=as.vector(bp), cex=cex[[genome]]/2, las=las[[genome]], line=1)
        legend("topright", c("All markers", "Non-overlapping markers"), fill=c("gray", "blue"), cex=1.5)
        }
    else # Assume thousands of scaffolds, do a sorted plot of number per Mbp of scaffold.
        {
        # Divide them into 50 bins and average the bins.
        binSize = as.integer(length(V)/50)
        Nbins = as.integer(length(V)/binSize)
        start.idxs = seq(1, by=binSize, length.out=Nbins)
        end.idxs = c(start.idxs[-1]-1, length(V))
        V = sapply(1:length(start.idxs), function(i) mean(V[start.idxs[i]:end.idxs[i]]))
        ylim = range(c(V, V2), na.rm=TRUE)
        ylim[2] = ylim[2] * 1.1
        V2[is.na(V2)] = 0
        V2 = sapply(1:length(start.idxs), function(i) mean(V2[start.idxs[i]:end.idxs[i]]))
        ord = order(V, decreasing=TRUE)
        V = V[ord]
        V2 = V2[ord]
        plot(1:length(V), V, type="l", ylim=ylim, col="black", cex.main=1.5, cex.axis=1.2, log="y",
            xlab="Number of FASTA Sequences", ylab="Density of IGG markers per 1 Mbp", cex.lab=2,
            main=paste("IGG marker density in genome '", genome, "'", sep=""))
        lines(1:length(V2), V2, col="blue")
        legend("topleft", c("All markers", "Non-overlapping markers"), lwd=2, col=c("black", "blue"), cex=1.5)
        }
    }

########################################
# Plot histogram of amplicon size differences.
########################################

ampDifs = integer()
for (genome in otherGenomeLtrs)
    ampDifs = c(ampDifs, dfNoOverlaps[,ampDifCol[genome]])
maxDif = max(abs(ampDifs))
xpretty = pretty(c(-maxDif, +maxDif), 20)
xlim = range(xpretty)
par(mai=c(1,1,1,5))
hist(ampDifs, xlim=xlim, breaks=50,
    xlab=paste("Amplicon length difference from genome", refGenomeLtr, sep=""),
    main="Distribution of amplicon length differences", cex.lab=1.5, cex.axis=1.5)

########################################
# Close pdf file.
########################################

dev.off()
catnow("Finished making plots of counts of number of good candidate markers, output file:\n", pdfCountsFile, "\n")

########################################
# Plot density of markers on chromosomes/scaffolds by drawing rectangles and plotting
# a bar on each one at each marker position, using black with an alpha that is low.
########################################

# Use the same scale for the y-axis of all plots.
maxSeqLenMbp = max(sapply(idlens, function(df) max(df$len)))
xlim = c(0, 1.2*maxSeqLenMbp)

# Set the color of the line segments drawn at each marker position.
markerColor = "black"
markerColor = col2rgb(markerColor)
markerColor = rgb(markerColor["red",], markerColor["green",], markerColor["blue",], alpha*255, maxColorValue=255)

# Make a separate .png plot file for each genome.
catnow("Making marker density plots:\n")
for (genome in genomeLtrs)
    {
    # Get column names for this genome.
    id.Col = idCol[genome]
    pos1.Col = ampPos1Col[genome]
    ids = rownames(idlens[[genome]])
    numIds = length(ids)

    # If there are zillions of IDs, as there would be with scaffolds, we don't
    # want to plot them all, so choose a limited number, and only the largest.
    maxIdsToPlot = 40
    if (numIds > maxIdsToPlot)
        {
        numIds = maxIdsToPlot
        theseIdLens = idlens[[genome]][ids, "len"]
        names(theseIdLens) = ids
        theseIdLens = sort(theseIdLens, decreasing=TRUE)
        ids = names(theseIdLens)[1:numIds]
        }

    # Title cex expansion factor needs to be smaller if the number of IDs is small.
    cex.main = (numIds-2)/4 + 2.5
    if (cex.main > 10)
        cex.main = 10

    # Make the png bitmap width be proportional to the number of sequences in this genome.
    pngDensityFile = paste(pngDensityFilePfx, "_", genome, ".plot.png", sep="")
    png(pngDensityFile, height=1200, width=500+100*numIds)
    par(mar=par("mar")+c(4,5,1,0))

    # Start the plot.
    xlim = c(0, numIds+1)
    ylim = c(0, maxSeqLenMbp)
    plot(NA, type="n", xlim=xlim, ylim=ylim, cex.main=cex.main, mgp=c(5,1,0), frame=FALSE,
        xlab=xlab[[genome]], ylab="Sequence position (Mbp)", cex.lab=6, axes=FALSE,
        main=paste("Non-overlapping IGG markers in genome '", genome, "'", sep=""))
    at = pretty(ylim)
    axis(2, at=at, line=-3, cex.axis=6)

    # Plot each sequence (chromosome/scaffold).
    for (i in 1:numIds)
        {
        id = ids[i]
        h = 0.4
        x = i - h
        len = idlens[[genome]][id, "len"]
        shortID = idlens[[genome]][id, "shortIDs"]
        mtext(shortID, 1, at=x, cex=cex[[genome]], las=las[[genome]], line=line[[genome]])
        rect(x-h, 0, x+h, len)
        pos = dfNoOverlaps[dfNoOverlaps[,id.Col] == id, pos1.Col]*1e-6
        if (length(pos) > 0)
            {
            h2 = h*0.9
            segments(x-h2, pos, x+h2, pos, col=markerColor)
            }
        }
    dev.off()
    catnow("Genome", genome, "density plot in file", pngDensityFile, "\n")
    }

catnow("Finished making chromosome density plots of good candidate markers\n")
}

################################################################################
# End of file.
################################################################################
