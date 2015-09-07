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

# cat() that immediately flushes to console.
catnow = function(...)
    {
    cat(...)
    flush.console()
    return(invisible(0))
    }

# Get arguments.
testing = 0
#testing = 1 # For testing only.  .test
#testing = 2 # For testing only.  .HP
{
if (testing == 0)
    args = commandArgs(TRUE)
else if (testing == 1)
    {
    args = c("~/Documents/UCDavis/BradyLab/Genomes/kmers/IGGPIPE", 2, 0.3,
        "outTestHP11/MarkersOverlapping_K11k2L100D10_2000A100_2000d10_100N2F0X20V3000W8M3G1.tsv",
        "outTestHP11/MarkersNonoverlapping_K11k2L100D10_2000A100_2000d10_100N2F0X20V3000W8M3G1.tsv",
        "outTestHP11/MarkerCounts_K11k2L100D10_2000A100_2000d10_100N2F0X20V3000W8M3G1",
        "outTestHP11/MarkerDensity_K11k2L100D10_2000A100_2000d10_100N2F0X20V3000W8M3G1",
        "outTestHP11/GenomeData/Genome_1.idlens", "outTestHP11/GenomeData/Genome_2.idlens")
    }
else if (testing == 2)
    {
    args = c("~/Documents/UCDavis/BradyLab/Genomes/kmers/IGGPIPE", 2, 0.25,
        "outHP14/MarkersOverlapping_K11k2L100D10_2000A100_2000d10_100N2F0X20V3000W8M3G1.tsv",
        "outHP14/MarkersNonoverlapping_K11k2L100D10_2000A100_2000d10_100N2F0X20V3000W8M3G1.tsv",
        "outHP14/MarkerCounts_K14k2L400D10_1500A400_1500d50_300N2F0X20V3000W8M3G1",
        "outHP14/MarkerDensity_K14k2L400D10_1500A400_1500d50_300N2F0X20V3000W8M3G1",
        "outHP14/GenomeData/Genome_1.idlens", "outHP14/GenomeData/Genome_2.idlens")
    }
else stop("Unknown value for 'testing'")
}

NexpectedMin = 9
if (length(args) < NexpectedMin)
    {
    usage = c(
        "Read a data frame of good candidate IGG markers that have passed all tests, and make",
        "marker bar plots and density plots showing numbers of markers on each chromosome.",
        "",
        "Usage: Rscript findPrimers.R <wd> <plotNDAmin> <alpha> <overlappingFile> <nonoverlappingFile> \\",
        "           <pdfCountsFilePfx> <pngDensityFilePfx> <idlensFile1> ...",
        "",
        "Arguments:",
        "   <wd>                 : Path of R working directory, specify other file paths relative to this.",
        "   <plotNDAmin>         : Minimum number of distinct amplicons for marker to be plotted.",
        "   <alpha>              : Value of alpha channel of color for plotting marker density lines.",
        "   <overlappingFile>    : Input file containing overlapping markers with e-PCR-checked primers.",
        "   <nonoverlappingFile> : Input file containing non-overlapping markers with e-PCR-checked primers.",
        "   <pdfCountsFilePfx>   : Prefix, including path, of output .pdf file to write, containing plots of",
        "                          number of markers per seq ID.  The suffix '.plot.pdf' is appended to this",
        "                          prefix to form the full file name.",
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

catnow("findPrimers.R arguments:\n")
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
# sequence IDs and their lengths.
idlens = list()
for (genome in genomeLtrs)
    {
    idlens[[genome]] = read.table(idlensFiles[genome], header=TRUE, row.names=1, sep="\t", stringsAsFactors=FALSE)
    idlens[[genome]]$len = idlens[[genome]]$len * 1e-6 # Convert to Mbp.
    }

########################################
# Remove markers not satisfying plotNDAmin.
########################################

dfMarkers = dfMarkers[dfMarkers$NDA >= plotNDAmin,]
dfNoOverlaps = dfNoOverlaps[dfNoOverlaps$NDA >= plotNDAmin,]

########################################
# Plot counts of markers on chromosomes as either a bar or line plot.
########################################

pdfCountsFile = paste(pdfCountsFilePfx, ".plot.pdf", sep="")
pdf(pdfCountsFile, height=8, width=11)
par(mar=par("mar")+c(4,2,0,0))
par(mgp=c(3,1,1))

for (genome in genomeLtrs)
    {
    id.Col = idCol[genome]
    pct.Col = pctCol[genome]
    pos1.Col = ampPos1Col[genome]
    V = tapply(1:nrow(dfMarkers), dfMarkers[,id.Col], function(ii)
        return(length(ii)/(idlens[[genome]][dfMarkers[ii[1], id.Col], "len"])))
    V2 = tapply(1:nrow(dfNoOverlaps), dfNoOverlaps[,id.Col], function(ii)
        return(length(ii)/(idlens[[genome]][dfNoOverlaps[ii[1], id.Col], "len"])))
    V[setdiff(rownames(idlens[[genome]]), names(V))] = 0
    V2[setdiff(rownames(idlens[[genome]]), names(V2))] = 0
    V = V[rownames(idlens[[genome]])]
    V2 = V2[rownames(idlens[[genome]])]
    # Make either a bar plot or a line plot, depending on number of chromosomes/scaffolds.
    if (length(V) <= 50) # Reasonable maximum number of chromosomes.
        {
        ylim = range(pretty(c(0, max(V)*1.1)))
        barplot(V, col="gray", las=2, cex.main=1.5, cex.names=1.5, ylim=ylim, yaxp=c(ylim, 10), tck=1,
            ylab="IGG markers per 1 Mbp", cex.lab=2, cex.axis=1.2,
            main=paste("Density of IGG markers in each FASTA sequence of genome '", genome, "'", sep=""))
        barplot(V2, width=0.5, space=c(1.4, 0.9), col="blue", add=TRUE, axes=FALSE, axisnames=FALSE, main="")
        legend("topright", c("All markers", "Non-overlapping markers"), fill=c("gray", "blue"), cex=1.5)
        }
    else # Assume thousands of scaffolds, do a sorted plot of number per Mbp of scaffold.
        {
        ord = order(V)
        V = V[ord]
        V2 = V2[ord]
        plot(1:length(V), V, type="l", col="black", cex.main=1.5, cex.axis=1.2,
            xlab="Number of FASTA Sequences", ylab="Density of IGG markers per 1 Mbp", cex.lab=2,
            main=paste("Density of IGG markers in FASTA sequences of genome '", genome, "'", sep=""))
        lines(1:length(V2), V2, col="blue")
        legend("topleft", c("All markers", "Non-overlapping markers"), lwd=2, col=c("black", "blue"), cex=1.5)
        }
    }

dev.off()
catnow("Finished making plots of counts of number of good candidate markers, output file:\n", pdfCountsFile, "\n")

########################################
# Plot density of markers on chromosomes/scaffolds by drawing rectangles and plotting
# a bar on each one at each marker position, using black with an alpha that is low.
########################################

# Use the same scale for the y-axis of all plots.
maxSeqLenMbp = max(unlist(idlens))
xlim = c(0, 1.2*maxSeqLenMbp)

# Set the color of the line segments drawn at each marker position.  Use an
# alpha value that is fairly small.
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

    # Title cex expansion factor needs to be smaller if the number of IDs is small.
    cex.main = (numIds-2)/4 + 1
    if (cex.main > 3)
        cex.main = 3

    # Make the png bitmap width be proportional to the number of sequences in this genome.
    pngDensityFile = paste(pngDensityFilePfx, "_", genome, ".plot.png", sep="")
    png(pngDensityFile, height=1200, width=500+100*numIds)
    par(mar=par("mar")+c(4,2,0,0))

    # Start the plot.
    xlim = c(0, numIds+1)
    ylim = c(0, maxSeqLenMbp)
    plot(NA, type="n", xlim=xlim, ylim=ylim, cex.main=cex.main,
        xlab="", ylab="Sequence position (Mbp)", cex.lab=3, axes=FALSE, frame=FALSE,
        main=paste("Non-overlapping IGG markers in each FASTA sequence of genome '", genome, "'", sep=""))
    at = pretty(ylim)
    axis(2, at=at, line=-2, cex.axis=3)

    # Plot each sequence (chromosome/scaffold).
    for (i in 1:numIds)
        {
        id = ids[i]
        h = 0.4
        x = i - h
        len = idlens[[genome]][id, "len"]
        mtext(id, 1, at=x, cex=2, las=2, line=-2)
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
