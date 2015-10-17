################################################################################
# See usage below for description.
# Author: Ted Toal
# Date: 2015
# Brady Lab, UC Davis
################################################################################

# Enclose everything in braces so stop statements will work correctly.
{

# Maximum number of points to plot on any plot.
maxPlotPoints = 20000

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
#testing = 2 # For testing only.  .HP
{
if (testing == 0)
    args = commandArgs(TRUE)
else if (testing == 1)
    {
    args = c("~/Documents/UCDavis/BradyLab/Genomes/IGGPIPE",
        "outTestHP11/MarkersNonoverlapping_K11k2L100D10_2000A100_2000d10_100N2F0X20V3000W8M3G1.indels.tsv",
        "outTestHP11/MarkersNonoverlapping_K11k2L100D10_2000A100_2000d10_100N2F0X20V3000W8M3G1.indels.pdf")
    }
else if (testing == 2)
    {
    args = c("~/Documents/UCDavis/BradyLab/Genomes/IGGPIPE",
        "outHP14/MarkersNonoverlapping_K14k2L400D10_1500A400_1500d50_300N2F0X20V3000W8M3G1.indels.tsv",
        "outHP14/MarkersNonoverlapping_K14k2L400D10_1500A400_1500d50_300N2F0X20V3000W8M3G1.indels.pdf")
    }
else stop("Unknown value for 'testing'")
}

Nexpected = 3
if (length(args) != Nexpected)
    {
    usage = c(
        "Read a data frame of Indel information and make plots of various aspects of the data.",
        "",
        "Usage: Rscript plotIndels.R <wd> <IndelsFile> <pdfPlotFile>",
        "",
        "Arguments:",
        "   <wd>                 : Path of R working directory, specify other file paths relative to this.",
        "   <IndelsFile>         : Name of file containing Indel position information.",
        "   <pdfPlotFile>        : Name of pdf file to be created."
        )
    for (S in usage)
        catnow(S, "\n", sep="")
    stop("Try again with correct number of arguments")
    }

catnow("plotIndels.R arguments:\n")
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
refGenomeLtr = genomeLtrs[1]
otherGenomeLtrs = genomeLtrs[-1]

# Get all pairs of genome letters.
genomePairs = t(combn(genomeLtrs, 2))
Npairs = nrow(genomePairs)

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
# Create pdf file.
########################################

pdf(pdfPlotFile, height=8, width=8)

########################################
# Make scatter plots (one per genome pair) where x-axis is amplicon or LCR length
# difference and y-axis is total number of Indels.
########################################

cat("Npairs = ", Npairs, "\n")
for (i in 1:Npairs)
    {
    pair = genomePairs[i,]
    g1 = pair[1]
    g2 = pair[2]
    delCol1 = delCol[g1]
    delCol2 = delCol[g2]
    ampliconLens = tapply(1:nrow(dfIndels), dfIndels$ID, function(ii)
        abs(sum(dfIndels[ii,delCol1])-sum(dfIndels[ii,delCol2])))
    numIndels = tapply(1:nrow(dfIndels), dfIndels$ID, length)
    numIndels = numIndels + rnorm(length(numIndels), 0, 0.1) # Add some jitter.
    xlim = range(pretty(c(0, ampliconLens)))
    ylim = range(pretty(c(0, numIndels)))
    N = length(ampliconLens)
    if (N > maxPlotPoints)
        {
        idxs = sample(1:N, maxPlotPoints)
        N = maxPlotPoints
        ampliconLens = ampliconLens[idxs]
        numIndels = numIndels[idxs]
        }
    plot(ampliconLens, numIndels, type="p", xlim=xlim, ylim=ylim, cex=0.5, pch=20,
        main=paste("Scatter plot of amplicon/LCR size diff. vs. number of indels, genomes ",
            g1, " and ", g2, sep=""),
        xlab="Amplicon size difference (bp)", ylab="Number of Indels")
    catnow("Plotted", N, "points in plot 1\n") 
    meanAmpliconLen = mean(ampliconLens)
    meanNumIndels = mean(numIndels)
    catnow("meanAmpliconLen =", meanAmpliconLen, " meanNumIndels =", meanNumIndels, "\n")
    points(meanAmpliconLen, meanNumIndels, pch=20, cex=2, col="tan")
    text(meanAmpliconLen, meanNumIndels, "Mean", cex=2, col="tan", pos=4)
    }

########################################
# Make scatter plots (one per genome pair) where x-axis is amplicon or LCR length
# difference and y-axis is Indel size.
########################################

cat("Npairs = ", Npairs, "\n")
for (i in 1:Npairs)
    {
    pair = genomePairs[i,]
    g1 = pair[1]
    g2 = pair[2]
    delCol1 = delCol[g1]
    delCol2 = delCol[g2]
    ampliconLens = unlist(tapply(1:nrow(dfIndels), dfIndels$ID, function(ii)
        rep(abs(sum(dfIndels[ii,delCol1])-sum(dfIndels[ii,delCol2])), length(ii))))
    indelLens = unlist(tapply(1:nrow(dfIndels), dfIndels$ID, function(ii)
        {
        V1 = dfIndels[ii, delCol1]
        V2 = dfIndels[ii, delCol2]
        V = c(V1, V2)
        V = V[V != 0]
        return(V)
        }))
    indelLens = indelLens + rnorm(length(indelLens), 0, 0.1) # Add some jitter.
    xlim = range(pretty(c(0, ampliconLens)))
    ylim = range(pretty(c(0, indelLens)))
    N = length(ampliconLens)
    if (N > maxPlotPoints)
        {
        idxs = sample(1:N, maxPlotPoints)
        N = maxPlotPoints
        ampliconLens = ampliconLens[idxs]
        indelLens = indelLens[idxs]
        }
    plot(ampliconLens, indelLens, type="p", xlim=xlim, ylim=ylim, cex=0.5, pch=20,
        main=paste("Scatter plot of amplicon/LCR size diff. vs. indels size, genomes ",
            g1, " and ", g2, sep=""),
        xlab="Amplicon size difference (bp)", ylab="Indels size (bp)")
    catnow("Plotted", N, "points in plot 2\n") 
    meanAmpliconLen = mean(ampliconLens)
    meanIndelsLens = mean(indelLens)
    catnow("meanAmpliconLen =", meanAmpliconLen, " meanIndelsLens =", meanIndelsLens, "\n")
    points(meanAmpliconLen, meanIndelsLens, pch=20, cex=2, col="tan")
    text(meanAmpliconLen, meanIndelsLens, "Mean", cex=2, col="tan", pos=4)
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
