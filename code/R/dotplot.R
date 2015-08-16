Issues:

1. Includes.
2. Libraries.
3. Args:
    - LCRs file
    - Genomes in it
    - idlens files
    - out file name
    - out image resolution
    - seqIds in each genome to be plotted
    - bkgd color
    - chr line color and width and spacing
    - chr text line and cex and color
    - points or line
    - point cex or line width
    - point or line color
4. Model like the other .R files.
5. Why not make args like mtext() etc. args, use type.convert(as.is=TRUE)

# Read a file of common unique k-mers between two genomes and plot them as an
# (x,y) point as a single point on a plot, making the plot axes extend over the
# entire genomic range.

R_INCLUDE_DIR = Sys.getenv("R_INCLUDE_DIR")
source(paste(R_INCLUDE_DIR, "Include_RootRsourceFile.R", sep=""))
source(Include_UtilityFunctions)

addThisLibrary("plotrix")
addThisLibrary("grid")

source(Include_TextFunctions)
source(Include_PlotFunctions)
source(Include_GenomeDb)
source(Include_solyGenome)
source(Include_sopeGenome)

setwd(paste(genomeDbDir, "kmers", sep=PATHSEP))

# Plot point for each kmer, or line segments for each LCR?
plotPoints = FALSE

inKmerFile = "ITAG2.3_PennV1/14/HeinzPennLCRs_K14Km4Lm400Dm20Dx1500.tsv"

chrFile1 = paste(ITAG2.3_dir, "ITAG2.3_genomic.idlens", sep=PATHSEP)
chrFile2 = paste(SopeV1genome_dir, "Sope.V1.idlens", sep=PATHSEP)

outFile = "ITAG2.3_PennV1/14/HeinzPennLCRs_K14Km4Lm400Dm20Dx1500.png"
outFile = "CheckChr12.png"
outFile.width.px = 2000
outFile.height.px = 2000

scale = 1e-6
col = ifelse(plotPoints, addAlpha("black", 0.2), "black")


# Genome 1.
ignoreChrs1 = "SL2.40ch00"
ignoreChrs1 = c("SL2.40ch00","SL2.40ch01","SL2.40ch02","SL2.40ch03","SL2.40ch04",
    "SL2.40ch05", "SL2.40ch06","SL2.40ch07","SL2.40ch08","SL2.40ch09","SL2.40ch10","SL2.40ch11")

chr1 = read.table(chrFile1, sep="\t", header=TRUE, stringsAsFactors=FALSE)
chr1 = chr1[!chr1$id %in% ignoreChrs1,]
chr1Lens = chr1$len*scale
names(chr1Lens) = chr1$id
Nchr1 = length(chr1Lens)
chr1Names = names(chr1Lens)

space.x = 5
chr1LensPlusSpace = chr1Lens + space.x
chr1LensPlusSpace[Nchr1] = chr1LensPlusSpace[Nchr1] - space.x
xmax = sum(chr1LensPlusSpace) - space.x

chr1Start = c(0, cumsum(chr1LensPlusSpace)[-Nchr1])
names(chr1Start) = chr1Names
chr1End = chr1Start + chr1Lens

xlim.start = chr1Start["SL2.40ch12"] # 0
xlim.width = chr1Lens["SL2.40ch12"] # xmax
xlim.end = xlim.start + xlim.width
xlim = c(xlim.start, xlim.end)

ignoreChrs1 = unique(c(ignoreChrs1, chr1Names[chr1Start > xlim[2] | chr1End < xlim[1]]))


# Genome 2.
ignoreChrs2 = "ch00"
ignoreChrs1 = c("ch00","ch01","ch02","ch03","ch04","ch05","ch06","ch07","ch08","ch09","ch10","ch11")

chr2 = read.table(chrFile2, sep="\t", header=TRUE, stringsAsFactors=FALSE)
chr2 = chr2[!chr2$id %in% ignoreChrs2,]
chr2Lens = chr2$len*scale
names(chr2Lens) = chr2$id
Nchr2 = length(chr2Lens)
chr2Names = names(chr2Lens)

space.y = 5
chr2LensPlusSpace = chr2Lens + space.y
chr2LensPlusSpace[Nchr2] = chr2LensPlusSpace[Nchr2] - space.y
ymax = sum(chr2LensPlusSpace) - space.y

chr2Start = c(0, cumsum(chr2LensPlusSpace)[-Nchr2])
names(chr2Start) = chr2Names
chr2End = chr2Start + chr2Lens

ylim.start = chr2Start["ch12"] # 0
ylim.width = chr2Lens["ch12"] # ymax
ylim.end = ylim.start + ylim.width
ylim = c(ylim.start, ylim.end)

ignoreChrs2 = unique(c(ignoreChrs2, chr2Names[chr2Start > ylim[2] | chr2End < ylim[1]]))


# Now plot the data.

# Read the kmer data file.
catnow("Reading file ", inKmerFile, "...", sep="")
df = read.table(inKmerFile, header=TRUE, sep="\t", row.names=1, stringsAsFactors=FALSE)
catnow("\n")
if (nrow(df) == 0) stop("Error reading ", inKmerFile)
saveDf = df
rownames(df) = NULL
df = df[, grepl("(\\.seqID|\\.pos|LCR)$", colnames(df)), drop=FALSE]
colnames(df) = c("id1", "x", "id2", "y", "LCR")

df = df[!(df$id1 %in% ignoreChrs1 | df$id2 %in% ignoreChrs2),]
df$x = chr1Start[df$id1] + df$x*scale
df$y = chr2Start[df$id2] + df$y*scale
df = df[df$x >= xlim[1] & df$x <= xlim[2] & df$y >= ylim[1] & df$y <= ylim[2],]

# This plot has millions of points, which in a pdf file are stored as individual
# vector commands, making the file huge and unloadable.  To get around this, we
# will plot to a .png file, which can be examined nicely in Preview.
png(outFile, width=outFile.width.px, height=outFile.height.px, antialias="none")
plot(NA, type="n", xlim=xlim, ylim=ylim, xlab="", ylab="", axes=FALSE, frame.plot=TRUE)

segments(chr1Start, 0, chr1Start+chr1Lens, 0, lwd=4, lend=1)
NN = 1:Nchr1
at = (chr1Start+chr1Start+chr1Lens)/2
ok = (at >= xlim[1] & at <= xlim[2])
at = at[ok]
NN = NN[ok]
if (length(NN) > 0)
    mtext(NN, 1, at=at, line=-1)

segments(0, chr2Start, 0, chr2Start+chr2Lens, lwd=4, lend=1)
NN = 1:Nchr2
at = (chr2Start+chr2Start+chr2Lens)/2
ok = (at >= ylim[1] & at <= ylim[2])
at = at[ok]
NN = NN[ok]
if (length(NN) > 0)
    mtext(NN, 2, at=at, line=-1)

if (nrow(df) > 0)
    {
    catnow("Plotting data...")
    if (plotPoints)
        points(df$x, df$y, cex=0.05, pch=20, col=col)
    else
        {
        df$LCR = paste("C", sub("^.*ch", "", df$id1), "_", df$LCR, sep="")
        x = tapply(1:nrow(df), df$LCR, function(ii) lines(df$x[ii], df$y[ii], lwd=1, col=col))
        }
    catnow("\n")
    }

dev.off()
