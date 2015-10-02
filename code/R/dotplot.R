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
#thisDir = "~/Documents/UCDavis/BradyLab/Genomes/kmers/IGGPIPE/code/R/" # For testing only.

# Source the necessary include files from the same directory containing this file.
source(paste(thisDir, "Include_Common.R", sep=""))

################################################################################
# Process program arguments.
################################################################################

# Get arguments.
testing = 0
#testing = 1 # For testing only.
#testing = 2 # For testing only.
#testing = 3 # For testing only.
{
if (testing == 0)
    args = commandArgs(TRUE)
else if (testing == 1)
    args = "dotplot.template"
else if (testing == 2)
    args = "dotplot/dotplot.wheatBD"
else if (testing == 3)
    args = "dotplot/dotplot.HsSl15_chr1"
else stop("Unknown value for 'testing'")
}

Nexpected = 1
if (length(args) != 1)
    {
    usage = c(
        "Read a parameter file that specifies parameters for a dot plot, including the",
        "name of a file of LCRs (a set of common unique k-mers between two genomes, each",
        "assigned to a locally conserved region), and plot them as a dot plot consisting",
        "of a series of line segments or dots, each one defined by the (x,y) points that",
        "are the positions of each common unique k-mer within the two genomes.  Genome 1",
        "is the x-axis and genome 2 is the y-axis.",
        "",
        "Usage: Rscript dotplot.R <paramFile>",
        "",
        "Arguments:",
        "   <paramFile> : Path of parameter file, see example parameter file dotplot.template."
        )
    for (S in usage)
        cat(S, "\n", sep="")
    stop("Try again with correct number of arguments")
    }

cat("dotplot.R arguments:\n")
paramFile = args[1]
cat("  paramFile: ", paramFile, "\n")
if (!file.exists(paramFile))
    stop("File doesn't exist: ", paramFile)

################################################################################
# Functions.
################################################################################

################################################################################
# Change colors to add alpha channel to them.  Alpha channel allows the color to be
# partially transparent.  When alpha=1, the color is solid; when 0, it is fully
# transparent; when 0.5, it is half-transparent.
#
# Arguments:
#   colors: vector of colors to which to add alpha channel.
#   alpha: vector of alpha channel values to use, recycled to length of "color" vector.
#
# Returns: vector of colors with alpha added.
################################################################################
addAlpha = function(colors, alpha)
    {
    alpha = rep(alpha, length.out=length(colors))
    col = col2rgb(colors)
    col = rgb(col["red",], col["green",], col["blue",], alpha*255, maxColorValue=255)
    names(col) = names(colors)
    return(col)
    }

################################################################################
# A problem is that canCoerce() will return TRUE when something actually is not
# coercible.  When you try to coerce it, it succeeds but with value NA and a warning
# message.  We will solve the problem by making our own "coerce" function that will
# muffle warnings and convert them to errors, and will catch errors and issue our own
# error.
#
# Arguments:
#   obj: object to be coerced to another type.
#   type: the class to which it is to be coerced.  If "integer", the type is first
#       coerced to "numeric", then if it isn't an integer numeric, an error is given,
#       after which it is coerced to "integer".
#   errfunc: error function to call when the object cannot be coerced.
#   ...: arguments to use to call errfunc() if there is an error.  If none, a default
#       error message is used as the argument to errfunc().
#   allowNA: TRUE to allow NAs in result, FALSE not to.
#   name: the name of the object, to be used in the default error message.  If NULL,
#       the calling argument name is used.
#
# Returns: obj coerced to type.
#
# Notes: This is a general-purpose function even though it is currently used only here.
################################################################################
coerceObj = function(obj, type, errfunc, ..., allowNA=FALSE, name=NULL)
    {
    if (!is(obj, type))
        {
        # We need to get 'name' ahead of time in case it is needed.
        if (is.null(name))
            name = deparse(substitute(obj))

        # Define our own local error function.
        L = list(...)
        ourError = function()
            {
            if (is.null(L) || length(L) == 0 || is.null(L[1]))
                errfunc(name, " must be of or coercible to type '", type, "'",
                    ifelse(!allowNA, " and not NA", ""), call.=FALSE)
            else
                {
                L[["call."]] = FALSE
                do.call(errfunc, L)
                }
            }

        # Function to try to convert obj to type using as().  Catch errors and call
        # our own error function instead.
        convertObj = function(obj, type)
            {
            x = tryCatch(as(obj, type), error=function(e) e)
            if (is(x, "error"))
                ourError()
            return(x)
            }

        # Use "numeric" in place of "integer".
        cvtTo = ifelse(type == "integer", "numeric", type)

        # Establish simpleWarning handler that converts a warning into an error,
        # and call convertObj() with that handler.
        withCallingHandlers(simpleWarning=function(w) ourError(),
            { obj = convertObj(obj, cvtTo) })

        # Check for NA if disallowed.
        if (!allowNA)
            if (any(is.na(obj)))
                ourError()

        # If original type was integer, make that conversion from numeric.
        if (type == "integer")
            {
            if (any(obj != as.integer(obj), na.rm=TRUE))
                ourError()
            obj = as.integer(obj)
            }
        }
    return(obj)
    }

################################################################################
# Read and verify parameter file, and convert parameters to their correct type.
################################################################################

S = readLines(paramFile)

# Remove comments.
hasQuote = grepl('^[^#]*"', S)
S[!hasQuote] = sub("#.*$", "", S[!hasQuote])
S[hasQuote] = sub('^([^#]*"[^"]*".*)#.*$', "\\1", S[hasQuote])

# Trim whitespace.
S = sub("^[ \t]+", "", S)
S = sub("[ \t]+$", "", S)

# Remove blank lines.
S = S[S != ""]

# Convert := to =.
S = sub(":=", "=", S)

# Remove whitespace around =.
S = sub("[ \t]+=[ \t]+", "=", S)

# Remove quote marks.
S = sub('="(.*)"', "=\\1", S)
S = sub("='(.*)'", "=\\1", S)

# Split on "="
S = strsplit(S, "=", fixed=TRUE)

# If any are missing the second element, set it empty.  Convert to a named vector.
S = sapply(S, function(V)
    {
    if (length(V) == 1)
        V[2] = ""
    t = V[2]
    names(t) = V[1]
    return(t)
    })

# Turn S into a list where the name is the parameter name, and the value is its
# value.  We use a list rather than a vector so we can change the types of elements
# to any class we want.
L = as.list(S)

# Expected parameters and their types.
expected = c(
    LCR.file="character",
    PLOT.genomes="character",
    PLOT.file="character",
    X.idlens="character",
    Y.idlens="character",
    PLOT.width="integer",
    PLOT.height="integer",
    PLOT.margin.top="numeric",
    PLOT.margin.bottom="numeric",
    PLOT.margin.left="numeric",
    PLOT.margin.right="numeric",
    BKGD.color="character",
    FGND.color="character",
    PLOT.which="character",
    POINT.size="numeric",
    LINE.width="numeric",
    PTLINE.color="character",
    PTLINE.alpha="numeric",
    X.grid.lines.per.seqid="integer",
    Y.grid.lines.per.seqid="integer",
    GRID.line.color="character",
    GRID.line.width="numeric",
    X.grid.label.position="numeric",
    Y.grid.label.position="numeric",
    GRID.label.color="character",
    GRID.label.size="numeric",
    X.plot.seqids="character",
    Y.plot.seqids="character",
    X.plot.start="numeric",
    Y.plot.start="numeric",
    X.plot.end="numeric",
    Y.plot.end="numeric",
    X.chr.color="character",
    Y.chr.color="character",
    X.chr.width="numeric",
    Y.chr.width="numeric",
    X.chr.spacing="numeric",
    Y.chr.spacing="numeric",
    PLOT.chr.label="logical",
    X.chr.label.position="numeric",
    Y.chr.label.position="numeric",
    X.chr.label.color="character",
    Y.chr.label.color="character",
    X.chr.label.size="numeric",
    Y.chr.label.size="numeric",
    X.chr.label.RE="character",
    X.chr.label.replace="character",
    Y.chr.label.RE="character",
    Y.chr.label.replace="character",
    X.axis.label.text="character",
    Y.axis.label.text="character",
    X.axis.label.position="numeric",
    Y.axis.label.position="numeric",
    X.axis.label.color="character",
    Y.axis.label.color="character",
    X.axis.label.size="numeric",
    Y.axis.label.size="numeric")

# Test for presence of all expected names and no unexpected names.
missing = setdiff(names(expected), names(L))
if (length(missing) > 0)
    cat("These parameters were expected but were not found:\n", paste("   ", missing, collapse="\n"))

unexpected = setdiff(names(L), names(expected))
if (length(unexpected) > 0)
    cat("These unknown parameters were not expected but were found:\n", paste("   ", unexpected, collapse="\n"))

if (length(missing) > 0 || length(unexpected) > 0)
    stop("Fix bad parameter file and re-run")

# Test for presence of all expected names and no unexpected names.
missing = setdiff(names(expected), names(S))

# Try to convert each value to its expected class.
for (n in names(L))
    L[[n]] = coerceObj(L[[n]], expected[n], errFunc=stop, allowNA=FALSE, name=n)
   
# Do additional testing and conversions on specific items.
if (!file.exists(L$LCR.file))
    stop("File ", L$LCR.file, " does not exist")

if (nchar(L$PLOT.genomes) != 2)
    stop("Must specify exactly two characters for PLOT.genomes, got: ", L$PLOT.genomes)

if (!file.exists(L$X.idlens))
    stop("File ", L$X.idlens, " does not exist")

if (!file.exists(L$Y.idlens))
    stop("File ", L$Y.idlens, " does not exist")

if (!L["PLOT.which"] %in% c("points", "lines"))
    stop("PLOT.which must be either 'points' or 'lines'")

L$X.plot.seqids = strsplit(gsub(" ", "", L$X.plot.seqids), ",", fixed=TRUE)[[1]]
L$Y.plot.seqids = strsplit(gsub(" ", "", L$Y.plot.seqids), ",", fixed=TRUE)[[1]]

################################################################################
# Read data files.
################################################################################

scaleToMbp = 1e-6

# Read the LCR data file.
catnow("Reading LCR file ", L$LCR.file, "...", sep="")
df = read.table(L$LCR.file, header=TRUE, sep="\t", row.names=1, stringsAsFactors=FALSE)
catnow("\n")
if (nrow(df) == 0) stop("Error reading ", L$LCR.file, ", no data found")
rownames(df) = NULL
# Drop unwanted columns.  We want X.seqID and X.pos columns where X is specified
# by PLOT.genomes.  Rename those columns "id.x", "x", "id.y", "y".
genomeLtrs = unlist(strsplit(L$PLOT.genomes, "", fixed=TRUE))
wantCols = c(paste(genomeLtrs, ".seqID", sep=""), paste(genomeLtrs, ".pos", sep=""), "LCR")
if (!all(wantCols %in% colnames(df)))
    stop("Expected columns '", paste(wantCols, collapse=","), "' but got columns '",
        paste(colnames(df), collapse=","), "'")
df = df[, wantCols]
colnames(df) = c("id.x", "id.y", "x", "y", "LCR")
catnow("Number of points in file: ", nrow(df), "\n")
# Convert the x and y positions to Mbp.
df$x = df$x*scaleToMbp
df$y = df$y*scaleToMbp

# Genome 1 idlens file.
genome1idlens = read.table(L$X.idlens, sep="\t", header=TRUE, stringsAsFactors=FALSE)
if (length(L$X.plot.seqids) != 0)
    genome1idlens = genome1idlens[genome1idlens$id %in% L$X.plot.seqids,]
if (nrow(genome1idlens) == 0)
    stop("No genome 1 sequence IDs to plot (names specified wrong?)")
L$X.plot.seqids = genome1idlens$id
genome1SeqLens = genome1idlens$len*scaleToMbp
names(genome1SeqLens) = genome1idlens$id
Ngenome1seqs = length(genome1SeqLens)
genome1SeqNames = names(genome1SeqLens)

# Genome 2 idlens file.
genome2idlens = read.table(L$Y.idlens, sep="\t", header=TRUE, stringsAsFactors=FALSE)
if (length(L$Y.plot.seqids) != 0)
    genome2idlens = genome2idlens[genome2idlens$id %in% L$Y.plot.seqids,]
if (nrow(genome2idlens) == 0)
    stop("No genome 2 sequence IDs to plot (names specified wrong?)")
L$Y.plot.seqids = genome2idlens$id
genome2SeqLens = genome2idlens$len*scaleToMbp
names(genome2SeqLens) = genome2idlens$id
Ngenome2seqs = length(genome2SeqLens)
genome2SeqNames = names(genome2SeqLens)

################################################################################
# Determine sequence positions along axes.
################################################################################

# Genome 1 x-axis span in Mbp, including spaces between chromosomes.
genome1SeqLensPlusSpace = genome1SeqLens + L$X.chr.spacing
genome1SeqLensPlusSpace[Ngenome1seqs] = genome1SeqLensPlusSpace[Ngenome1seqs] - L$X.chr.spacing

# Genome 2 y-axis span in Mbp, including spaces between chromosomes.
genome2SeqLensPlusSpace = genome2SeqLens + L$Y.chr.spacing
genome2SeqLensPlusSpace[Ngenome2seqs] = genome2SeqLensPlusSpace[Ngenome2seqs] - L$Y.chr.spacing

# Genome 1 x-axis starting and eding position of each plotted sequence.
genome1SeqStart = c(0, cumsum(genome1SeqLensPlusSpace)[-Ngenome1seqs])
names(genome1SeqStart) = genome1SeqNames
genome1SeqEnd = genome1SeqStart + genome1SeqLens

# Genome 2 y-axis starting and ending position of each plotted sequence.
genome2SeqStart = c(0, cumsum(genome2SeqLensPlusSpace)[-Ngenome2seqs])
names(genome2SeqStart) = genome2SeqNames
genome2SeqEnd = genome2SeqStart + genome2SeqLens

# If just one sequence is plotted and only a portion of the sequence is to be
# plotted, adjust the sequence lengths, axis spans and start/end positions.

# x-axis.
if (Ngenome1seqs == 1 && L$X.plot.start > 0 && L$X.plot.end > L$X.plot.start && L$X.plot.start < genome1SeqLens)
    {
    if (L$X.plot.end > genome1SeqLens)
        L$X.plot.end = genome1SeqLens
    genome1SeqLens[1] = L$X.plot.end - L$X.plot.start
    genome1SeqLensPlusSpace[1] = genome1SeqLens[1]
    genome1SeqStart[1] = L$X.plot.start
    genome1SeqEnd[1] = L$X.plot.end
    }

# y-axis.
if (Ngenome2seqs == 1 && L$Y.plot.start > 0 && L$Y.plot.end > L$Y.plot.start && L$Y.plot.start < genome2SeqLens)
    {
    if (L$Y.plot.end > genome2SeqLens)
        L$Y.plot.end = genome2SeqLens
    genome2SeqLens[1] = L$Y.plot.end - L$Y.plot.start
    genome2SeqLensPlusSpace[1] = genome2SeqLens[1]
    genome2SeqStart[1] = L$Y.plot.start
    genome2SeqEnd[1] = L$Y.plot.end
    }

# x-axis limits.
xlim.start = min(genome1SeqStart)
xlim.end = max(genome1SeqEnd)
xlim = c(xlim.start, xlim.end)

# y-axis limits.
ylim.start = min(genome2SeqStart)
ylim.end = max(genome2SeqEnd)
ylim = c(ylim.start, ylim.end)

################################################################################
# Map the position data to the axes.
################################################################################

# Get only the k-mer data we need to plot.
df = df[df$id.x %in% L$X.plot.seqids & df$id.y %in% L$Y.plot.seqids,]

# Map the data x and y positions from sequence-relative to axis-origin-relative.
df$x = df$x + genome1SeqStart[df$id.x] - xlim[1]
df$y = df$y + genome2SeqStart[df$id.y] - ylim[1]

# Remove any data not strictly within plotting limits.
df = df[df$x >= xlim[1] & df$x <= xlim[2] & df$y >= ylim[1] & df$y <= ylim[2],]
catnow("Number of points to plot: ", nrow(df), "\n")

################################################################################
# Now plot the data.
################################################################################

# The plot may have millions of points, which in a pdf file are stored as individual
# vector commands, making the file huge and unloadable.  To get around this, we
# will plot to a .png file, which can be examined nicely in Preview.
png(L$PLOT.file, width=L$PLOT.width, height=L$PLOT.height, antialias="none", bg=L$BKGD.color)

# Set margins and default colors.
par(mai=c(L$PLOT.margin.bottom, L$PLOT.margin.left, L$PLOT.margin.top, L$PLOT.margin.right))
par(bg=L$BKGD.color)
par(col=L$FGND.color)
par(col.axis=L$FGND.color)
par(col.lab=L$FGND.color)
par(col.main=L$FGND.color)
par(col.sub=L$FGND.color)

# Open the plot page.
plot(NA, type="n", xlim=xlim, ylim=ylim, xlab="", ylab="", axes=FALSE, frame.plot=FALSE)

# Plot the chromosomes.

# x-axis
segments(genome1SeqStart, ylim[1], genome1SeqEnd, ylim[1], lwd=L$X.chr.width, lend=1, col=L$X.chr.color)
# y-axis
segments(xlim[1], genome2SeqStart, xlim[1], genome2SeqEnd, lwd=L$Y.chr.width, lend=1, col=L$Y.chr.color)

# Plot the chromosome labels.

if (L$PLOT.chr.label)
    {
    # x-axis
    txt = genome1SeqNames
    if (L$X.chr.label.RE != "")
        txt = sub(L$X.chr.label.RE, L$X.chr.label.replace, txt)
    at = (genome1SeqStart+genome1SeqEnd)/2
    ok = (at >= xlim[1] & at <= xlim[2])
    at = at[ok]
    txt = txt[ok]
    if (length(txt) > 0)
        mtext(txt, 1, at=at, line=L$X.chr.label.position, col=L$X.chr.label.color, cex=L$X.chr.label.size)

    # y-axis
    txt = genome2SeqNames
    if (L$Y.chr.label.RE != "")
        txt = sub(L$Y.chr.label.RE, L$Y.chr.label.replace, txt)
    at = (genome2SeqStart+genome2SeqEnd)/2
    ok = (at >= ylim[1] & at <= ylim[2])
    at = at[ok]
    txt = txt[ok]
    if (length(txt) > 0)
        mtext(txt, 2, at=at, line=L$Y.chr.label.position, col=L$Y.chr.label.color, cex=L$Y.chr.label.size)
    }

# Plot the axis labels.

# x-axis
if (L$X.axis.label.text != "")
    mtext(L$X.axis.label.text, 1, line=L$X.axis.label.position, col=L$X.axis.label.color, cex=L$X.axis.label.size)
# y-axis
if (L$Y.axis.label.text != "")
    mtext(L$Y.axis.label.text, 2, line=L$Y.axis.label.position, col=L$Y.axis.label.color, cex=L$Y.axis.label.size)

# Plot the grid if requested.

# x-axis
if (L$X.grid.lines.per.seqid > 0)
    {
    # Vertical lines.
    Nlines = L$X.grid.lines.per.seqid
    col = L$GRID.line.color
    lwd = L$GRID.line.width

    if (Nlines == 1)
        x.offset = genome1SeqLens/2
    else
        {
        spacing = genome1SeqLens/(Nlines-1)
        x.offset = sapply(1:Ngenome1seqs, function(i) seq(0, by=spacing[i], length.out=Nlines), simplify=FALSE)
        }
    x = sapply(1:Ngenome1seqs, function(i) genome1SeqStart[i] + x.offset[[i]])
    x = unlist(x)
    segments(x, ylim[1], x, ylim[2], col=col, lwd=lwd)

    # Grid labels.  Omit label 0 if more than one sequence ID, to reduce crowding.
    if (L$GRID.label.size > 0)
        {
        x.pos = signif(xlim[1]+unlist(x.offset), 3) # xlim[1] for case where part of one chr is plotted.
        if (Ngenome1seqs > 1)
            {
            x = x[x.pos != 0]
            x.pos = x.pos[x.pos != 0]
            }
        mtext(x.pos, 1, las=2, line=L$X.grid.label.position, at=x, col=L$GRID.label.color, cex=L$GRID.label.size)
        }
    }

# y-axis
if (L$Y.grid.lines.per.seqid > 0)
    {
    # Horizontal lines.
    Nlines = L$Y.grid.lines.per.seqid
    col = L$GRID.line.color
    lwd = L$GRID.line.width

    if (Nlines == 1)
        y.offset = genome2SeqLens/2
    else
        {
        spacing = genome2SeqLens/(Nlines-1)
        y.offset = sapply(1:Ngenome2seqs, function(i) seq(0, by=spacing[i], length.out=Nlines), simplify=FALSE)
        }
    y = sapply(1:Ngenome2seqs, function(i) genome2SeqStart[i] + y.offset[[i]])
    y = unlist(y)
    segments(xlim[1], y, xlim[2], y, col=col, lwd=lwd)

    # Grid labels.  Omit label 0 if more than one sequence ID, to reduce crowding.
    if (L$GRID.label.size > 0)
        {
        y.pos = signif(ylim[1]+unlist(y.offset), 3) # ylim[1] for case where part of one chr is plotted.
        if (Ngenome2seqs > 1)
            {
            y = y[y.pos != 0]
            y.pos = y.pos[y.pos != 0]
            }
        mtext(y.pos, 2, las=2, line=L$Y.grid.label.position, at=y, col=L$GRID.label.color, cex=L$GRID.label.size)
        }
    }

# Finally, plot the dots.

if (nrow(df) > 0)
    {
    catnow("Plotting", nrow(df), "dots...")
    col = addAlpha(L$PTLINE.color, L$PTLINE.alpha)

    if (L$PLOT.which == "points")
        points(df$x, df$y, cex=L$POINT.SIZE, pch=20, col=col)
    else
        x = tapply(1:nrow(df), list(df$id.x, df$LCR), function(ii) lines(df$x[ii], df$y[ii], lwd=L$LINE.width, col=col))
    catnow("\n")
    }

dev.off()
cat("Finished dot plot.\n")
}

################################################################################
# End of file.
################################################################################
