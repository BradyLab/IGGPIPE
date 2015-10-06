################################################################################
# Code common to more than one R file in the IGGPIPE pipeline, sourced by those
# files.
# Author: Ted Toal
# Date: 2015
# Brady Lab, UC Davis
################################################################################

# For some reason this is needed when using RScript, sometimes.
library("methods")

################################################################################
# cat() that immediately flushes to console.
################################################################################
catnow = function(...)
    {
    cat(...)
    flush.console()
    return(invisible(0))
    }

################################################################################
# Crude function to print an R object's value, for simple debugging.
################################################################################
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

################################################################################
# Crude function to display info when investigating, for simple debugging.
################################################################################
inv = function(a, title="") { if (investigate) objPrint(a, title) }

################################################################################
# An alternative to rbind() that is far, far faster.
#
# Arguments:
#   df1: First data frame to be bound together row-wise, or a list previously
#       returned by this function, or NULL to indicate an empty data frame.
#   df2: Second data frame to be bound together row-wise, or a list previously
#       returned by this function, or NULL to indicate an empty data frame.
#
# Returns: a list containing the data in df1 and df2.  To convert this list to
#   a data frame, use rbind.fast.finish().  The returned list has an attribute
#   named "nrow", containing the number of rows in total.  Do not add to or
#   subtract from the returned lists yourself.
#
# Note: When df1 or df2 is a data frame, if the first row name is not composed
# entirely of digits, the row names are saved in sublist element ".rownames",
# and are also accumulated when each new data frame is added.  Both data frames
# must have row names if they are present at all, and the row names must all be
# unique.  No checking is done here for that.
#
# Example of speed:
#
# N = 10000
# Na = 1:N
# Nb = (N+1):(2*N)
#
# f1 = function()
#   {
#   df1 = NULL
#   for (i in 1:100)
#       df1 = rbind(df1, data.frame(a=Na, b=Nb, row.names=paste("A", i, "_", Na, sep="")))
#   return(df1)
#   }
#
# f2 = function()
#   {
#   df2 = NULL
#   for (i in 1:100)
#       df2 = rbind.fast(df2, data.frame(a=Na, b=Nb, row.names=paste("A", i, "_", Na, sep="")))
#   df2 = rbind.fast.finish(df2)
#   return(df2)
#   }
#
# system.time(df1 <- f1())
#    user  system elapsed 
#   7.901   3.395  11.460 
# system.time(df2 <- f2())
#    user  system elapsed 
#   0.056   0.004   0.070 
# identical(df1, df2)
# [1] TRUE
################################################################################
rbind.fast = function(df1, df2)
    {
    if (!(is.null(df1) || is.data.frame(df1) || (is.list(df1) && is.list(df1[[1]]))) ||
        !(is.null(df2) || is.data.frame(df2) || (is.list(df2) && is.list(df2[[1]]))))
        stop("rbind.fast only works with NULL or data frame or list of lists arguments")

    # Convert df1 from data frame to list if not already done.  Each list element
    # is itself a list.
    if (is.data.frame(df1))
        {
        nr = nrow(df1)
        R = rownames(df1)
        df1 = as.list(df1)
        for (n in names(df1))
            df1[[n]] = list(df1[[n]])
        if (length(R) > 0 && !grepl("^[0-9]+$", R[1]))
            df1$.rownames = R
        attr(df1, "nrow") = nr
        }

    # Convert df2 from data frame to list if not already done.  Each list element
    # is itself a list.
    if (is.data.frame(df2))
        {
        nr = nrow(df2)
        R = rownames(df2)
        df2 = as.list(df2)
        for (n in names(df2))
            df2[[n]] = list(df2[[n]])
        if (length(R) > 0 && !grepl("^[0-9]+$", R[1]))
            df2$.rownames = R
        attr(df2, "nrow") = nr
        }

    # Handle NULL.
    if (is.null(df1))
        return(df2)
    if (is.null(df2))
        return(df1)

    # rbind df2 to df1 by extending each sublist of df1 by the corresponding
    # sublist of df2.
    if (any(names(df1) == ".rownames") != any(names(df2) == ".rownames"))
        stop("rbind.fast: one but not both data frames have row names")
    for (n in names(df2)) # Use df2 in case df1 is NULL
        df1[[n]] = c(df1[[n]], df2[[n]])
    attr(df1, "nrow") = attr(df1, "nrow") + attr(df2, "nrow")
    return(df1)
    }

################################################################################
# Convert the list returned by rbind.fast back into the appropriate data frame.
#
# Arguments:
#   df: List returned by rbind.fast(), or a data frame.
#   stringsAsFactors: passed to as.data.frame().
#
# Returns: a data frame corresponding to df.  If df is already a data frame, it
# is returned, else the list is collapsed back into a data frame.
#
# Note: if df is a list and has element ".rownames", these are made into row
# names of the returned data frame.
#
# Note: see rbind.fast().
################################################################################
rbind.fast.finish = function(df, stringsAsFactors=FALSE)
    {
    if (is.list(df))
        {
        R = NULL
        if (any(names(df) == ".rownames"))
            {
            R = df$.rownames
            df$.rownames = NULL
            }
        for (n in names(df))
            df[[n]] = unlist(df[[n]], use.names=FALSE)
        df = as.data.frame(df, row.names=R, stringsAsFactors=FALSE)
        }
    return(df)   
    }

################################################################################
# Does do.call(rbind, <data frame list>) using rbind.fast to vastly speed it up.
#
# Arguments:
#   L: List of data frames to be rbound.
#   stringsAsFactors: passed to as.data.frame().
#
# Returns: a data frame containing all data frames in list L, bound together
# one below the other, i.e. the number of columns stays constant.
################################################################################
do.call.rbind.fast = function(L, stringsAsFactors=FALSE)
    {
    dfAll = NULL
    for (df in L)
        dfAll = rbind.fast(dfAll, df)
    return(rbind.fast.finish(dfAll, stringsAsFactors))
    }

################################################################################
# Like pretty() without all the extra arguments, except that this finds a
# sequence of "round" values appropriate when the axis is logarithmic.  The
# returned sequence includes powers of 10, and optionally, 2 * powers of 10,
# and 5 * powers of 10.
#
# Arguments:
#   x: an object coercible to numeric().
#   include: one of the values "X1", "X12", "X15", or "X125" indicating inclusion
#       of powers of 10 only (X1), 1* and 2* powers of 10 (X12), 1* and 5* (X15),
#       or 1*, 2*, and 5* (X125).
#   P10.StartEnd: TRUE if sequence must start and end with a power of 10, FALSE
#       if it may start or end with 2* or 5* a power of 10 (depending of course
#       on 'include').
#
# Returns: the computed sequence.
################################################################################
pretty.log = function(x, include="X125", P10.StartEnd=FALSE)
    {
    x = as.numeric(x)
    r = range(x)

    # Compute first element of the sequence.
    x1 = log10(r[1])
    x1 = ifelse(x1 < 0, as.integer(x1)-1, as.integer(x1))
    x1 = 10^x1
    # Raise first sequence element by *5 if possible, else *2 if possible.
    cur = "atX1"
    if (!P10.StartEnd && (include %in% c("X15", "X125")) && (x1*5 <= r[1]))
        {
        x1 = x1*5
        cur = "atX5"
        }
    else if (!P10.StartEnd && (include %in% c("X12", "X125")) && (x1*2 <= r[1]))
        {
        x1 = x1*2
        cur = "atX2"
        }

    # Compute last element of the sequence.
    x2 = log10(r[2])
    x2 = ifelse(x2 < 0, as.integer(x2), as.integer(x2)+1)
    x2 = 10^x2
    # Lower last sequence element by *2/10 if possible, else *5/10 if possible.
    if (!P10.StartEnd && (include %in% c("X12", "X125")) && (x2*2/10 >= r[2]))
        x2 = x2*2/10
    else if (!P10.StartEnd && (include %in% c("X15", "X125")) && (x2*5/10 >= r[2]))
        x2 = x2*5/10

    # Value to assign to 'cur' when moving from one sequence element to the next,
    # depending on value of 'include' argument as well as value of 'cur'.
    nextCur = list(X1=c(atX1="atX1"),
                   X12=c(atX1="atX2", atX2="atX1"),
                   X15=c(atX1="atX5", atX5="atX1"),
                   X125=c(atX1="atX2", atX2="atX5", atX5="atX1"))

    # Value to multiply current sequence member by, to get next sequence member,
    # depending on value of 'include' argument as well as value of 'cur'.
    multFactor = list(X1=c(atX1=10),
                      X12=c(atX1=2, atX2=5),
                      X15=c(atX1=5, atX5=2),
                      X125=c(atX1=2, atX2=5/2, atX5=2))

    # Compute the full sequence.
    seq = x1
    xt = x1
    while (xt <= x2)
        {
        xt = xt * multFactor[[include]][cur]
        cur = nextCur[[include]][cur]
        seq = c(seq, xt)
        }
    names(seq) = NULL
    return(seq)
    }

################################################################################
# Conservative 8-color palette adapted for color blindness, with first color = "black".
#
# Wong, Bang. "Points of view: Color blindness." nature methods 8.6 (2011): 441-441.
################################################################################
colorBlind.8 = c(black="#000000", orange="#E69F00", skyblue="#56B4E9",
    bluegreen="#009E73", yellow="#F0E442", blue="#0072B2", vermillion="#D55E00",
    reddishpurple="#CC79A7")
# plot(NA, xlim=0:1, ylim=0:1, axes=FALSE, xlab="", ylab="")
# rect((0:7)/8, 0, (1:8)/8, 1, border=NA, col=colorBlind.8)

#######################################################################################
# Sample standard error of mean.
#######################################################################################
sem = function(V) sd(V)/sqrt(length(V))

################################################################################
# End of file.
################################################################################
