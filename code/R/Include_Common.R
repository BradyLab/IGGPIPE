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
# End of file.
################################################################################
