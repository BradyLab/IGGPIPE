#######################################################################################
# This file contains R definitions and functions for merging data frame data using
# positional information in the data frames.
# Author: Ted Toal
# Lab: Brady
# 2015
#######################################################################################

#######################################################################################
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
#######################################################################################
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

#######################################################################################
# Find the rows of df1 and df2 such that df2's position is within df1's contig.
#
# Arguments:
#   df1: data frame with columns start, end, and idx.  start and end define the starting
#       and ending positions of a contig, and idx is returned from each df1 row of all
#       (df1, df2) matching pairs of rows.
#   df2: data frame with columns like df1.  idx is returned from each df2 row of all
#       (df1, df2) matching pairs of rows.
#
# Returns: a vector with an even length, composed of consecutive pairs (x, y), where,
# for each row of df1 and row of df2 such that df1[["start"]] <= df2[["pos"]] <= df1[["end"]],
# x is df1[["idx"]] and y is df2[["idx"]].
#
# This operates in time proportional to nrow(df1)*nrow(df2) and so can be very
# slow for large data frames.  Use findContainsIdxs.recursive instead.
#######################################################################################
findContainsIdxs = function(df1, df2)
    {
    # cat("df2$idx[1]=", df2$idx[1], "\n")
    # This can be done using an "apply" function on either df1 or df2.  Do it on
    # the smaller of the two.
    if (nrow(df1) <= nrow(df2))
        {
        V = unlist(sapply(1:nrow(df1), function(i)
            {
            y = df2[["idx"]][df2[["pos"]] >= df1[["start"]][i] & df2[["pos"]] <= df1[["end"]][i]]
            x = rep(df1[["idx"]][i], length(y))
            return(c(rbind(x,y)))
            }))
        }
    else
        {
        V = unlist(sapply(1:nrow(df2), function(i)
            {
            x = df1[["idx"]][df2[["pos"]][i] >= df1[["start"]] & df2[["pos"]][i] <= df1[["end"]]]
            y = rep(df2[["idx"]][i], length(x))
            return(c(rbind(x,y)))
            }))
        }
    return(V)
    }

#######################################################################################
# Perform the same thing as findContainsIdxs(), but subdivide the problem to vastly
# increase efficiency (provided that the contigs in df1 generally do not all overlap
# each other).
#
# Arguments:
#   df1, df2: same as findContainsIdxs() arguments.
#   start, end: two positions between which all df1 and df2 positions lie.
#
# Returns: a vector with an even length, composed of consecutive pairs (x, y), where,
# for each row of df1 and row of df2 such that df1[["start"]] <= df2[["pos"]] <= df1[["end"]],
# x is df1[["idx"]] and y is df2[["idx"]].  To convert the return vector V into a 2-column matrix
# of index pairs, do this:
#   mtx = matrix(V, ncol=2, byrow=TRUE)
# The first column is indexes into df1 and the second column is indexes into df2.
#######################################################################################
findContainsIdxs.recursive = function(df1, df2, start, end)
    {
    # cat("df1$idx[1]=", df1$idx[1], " nrow(df1)=", nrow(df1),
    #    " df2$idx[1]=", df2$idx[1], "nrow(df2)=", nrow(df2), " start=", start, " end=", end, "\n")
    # If df1 and/or df2 and/or start..end is small enough, just use findContainsIdxs()
    # to do the job.
    N = end - start + 1
    if (nrow(df1) <= 10 || nrow(df2) <= 10 || N <= 10 || (nrow(df1) <= 100 && nrow(df2) <= 100))
        return(findContainsIdxs(df1, df2))

    # Otherwise, the approach is this: split the start..end interval into two equal
    # halves A and B, and split df1 into three sets: those in A, those in B, and
    # those which start in A and end in B (set AB).  Also split df2 into two sets,
    # those in A and those in B.  Then, test the three smaller sets, A, B, and AB,
    # since it will never be the case that a df1 contig in one set contains the df2
    # position in a different set.
    M = start + N/2
    df1.A = df1[df1[["end"]] <= M,]
    df1.B = df1[df1[["start"]] > M,]
    df1.AB = df1[df1[["start"]] <= M & df1[["end"]] > M,]
    df2.A = df2[df2[["pos"]] <= M,]
    df2.B = df2[df2[["pos"]] > M,]

    # Now test each of the sets A, B, and AB, provided neither df1 nor df2 is empty.
    V = NULL
    if (nrow(df1.A) > 0 && nrow(df2.A) > 0)
        V = c(V, findContainsIdxs.recursive(df1.A, df2.A, start, M))
    if (nrow(df1.B) > 0 && nrow(df2.B) > 0)
        V = c(V, findContainsIdxs.recursive(df1.B, df2.B, M, end))

    # Set AB has to be tested with findContainsIdx().  It needs to be tested with
    # all df2 rows, except that we might get an improvement in some cases by
    # eliminating obviously non-intersecting rows.
    if (nrow(df1.AB) > 0)
        {
        minStart = min(df1.AB[["start"]])
        maxEnd = max(df1.AB[["end"]])
        df2 = df2[df2[["pos"]] >= minStart & df2[["pos"]] <= maxEnd,]
        if (nrow(df2) > 0)
            V = c(V, findContainsIdxs(df1.AB, df2))
        }
    return(V)
    }

#######################################################################################
# Helper function for mergeOnMatches() below.  Find the rows of df1 and df2 such that
# df2's position is within df1's contig.
#
# Arguments:
#   df1: data frame with columns start and end that define the starting and ending
#       positions of a contig.
#   df2: data frame with columns like df1.
#
# Returns: a matrix with 2 columns.  The first column is row numbers of df1 and the
# second column is row numbers of df2.  Each row of the matrix identifies a pair of
# rows (in df1 and df2) that match.
#######################################################################################
findContainsIdxs.rows = function(df1, df2)
    {
    df1[["idx"]] = 1:nrow(df1)
    df2[["idx"]] = 1:nrow(df2)
    V = findContainsIdxs.recursive(df1, df2,
        min(df1[["start"]], df2[["pos"]]), max(df1[["end"]], df2[["pos"]]))
    mtx = matrix(V, ncol=2, byrow=TRUE)
    return(mtx)
    }

#######################################################################################
# Helper function for mergeOnMatches() below.  Find matches of positions in t.position
# and s.position according to the parameters in dist.
#
# Arguments:
#   t.position: data frame with columns start and end that define starting and ending
#       positions.
#   s.position: data frame with columns like t.position.
#   dist: dist argument of mergeOnMatches().
#
# Returns: a matrix with 2 columns.  The first column is indexes of t.position rows
# and the second column is indexes of s.position rows.  Each row of the matrix
# identifies a pair of rows (in t.position and s.position) that match.
#######################################################################################
getMatchIdxs = function(t.position, s.position, dist)
    {
    # Combine rows of indexes in two matrices of indexes to eliminate duplicates.
    # Return a matrix of unique indexes.
    getUniqueIdxs = function(idxs1, idxs2)
        {
        x = union(paste(idxs1[,1], idxs1[,2], sep="_"), paste(idxs2[,1], idxs2[,2], sep="_"))
        idxs = matrix(as.integer(unlist(strsplit(x, "_", fixed=TRUE))), ncol=2, byrow=TRUE)
        return(idxs)
        }

    # Find rows containing identical pairs of indexes in two matrices of indexes.
    # Return a matrix of the common indexes.
    getCommonIdxs = function(idxs1, idxs2)
        {
        x = intersect(paste(idxs1[,1], idxs1[,2], sep="_"), paste(idxs2[,1], idxs2[,2], sep="_"))
        idxs = matrix(as.integer(unlist(strsplit(x, "_", fixed=TRUE))), ncol=2, byrow=TRUE)
        return(idxs)
        }

    if (nrow(t.position) == 0) stop("getMatchIdxs: software error, t.position is empty")
    if (nrow(s.position) == 0) stop("getMatchIdxs: software error, s.position is empty")

    # Method depends on dist[["method"]].

    # "OVERLAP": Match only if s.df and t.df overlap by at least one nucleotide.
    #   If both are SNPs, the positions must match.
    #   0: match only if s.df and t.df overlap by at least one nucleotide.  If both
    #       are SNPs, the positions must match.
    if (dist[["method"]] == "OVERLAP")
        {
        # To test for overlap, note that when an overlap exists, one of the four endpoints
        # (start and end for t.position, start and end for s.position) must lie between
        # the two endpoints of the other position.  Test all four situations.
        s.position[["pos"]] = s.position[["start"]]
        idxs = findContainsIdxs.rows(t.position, s.position)

        s.position[["pos"]] = s.position[["end"]]
        idxs2 = findContainsIdxs.rows(t.position, s.position)
        idxs = getUniqueIdxs(idxs, idxs2)

        t.position[["pos"]] = t.position[["start"]]
        idxs2 = findContainsIdxs.rows(s.position, t.position)
        # We must swap the two idxs2 columns so that column 1 is for t.df.
        idxs = getUniqueIdxs(idxs, idxs2[,2:1])

        t.position[["pos"]] = t.position[["end"]]
        idxs2 = findContainsIdxs.rows(s.position, t.position)
        # We must swap the two idxs2 columns so that column 1 is for t.df.
        idxs = getUniqueIdxs(idxs, idxs2[,2:1])
        }

    # "s.TINY": The positions in s.df are either SNPs or very small contigs that
    #   are much smaller that the contigs in t.df, and matching requires that
    #   the t.df contig completely encompasses the s.df SNP or contig.
    else if (dist[["method"]] == "s.TINY")
        {
        s.position[["pos"]] = s.position[["start"]]
        idxs = findContainsIdxs.rows(t.position, s.position)
        if (any(s.position[["end"]] != s.position[["start"]]))
            {
            s.position[["pos"]] = s.position[["end"]]
            idxs2 = findContainsIdxs.rows(t.position, s.position)
            idxs = getCommonIdxs(idxs, idxs2)
            }
        }

    # "t.TINY": The opposite of s.TINY, swap s.df and t.df roles.
    else if (dist[["method"]] == "t.TINY")
        {
        t.position[["pos"]] = t.position[["start"]]
        idxs = findContainsIdxs.rows(s.position, t.position)
        if (any(t.position[["end"]] != t.position[["start"]]))
            {
            t.position[["pos"]] = t.position[["end"]]
            idxs2 = findContainsIdxs.rows(s.position, t.position)
            idxs = getCommonIdxs(idxs, idxs2)
            }
        # We must swap the two idxs columns so that column 1 is for t.df.
        idxs = idxs[,2:1]
        }

    # "x.NEAR", x = s/t, y = t/s: The positions in x.df are not large compared to
    #   those in y.df, either one may be a SNP or a CONTIG, and matching requires
    #   that the two are near to one another to the degree specified by the members
    #   "closest" and "start.up", "start.down", "end.up", "end.down".

    else # (dist[["method"]] == "s.NEAR" || dist[["method"]] == "t.NEAR")
        {
        # For t.NEAR we will swap s.position and t.position, and at the end, swap
        # the two columns of idxs.
        if (dist[["method"]] == "t.NEAR")
            {
            tmp = t.position
            t.position = s.position
            s.position = tmp
            }

        # Start by getting all indexes that satisfy start.up/start.down/end.up/end.down.

        # Here, we are using x=s, i.e. s.NEAR:
        # start.up: s.start must be no less than t.start-start.up
        # start.down: s.start must be no more than t.end+start.down
        # end.up: s.end must be no less than t.start-end.up
        # end.down: s.end must be no more than t.end+end.down

        dontCareUpDist = max(abs(diff(s.position[order(s.position[["start"]]), "start"])))
        dontCareDownDist = max(abs(diff(s.position[order(s.position[["end"]]), "end"])))
        #cat("dontCareUpDist=", dontCareUpDist, "\n")
        #cat("dontCareDownDist=", dontCareDownDist, "\n")

        s.position[["pos"]] = s.position[["start"]]
        tmp.position = t.position
        maxDist = dist[["start.up"]]
        if (is.na(maxDist))
            maxDist = dontCareUpDist
        tmp.position[["start"]] = tmp.position[["start"]] - maxDist
        maxDist = dist[["start.down"]]
        if (is.na(maxDist))
            maxDist = dontCareDownDist
        tmp.position[["end"]] = tmp.position[["end"]] + maxDist
        idxs = findContainsIdxs.rows(tmp.position, s.position)

        s.position[["pos"]] = s.position[["end"]]
        tmp.position = t.position
        maxDist = dist[["end.up"]]
        if (is.na(maxDist))
            maxDist = dontCareUpDist
        tmp.position[["start"]] = tmp.position[["start"]] - maxDist
        maxDist = dist[["end.down"]]
        if (is.na(maxDist))
            maxDist = dontCareDownDist
        tmp.position[["end"]] = tmp.position[["end"]] + maxDist
        idxs2 = findContainsIdxs.rows(tmp.position, s.position)
        idxs = getCommonIdxs(idxs, idxs2)

        # Now test dist[["closest"]] to decide how to further break down the method.

        # 0: ALL contigs satisfying the four position limits are taken as matches.
        #   else (dist[["closest"]] == 0).  This result is already computed in "idxs".
        # 1: Like "OVERLAP" but when there is no overlap, only the NEAREST to each t.df
        #   row of all non-overlapping matches that satisfy the four position limits
        #   is taken as a match.
        # 2: like 1, but allows one match upstream of t.df and a second downstream,
        #   in both cases the NEAREST one.

        # We already have the indexes in idxs for closest = 0.  For closest = 1 or 2,
        # processing is almost identical.
        if (dist[["closest"]] == 1 || dist[["closest"]] == 2)
            {
            # We must exclude from idxs those that are non-overlapping and are not
            # the nearest of the non-overlapping ones.

            # Compute distance upstream and downstream (of t.df) for each idxs row.
            dist.s.upstreamOf.t = t.position[idxs[,1], "start"] - s.position[idxs[,2], "end"]
            dist.s.downstreamOf.t = s.position[idxs[,2], "start"] - t.position[idxs[,1], "end"]

            # Determine which idxs rows have overlaps.
            overlaps = (dist.s.upstreamOf.t <= 0) & (dist.s.downstreamOf.t <= 0)

            # The remaining rows are cases of near but not overlapping.  We need the set
            # of non-overlapping idxs rows.  However, we want to exclude rows whose [,1]
            # index is present in idxs[overlaps,1], because when there IS an overlap, we
            # do NOT include ANY nearest one.
            no.overlaps = !overlaps & !(idxs[,1] %in% idxs[overlaps,1])

            # Split idxs into two sets, overlapping, and strictly non-overlapping.
            idxs.overlap = idxs[overlaps,]
            idxs.no.overlap = idxs[no.overlaps,]

            # Get upstream and downstream distances for the non-overlapping set.
            dist.s.upstreamOf.t = dist.s.upstreamOf.t[no.overlaps]
            dist.s.downstreamOf.t = dist.s.downstreamOf.t[no.overlaps]

            # If we have any non-overlapping indexes, we must handle them according to "closest".
            if (nrow(idxs.no.overlap) > 0)
                {
                # Handle the two "closest" cases separately from here.
                if (dist[["closest"]] == 1)
                    {
                    # We need the distance away, which for each index row will be either
                    # the dist.s.upstreamOf.t or dist.s.downstreamOf.t value, whichever one is
                    # not negative.
                    dist.t = pmax(dist.s.upstreamOf.t, dist.s.downstreamOf.t)
                    L = tapply(1:nrow(idxs.no.overlap), idxs.no.overlap[,1], function(ii)
                        {
                        # Get value of ii (idxs.no.overlap row number) of the row that has the
                        # smallest value in dist.t.
                        i = ii[which.min(dist.t[ii])]
                        return(idxs.no.overlap[i,])
                        })
                    }
                else # (dist[["closest"]] == 2)
                    {
                    # Split the data into upstream and downstream data.
                    upstream = (dist.s.upstreamOf.t > 0)
                    idxs.upstream = idxs.no.overlap[upstream,]
                    idxs.downstream = idxs.no.overlap[!upstream,]
                    dist.s.upstreamOf.t = dist.s.upstreamOf.t[upstream]
                    dist.s.downstreamOf.t = dist.s.downstreamOf.t[!upstream]

                    # Get separate lists of nearest upstream and nearest downstream index pairs.
                    L.upstream = tapply(1:nrow(idxs.upstream), idxs.upstream[,1], function(ii)
                        {
                        i = ii[which.min(dist.s.upstreamOf.t[ii])]
                        return(idxs.upstream[i,])
                        })

                    L.downstream = tapply(1:nrow(idxs.downstream), idxs.downstream[,1], function(ii)
                        {
                        i = ii[which.min(dist.s.downstreamOf.t[ii])]
                        return(idxs.downstream[i,])
                        })

                    L = c(L.upstream, L.downstream)
                    }
                # Convert list L of one-row data frames to a single data frame
                # which is the new idxs.no.overlap.
                idxs.no.overlap = do.call(rbind, L)
                }

            # Bind idxs.overlap and idxs.no.overlap to form new idxs matrix.
            idxs = rbind(idxs.overlap, idxs.no.overlap)
            }

        # For t.NEAR we now must swap the two columns of idxs.
        if (dist[["method"]] == "t.NEAR")
            idxs = idxs[,2:1]
        }

    # Sort idxs by column 1 then column 2.
    idxs = idxs[order(idxs[,1], idxs[,2]),]
    return(idxs)
    }

#######################################################################################
# Helper function for formatData() below.  This generates strings for the format codes
# "#t", "#s", "%t", "%s".
#
# Arguments:
#   code: one of "#t", "#s", "%t", "%s".
#   t.position: see formatData() t.position argument.
#   s.position: see formatData() s.position argument.
#   idxs: a 2-column matrix of row indexes into t.df (column 1) and s.df (column 2) of
#       the matching pairs.
#
# Returns:
#   A vector of data formatted according to "code", and of length equal to nrow(idxs),
#   to become part of a new t.df column.
#######################################################################################
positionString = function(code, t.position, s.position, idxs)
    {
    #       {#t} : bp position of s.start in t.df contig, as: -#, +#, or @#, see below.
    #       {#s} : bp position of t.start in s.df contig, as: -#, +#, or @#, see below.
    #       {%t} : percent position of s.start in t.df contig, as: -#%, +#%, or @#%, see below.
    #       {%s} : percent position of t.start in s.df contig, as: -#%, +#%, or @#%, see below.
    #   A start position may lie UPSTREAM, WITHIN, or DOWNSTREAM of a contig, and the prefix
    #   characters "-", "@", and "+", respectively, are used in the four distance/percent
    #   specifiers above to indicate which one is the case.

    # Make it so it is always like {#t} or {%t}.  Must swap idxs columns also.
    if (code == "#s" || code == "%s")
        {
        t = t.position
        t.position = s.position
        s.position = t
        rm(t)
        idxs = idxs[,2:1]
        }

    # Figure out {#t}.
    len = (t.position[["end"]] - t.position[["start"]] + 1)[idxs[,1]]
    amount = (s.position[["start"]][idxs[,2]] - t.position[["start"]][idxs[,1]])
    isUpstream = (amount < 0)
    isDownstream = (amount >= len)
    isWithin = (!isUpstream & !isDownstream)
    amount[isDownstream] = amount[isDownstream] - len[isDownstream] + 1

    # Figure out {%t}.
    if (code == "%t" || code == "%s")
        amount = as.integer(100*amount/len)

    # Make character strings.
    S = as.character(amount)
    # S[isUpstream] = paste("-", S[isUpstream], sep="") # Already has a "-" sign.
    S[isDownstream] = paste("+", S[isDownstream], sep="")
    S[isWithin] = paste("@", S[isWithin], sep="")
    if (code == "%t" || code == "%s")
        S = paste(S, "%", sep="")
    return(S)
    }

#######################################################################################
# Helper function for mergeOnMatches() below.  This parses cols[[i]][["format"]] strings
# and checks their validity, generating an error if invalid.  Optionally it creates the
# new column values using the format string and match data.
#
# Arguments:
#   format: a cols[[i]][["format"]] string.
#   t.colnames: the column names of the t.df data frame.
#   s.colnames: the column names of the s.df data frame.
#   t.df: NULL for parsing only, else the t.df data frame.
#   s.df: NULL for parsing only, else the s.df data frame.
#   idxs: NULL for parsing only, else a 2-column matrix of row indexes into t.df
#       (column 1) and s.df (column 2) of the matching pairs.
#   t.position: NULL for parsing only, else data frame with same number of rows as t.df,
#       and columns "start" and "end" giving positions of t.df row contigs or SNPs.
#   s.position: NULL for parsing only, else data frame with same number of rows as s.df,
#       and columns "start" and "end" giving positions of s.df row contigs or SNPs.
#
# Returns:
#   If parsing only, NULL is returned invisibly.  Otherwise, a vector of
#   data formatted according to "format", and of length equal to nrow(idxs).
#######################################################################################
formatData = function(format, t.colnames, s.colnames, t.df=NULL, s.df=NULL, idxs=NULL,
    t.position=NULL, s.position=NULL)
    {
    error = function(...) stop("formatData: ", ..., call.=FALSE)

    chk.col = function(col, st)
        {
        if (st == "s" && !any(col == s.colnames))
            error("'", col, "' is not a column name of s.df")
        if (st == "t" && !any(col == t.colnames))
            error("'", col, "' is not a column name of t.df")
        }

    # Function to try to do sub() with given arguments, but catch errors and
    # warnings and call our own error function with ... instead.
    sub.catch = function(RE, RE.replace, x, ...)
        {
        # Establish simpleWarning handler that converts a warning into an error,
        # and use tryCatch to catch errors.
        withCallingHandlers(simpleWarning=function(w) error(...),
            {
            x = tryCatch(sub(RE, RE.replace, x), error=function(e) e, warning=function(w) w)
            if (is(x, "error") || is(x, "warning"))
                error(...)
            })
        return(x)
        }

    # If generating formatted strings and there are no indexes, return an empty vector.
    generate = !is.null(t.df)
    if (generate && nrow(idxs) == 0)
        return(character())

    # Split on "{" or "}".
    V = unlist(strsplit(format, "[{}]"))
    if (length(V) %% 2 == 1)
        V = c(V, "")
    if (length(V) == 0) stop("formatData: software error, V is empty")
    V = matrix(V, ncol=2, byrow=TRUE, dimnames=list(NULL, c("verbatim", "format")))

    # Test each format code for validity.  If requested, generate formatted strings.
    fmtstrs = NULL
    if (generate)
        fmtstrs = rep("", nrow(idxs))
    for (i in 1:nrow(V))
        {
        code = V[i,2]

        # Add the verbatim text.
        if (generate)
            fmtstrs = paste(fmtstrs, V[i, "verbatim"], sep="")

        # Get first character of the code after "{"
        firstChar = substring(code, 1, 1)

        # Process the code.

        # {lb} : a verbatim left brace
        if (code == "lb")
            {
            if (generate)
                fmtstrs = paste(fmtstrs, "{", sep="")
            }

        # {rb} : a verbatim right brace
        else if (code == "rb")
            {
            if (generate)
                fmtstrs = paste(fmtstrs, "}", sep="")
            }

        # {#t}, {#s}, {%t}, {%s}
        else if (code == "#t" || code == "#s" || code == "%t" || code == "%s")
            {
            if (generate)
                fmtstrs = paste(fmtstrs, positionString(code, t.position, s.position, idxs), sep="")
            }

        # {+col} : the value in s.df column "col" of the matching row
        else if (firstChar == "+")
            {
            col = sub("^\\+(.*)$", "\\1", code)
            chk.col(col, "s")
            if (generate)
                fmtstrs = paste(fmtstrs, s.df[,col][idxs[,2]], sep="")
            }

        # {*s*col*val*dgts} and {*t*col*val*dgts}
        else if (firstChar == "*")
            {
            s = unlist(strsplit(code, "*", fixed=TRUE))
            if (length(s) != 5)
                error("expected {*s|t*col*val*dgts} but got: ", code)
            if (s[2] != "s" && s[2] != "t")
                error("expected {*s... or {*t... but got: ", code)
            col = s[3]
            chk.col(col, s[2])
            val = coerceObj(s[4], "numeric",
                error, "expected numeric val in {*s|t*col*val*dgts} but got ", s[4])
            dgts = coerceObj(s[5], "integer",
                error, "expected integer dgts in {*s|t*col*val*dgts} but got ", s[5])
            if (generate)
                {
                if (s[2] == "s")
                    fmtstrs = paste(fmtstrs, round(s.df[,col][idxs[,2]]*val, dgts), sep="")
                else
                    fmtstrs = paste(fmtstrs, round(t.df[,col][idxs[,1]]*val, dgts), sep="")
                }
            }

        # {/s/col/RE/RE.replace} and {/t/col/RE/RE.replace} :
        # regular expression search/replace of s.df column "col"
        else if (firstChar == "/")
            {
            s = unlist(strsplit(code, "/", fixed=TRUE))
            if (length(s) == 4)
                s[5] = ""
            if (length(s) != 5)
                error("expected {/s|t/col/RE/RE.replace} but got: ", code)
            if (s[2] != "s" && s[2] != "t")
                error("expected {*s... or {*t... but got: ", code)
            col = s[3]
            chk.col(col, s[2])
            RE = s[4]
            RE.replace = s[5]
            if (generate)
                {
                if (s[2] == "s")
                    S = sub.catch(RE, RE.replace, s.df[,col][idxs[,2]])
                else
                    S = sub.catch(RE, RE.replace, t.df[,col][idxs[,1]])
                fmtstrs = paste(fmtstrs, S, sep="")
                }
            }

        # else if not empty, it is an invalid code.
        else if (code != "")
            {
            error("invalid format code: '", code, "'")
            }
        }

    # Return result.
    if (!generate)
        return(invisible(NULL))
    return(fmtstrs)
    }

#######################################################################################
# Find position matches between rows of two different data frames, where each data
# frame contains columns identifying a position.  Copy specified column data from the
# matching rows of one data frame into the corresponding matching rows of the other.
#
# Arguments:
#   t.df: target data frame to which columns will be added and containing positions.
#   s.df: source data frame containing data to be added to t.df and containing positions.
#   t.pos: list giving column names in t.df that define a position.  Members:
#           "start" (required) : start position or main position.
#           "id", "len", "end" (optional) : sequence ID, length, and end position in bp.
#   s.pos: like t.pos, for s.df.
#   dist: list defining how to find position matches between t.df and s.df.  Members:
#           "method" (required) : one of "OVERLAP", "s.TINY", "t.TINY", "s.NEAR", "t.NEAR"
#           "closest" (optional) : for s/t.NEAR, 0 (all), 1 (nearest 1), or 2 (nearest 2)
#           "start.up", "start.down", "end.up", "end.down" (optional) : for x.NEAR, x = s/t:
#               start.up: x.start must be no less than y.start-start.up
#               start.down: x.start must be no more than y.end+start.down
#               end.up: x.end must be no less than y.start-end.up
#               end.down: x.end must be no more than y.end+end.down
#   cols: list defining the column names in s.df to be copied to t.df, the format of
#       the column data that is copied, and the column names in t.df to which the data
#       is copied.  Members are sublists, one per column to be added to t.df.  Members:
#           "col" : name of column to add
#           "before" : column to add it before, "" to add to end
#           "format" : constructs are: {lb}, {rb}, {+col}, {#t}, {#s}, {%t}, {%s},
#               {*s*col*val*dgts}, {*t*col*val*dgts}, {/s/col/RE/RE.replace}, {/t/col/RE/RE.replace}
#           "maxMatch" : maximum number of matches per t.df row, 0 for no limit, default 0.
#           "join" : "YES" to join all match strings for the t.df row into one column,
#               "NO" to put them in separate columns with a number appended to column name
#           "joinStart", "joinSep", "joinEnd" : strings to separate joined match strings.
#
#   *** See Details below for complete information. ***
#
# Returns: modified dft data frame.
#
# Details:
#
# This description is written assuming that the position information is nucleotide
# genomic positions, but the code doesn't care what it is as long as it is integers
# (and an optional arbitrary ID such as a chromosome name).
#
# The position information in t.df and s.df may be the position of either a CONTIG
# consisting of one or more nucleotides, or a SINGLE NUCLEOTIDE (e.g. a SNP position).
# For convenience below, we will use the term "SNP" to mean a "single nucleotide
# position", so we take "P" to mean "position", not "polymorphism".  The t.df data
# frame may use one form (CONTIG or SNP) and s.df may use the same or opposite form.
#
# The arguments "t.pos" and "s.pos" are lists that describe where to find the position
# information in t.df and s.df. They must have AT LEAST a member named "start", and may
# also have other members.  The possible list members and their meanings are:
#       id: (optional) NAME OF THE COLUMN containing a POSITION ID, typically a sequence
#           ID such as a chromosome ID, applying to that row of the data frame.  If 'id'
#           is specified, it must be specified in both t.pos and s.pos, and only rows
#           containing the same ID will match.
#       start: (required) NAME OF THE COLUMN containing a number that is the START
#           POSITION of the CONTIG or the POSITION of the SINGLE NUCLEOTIDE designated
#           by each row of the data frame.
#       len: (optional) NAME OF THE COLUMN containing the NUMBER OF NUCLEOTIDE BASE PAIRS
#           of the CONTIG.
#       end: (optional) NAME OF THE COLUMN containing a number that is the END POSITION
#           of the CONTIG.  Both "len" and "end" should not be specified, and if they
#           are, "len" takes precedence.  If neither "len" nor "end" are specified, the
#           length is taken as 1 bp for all data frame rows.
# Example: t.pos=list(id="chr", start="pos", len="length")
#
# The argument "dist" is a list that describes how to compare the position information
# in t.df and s.df to find a match.  It must have AT LEAST a member named "method", and
# may also have other members.  The possible list members and their meanings are:
#       method: one of the following values indicating the general method of matching:
#           "OVERLAP": Match only if s.df and t.df overlap by at least one nucleotide.
#               If both are SNPs, the positions must match.
#           "s.TINY": The positions in s.df are either SNPs or very small contigs that
#               are much smaller that the contigs in t.df, and matching requires that
#               the t.df contig completely encompasses the s.df SNP or contig.
#           "t.TINY": The opposite of s.TINY, swap s.df and t.df roles.
#           "s.NEAR": The positions in s.df are not large compared to those in t.df and
#               either one may be a SNP or a CONTIG, and matching requires that the two
#               are near to one another to the degree specified by the members "closest"
#               and "start.up", "start.down", "end.up", "end.down".
#           "t.NEAR": The opposite of t.NEAR, swap s.df and t.df roles.
#       closest: this is only applicable when method is s/t.NEAR, and it is optional and
#               if not specified is taken as 0.  It can have one of these values:
#           0: ALL contigs satisfying the four position limits below are taken as matches.
#           1: Like "OVERLAP" but when there is no overlap, only the NEAREST to each t.df
#               (s.NEAR) or s.df (t.NEAR) row of all non-overlapping matches that satisfy
#               the four position limits below is taken as a match.
#           2: like 1, but allows one match upstream of s/t.df and a second downstream,
#               in both cases the NEAREST one.
#       start.up, start.down, end.up, end.down: these are only applicable when 'method'
#           is s/t.NEAR.  Each is a signed integer that may be positive OR NEGATIVE OR NA.
#           Each is optional and taken as NA if unspecified.  All four function similarly,
#           being used as limits to how far the start and end of the s.df (s.NEAR) or t.df
#           (t.NEAR) position can be from the start and end of the t.df (s.NEAR) or s.df
#           (t.NEAR) position for a match to occur.  (For SNPs, the start and end positions
#           are equal).
#
#           For the following description, we use the terms s.start, s.end, t.start, and
#           t.end to designate the start and end positions of the s.df and t.df SNPs or
#           CONTIGs.  The four values function as follows for method x.NEAR (x = s or t,
#           y = the opposite, t or s):
#               start.up: x.start must be no less than y.start-start.up
#               start.down: x.start must be no more than y.end+start.down
#               end.up: x.end must be no less than y.start-end.up
#               end.down: x.end must be no more than y.end+end.down
#           If any of the values is NA, the test above uses the largest adjacent
#           x.start-to-x.start or x.end-to-x.end distance 
# Example: dist=list(method="s.NEAR", start.up=NA, start.down=1000, end.up=1000, end.down=NA, closest=2)
#
# The argument "cols" is a list that defines the column names in s.df to be copied to t.df,
# the way to format the data that is copied, and the column names in t.df to which the data
# is copied.  Each element of the "cols" list is a sublist and it creates one new column in
# t.df.  Each sublist must have these members:
#       col: name of the column to be created in t.df
#       before: name of column in t.df before which to put the new column 'col', or "" to
#           make it the last column of t.df.
#       format: character string defining the format of the data placed in the new column.
#           Each character is copied verbatim to the new column except for values
#           surrounded by {} braces, and these are defined as follows:
#               {lb} : a verbatim left brace
#               {rb} : a verbatim right brace
#               {+col} : the value in s.df column "col" of the matching row
#               {#t} : bp position of s.start in t.df contig, as: -#, +#, or @#, see below.
#               {#s} : bp position of t.start in s.df contig, as: -#, +#, or @#, see below.
#               {%t} : percent position of s.start in t.df contig, as: -#%, +#%, or @#%, see below.
#               {%s} : percent position of t.start in s.df contig, as: -#%, +#%, or @#%, see below.
#               {*s*col*val*dgts} : multiply the value in s.df column "col" by "val" and round to
#                   "dgts" digits after decimal (e.g. {*s*start*1e-6*0}).
#               {*t*col*val*dgts} : likewise for t.df
#               {/s/col/RE/RE.replace} : regular expression search/replace of s.df column "col"
#                   (e.g. {/s/chr/SL2.40ch0?/})
#               {/t/col/RE/RE.replace} : likewise for t.df
#           A start position may lie UPSTREAM, WITHIN, or DOWNSTREAM of a contig, and the prefix characters
#           "-", "@", and "+", respectively, are used in the four distance/percent specifiers above to
#           indicate which one is the case.
#       maxMatch: optional integer specifying the maximum number of matches of multiple s.df rows to a single
#           t.df row, whose data is to be copied to "col".  If 0 or if not specified, there is no limit, all
#           matches are copied.  Otherwise, if there are more matches than this integer, the remaining matches
#           are ignored.
#       join: optional character string that is "YES" to join together the strings that are generated from the
#           'format' member for each match, using the members 'joinStart', 'joinSep', and 'joinEnd' below as
#           separators between the joined strings, or is "NO" to not join the strings but instead put each one
#           in a separate column of t.df by appending numbers 1,2,...maxMatch to the column name specified by
#           'col' (unless "maxMatch" is 1).  If not specified, it defaults to "YES".
#       joinStart: optional character string to prepend to the beginning of the joined strings when 'join' is
#           "YES".  Defaults to "" (nothing is prepended).
#       joinSep: optional character string to separate each string generated from each match via the 'format'
#           member when 'join' is "YES".  Defaults to "," (comma separator).
#       joinEnd: optional character string to append to the end of the joined strings when 'join' is "YES".
#           YES.  Defaults to "" (nothing is appended).
# Example: cols=list(col="genes", before="", format="{+gene_id}({#.in.s})")
#######################################################################################
mergeOnMatches = function(t.df, s.df, t.pos, s.pos, dist, cols)
    {

    # Check arguments.  Try to coerce args of the wrong type to the right type.
    
    # Invoke stop with message given by ... or list L (if not NULL), including
    # name of this function in the message.
    error = function(..., L=NULL)
        {
        if (is.null(L) || length(L) == 0 || is.null(L[1]))
            L = list(...)
        L[["call."]] = FALSE
        L = c("mergeOnMatches: ", L)
        do.call(stop, L)
        }

    # Call coerceObj(), then make sure 'obj' is of length 1, and if not, stop
    # stop with an appropriate message including the calling name of 'obj', or
    # use ... as the stop message if ... is specified.  If 'name' is specified,
    # use that as the name of 'obj' in the stop message.
    coerceJustOne = function(obj, type, ..., allowNA=FALSE, name=NULL)
        {
        obj = coerceObj(obj, type, error, ..., allowNA=allowNA, name=name)
        if (is.null(name))
            name = deparse(substitute(obj))
        if (length(obj) != 1)
            error(name, " must be of length 1", L=list(...))
        return(obj)
        }
                    
    # Call coerceJustOne() with type=class(V), then make sure 'obj' is one of
    # the elements in V, and if not, stop with an appropriate message including
    # the calling name of 'obj', or use ... as the stop message if ... is specified.
    # If 'name' is specified, use that as the name of 'obj' in the stop message.
    coerceJustOneOf = function(obj, V, ..., name=NULL)
        {
        if (is.null(name))
            name = deparse(substitute(obj))
        obj = coerceJustOne(obj, class(V), ..., allowNA=FALSE, name=name)
        if (!obj %in% V)
            error(name, " must be one of ", paste(V, collapse=", "), L=list(...))
        return(obj)
        }
                    
    # Make sure 'obj' has a member named 'mem' and if not, stop with an
    # appropriate message including the calling name of 'obj'.
    requireMember = function(obj, mem)
        {
        if (is.null(obj[[mem]]))
            error(deparse(substitute(obj)), " must have member '", mem, "'")
        }

    # Make sure member 'mem' of 'obj' is a single character string that is a
    # member of vector 'V', where 'Vname' is a name descriptive of 'V' (for
    # error messages that say "must be a <Vname"), and if not, stop with in
    # appropriate message that includes the calling name of obj.
    requireNonNAcolName = function(obj, mem, V, Vname)
        {
        objName = deparse(substitute(obj))
        if (!is.character(obj[[mem]]) || length(obj[[mem]]) != 1)
            error(objName, "[[", mem, "]] must be a single character string")
        if (is.na(obj[[mem]]))
            error(objName, "[[", mem, "]] must be not be NA")
        if (!obj[[mem]] %in% V)
            error(objName, "[[", mem, "]] must be a ", Vname)
        }

    # Coerce vector 'V' to an integer and make sure it has no NAs.  If this fails,
    # stop with an appropriate error message that includes 'Vname' in the form
    # "<Vname> must not be NA".  Return the coerced V.
    coerceIntegerVectorNoNAs = function(V, Vname)
        {
        # Coerce to type 'int' and from there to 'integer', to catch truncation of digits after decimal point.
        V = coerceObj(V, allowNA=FALSE, "integer", error, Vname, " must be of or coercible to type integer")
        return(V)
        }

    # Coerce vector 'V', which is column 'col' of a data frame whose name is
    # 'dfName', to an integer and make sure it has no NAs.  If this fails, stop
    # with an appropriate error message that includes "Column <colName> of <dfName>...".
    # Return the coerced V.
    coerceIntegerColumnNoNAs = function(V, col, dfName)
        {
        return(coerceIntegerVectorNoNAs(V, paste("column", col, "of", dfName)))
        }

    # t.df and s.df
    t.df = coerceObj(t.df, "data.frame", error)
    s.df = coerceObj(s.df, "data.frame", error)
    if (nrow(t.df) == 0 || ncol(t.df) == 0) error("t.df is empty")
    if (nrow(s.df) == 0 || ncol(s.df) == 0) error("s.df is empty")

    # t.pos and s.pos
    t.pos = coerceObj(t.pos, "list", error)
    s.pos = coerceObj(s.pos, "list", error)
    requireMember(t.pos, "start")
    requireMember(s.pos, "start")

    mems.cols = c("start", "id", "len", "end")
    mems = intersect(names(t.pos), mems.cols)
    for (mem in mems)
        {
        requireNonNAcolName(t.pos, mem, colnames(t.df), "column name of t.df")
        col = t.pos[[mem]]
        if (mem != "id")
            t.df[,col] = coerceIntegerColumnNoNAs(t.df[,col], col, "t.df")
        }
    if (!is.null(t.pos[["len"]]) && any(t.df[,t.pos[["len"]]] <= 0))
        error("One or more in column t.df$", t.pos[["len"]], " is <= 0")
    if (!is.null(t.pos[["end"]]) && any(t.df[,t.pos[["start"]]] > t.df[,t.pos[["end"]]]))
        error("One or more in column t.df$", t.pos[["start"]], " is greater than corresponding t.df$", t.pos[["end"]])

    mems = intersect(names(s.pos), mems.cols)
    for (mem in mems)
        {
        requireNonNAcolName(s.pos, mem, colnames(s.df), "column name of s.df")
        col = s.pos[[mem]]
        if (mem != "id")
            s.df[,col] = coerceIntegerColumnNoNAs(s.df[,col], col, "s.df")
        }
    if (!is.null(s.pos[["len"]]) && any(s.df[,s.pos[["len"]]] <= 0))
        error("One or more in column s.df$", s.pos[["len"]], " is <= 0")
    if (!is.null(s.pos[["end"]]) && any(s.df[,s.pos[["start"]]] > s.df[,s.pos[["end"]]]))
        error("One or more in column s.df$", s.pos[["start"]], " is greater than corresponding s.df$", s.pos[["end"]])

    # dist
    dist = coerceObj(dist, "list", error)
    requireMember(dist, "method")
    dist[["method"]] = coerceJustOneOf(dist[["method"]], c("OVERLAP", "s.TINY", "t.TINY", "s.NEAR", "t.NEAR"))
    if (dist[["method"]] == "s.NEAR" || dist[["method"]] == "t.NEAR")
        for (S in c("closest", "start.up", "start.down", "end.up", "end.down"))
            {
            if (is.null(dist[[S]]))
                {
                if (S == "closest")
                    dist[[S]] = 0
                else
                    dist[[S]] = NA
                }
            if (S == "closest")
                dist[["closest"]] = coerceJustOneOf(dist[["closest"]], c(0, 1, 2))
            else
                dist[[S]] = coerceJustOne(dist[[S]], "integer", allowNA=TRUE, name=paste('dist[["', S, '"]]', sep=""))
            }

    # cols
    cols = coerceObj(cols, "list", error)
    if (length(cols) == 0) error("cols must have at least one element")
    for (i in 1:length(cols))
        {
        L = cols[[i]]
        if (!is.list(L)) error("Each cols element must be a sub-list")
        if (is.null(L[["maxMatch"]]))
            L[["maxMatch"]] = 0
        # This does not work.  L[["x"]] is equivalent to L[["x", exact=FALSE]]  !!!!!   Better not use $!
        #if (is.null(L[["join"]]))
        #    L[["join"]] = "YES"
        if (is.null(L[["join"]]))
            L[["join"]] = "YES"
        if (is.null(L[["joinStart"]]))
            L[["joinStart"]] = ""
        if (is.null(L[["joinSep"]]))
            L[["joinSep"]] = ","
        if (is.null(L[["joinEnd"]]))
            L[["joinEnd"]] = ""
            
        for (S in c("col", "before", "format", "maxMatch", "join", "joinStart", "joinSep", "joinEnd"))
            {
            if (is.null(L[[S]]))
                error("Each cols sublist must have a member named '", S, "'")
            req.type = ifelse(S != "maxMatch", "character", "integer")
            L[[S]] = coerceJustOne(L[[S]], req.type, allowNA=FALSE, name=paste("cols[[", i, "]][['", S, "']]", sep=""))
            if (S == "col" && L[["col"]] %in% colnames(t.df))
                error("cols[[", i, "]][['col']] must NOT be a column name of t.df")
            if (S == "before" && L[["before"]] != "" && !L[["before"]] %in% colnames(t.df))
                error("cols[[", i, "]][['before']] must be an empty string or a column name of t.df")
            if (S == "join" && !L[["join"]] %in% c("YES", "NO"))
                error("cols[[", i, "]][['join']] must be 'YES' or 'NO'")
            if (S == "maxMatch" && L[["maxMatch"]] < 0)
                error("cols[[", i, "]][['maxMatch']] must be >= 0")
            }
        cols[[i]] = L

        # Check the format code without generating result strings.
        formatData(L[["format"]], colnames(t.df), colnames(s.df))
        }

    # Now do it.  Start by finding matches.

    # First split up t.df and s.df by ID, if an ID column was specified.  Lists
    # L.t.df and L.s.df contain subsets of t.df and s.df for each ID, and the
    # list index is the ID.  If no ID column was specified, L.t.df and L.s.df
    # contain only one member each, the entire t.df and s.df.
    {
    L.t.df = list()
    L.s.df = list()
    if (is.null(t.pos[["id"]]))
        {
        L.t.df[[1]] = t.df
        L.s.df[[1]] = s.df
        }
    else
        {
        # Matches can only occur with IDs that are present in BOTH t.df and s.df.
        IDs = intersect(unique(t.df[,t.pos[["id"]]]), unique(s.df[,s.pos[["id"]]]))
        if (length(IDs) == 0)
            error("There are no IDs in common between t.df[,'", t.pos[["id"]], "'] and s.df[,'", s.pos[["id"]], "']")
        for (ID in IDs)
            {
            L.t.df[[ID]] = t.df[t.df[,t.pos[["id"]]] == ID,]
            L.s.df[[ID]] = s.df[s.df[,s.pos[["id"]]] == ID,]
            }
        }
    }

    # Process each member of L.t.df and L.s.df, adding columns to the L.t.df data frames.
    # Note: it is ok to reuse t.df and s.df here, we don't need the data frames any more
    # because they are in L.t.df and L.s.df.
    anyMatchesFound = FALSE
    for (i in 1:length(L.t.df))
        {
        t.df = L.t.df[[i]]
        s.df = L.s.df[[i]]

        # Make data frames t.position and s.position with columns "start" and "end".
        {
        t.position = data.frame(start=t.df[,t.pos[["start"]]], end=NA)
        if (!is.null(t.pos[["len"]]))
            t.position[["end"]] = t.position[["start"]] + t.df[,t.pos[["len"]]] - 1
        else if (!is.null(t.pos[["end"]]))
            t.position[["end"]] = t.df[,t.pos[["end"]]]
        else
            t.position[["end"]] = t.position[["start"]]

        s.position = data.frame(start=s.df[,s.pos[["start"]]], end=NA)
        if (!is.null(s.pos[["len"]]))
            s.position[["end"]] = s.position[["start"]] + s.df[,s.pos[["len"]]] - 1
        else if (!is.null(s.pos[["end"]]))
            s.position[["end"]] = s.df[,s.pos[["end"]]]
        else
            s.position[["end"]] = s.position[["start"]]
        }

        # Search for matches based on dist[["method"]].  Matrix idxs holds indexes
        # of t.df (idxs column 1) rows and s.df (idxs column 2) rows that match.
        idxs = getMatchIdxs(t.position, s.position, dist)

        # With the matching indexes at hand, we can now add columns to t.df under
        # control of 'cols' (if the idxs data frame is not empty).
        if (nrow(idxs) > 0)
            {
            anyMatchesFound = TRUE
            for (j in 1:length(cols))
                {
                L = cols[[j]]
                V = formatData(L[["format"]], colnames(t.df), colnames(s.df), t.df, s.df,
                    idxs, t.position, s.position)

                # Collect together the formatted strings for individual rows of t.df
                # as specified by idxs[,1].  Discard any beyond L[["maxMatch"]].  Join the
                # strings if L[["join"]] is "YES".
                colStrs = tapply(V, idxs[,1], function(S)
                    {
                    N = length(S)
                    if (L[["maxMatch"]] > 0 && N > L[["maxMatch"]])
                        S = S[1:L[["maxMatch"]]]
                    if (L[["join"]] == "YES")
                        S = paste(L[["joinStart"]], paste(S, collapse=L[["joinSep"]]), L[["joinEnd"]], sep="")
                    return(S)
                    }, simplify=FALSE)
                colStrsIdxs = as.integer(names(colStrs))

                # The elements of colStrs generally have varying numbers of strings
                # in them, in which colStrs will be a list.  If L[["join"]] is "YES", it
                # will be a 1-dimensional array of character strings.  If it is a list,
                # add empty strings to elements so that all elements are the same length.
                # This will turn colStrs into a matrix, but we need to transpose it to
                # get the right configuration.
                numNewCols = 1
                if (is.list(colStrs))
                    {
                    numNewCols = max(sapply(colStrs, length)) # Can't be more than L[["maxMatch"]].
                    colStrs = sapply(colStrs, function(V) c(V, rep("", numNewCols-length(V))))
                    if (is.matrix(colStrs))
                        colStrs = t(colStrs)
                    }

                # Create a data frame with the new column(s).  If L[["join"]] is "YES" or if
                # there is just one column (L[["maxMatch"]] = 1), the data frame will have
                # one column, else it will (most likely) have L[["maxMatch"]] columns.  The
                # actual number of columns is in numNewCols, computed above.
                newColNames = L[["col"]]
                if (numNewCols > 1)
                    newColNames = paste(newColNames, 1:numNewCols, sep="")
                newCols = matrix("", ncol=numNewCols, nrow=nrow(t.df),
                    dimnames=list(NULL, newColNames))
                rownames(newCols) = rownames(t.df)
                newCols[colStrsIdxs,] = colStrs
                newCols = as.data.frame(newCols, stringsAsFactors=FALSE)

                # Insert the new columns into t.df at the appropriate position and
                # we're done.
                if (L[["before"]] == "")
                    t.df = cbind(t.df, newCols)
                else if (L[["before"]] == colnames(t.df)[1])
                    t.df = cbind(newCols, t.df)
                else
                    {
                    k = which(L[["before"]] == colnames(t.df))[1]
                    t.df = cbind(t.df[,1:(k-1),drop=FALSE], newCols, t.df[, k:ncol(t.df),drop=FALSE])
                    }
                }
            }

        # Save modified t.df in L.t.df.
        L.t.df[[i]] = t.df
        }

    # There had better have been at least some matches found!
    if (!anyMatchesFound)
        error("No matches found")

    # Join the separate t.df data frames back together and return that data frame.
    # This is not as trivial as a do.call(rbind, L.t.df) because different members
    # of L.t.df might have different numbers of columns because some might have not
    # had ANY matching indexes, and because maxMatches may be non-zero and some
    # members might have had fewer than maxMatches matches.
    t.df = NULL
    colnamesMismatch = function() error("Program error: mismatch of column names: ",
                    paste(colnames(t.df), collapse=","), paste(colnames(df2), collapse=","))

    for (df2 in L.t.df)
        {
        if (is.null(t.df))
            t.df = df2

        # It should not be possible for both t.df and df2 to have columns the other
        # does not have.  The same column names are added in the same order whenever
        # columns are added.
        else if (ncol(t.df) == ncol(df2))
            {
            if (!all(colnames(t.df) == colnames(df2)))
                colnamesMismatch()
            t.df = rbind(t.df, df2)
            }

        # Add empty columns ("" strings) to make both t.df and df2 have the same columns.
        else if (ncol(t.df) > ncol(df2))
            {
            if (!all(colnames(df2) %in% colnames(t.df)))
                colnamesMismatch()
            newCols = setdiff(colnames(t.df), colnames(df2))
            df2[,newCols] = ""
            df2 = df2[, colnames(t.df), drop=FALSE]
            t.df = rbind(t.df, df2)
            }
        else # (ncol(df2) > ncol(t.df))
            {
            if (!all(colnames(t.df) %in% colnames(df2)))
                colnamesMismatch()
            newCols = setdiff(colnames(df2), colnames(t.df))
            t.df[,newCols] = ""
            t.df = t.df[, colnames(df2), drop=FALSE]
            t.df = rbind(t.df, df2)
            }
        }

    return(t.df)
    }

#######################################################################################
# End of file.
#######################################################################################
