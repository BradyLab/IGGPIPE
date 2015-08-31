#######################################################################################
# This file contains R definitions and functions for merging data frame data using
# positional information in the data frames.
# Author: Ted Toal
# Date: 2015
# Brady Lab, UC Davis
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
# Helper function for mergeOnMatches() below.  Find matches of positions in T.position
# and S.position according to the parameters in "match".
#
# Arguments:
#   T.position: data frame with columns start and end that define starting and ending
#       positions.
#   S.position: data frame with columns like T.position.
#   match: "match" argument of mergeOnMatches().
#
# Returns: a matrix with 2 columns.  The first column is indexes of T.position rows
# and the second column is indexes of S.position rows.  Each row of the matrix
# identifies a pair of rows (in T.position and S.position) that match.
#######################################################################################
getMatchIdxs = function(T.position, S.position, match)
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

    if (nrow(T.position) == 0) stop("getMatchIdxs: software error, T.position is empty")
    if (nrow(S.position) == 0) stop("getMatchIdxs: software error, S.position is empty")

    # Method depends on match[["method"]].

    # "OVERLAP": Match only if S.df and T.df overlap by at least one nucleotide.
    #   If both are SNPs, the positions must match.
    #   0: match only if S.df and T.df overlap by at least one nucleotide.  If both
    #       are SNPs, the positions must match.
    if (match[["method"]] == "OVERLAP")
        {
        # To test for overlap, note that when an overlap exists, one of the four endpoints
        # (start and end for T.position, start and end for S.position) must lie between
        # the two endpoints of the other position.  Test all four situations.
        S.position[["pos"]] = S.position[["start"]]
        idxs = findContainsIdxs.rows(T.position, S.position)

        S.position[["pos"]] = S.position[["end"]]
        idxs2 = findContainsIdxs.rows(T.position, S.position)
        idxs = getUniqueIdxs(idxs, idxs2)

        T.position[["pos"]] = T.position[["start"]]
        idxs2 = findContainsIdxs.rows(S.position, T.position)
        # We must swap the two idxs2 columns so that column 1 is for T.df.
        idxs = getUniqueIdxs(idxs, idxs2[,2:1])

        T.position[["pos"]] = T.position[["end"]]
        idxs2 = findContainsIdxs.rows(S.position, T.position)
        # We must swap the two idxs2 columns so that column 1 is for T.df.
        idxs = getUniqueIdxs(idxs, idxs2[,2:1])
        }

    # "S.TINY": The positions in S.df are either SNPs or very small contigs that
    #   are much smaller that the contigs in T.df, and matching requires that
    #   the T.df contig completely encompasses the S.df SNP or contig.
    else if (match[["method"]] == "S.TINY")
        {
        S.position[["pos"]] = S.position[["start"]]
        idxs = findContainsIdxs.rows(T.position, S.position)
        if (any(S.position[["end"]] != S.position[["start"]]))
            {
            S.position[["pos"]] = S.position[["end"]]
            idxs2 = findContainsIdxs.rows(T.position, S.position)
            idxs = getCommonIdxs(idxs, idxs2)
            }
        }

    # "T.TINY": The opposite of S.TINY, swap S.df and T.df roles.
    else if (match[["method"]] == "T.TINY")
        {
        T.position[["pos"]] = T.position[["start"]]
        idxs = findContainsIdxs.rows(S.position, T.position)
        if (any(T.position[["end"]] != T.position[["start"]]))
            {
            T.position[["pos"]] = T.position[["end"]]
            idxs2 = findContainsIdxs.rows(S.position, T.position)
            idxs = getCommonIdxs(idxs, idxs2)
            }
        # We must swap the two idxs columns so that column 1 is for T.df.
        idxs = idxs[,2:1]
        }

    # "x.NEAR", x = S/T, y = T/S: The positions in x.df are not large compared to
    #   those in y.df, either one may be a SNP or a CONTIG, and matching requires
    #   that the two are near to one another to the degree specified by the members
    #   "closest" and "start.up", "start.down", "end.up", "end.down".

    else # (match[["method"]] == "S.NEAR" || match[["method"]] == "T.NEAR")
        {
        # For T.NEAR we will swap S.position and T.position, and at the end, swap
        # the two columns of idxs.
        if (match[["method"]] == "T.NEAR")
            {
            tmp = T.position
            T.position = S.position
            S.position = tmp
            }

        # Start by getting all indexes that satisfy start.up/start.down/end.up/end.down.

        # Here, we are using x=s, i.e. S.NEAR:
        # start.up: S.start must be no less than T.start-start.up
        # start.down: S.start must be no more than T.end+start.down
        # end.up: S.end must be no less than T.start-end.up
        # end.down: S.end must be no more than T.end+end.down

        dontCareUpDist = max(abs(diff(S.position[order(S.position[["start"]]), "start"])))
        dontCareDownDist = max(abs(diff(S.position[order(S.position[["end"]]), "end"])))
        #cat("dontCareUpDist=", dontCareUpDist, "\n")
        #cat("dontCareDownDist=", dontCareDownDist, "\n")

        S.position[["pos"]] = S.position[["start"]]
        tmp.position = T.position
        maxDist = match[["start.up"]]
        if (is.na(maxDist))
            maxDist = dontCareUpDist
        tmp.position[["start"]] = tmp.position[["start"]] - maxDist
        maxDist = match[["start.down"]]
        if (is.na(maxDist))
            maxDist = dontCareDownDist
        tmp.position[["end"]] = tmp.position[["end"]] + maxDist
        idxs = findContainsIdxs.rows(tmp.position, S.position)

        S.position[["pos"]] = S.position[["end"]]
        tmp.position = T.position
        maxDist = match[["end.up"]]
        if (is.na(maxDist))
            maxDist = dontCareUpDist
        tmp.position[["start"]] = tmp.position[["start"]] - maxDist
        maxDist = match[["end.down"]]
        if (is.na(maxDist))
            maxDist = dontCareDownDist
        tmp.position[["end"]] = tmp.position[["end"]] + maxDist
        idxs2 = findContainsIdxs.rows(tmp.position, S.position)
        idxs = getCommonIdxs(idxs, idxs2)

        # Now test match[["closest"]] to decide how to further break down the method.

        # 0: ALL contigs satisfying the four position limits are taken as matches.
        #   else (match[["closest"]] == 0).  This result is already computed in "idxs".
        # 1: Like "OVERLAP" but when there is no overlap, only the NEAREST to each T.df
        #   row of all non-overlapping matches that satisfy the four position limits
        #   is taken as a match.
        # 2: like 1, but allows one match upstream of T.df and a second downstream,
        #   in both cases the NEAREST one.

        # We already have the indexes in idxs for closest = 0.  For closest = 1 or 2,
        # processing is almost identical.
        if (match[["closest"]] == 1 || match[["closest"]] == 2)
            {
            # We must exclude from idxs those that are non-overlapping and are not
            # the nearest of the non-overlapping ones.

            # Compute distance upstream and downstream (of T.df) for each idxs row.
            dist.S.upstreamOf.T = T.position[idxs[,1], "start"] - S.position[idxs[,2], "end"]
            dist.S.downstreamOf.T = S.position[idxs[,2], "start"] - T.position[idxs[,1], "end"]

            # Determine which idxs rows have overlaps.
            overlaps = (dist.S.upstreamOf.T <= 0) & (dist.S.downstreamOf.T <= 0)

            # The remaining rows are cases of near but not overlapping.  We need the set
            # of non-overlapping idxs rows.  However, we want to exclude rows whose [,1]
            # index is present in idxs[overlaps,1], because when there IS an overlap, we
            # do NOT include ANY nearest one.
            no.overlaps = !overlaps & !(idxs[,1] %in% idxs[overlaps,1])

            # Split idxs into two sets, overlapping, and strictly non-overlapping.
            idxs.overlap = idxs[overlaps,]
            idxs.no.overlap = idxs[no.overlaps,]

            # Get upstream and downstream distances for the non-overlapping set.
            dist.S.upstreamOf.T = dist.S.upstreamOf.T[no.overlaps]
            dist.S.downstreamOf.T = dist.S.downstreamOf.T[no.overlaps]

            # If we have any non-overlapping indexes, we must handle them according to "closest".
            if (nrow(idxs.no.overlap) > 0)
                {
                # Handle the two "closest" cases separately from here.
                if (match[["closest"]] == 1)
                    {
                    # We need the distance away, which for each index row will be either
                    # the dist.S.upstreamOf.T or dist.S.downstreamOf.T value, whichever one is
                    # not negative.
                    dist.T = pmax(dist.S.upstreamOf.T, dist.S.downstreamOf.T)
                    L = tapply(1:nrow(idxs.no.overlap), idxs.no.overlap[,1], function(ii)
                        {
                        # Get value of ii (idxs.no.overlap row number) of the row that has the
                        # smallest value in dist.T.
                        i = ii[which.min(dist.T[ii])]
                        return(idxs.no.overlap[i,])
                        })
                    }
                else # (match[["closest"]] == 2)
                    {
                    # Split the data into upstream and downstream data.
                    upstream = (dist.S.upstreamOf.T > 0)
                    idxs.upstream = idxs.no.overlap[upstream,]
                    idxs.downstream = idxs.no.overlap[!upstream,]
                    dist.S.upstreamOf.T = dist.S.upstreamOf.T[upstream]
                    dist.S.downstreamOf.T = dist.S.downstreamOf.T[!upstream]

                    # Get separate lists of nearest upstream and nearest downstream index pairs.
                    L.upstream = tapply(1:nrow(idxs.upstream), idxs.upstream[,1], function(ii)
                        {
                        i = ii[which.min(dist.S.upstreamOf.T[ii])]
                        return(idxs.upstream[i,])
                        })

                    L.downstream = tapply(1:nrow(idxs.downstream), idxs.downstream[,1], function(ii)
                        {
                        i = ii[which.min(dist.S.downstreamOf.T[ii])]
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

        # For T.NEAR we now must swap the two columns of idxs.
        if (match[["method"]] == "T.NEAR")
            idxs = idxs[,2:1]
        }

    # Sort idxs by column 1 then column 2.
    idxs = idxs[order(idxs[,1], idxs[,2]),]
    return(idxs)
    }

#######################################################################################
# Helper function for formatData() below.  This generates strings for the format codes
# "#T", "#S", "%T", "%S".
#
# Arguments:
#   code: one of "#T", "#S", "%T", "%S".
#   T.position: see formatData() T.position argument.
#   S.position: see formatData() S.position argument.
#   idxs: a 2-column matrix of row indexes into T.df (column 1) and S.df (column 2) of
#       the matching pairs.
#
# Returns:
#   A vector of data formatted according to "code", and of length equal to nrow(idxs),
#   to become part of a new T.df column.
#######################################################################################
positionString = function(code, T.position, S.position, idxs)
    {
    #       {#T} : bp position of S.start in T.df contig, as: -#, +#, or @#, see below.
    #       {#S} : bp position of T.start in S.df contig, as: -#, +#, or @#, see below.
    #       {%T} : percent position of S.start in T.df contig, as: -#%, +#%, or @#%, see below.
    #       {%S} : percent position of T.start in S.df contig, as: -#%, +#%, or @#%, see below.
    #   A start position may lie UPSTREAM, WITHIN, or DOWNSTREAM of a contig, and the prefix
    #   characters "-", "@", and "+", respectively, are used in the four distance/percent
    #   specifiers above to indicate which one is the case.

    # Make it so it is always like {#T} or {%T}.  Must swap idxs columns also.
    if (code == "#S" || code == "%S")
        {
        t = T.position
        T.position = S.position
        S.position = t
        rm(t)
        idxs = idxs[,2:1]
        }

    # Figure out {#T}.
    len = (T.position[["end"]] - T.position[["start"]] + 1)[idxs[,1]]
    amount = (S.position[["start"]][idxs[,2]] - T.position[["start"]][idxs[,1]])
    isUpstream = (amount < 0)
    isDownstream = (amount >= len)
    isWithin = (!isUpstream & !isDownstream)
    amount[isDownstream] = amount[isDownstream] - len[isDownstream] + 1

    # Figure out {%T}.
    if (code == "%T" || code == "%S")
        amount = as.integer(100*amount/len)

    # Make character strings.
    S = as.character(amount)
    # S[isUpstream] = paste("-", S[isUpstream], sep="") # Already has a "-" sign.
    S[isDownstream] = paste("+", S[isDownstream], sep="")
    S[isWithin] = paste("@", S[isWithin], sep="")
    if (code == "%T" || code == "%S")
        S = paste(S, "%", sep="")
    return(S)
    }

#######################################################################################
# Helper function for mergeOnMatches() below.  This parses mergeCols[[i]][["format"]]
# strings and checks their validity, generating an error if invalid.  Optionally it
# creates the new column values using the format string and match data.
#
# Arguments:
#   format: a mergeCols[[i]][["format"]] string.
#   T.colnames: the column names of the T.df data frame.
#   S.colnames: the column names of the S.df data frame.
#   T.df: NULL for parsing only, else the T.df data frame.
#   S.df: NULL for parsing only, else the S.df data frame.
#   idxs: NULL for parsing only, else a 2-column matrix of row indexes into T.df
#       (column 1) and S.df (column 2) of the matching pairs.
#   T.position: NULL for parsing only, else data frame with same number of rows as T.df,
#       and columns "start" and "end" giving positions of T.df row contigs or SNPs.
#   S.position: NULL for parsing only, else data frame with same number of rows as S.df,
#       and columns "start" and "end" giving positions of S.df row contigs or SNPs.
#
# Returns:
#   If parsing only, NULL is returned invisibly.  Otherwise, a vector of
#   data formatted according to "format", and of length equal to nrow(idxs).
#######################################################################################
formatData = function(format, T.colnames, S.colnames, T.df=NULL, S.df=NULL, idxs=NULL,
    T.position=NULL, S.position=NULL)
    {
    error = function(...) stop("formatData: ", ..., call.=FALSE)

    chk.col = function(col, st)
        {
        if (st == "S" && !any(col == S.colnames))
            error("'", col, "' is not a column name of S.df")
        if (st == "T" && !any(col == T.colnames))
            error("'", col, "' is not a column name of T.df")
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
    generate = !is.null(T.df)
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

        # {#T}, {#S}, {%T}, {%S}
        else if (code == "#T" || code == "#S" || code == "%T" || code == "%S")
            {
            if (generate)
                fmtstrs = paste(fmtstrs, positionString(code, T.position, S.position, idxs), sep="")
            }

        # {+col} : the value in S.df column "col" of the matching row
        else if (firstChar == "+")
            {
            col = sub("^\\+(.*)$", "\\1", code)
            chk.col(col, "S")
            if (generate)
                fmtstrs = paste(fmtstrs, S.df[,col][idxs[,2]], sep="")
            }

        # {*S*col*val*dgts} and {*T*col*val*dgts}
        else if (firstChar == "*")
            {
            s = unlist(strsplit(code, "*", fixed=TRUE))
            if (length(s) != 5)
                error("expected {*S|T*col*val*dgts} but got: ", code)
            if (s[2] != "S" && s[2] != "T")
                error("expected {*S... or {*T... but got: ", code)
            col = s[3]
            chk.col(col, s[2])
            val = coerceObj(s[4], "numeric",
                error, "expected numeric val in {*S|T*col*val*dgts} but got ", s[4])
            dgts = coerceObj(s[5], "integer",
                error, "expected integer dgts in {*S|T*col*val*dgts} but got ", s[5])
            if (generate)
                {
                if (s[2] == "S")
                    fmtstrs = paste(fmtstrs, round(S.df[,col][idxs[,2]]*val, dgts), sep="")
                else
                    fmtstrs = paste(fmtstrs, round(T.df[,col][idxs[,1]]*val, dgts), sep="")
                }
            }

        # {/S/col/RE/RE.replace} and {/T/col/RE/RE.replace} :
        # regular expression search/replace of S.df column "col"
        else if (firstChar == "/")
            {
            s = unlist(strsplit(code, "/", fixed=TRUE))
            if (length(s) == 4)
                s[5] = ""
            if (length(s) != 5)
                error("expected {/S|T/col/RE/RE.replace} but got: ", code)
            if (s[2] != "S" && s[2] != "T")
                error("expected {*S... or {*T... but got: ", code)
            col = s[3]
            chk.col(col, s[2])
            RE = s[4]
            RE.replace = s[5]
            if (generate)
                {
                if (s[2] == "S")
                    S = sub.catch(RE, RE.replace, S.df[,col][idxs[,2]])
                else
                    S = sub.catch(RE, RE.replace, T.df[,col][idxs[,1]])
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
#   T.df: target data frame to which columns will be added and containing positions.
#   S.df: source data frame containing data to be added to T.df and containing positions.
#   T.pos: list giving column names in T.df that define a position.  Members:
#           "start" (required) : start position or main position.
#           "id", "len", "end" (optional) : sequence ID, length, and end position in bp.
#   S.pos: like T.pos, for S.df.
#   match: list defining how to find position matches between T.df and S.df.  Members:
#           "method" (required) : one of "OVERLAP", "S.TINY", "T.TINY", "S.NEAR", "T.NEAR"
#           "closest" (optional) : for s/T.NEAR, 0 (all), 1 (nearest 1), or 2 (nearest 2)
#           "start.up", "start.down", "end.up", "end.down" (optional) : for x.NEAR, x = S/T:
#               start.up: x.start must be no less than y.start-start.up
#               start.down: x.start must be no more than y.end+start.down
#               end.up: x.end must be no less than y.start-end.up
#               end.down: x.end must be no more than y.end+end.down
#   mergeCols: list defining the column names in S.df to be copied to T.df, the format of
#       the column data that is copied, and the column names in T.df to which the data
#       is copied.  Members are sublists, one per column to be added to T.df.  Members:
#           "col" : name of column to add
#           "before" : column to add it before, "" to add to end
#           "format" : constructs are: {lb}, {rb}, {+col}, {#T}, {#S}, {%T}, {%S},
#               {*S*col*val*dgts}, {*T*col*val*dgts}, {/S/col/RE/RE.replace}, {/T/col/RE/RE.replace}
#           "maxMatch" : maximum number of matches per T.df row, 0 for no limit, default 0.
#           "join" : "TRUE" to join all match strings for the T.df row into one column,
#               "FALSE" to put them in separate columns with a number appended to column name
#           "joinPfx", "joinSep", "joinSfx" : strings to separate joined match strings.
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
# The position information in T.df and S.df may be the position of either a CONTIG
# consisting of one or more nucleotides, or a SINGLE NUCLEOTIDE (e.g. a SNP position).
# For convenience below, we will use the term "SNP" to mean a "single nucleotide
# position", so we take "P" to mean "position", not "polymorphism".  The T.df data
# frame may use one form (CONTIG or SNP) and S.df may use the same or opposite form.
#
# The arguments "T.pos" and "S.pos" are lists that describe where to find the position
# information in T.df and S.df. They must have AT LEAST a member named "start", and may
# also have other members.  The possible list members and their meanings are:
#       id: (optional) NAME OF THE COLUMN containing a POSITION ID, typically a sequence
#           ID such as a chromosome ID, applying to that row of the data frame.  If 'id'
#           is specified, it must be specified in both T.pos and S.pos, and only rows
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
# Example: T.pos=list(id="chr", start="pos", len="length")
#
# The argument "match" is a list that describes how to compare the position information
# in T.df and S.df to find a match.  It must have AT LEAST a member named "method", and
# may also have other members.  The possible list members and their meanings are:
#       method: one of the following values indicating the general method of matching:
#           "OVERLAP": Match only if S.df and T.df overlap by at least one nucleotide.
#               If both are SNPs, the positions must match.
#           "S.TINY": The positions in S.df are either SNPs or very small contigs that
#               are much smaller that the contigs in T.df, and matching requires that
#               the T.df contig completely encompasses the S.df SNP or contig.
#           "T.TINY": The opposite of S.TINY, swap S.df and T.df roles.
#           "S.NEAR": The positions in S.df are not large compared to those in T.df and
#               either one may be a SNP or a CONTIG, and matching requires that the two
#               are near to one another to the degree specified by the members "closest"
#               and "start.up", "start.down", "end.up", "end.down".
#           "T.NEAR": The opposite of T.NEAR, swap S.df and T.df roles.
#       closest: this is only applicable when method is s/T.NEAR, and it is optional and
#               if not specified is taken as 0.  It can have one of these values:
#           0: ALL contigs satisfying the four position limits below are taken as matches.
#           1: Like "OVERLAP" but when there is no overlap, only the NEAREST to each T.df
#               (S.NEAR) or S.df (T.NEAR) row of all non-overlapping matches that satisfy
#               the four position limits below is taken as a match.
#           2: like 1, but allows one match upstream of s/T.df and a second downstream,
#               in both cases the NEAREST one.
#       start.up, start.down, end.up, end.down: these are only applicable when 'method'
#           is s/T.NEAR.  Each is a signed integer that may be positive OR NEGATIVE OR NA.
#           Each is optional and taken as NA if unspecified.  All four function similarly,
#           being used as limits to how far the start and end of the S.df (S.NEAR) or T.df
#           (T.NEAR) position can be from the start and end of the T.df (S.NEAR) or S.df
#           (T.NEAR) position for a match to occur.  (For SNPs, the start and end positions
#           are equal).
#
#           For the following description, we use the terms S.start, S.end, T.start, and
#           T.end to designate the start and end positions of the S.df and T.df SNPs or
#           CONTIGs.  The four values function as follows for method x.NEAR (x = s or t,
#           y = the opposite, t or s):
#               start.up: x.start must be no less than y.start-start.up
#               start.down: x.start must be no more than y.end+start.down
#               end.up: x.end must be no less than y.start-end.up
#               end.down: x.end must be no more than y.end+end.down
#           If any of the values is NA, the test above uses the largest adjacent
#           x.start-to-x.start or x.end-to-x.end distance 
# Example: match=list(method="S.NEAR", start.up=NA, start.down=1000, end.up=1000, end.down=NA, closest=2)
#
# The argument "mergeCols" is a list that defines the column names in S.df to be copied to T.df,
# the way to format the data that is copied, and the column names in T.df to which the data
# is copied.  Each element of the "mergeCols" list is a sublist and it creates one new column in
# T.df.  Each sublist must have these members:
#       col: name of the column to be created in T.df
#       before: name of column in T.df before which to put the new column 'col'.  If "" or
#           not specified, the new column is placed as the last column of T.df.
#       format: character string defining the format of the data placed in the new column.
#           Each character is copied verbatim to the new column except for values
#           surrounded by {} braces, and these are defined as follows:
#               {lb} : a verbatim left brace
#               {rb} : a verbatim right brace
#               {+col} : the value in S.df column "col" of the matching row
#               {#T} : bp position of S.start in T.df contig, as: -#, +#, or @#, see below.
#               {#S} : bp position of T.start in S.df contig, as: -#, +#, or @#, see below.
#               {%T} : percent position of S.start in T.df contig, as: -#%, +#%, or @#%, see below.
#               {%S} : percent position of T.start in S.df contig, as: -#%, +#%, or @#%, see below.
#               {*S*col*val*dgts} : multiply the value in S.df column "col" by "val" and round to
#                   "dgts" digits after decimal (e.g. {*S*start*1e-6*0}).
#               {*T*col*val*dgts} : likewise for T.df
#               {/S/col/RE/RE.replace} : regular expression search/replace of S.df column "col"
#                   (e.g. {/S/chr/SL2.40ch0?/})
#               {/T/col/RE/RE.replace} : likewise for T.df
#           A start position may lie UPSTREAM, WITHIN, or DOWNSTREAM of a contig, and the prefix characters
#           "-", "@", and "+", respectively, are used in the four distance/percent specifiers above to
#           indicate which one is the case.
#       maxMatch: optional integer specifying the maximum number of matches of multiple S.df rows to a single
#           T.df row, whose data is to be copied to "col".  If 0 or if not specified, there is no limit, all
#           matches are copied.  Otherwise, if there are more matches than this integer, the remaining matches
#           are ignored.
#       join: optional character string that is "TRUE" to join together the strings that are generated from the
#           'format' member for each match, using the members 'joinPfx', 'joinSep', and 'joinSfx' below as
#           separators between the joined strings, or is "FALSE" to not join the strings but instead put each one
#           in a separate column of T.df by appending numbers 1,2,...maxMatch to the column name specified by
#           'col' (unless "maxMatch" is 1).  If not specified, it defaults to "TRUE".
#       joinPfx: optional character string to prepend to the beginning of the joined strings when 'join' is
#           "TRUE".  Defaults to "" (nothing is prepended).
#       joinSep: optional character string to separate each string generated from each match via the 'format'
#           member when 'join' is "TRUE".  Defaults to "," (comma separator).
#       joinSfx: optional character string to append to the end of the joined strings when 'join' is "TRUE".
#           Defaults to "" (nothing is appended).
# Example: mergeCols=list(col="genes", before="", format="{+gene_id}({#S})")
#######################################################################################
mergeOnMatches = function(T.df, S.df, T.pos, S.pos, match, mergeCols)
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

    # T.df and S.df
    T.df = coerceObj(T.df, "data.frame", error)
    S.df = coerceObj(S.df, "data.frame", error)
    if (nrow(T.df) == 0 || ncol(T.df) == 0) error("T.df is empty")
    if (nrow(S.df) == 0 || ncol(S.df) == 0) error("S.df is empty")

    # T.pos and S.pos
    T.pos = coerceObj(T.pos, "list", error)
    S.pos = coerceObj(S.pos, "list", error)
    requireMember(T.pos, "start")
    requireMember(S.pos, "start")

    mems.posCols = c("start", "id", "len", "end")
    mems = intersect(names(T.pos), mems.posCols)
    for (mem in mems)
        {
        requireNonNAcolName(T.pos, mem, colnames(T.df), "column name of T.df")
        col = T.pos[[mem]]
        if (mem != "id")
            T.df[,col] = coerceIntegerColumnNoNAs(T.df[,col], col, "T.df")
        }
    if (!is.null(T.pos[["len"]]) && any(T.df[,T.pos[["len"]]] <= 0))
        error("One or more in column T.df$", T.pos[["len"]], " is <= 0")
    if (!is.null(T.pos[["end"]]) && any(T.df[,T.pos[["start"]]] > T.df[,T.pos[["end"]]]))
        error("One or more in column T.df$", T.pos[["start"]], " is greater than corresponding T.df$", T.pos[["end"]])

    mems = intersect(names(S.pos), mems.posCols)
    for (mem in mems)
        {
        requireNonNAcolName(S.pos, mem, colnames(S.df), "column name of S.df")
        col = S.pos[[mem]]
        if (mem != "id")
            S.df[,col] = coerceIntegerColumnNoNAs(S.df[,col], col, "S.df")
        }
    if (!is.null(S.pos[["len"]]) && any(S.df[,S.pos[["len"]]] <= 0))
        error("One or more in column S.df$", S.pos[["len"]], " is <= 0")
    if (!is.null(S.pos[["end"]]) && any(S.df[,S.pos[["start"]]] > S.df[,S.pos[["end"]]]))
        error("One or more in column S.df$", S.pos[["start"]], " is greater than corresponding S.df$", S.pos[["end"]])

    # match
    match = coerceObj(match, "list", error)
    requireMember(match, "method")
    match[["method"]] = coerceJustOneOf(match[["method"]], c("OVERLAP", "S.TINY", "T.TINY", "S.NEAR", "T.NEAR"))
    if (match[["method"]] == "S.NEAR" || match[["method"]] == "T.NEAR")
        for (S in c("closest", "start.up", "start.down", "end.up", "end.down"))
            {
            if (is.null(match[[S]]))
                {
                if (S == "closest")
                    match[[S]] = 0
                else
                    match[[S]] = NA
                }
            if (S == "closest")
                match[["closest"]] = coerceJustOneOf(match[["closest"]], c(0, 1, 2))
            else
                match[[S]] = coerceJustOne(match[[S]], "integer", allowNA=TRUE, name=paste('match[["', S, '"]]', sep=""))
            }

    # mergeCols
    mergeCols = coerceObj(mergeCols, "list", error)
    if (length(mergeCols) == 0) error("mergeCols must have at least one element")
    for (i in 1:length(mergeCols))
        {
        L = mergeCols[[i]]
        if (!is.list(L)) error("Each mergeCols element must be a sub-list")
        if (is.null(L[["before"]]))
            L[["before"]] = ""
        if (is.null(L[["maxMatch"]]))
            L[["maxMatch"]] = 0
        # This does not work.  L[["x"]] is equivalent to L[["x", exact=FALSE]]  !!!!!   Better not use $!
        #if (is.null(L[["join"]]))
        #    L[["join"]] = TRUE
        if (is.null(L[["join"]]))
            L[["join"]] = TRUE
        if (is.null(L[["joinPfx"]]))
            L[["joinPfx"]] = ""
        if (is.null(L[["joinSep"]]))
            L[["joinSep"]] = ","
        if (is.null(L[["joinSfx"]]))
            L[["joinSfx"]] = ""
            
        for (S in c("col", "before", "format", "maxMatch", "join", "joinPfx", "joinSep", "joinSfx"))
            {
            if (is.null(L[[S]]))
                error("Each mergeCols sublist must have a member named '", S, "'")
            req.type = ifelse(S != "maxMatch", "character", "integer")
            if (S == "join")
                req.type = "logical"
            L[[S]] = coerceJustOne(L[[S]], req.type, allowNA=FALSE,
                name=paste("mergeCols[[", i, "]][['", S, "']]", sep=""))
            if (S == "col" && L[["col"]] %in% colnames(T.df))
                error("mergeCols[[", i, "]][['col']] must NOT be a column name of T.df")
            if (S == "before" && L[["before"]] != "" && !L[["before"]] %in% colnames(T.df))
                error("mergeCols[[", i, "]][['before']] must be an empty string or a column name of T.df")
            if (S == "join" && !is.logical(L[["join"]]))
                error("mergeCols[[", i, "]][['join']] must be TRUE of FALSE")
            if (S == "maxMatch" && L[["maxMatch"]] < 0)
                error("mergeCols[[", i, "]][['maxMatch']] must be >= 0")
            }
        mergeCols[[i]] = L

        # Check the format code without generating result strings.
        formatData(L[["format"]], colnames(T.df), colnames(S.df))
        }

    # Now do it.  Start by finding matches.

    # First split up T.df and S.df by ID, if an ID column was specified.  Lists
    # L.T.df and L.S.df contain subsets of T.df and S.df for each ID, and the
    # list index is the ID.  If no ID column was specified, L.T.df and L.S.df
    # contain only one member each, the entire T.df and S.df.
    {
    L.T.df = list()
    L.S.df = list()
    if (is.null(T.pos[["id"]]))
        {
        L.T.df[[1]] = T.df
        L.S.df[[1]] = S.df
        }
    else
        {
        # Matches can only occur with IDs that are present in BOTH T.df and S.df.
        IDs = intersect(unique(T.df[,T.pos[["id"]]]), unique(S.df[,S.pos[["id"]]]))
        if (length(IDs) == 0)
            error("There are no IDs in common between T.df[,'", T.pos[["id"]], "'] and S.df[,'", S.pos[["id"]], "']")
        for (ID in IDs)
            {
            L.T.df[[ID]] = T.df[T.df[,T.pos[["id"]]] == ID,]
            L.S.df[[ID]] = S.df[S.df[,S.pos[["id"]]] == ID,]
            }
        }
    }

    # Process each member of L.T.df and L.S.df, adding columns to the L.T.df data frames.
    # Note: it is ok to reuse T.df and S.df here, we don't need the data frames any more
    # because they are in L.T.df and L.S.df.
    anyMatchesFound = FALSE
    for (i in 1:length(L.T.df))
        {
        T.df = L.T.df[[i]]
        S.df = L.S.df[[i]]

        # Make data frames T.position and S.position with columns "start" and "end".
        {
        T.position = data.frame(start=T.df[,T.pos[["start"]]], end=NA)
        if (!is.null(T.pos[["len"]]))
            T.position[["end"]] = T.position[["start"]] + T.df[,T.pos[["len"]]] - 1
        else if (!is.null(T.pos[["end"]]))
            T.position[["end"]] = T.df[,T.pos[["end"]]]
        else
            T.position[["end"]] = T.position[["start"]]

        S.position = data.frame(start=S.df[,S.pos[["start"]]], end=NA)
        if (!is.null(S.pos[["len"]]))
            S.position[["end"]] = S.position[["start"]] + S.df[,S.pos[["len"]]] - 1
        else if (!is.null(S.pos[["end"]]))
            S.position[["end"]] = S.df[,S.pos[["end"]]]
        else
            S.position[["end"]] = S.position[["start"]]
        }

        # Search for matches based on match[["method"]].  Matrix idxs holds indexes
        # of T.df (idxs column 1) rows and S.df (idxs column 2) rows that match.
        idxs = getMatchIdxs(T.position, S.position, match)

        # With the matching indexes at hand, we can now add columns to T.df under
        # control of 'mergeCols' (if the idxs data frame is not empty).
        if (nrow(idxs) > 0)
            {
            anyMatchesFound = TRUE
            for (j in 1:length(mergeCols))
                {
                L = mergeCols[[j]]
                V = formatData(L[["format"]], colnames(T.df), colnames(S.df), T.df, S.df,
                    idxs, T.position, S.position)

                # Collect together the formatted strings for individual rows of T.df
                # as specified by idxs[,1].  Discard any beyond L[["maxMatch"]].  Join the
                # strings if L[["join"]] is TRUE.
                colStrs = tapply(V, idxs[,1], function(S)
                    {
                    N = length(S)
                    if (L[["maxMatch"]] > 0 && N > L[["maxMatch"]])
                        S = S[1:L[["maxMatch"]]]
                    if (L[["join"]])
                        S = paste(L[["joinPfx"]], paste(S, collapse=L[["joinSep"]]), L[["joinSfx"]], sep="")
                    return(S)
                    }, simplify=FALSE)
                colStrsIdxs = as.integer(names(colStrs))

                # The elements of colStrs generally have varying numbers of strings
                # in them, in which colStrs will be a list.  If L[["join"]] is TRUE, it
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

                # Create a data frame with the new column(s).  If L[["join"]] is TRUE or if
                # there is just one column (L[["maxMatch"]] = 1), the data frame will have
                # one column, else it will (most likely) have L[["maxMatch"]] columns.  The
                # actual number of columns is in numNewCols, computed above.
                newColNames = L[["col"]]
                if (numNewCols > 1)
                    newColNames = paste(newColNames, 1:numNewCols, sep="")
                newCols = matrix("", ncol=numNewCols, nrow=nrow(T.df),
                    dimnames=list(NULL, newColNames))
                rownames(newCols) = rownames(T.df)
                newCols[colStrsIdxs,] = colStrs
                newCols = as.data.frame(newCols, stringsAsFactors=FALSE)

                # Insert the new columns into T.df at the appropriate position and
                # we're done.
                if (L[["before"]] == "")
                    T.df = cbind(T.df, newCols)
                else if (L[["before"]] == colnames(T.df)[1])
                    T.df = cbind(newCols, T.df)
                else
                    {
                    k = which(L[["before"]] == colnames(T.df))[1]
                    T.df = cbind(T.df[,1:(k-1),drop=FALSE], newCols, T.df[, k:ncol(T.df),drop=FALSE])
                    }
                }
            }

        # Save modified T.df in L.T.df.
        L.T.df[[i]] = T.df
        }

    # There had better have been at least some matches found!
    if (!anyMatchesFound)
        error("No matches found")

    # Join the separate T.df data frames back together and return that data frame.
    # This is not as trivial as a do.call(rbind, L.T.df) because different members
    # of L.T.df might have different numbers of columns because some might have not
    # had ANY matching indexes, and because maxMatches may be non-zero and some
    # members might have had fewer than maxMatches matches.
    T.df = NULL
    colnamesMismatch = function() error("Program error: mismatch of column names: ",
                    paste(colnames(T.df), collapse=","), paste(colnames(df2), collapse=","))

    for (df2 in L.T.df)
        {
        if (is.null(T.df))
            T.df = df2

        # It should not be possible for both T.df and df2 to have columns the other
        # does not have.  The same column names are added in the same order whenever
        # columns are added.
        else if (ncol(T.df) == ncol(df2))
            {
            if (!all(colnames(T.df) == colnames(df2)))
                colnamesMismatch()
            T.df = rbind(T.df, df2)
            }

        # Add empty columns ("" strings) to make both T.df and df2 have the same columns.
        else if (ncol(T.df) > ncol(df2))
            {
            if (!all(colnames(df2) %in% colnames(T.df)))
                colnamesMismatch()
            newCols = setdiff(colnames(T.df), colnames(df2))
            df2[,newCols] = ""
            df2 = df2[, colnames(T.df), drop=FALSE]
            T.df = rbind(T.df, df2)
            }
        else # (ncol(df2) > ncol(T.df))
            {
            if (!all(colnames(T.df) %in% colnames(df2)))
                colnamesMismatch()
            newCols = setdiff(colnames(df2), colnames(T.df))
            T.df[,newCols] = ""
            T.df = T.df[, colnames(df2), drop=FALSE]
            T.df = rbind(T.df, df2)
            }
        }

    return(T.df)
    }

#######################################################################################
# End of file.
#######################################################################################
