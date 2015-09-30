################################################################################
# This file contains R definitions and functions for removing overlapping rows
# of a data frame that contains line segment position data, leaving only rows
# that do not overlap one another.  This is to be sourced by R code that needs
# to use it.
# Author: Ted Toal
# Date: 2015
# Brady Lab, UC Davis
################################################################################

################################################################################
# Remove overlapping rows of a data frame whose rows contain start and end
# positions of line segments.
#
# Arguments:
#   df: data frame to be processed.
#   minmax: "MIN" to prefer retention of shorter segments when removing segments
#       that overlap, "MAX" to prefer retention of longer segments.
#   segName: a name by which to refer to a line segment, e.g. "indel" or "marker".
#   nameSets: a name by which to refer to an element of setNames, e.g. "genome".
#   setNames: vector of names for the set of line segments contained in one row
#       of df.  The length is equal to the number of different line segments in
#       a row, each of which, will independently of the others, be analyzed for
#       overlapping segments.  For example, these names might be genome names,
#       if each row of df contains a line segment for each of two or more genomes.
#   idCols: columns of df that contain IDs that break the line segments up into
#       different subsets, each subset of which is to be processed independently.
#       For example, these columns might be chromosome or scaffold IDs.  The
#       length of this must equal the length of setNames.
#   pos1Cols: columns of df that contain the position of one endpoint of each
#       line segment.  The length of this must equal the length of setNames.
#       This position can be either the start or end position.
#   pos2Cols: columns of df that contain the position of the OTHER endpoint of
#       each line segment.  The length of this must equal the length of setNames.
#       This position is either the start or end position, opposite of pos1Col.
#   verbose: TRUE to give progress reports.
#
# Returns: df with overlapping rows removed, and sorted by position within the
# first set of segments (i.e. setNames[1]).
################################################################################
removeOverlappingRows = function(df, segName, nameSets, setNames, idCols, pos1Cols, pos2Cols, verbose=TRUE)
    {
    if (length(idCols) != length(setNames))
        stop("removeOverlappingRows: idCols length must be same as setNames length")
    if (length(pos1Cols) != length(setNames))
        stop("removeOverlappingRows: pos1Cols length must be same as setNames length")
    if (length(pos2Cols) != length(setNames))
        stop("removeOverlappingRows: pos2Cols length must be same as setNames length")

    names(idCols) = setNames
    names(pos1Cols) = setNames
    names(pos2Cols) = setNames

    # Test each setName, one by one, for overlaps, and remove segments to get rid of them.
    #gcinfo(TRUE)
    for (setName in setNames)
        {
        if (verbose)
            catnow("Remove ", segName, " overlaps for ", nameSets, ": ", setName, "\n", sep="")

        #gc()

        id.Col = idCols[setName]
        pos1.Col = pos1Cols[setName]
        pos2.Col = pos2Cols[setName]

        # The idea behind the method is: if we make a vector of all positions,
        # whether start or end (but keeping track of which are start and which are
        # end), and sort them by position, then do a cumulative sum (i.e. cumsum())
        # of point == start point flags, and another cumsum of point == end point
        # flags, then any endpoint whose cumsums are identical MUST be an endpoint
        # preceding a break where there are no overlaps.  This is because the cumsum
        # of endpoints gives the number of endpoints at and left of a given endpoint,
        # and the number of start points left of it must be equal to that value, or
        # greater.  The only way it can be greater is if a segment starts left of it
        # and ends right of it, which means it is not a break position.
        #
        # Having identified break positions, the diff() of the cumsum values at the
        # start points gives the number of line segments between each set of breaks.
        # If that value is one, there are no overlaps to be handled between those
        # breaks, else there are overlaps.
        #
        # For each region of overlaps, find the segment that is either smallest (MIN)
        # or largest (MAX), and remove all other segments in the region that overlap
        # that segment.  Then, repeat everything again, until there are no overlaps.
        #
        # It could be a problem if a start point is equal to an end point, since when
        # we sort by position, the start point may end up before or after the end point.
        # To force start points to sort AFTER end points, we will sort using the "start"
        # column as the third sort key (after "id" and "pos").  This will result in
        # us considering two segments that have an end point of one in common with the
        # start point of the other to be non-overlapping.

        # Start by preparing a data frame dfP containing a subset of the data for this
        # setName, and with two rows for each row of df, one row containing the
        # start position and one row containing the end position.  The data frame has
        # these columns: id, pos, start, and idx.  The "pos" column has either a
        # start or end position, where start position is always the lower-numbered
        # position of the two for any given segment.  The "start" column is TRUE if
        #  "pos" is a start position, FALSE if an end position.  The "idx" column has
        # the row number in df of the data for that dfP row.  Additional working
        # columns may be added below as we make use of dfP, but these are recomputed
        # each time through the loop below.
        # Also, we will use a factor for the id, since this array could be huge and
        # we want to save space (should have used factors for lots of things but just
        # haven't gotten into it yet).
        # We build dfP by putting all the df start positions in the first N rows
        # and all the end positions in the second N rows.  As overlapping segments
        # are identified, their row index in df is added to vector "idxsToRemove", and
        # they are removed from dfP, but df remains unmodified so that the "idx"
        # column values in dfP will remain valid.  We will also remove from dfP any
        # segments that we know don't overlap anything else.  So, dfP contains the
        # start/end points and indexes of segments that MIGHT overlap and are being
        # tested for overlap.  When finished testing all segments for this setName,
        # df is modified by removing the rows given by "idxsToRemove".

        # Add columns 'start', 'end', and 'len' to df, containing the start and
        # end positions and length of each row's segment in setName 'setName'.  The
        # 'start' and 'end' values differ from the pos1.Col and pos2.Col values
        # because start < end, whereas it is possible that pos1.Col > pos2.Col.
        df$start = df[, pos1.Col]
        df$end = df[, pos2.Col]
        swapIt = (df$start > df$end)
        df$start[swapIt] = df[swapIt, pos2.Col]
        df$end[swapIt] = df[swapIt, pos1.Col]
        df$len = df$end - df$start + 1

        # Create the dfP data frame.
        N = nrow(df)
        dfP = data.frame(
            id=factor(c(as.character(df[,id.Col]), as.character(df[,id.Col])), ordered=TRUE),
            pos=c(df$start, df$end), start=c(rep(TRUE, N), rep(FALSE, N)),
            idx=c(1:N,1:N), stringsAsFactors=TRUE)

        # Create empty 'idxsToRemove' vector.
        idxsToRemove = integer()

        # Sort by position, including "start" column as final key.  Once sorted, it
        # never needs to be sorted again.
        dfP = dfP[order(dfP$id, dfP$pos, dfP$start),]
        rownames(dfP) = NULL

        # Loop until no more segments are found to overlap in this setName.
        # Use dfP as described in comments above, removing elements from dpP and
        # adding the index of elements that are to be removed from df to idxsToRemove.
        # We know we are done when dfP becomes empty.
        while (nrow(dfP) > 0)
            {
            if (verbose)
                catnow(segName, "elements remaining to be tested:", nrow(dfP), "\n")

            # Compute cumulative sum of "start" flag.  The cumulative sum of !start
            # (i.e. of endpoint flags) is simply 1:nrow(dfP) - cumsumStart.
            cumsumStart = cumsum(dfP$start)

            # Get logical vector for dfP rows which are end points that precede a
            # "break", i.e. no segments overlap the end point.  These are the points
            # where the number of end points at or left of this end point is equal
            # to the number of start points left of it.
            breaksAt = (!dfP$start & (cumsumStart == (1:nrow(dfP)) - cumsumStart))

            # Note that since we sorted by dfP$id first, and since every row of df
            # produced two rows of dfP (one start and one end), it MUST be true that
            # at the endpoint of the last segment of each chromosome, a break occurs.
            # Double-check this during debugging.
            #lastChrEndpoint = (c(dfP$id[-1] != dfP$id[-nrow(dfP)], TRUE))
            #sum(lastChrEndpoint)
            #all((lastChrEndpoint & breaksAt) == lastChrEndpoint)
            #rm(lastChrEndpoint)

            # Get logical vector of dfP rows which do not overlap any other rows, i.e.
            # those rows that are both preceded and followed by a break.
            justOne = (diff(c(0, cumsumStart[breaksAt])) == 1) # Flags whether group just before each break in breaksAt has 1, or more, segments.
            if (sum(justOne) > 0)
                {
                # Every group that consists of just a SINGLE segment is a non-overlapping
                # segment.  Remove those non-overlapping rows from dfP.  We want breaksAt
                # to still be valid, so remove them from it also.  First, get the indexes
                # of the start and end points in dfP.
                dfP.idx.1.end = which(breaksAt)[justOne]
                dfP.idx.1.start = dfP.idx.1.end-1

                # We need to remove TWO ENTRIES from both dfP and breaksAt for each TRUE
                # member of justOne, because the TRUE member is for the END POINT, but the
                # START POINT should immediately precede it.  Test that the start point
                # DOES immediately precede the end point, during debugging.
                #if (!all(dfP[dfP.idx.1.start, "start"] == TRUE))
                #   stop("Point immediately preceding end point of non-overlapping segments is not always a start point")
                #if (!all(dfP[dfP.idx.1.start, "idx"] == dfP[dfP.idx.1.end, "idx"]))
                #   stop("Start point doesn't immediately precede end point of non-overlapping segments")
                # Perfect!

                # Ok, remove the two entries from dfP and breaksAt.
                dfP = dfP[-c(dfP.idx.1.start, dfP.idx.1.end),]
                breaksAt = breaksAt[-c(dfP.idx.1.start, dfP.idx.1.end)]
                rm(dfP.idx.1.start, dfP.idx.1.end)
                }
            rm(justOne, cumsumStart)      

            # If no more segments to check for overlap, break out of loop.
            if (nrow(dfP) == 0)
                break

            # Define endsAt to contain indexes into dfP of last end point of each group.
            endsAt = which(breaksAt)
            rm(breaksAt)

            # Define startsAt to contain indexes into dfP of first start point of
            # each group.
            startsAt = c(1, 1+endsAt[-length(endsAt)])

            # Find the shortest or longest segment within each overlap group and
            # put its df index into dfP$minmax.
            {
            if (minMax == "MIN")
                {
                dfP$minmax = unlist(sapply(1:length(startsAt), function(i)
                    {
                    ii = startsAt[i]:endsAt[i]
                    idxs = dfP[ii, "idx"][dfP[ii, "start"]]
                    keep.df.idx = idxs[which.min(df[idxs, "len"])[1]] # ***** MIN *****
                    return(rep(keep.df.idx, length(ii)))
                    }, simplify=FALSE), use.names=FALSE)
                }
            else
                {
                dfP$minmax = unlist(sapply(1:length(startsAt), function(i)
                    {
                    ii = startsAt[i]:endsAt[i]
                    idxs = dfP[ii, "idx"][dfP[ii, "start"]]
                    keep.df.idx = idxs[which.max(df[idxs, "len"])[1]] # ***** MAX *****
                    return(rep(keep.df.idx, length(ii)))
                    }, simplify=FALSE), use.names=FALSE)
                }
            }

            # Set dfP$overlap TRUE if the df segment indexed by dfP$idx overlaps
            # the df segment indexed by dfP$minmax.
            dfP$overlap = !(df$end[dfP$idx] <= df$start[dfP$minmax]) &
                    !(df$start[dfP$idx] >= df$end[dfP$minmax])

            # Get the set of df indexes to be discarded because they overlap.
            # This includes all the overlaps identified above EXCEPT the minmax
            # segments, which overlap themselves and so have the overlap flag set.
            # Add these indexes to the idxsToRemove vector and remove them from dfP.
            df.idxsToDiscard = setdiff(dfP$idx[dfP$overlap], dfP$minmax)
            idxsToRemove = c(idxsToRemove, df.idxsToDiscard)
            dfP = dfP[!dfP$idx %in% df.idxsToDiscard,]
            if (verbose)
                {
                cat("  Total number of", segName, "elements deleted:", length(idxsToRemove), "\n")
                cat("  Total number of", segName, "elements remaining:", nrow(df)-length(idxsToRemove), "\n")
                }
            }

        # Remove the df rows indexed by 'idxsToRemove'.
        if (length(idxsToRemove) > 0)
            df = df[-idxsToRemove,]
        if (verbose)
            catnow(" Finished ", nameSets, " ", setName, ", ", nrow(df), " ", segName, " elements remaining.\n", sep="")
        if (nrow(df) == 0)
            stop("There are no non-overlapping ", segName, " elements.")
        }
    # Get rid of the extra columns we added to df.
    df = df[, !colnames(df) %in% c("start", "end", "len")]
    if (verbose)
        catnow("Finished all ", nameSets, "s, have a total of ", nrow(df), " non-overlapping ", segName, " elements\n", sep="")

    # Put the data frame in order by reference setName position.
    if (verbose)
        catnow("Sorting by reference", nameSets, "position...")
    df = df[order(df[, idCols[1]], df[, pos1Cols[1]]),]
    rownames(df) = NULL
    if (verbose)
        catnow("\n")
    return(df)
    }

################################################################################
# End of file.
################################################################################
