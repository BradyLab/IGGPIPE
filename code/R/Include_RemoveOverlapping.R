################################################################################
# This file contains R definitions and functions for removing overlapping rows
# of a data frame (that containing genome position data), leaving only rows that
# do not overlap one another.  This is to be sourced by R code that needs to
# use it.
# Author: Ted Toal
# Date: 2015
# Brady Lab, UC Davis
################################################################################

################################################################################
# Remove overlapping rows of a data frame.
#
# Arguments:
#   df: data frame to be processed.
#
# Returns: data frame containing the data, with columns "seqname", "source", "feature",
#   "start", "end", "score", "strand", "frame", "attributes".
################################################################################


catnow("Getting non-overlapping Indel Groups...")
inv(nrow(dfIGs), "\nNumber of Indel Groups including overlapping ones")

# We will use the same data frame that holds the overlapping Indel Groups, dfIGs,
# to hold the non-overlapping InDel Groups, since there might be a LOT of them
# (so copying the data frame would be costly) and we need to start out this algorithm
# with the data frame containing the overlapping InDel Groups, which dfIGs already
# does.  The data frame will be transformed into one containing only non-overlapping
# Indel Groups.  We no longer need the overlapping groups (except for this).

# Test each genome, one by one, for Indel Group overlaps, and remove Indel Groups
# to get rid of them.
for (genome in genomes)
    {
    inv(genome, "Remove overlaps")
    id.Col = idCol.I[genome]
    pos1.Col = pos1Col.I[genome]
    pos2.Col = pos2Col.I[genome]

    # The method will be:
    #   1. Find, for each Indel Group X with end position X.E, the row index of
    #       the Indel Group Y with starting position Y.S such that X.E > Y.S.
    #       (We count X.E == Y.S as non-overlapping, i.e. Indel Group AB does
    #       not overlap Indel Group BC).  Then, set new column "overlap" TRUE
    #       if that row index is not equal to the row index of Indel Group X.
    #   2. Find the start and end index of each continuous set of TRUE overlap
    #       values, those being the start and end of a group of Indel Groups
    #       that all overlap one another in some manner.
    #   3. For each such start and end index, find the index of the Indel Group
    #       in that row index range which has the smallest (MIN) or largest (MAX)
    #       length.
    #   4. Again for each such start and end index, remove all Indel Groups in
    #       that index range which overlap with the one with the shortest or
    #       longest length.
    #   5. Repeat steps 1-4 until there are no more overlaps.

    # Loop until no more Indel Groups are found to overlap in this genome.
    while (TRUE)
        {
        inv(nrow(dfIGs), "Loop with # Indel Groups remaining")

        # Copy start or end position, whichever is smaller (it varies depending on strand)
        # into new column "start", and vice-versa, copy the larger position into new
        # column "end".
        dfIGs$start = dfIGs[, pos1.Col]
        dfIGs$end = dfIGs[, pos2.Col]
        pos2IsSmaller = (dfIGs[,pos1.Col] > dfIGs[,pos2.Col])
        dfIGs$start[pos2IsSmaller] = dfIGs[pos2IsSmaller, pos2.Col]
        dfIGs$end[pos2IsSmaller] = dfIGs[pos2IsSmaller, pos1.Col]

        # Add new column "len" equal to length of the Indel Group segments.
        dfIGs$len = dfIGs$end - dfIGs$start + 1

        # Sort by ID and start position.
        dfIGs = dfIGs[order(dfIGs[, id.Col], dfIGs$start),]
        N = nrow(dfIGs)

        # Find index of Indel Group through which each Indel Group overlaps and
        # put it in column thruIdx.
        dfIGs$thruIdx = NA
        for (id in unique(dfIGs[, id.Col]))
            {
            thisId = (dfIGs[, id.Col] == id)
            dfIGs$thruIdx[thisId] = match(TRUE, thisId) - 1 +
                findInterval(dfIGs$end[thisId], dfIGs$start[thisId]+1)
            }

        # Set thruIdx of Indel Groups that do not overlap even the next Indel
        # Group to 0.
        dfIGs$thruIdx[dfIGs$thruIdx == 1:nrow(dfIGs)] = 0

        # Get the set of indexes of all Indel Groups which overlap at least one
        # other Indel Group.
        overlapIdxs = sapply(1:nrow(dfIGs), function(i)
            {
            if (dfIGs$thruIdx[i] == 0)
                return(0)
            return(i:dfIGs$thruIdx[i])
            })
        overlapIdxs = unique(unlist(overlapIdxs))
        # Remove index 0, which comes from non-overlapping Indel Groups.
        overlapIdxs = overlapIdxs[overlapIdxs != 0]

        # If no Indel Groups overlap, break out of loop.
        inv(length(overlapIdxs), "Number of overlapping Indel Groups")
        if (length(overlapIdxs) == 0)
            break

        # Set column "overlap" TRUE for each of those overlapping Indel Groups.
        dfIGs$overlap = FALSE
        dfIGs$overlap[overlapIdxs] = TRUE

        # Get the start and end index of each group of mutually overlapping Indel Groups.
        startOverlap = which(!c(FALSE, dfIGs$overlap[-N]) & dfIGs$overlap)
        endOverlap = which(dfIGs$overlap & !c(dfIGs$overlap[-1], FALSE))
        if (length(startOverlap) != length(endOverlap)) stop("Expected equal start/end overlap vectors")

        # Find the index within each group of that Indel Group with the shortest
        # or longest length.  Then get the indexes within the group of Indel Groups
        # that overlap that shortest/longest Indel Group.
        d = dfIGs[, c("len", "start", "end")] # Smaller data frame to work with.
        if (minMax == "MIN")
            idxsToRemove = sapply(1:length(startOverlap), function(i)
                {
                idxs = startOverlap[i]:endOverlap[i]
                dft = d[idxs,,drop=FALSE]
                idx.SL = which.min(dft$len)[1] # ***** MIN *****
                dft.SL = dft[idx.SL,,drop=FALSE]
                dft = dft[-idx.SL,,drop=FALSE]
                idxs = idxs[-idx.SL]
                return(idxs[!(dft$end <= dft.SL$start) & !(dft$start >= dft.SL$end)])
                })
        else
            idxsToRemove = sapply(1:length(startOverlap), function(i)
                {
                idxs = startOverlap[i]:endOverlap[i]
                dft = d[idxs,,drop=FALSE]
                idx.SL = which.max(dft$len)[1] # ***** MAX *****
                dft.SL = dft[idx.SL,,drop=FALSE]
                dft = dft[-idx.SL,,drop=FALSE]
                idxs = idxs[-idx.SL]
                return(idxs[!(dft$end <= dft.SL$start) & !(dft$start >= dft.SL$end)])
                })
        idxsToRemove = sapply(1:length(startOverlap), function(i)
            {
            idxs = startOverlap[i]:endOverlap[i]
            len.SL = ifelse(minMax == "MIN", min(dfIGs$len[idxs]), max(dfIGs$len[idxs]))
            idx.SL = idxs[dfIGs$len[idxs] == len.SL][1] # If more than one, pick the first.
            return(idxs[!(dfIGs$end[idxs] <= dfIGs$start[idx.SL]) &
                !(dfIGs$start[idxs] >= dfIGs$end[idx.SL]) & idxs != idx.SL])
            })
        idxsToRemove = unlist(idxsToRemove)

        # Remove those Indel Groups.
        dfIGs = dfIGs[-idxsToRemove,]
        inv(length(idxsToRemove), "Number of Indel Groups deleted")
        inv(nrow(dfIGs), "Number of Indel Groups remaining")
        }
    }
dfIGs = dfIGs[, !colnames(dfIGs) %in% c("start", "end", "len", "thruIdx", "overlap")]
inv(nrow(dfIGs), "Number of Indel Groups with overlapping Indel Groups removed")
catnow("\n")

# Put the data frame in order by reference genome position.
catnow("Sorting by reference genome position...")
dfIGs = dfIGs[order(dfIGs[, idCol.I[refGenome]], dfIGs[, pos1Col.I[refGenome]]),]
rownames(dfIGs) = NULL
if (nrow(dfIGs) == 0)
    stop("There are no non-overlapping Indel Groups.")
catnow("\n")


################################################################################
# End of file.
################################################################################
