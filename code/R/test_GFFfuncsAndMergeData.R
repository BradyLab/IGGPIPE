################################################################################
# Perform extensive testing of functions in Include_GFFfuncs.R and
# Include_MergeDataUsingPosition.R, using test data in folder test_GFFfuncsAndMergeData.
#
# This requires the person running this to know what is being done and examine output
# for correctness.
#
# Author: Ted Toal
# Date: 2015
# Brady Lab, UC Davis
################################################################################

# Set working directory.
setwd("~/Documents/UCDavis/BradyLab/Genomes/IGGPIPE/code/R")

# Source the files we are testing.
source("Include_Common.R")
source("Include_GFFfuncs.R")
source("Include_MergeDataUsingPosition.R")

# Read a test file data: a file of markers generated between Soly and Sope on
# the first 14 Mbp of ch01 and ch02.
dfMarkers = read.table("../../goodTest/MarkersOverlapping_K11k2L100D10_2000A100_2000d10_100N2F0X20V3000W8M3G1.tsv",
    header=TRUE, sep="\t", stringsAsFactors=FALSE)
head(dfMarkers)

# Read another test data file: IL introgression positions.
dfILs = read.table("test_GFFfuncsAndMergeData/ILintrogressions.tsv", header=TRUE, sep="\t", stringsAsFactors=FALSE)
dfILs$id = sprintf("SL2.50ch%02d", dfILs$chr)

################################################################################
# Testing of Include_GFFfuncs.R.
################################################################################

# Read additional test file data: a .gff3 and .gtf gene model file for Soly, truncated
# to first 14 Mbp of ch 01 and ch 02.
dfGff = readFile_GFF3("test_GFFfuncsAndMergeData/ITAG2.4_test.gff3")
head(dfGff)
dfGtf = readFile_GTF("test_GFFfuncsAndMergeData/ITAG2.4_test.gtf")
head(dfGtf)

getFeatures(dfGff)
getFeatures(dfGtf)

df1 = convertAttrsToCols(dfGff)
df2 = convertAttrsToCols(dfGtf)
dim(df1)
dim(df2)
df1[1,]
df2[1,]
df1[10,]
df2[10,]
df1a = convertColsToAttrs(df1, colnames(df1)[9:ncol(df1)])
df2a = convertColsToAttrs(df2, colnames(df2)[9:ncol(df2)])
df1a[1,]
df2a[1,]
df1a[10,]
df2a[10,]

dfx = convertAttrsToCols(dfGff, missingAttrVals="missing")
head(dfx)

dfx = convertAttrsToCols(dfGff, includeAttrs=c("Ontology_term", "nb_exon"),
    missingAttrVals=c(NA, "R.NA"))
head(dfx)

dfx = convertAttrsToCols(dfGff, includeAttrs=c("Ontology_term", "nb_exon"),
    missingAttrVals=c(NA, "R.NA"), newAttrCols=c("term", "exon"))
head(dfx)

dfx = convertAttrsToCols(dfGff, excludeAttrs=c("Ontology_term", "nb_exon"), removeAttrCol=FALSE)
head(dfx)


dfx = convertColsToAttrs(df1, "Alias", newAttrNames="synonym", noAttrValues="(None)",
    merge=FALSE, remove=FALSE)
head(dfx)


dir.create("testOutput", showWarnings=FALSE)
writeFile_GFF3_GTF(df1a, "testOutput/gff1.gff3")
writeFile_GFF3_GTF(df2a, "testOutput/gff1.gtf")

df1b = readFile_GFF3("testOutput/gff1.gff3")
df2b = readFile_GTF("testOutput/gff1.gtf")
df1b = convertAttrsToCols(df1b)
df2b = convertAttrsToCols(df2b)
identical(df1, df1b)
identical(df2, df2b)

df1c = clean_GFF3(dfGtf)
df2c = clean_GTF(dfGff)
df1c[1,]
df2c[1,]
identical(df1c, dfGff)
identical(df2c, dfGtf)

################################################################################
# Testing of Include_MergeDataUsingPosition.R.
################################################################################

# First test.
T.df = dfMarkers
S.df = dfILs
T.pos = list(id="Hid", start="HampPos1", end="HampPos2")
S.pos = list(id="id", start="start_right", end="end_left")
dist = list(method="OVERLAP")
cols = list(list(col="ILs", before="prmSeqL", maxMatch=0, join=TRUE, joinSep=";", format="{+IL_segment}({%S})"))
df = mergeOnMatches(T.df, S.df, T.pos, S.pos, dist, cols)
df[1:10,1:11]
df[1500:1510,1:11]
df[1600:1610,1:11]
tail(df[,1:11], 20)

#############################################
# Check cols arg
#############################################

# Check that we can insert as last column.
cols[[1]]$before = ""
df = mergeOnMatches(T.df, S.df, T.pos, S.pos, dist, cols)
df[1,]

# Check that we can insert just before last column.
cols[[1]]$before = "Pseq2"
df = mergeOnMatches(T.df, S.df, T.pos, S.pos, dist, cols)
df[1,]

# Check that we can insert as first column.
cols[[1]]$before = "NDA"
df = mergeOnMatches(T.df, S.df, T.pos, S.pos, dist, cols)
df[1,]

# Check that we can insert as second column.
cols[[1]]$before = "Hid"
df = mergeOnMatches(T.df, S.df, T.pos, S.pos, dist, cols)
df[1,]

# Restore "before" column.
cols[[1]]$before = "prmSeqL"

#############################################
# Check dist$methods arg
#############################################

dist$method = "OVERLAP"
df = mergeOnMatches(T.df, S.df, T.pos, S.pos, dist, cols)
dim(df)
sum(df$ILs != "") # No surprise, every amplicon is in at least one introgression.

dist$method = "S.TINY" # ILs are not tiny.  Should get empty result = error.
df = mergeOnMatches(T.df, S.df, T.pos, S.pos, dist, cols) # Yup.

dist$method = "T.TINY"
df = mergeOnMatches(T.df, S.df, T.pos, S.pos, dist, cols)
sum(df$ILs != "") # No surprise, all amplicons are INSIDE the introgressions, which are huge.

dist = list(method="S.NEAR", closest=0, start.up=0, end.up=0, start.down=0, end.down=0)
# ILs are huge compared to markers.  Will probably get empty result = error.
df = mergeOnMatches(T.df, S.df, T.pos, S.pos, dist, cols) # Yup.

dist$method = "T.NEAR"
df = mergeOnMatches(T.df, S.df, T.pos, S.pos, dist, cols)
sum(df$ILs != "") # No surprise, all amplicons are INSIDE the introgressions, which are huge.

# We will re-test dist$method with genes.  But now, let's test various output formats.

# Try inserting two different columns of S.df.
cols = list(
    list(col="ILs", before="prmSeqL", maxMatch=0, join=TRUE, joinSep=";", format="{+IL_segment}({%S})"),
    list(col="ILstart", before="prmSeqL", maxMatch=0, join=TRUE, joinSep=";", format="{+start_right}")
    )
df = mergeOnMatches(T.df, S.df, T.pos, S.pos, dist, cols)
df[1,]

# Try not joining columns.
dist$method = "OVERLAP"
cols = list(list(col="ILs", before="prmSeqL", maxMatch=0, join="NO", format="{+IL_segment}({%S})"))
df = mergeOnMatches(T.df, S.df, T.pos, S.pos, dist, cols)
df[1,]
unique(df$ILs1)
unique(df$ILs2)
unique(df$ILs3)
unique(df$ILs4)

# Test smaller maxMatch.
cols[[1]]$maxMatch=3
df = mergeOnMatches(T.df, S.df, T.pos, S.pos, dist, cols)
df[1,]

# Join columns with start and end strings.
cols = list(list(col="ILs", before="prmSeqL", maxMatch=0, join=TRUE, joinSep=";", joinPfx="[", joinSfx="]",
    format="{+IL_segment}({%S})"))
df = mergeOnMatches(T.df, S.df, T.pos, S.pos, dist, cols)
df[1,]

#############################################
# Test format options
#############################################

# Test {lb} and {rb}
cols = list(list(col="ILs", before="prmSeqL", maxMatch=0, join=TRUE, joinSep=";", format="{+IL_segment}({%S})"))
cols[[1]]$format = "{+IL_segment}{lb}{%S}{rb}"
df = mergeOnMatches(T.df, S.df, T.pos, S.pos, dist, cols)
df[100,]
dfILs[dfILs$IL_segment == "IL1-1-3",]
100*(260327-20659)/(4416584-20659)

# Test {#S}
cols[[1]]$format = "{+IL_segment}({#S})"
df = mergeOnMatches(T.df, S.df, T.pos, S.pos, dist, cols)
df[100,]
(260327-20659)

# Test {%T}
cols[[1]]$format = "{+IL_segment}({%T})"
df = mergeOnMatches(T.df, S.df, T.pos, S.pos, dist, cols)
df[100,]
100*(20659-260327)/(261232-260327)

# Test {#T}
cols[[1]]$format = "{+IL_segment}({#T})"
df = mergeOnMatches(T.df, S.df, T.pos, S.pos, dist, cols)
df[100,]
(20659-260327)

# Test {*S*col*val*dgts}
cols[[1]]$format = "{+IL_segment}[{+start_right}]({*S*start_right*1e-3*0}Kbp)"
df = mergeOnMatches(T.df, S.df, T.pos, S.pos, dist, cols)
df[100,]

# Test {*T*col*val*dgts}
cols[[1]]$format = "{+IL_segment}({*T*PampPos1*1e-3*1}Kbp)"
df = mergeOnMatches(T.df, S.df, T.pos, S.pos, dist, cols)
df[100,]

# Test {/S/col/RE/RE.replace}
cols[[1]]$format = "{+IL_segment}({/S/id/SL2.50ch0?/})"
df = mergeOnMatches(T.df, S.df, T.pos, S.pos, dist, cols)
df[100,]

cols[[1]]$format = "{+IL_segment}({/S/id/SL2.50ch0?/chrom})"
df = mergeOnMatches(T.df, S.df, T.pos, S.pos, dist, cols)
df[100,]

#############################################
# Test various bad arguments.
#############################################

df = mergeOnMatches(S.df=S.df, T.pos=T.pos, S.pos=S.pos, dist=dist, cols=cols)
df = mergeOnMatches(T.df=T.df, T.pos=T.pos, S.pos=S.pos, dist=dist, cols=cols)
df = mergeOnMatches(T.df=T.df, S.df=S.df, S.pos=S.pos, dist=dist, cols=cols)
df = mergeOnMatches(T.df=T.df, S.df=S.df, T.pos=T.pos, dist=dist, cols=cols)
df = mergeOnMatches(T.df=T.df, S.df=S.df, T.pos=T.pos, S.pos=S.pos, cols=cols)
df = mergeOnMatches(T.df=T.df, S.df=S.df, T.pos=T.pos, S.pos=S.pos, dist=dist)
df = mergeOnMatches(NULL, S.df, T.pos, S.pos, dist, cols)
df = mergeOnMatches(NA, S.df, T.pos, S.pos, dist, cols)
df = mergeOnMatches("A", S.df, T.pos, S.pos, dist, cols)
df = mergeOnMatches(1:10, S.df, T.pos, S.pos, dist, cols)
df = mergeOnMatches(T.df, NULL, T.pos, S.pos, dist, cols)
df = mergeOnMatches(T.df, NA, T.pos, S.pos, dist, cols)
df = mergeOnMatches(T.df, "A", T.pos, S.pos, dist, cols)
df = mergeOnMatches(T.df, 1:10, T.pos, S.pos, dist, cols)
df = mergeOnMatches(T.df, S.df, NULL, S.pos, dist, cols)
df = mergeOnMatches(T.df, S.df, NA, S.pos, dist, cols)
df = mergeOnMatches(T.df, S.df, list(end=0), S.pos, dist, cols)
df = mergeOnMatches(T.df, S.df, list(start="id"), S.pos, dist, cols)
df = mergeOnMatches(T.df, S.df, list(start="kmer1"), S.pos, dist, cols)
df = mergeOnMatches(T.df, S.df, list(start="Hpct"), S.pos, dist, cols)
df = mergeOnMatches(T.df, S.df, list(end="HampPos1", start="HampPos2"), S.pos, dist, cols)
df = mergeOnMatches(T.df, S.df, list(start="HampPos2", len="kmer1offset"), S.pos, dist, cols)
df = mergeOnMatches(T.df, S.df, T.pos, NULL, dist, cols)
df = mergeOnMatches(T.df, S.df, T.pos, NA, dist, cols)
df = mergeOnMatches(T.df, S.df, T.pos, list(end=0), dist, cols)
df = mergeOnMatches(T.df, S.df, T.pos, list(start="ID"), dist, cols)
df = mergeOnMatches(T.df, S.df, T.pos, list(start="id"), dist, cols)
df = mergeOnMatches(T.df, S.df, T.pos, S.pos, NULL, cols)
df = mergeOnMatches(T.df, S.df, T.pos, S.pos, NA, cols)
df = mergeOnMatches(T.df, S.df, T.pos, S.pos, list(ted=0), cols)
df = mergeOnMatches(T.df, S.df, T.pos, S.pos, list(method=0), cols)
df = mergeOnMatches(T.df, S.df, T.pos, S.pos, list(method="overlap"), cols)
df = mergeOnMatches(T.df, S.df, T.pos, S.pos, list(method="S.NEAR", closest="A"), cols)
df = mergeOnMatches(T.df, S.df, T.pos, S.pos, list(method="S.NEAR", closest=3), cols)
df = mergeOnMatches(T.df, S.df, T.pos, S.pos, list(method="S.NEAR", closest=NA), cols)
df = mergeOnMatches(T.df, S.df, T.pos, S.pos, list(method="S.NEAR", closest=1, start.up="A"), cols)
df = mergeOnMatches(T.df, S.df, T.pos, S.pos, list(method="S.NEAR", closest=1, start.up=7.5), cols)
df = mergeOnMatches(T.df, S.df, T.pos, S.pos, list(method="S.NEAR", closest=1, start.up=NA), cols)
df = mergeOnMatches(T.df, S.df, T.pos, S.pos, list(method="S.NEAR", closest=1, start.down=7.5), cols)
df = mergeOnMatches(T.df, S.df, T.pos, S.pos, list(method="S.NEAR", closest=1, end.up=7.5), cols)
df = mergeOnMatches(T.df, S.df, T.pos, S.pos, list(method="S.NEAR", closest=1, end.down=7.5), cols)
df = mergeOnMatches(T.df, S.df, T.pos, S.pos, dist)
df = mergeOnMatches(T.df, S.df, T.pos, S.pos, dist, NULL)
df = mergeOnMatches(T.df, S.df, T.pos, S.pos, dist, NA)
df = mergeOnMatches(T.df, S.df, T.pos, S.pos, dist, 3)
df = mergeOnMatches(T.df, S.df, T.pos, S.pos, dist, list())
df = mergeOnMatches(T.df, S.df, T.pos, S.pos, dist, list(5))
df = mergeOnMatches(T.df, S.df, T.pos, S.pos, dist, list(list()))
df = mergeOnMatches(T.df, S.df, T.pos, S.pos, dist, list(list(col=NA)))
df = mergeOnMatches(T.df, S.df, T.pos, S.pos, dist, list(list(col="NDA")))
df = mergeOnMatches(T.df, S.df, T.pos, S.pos, dist, list(list(col="X")))
df = mergeOnMatches(T.df, S.df, T.pos, S.pos, dist, list(list(col="X", before="X")))
df = mergeOnMatches(T.df, S.df, T.pos, S.pos, dist, list(list(col="X", before="")))
df = mergeOnMatches(T.df, S.df, T.pos, S.pos, dist, list(list(col="X", before="", format=NA)))
df = mergeOnMatches(T.df, S.df, T.pos, S.pos, dist, list(list(col="X", before="", format="{junk}")))
df = mergeOnMatches(T.df, S.df, T.pos, S.pos, dist, list(list(col="X", before="", format="X", maxMatch="A")))
df = mergeOnMatches(T.df, S.df, T.pos, S.pos, dist, list(list(col="X", before="", format="X", maxMatch=NA)))
df = mergeOnMatches(T.df, S.df, T.pos, S.pos, dist, list(list(col="X", before="", format="X", maxMatch=-1)))
df = mergeOnMatches(T.df, S.df, T.pos, S.pos, dist, list(list(col="X", before="", format="X", joinPfx=T.df)))
df = mergeOnMatches(T.df, S.df, T.pos, S.pos, dist, list(list(col="X", before="", format="X", joinSfx=T.df)))

#############################################
# Test using a GFF3 gene file with various dist$methods.
#############################################

# Read .gff3 file of genes in same chromosomes as the sequence used to create the markers in dfMarkers.
dfGff = readFile_GFF3("test_GFFfuncsAndMergeData/ITAG2.4_test.gff3")
head(dfGff)
dfGff = clean_GFF3(dfGff)
dfG = selectFeatures(dfGff, "gene")
dfG = convertAttrsToCols(dfG, includeAttrs="Alias")
dfG = dfG[, colnames(dfG) %in% c("Alias", "seqname", "start", "end", "strand")]
colnames(dfG) = c("id", "start", "end", "strand", "gene")
dfG = dfG[, c("gene", "id", "start", "end", "strand")]

# Now annotate markers with genes using different methods.

# method = OVERLAP
T.df = dfMarkers
S.df = dfG
T.pos = list(id="Hid", start="HampPos1", end="HampPos2")
S.pos = list(id="id", start="start", end="end")
dist = list(method="OVERLAP")
cols = list(list(col="genes", before="prmSeqL", maxMatch=0, join=TRUE, joinSep=";", format="{+gene}({#S})"))
df = mergeOnMatches(T.df, S.df, T.pos, S.pos, dist, cols)
sum(df$genes != "")
df[1:10,1:11]
df[1500:1510,1:11]
df[1600:1610,1:11]
tail(df[,1:11], 20)
df[1,]
head(dfG)
# Correct

# Try #T in place of #S
dist = list(method="OVERLAP")
cols = list(list(col="genes", before="prmSeqL", maxMatch=0, join=TRUE, joinSep=";", format="{+gene}({#T})"))
df = mergeOnMatches(T.df, S.df, T.pos, S.pos, dist, cols)
sum(df$genes != "")
df[1:10,1:11]
df[1500:1510,1:11]
df[1600:1610,1:11]
tail(df[,1:11], 20)
df[1,]
head(dfG)
# Correct
df.overlap = df[df$genes != "",]
dim(df.overlap)

# method = S.TINY : amplicon must contain entire gene
dist = list(method="S.TINY")
df = mergeOnMatches(T.df, S.df, T.pos, S.pos, dist, cols)
sum(df$genes != "") # Not very many (8), sounds good.
df[df$genes != "",1:11]
df[1488,]
dfG[dfG$gene == "Solyc01g009550",]
# Correct

# method = T.TINY : gene must contain entire amplicon
dist = list(method="T.TINY")
df = mergeOnMatches(T.df, S.df, T.pos, S.pos, dist, cols)
sum(df$genes != "")
dim(df) # Interesting, about 2/3 of the amplicons ARE within genes!!!
df[1:10,1:11]
df[1500:1510,1:11]
df[1600:1610,1:11]
tail(df[,1:11], 20)
df[1956,]
dfG[dfG$gene == "Solyc02g005520",]
# Correct

#############################################
# method = S.NEAR, closest = 0, distance 100bp : gene must be within 100 bp of amplicon (or contain it)
#############################################

#               start.up: G.start must be no less than M.start-twice max gene start-start distance
#               start.down: G.start must be no more than M.end+100 bp
#               end.up: G.end must be no less than M.start-100 bp
#               end.down: G.end must be no more than M.end+twice max gene end-end distance

dist = list(method="S.NEAR", closest=0, start.down=100, end.up=100)
df = mergeOnMatches(T.df, S.df, T.pos, S.pos, dist, cols)
sum(df$genes != "") # 1671
dim(df.overlap) # 1647
df.S.NEAR = df[df$genes != "",]
# Some of the markers overlapping genes should now also be near another gene, so
# the marker's gene list should have an additional gene in it.  Let's use HampPos1_HampPos2
# as row names for df.overlap and df.S.NEAR, which will help us.
rownames(df.overlap) = paste(df.overlap$HampPos1, df.overlap$HampPos2, sep="_")
rownames(df.S.NEAR) = paste(df.S.NEAR$HampPos1, df.S.NEAR$HampPos2, sep="_")
# Get annotated markers common between OVERLAP and S.NEAR, and those specific to each.
common = intersect(rownames(df.overlap), rownames(df.S.NEAR))
only.OVERLAP = setdiff(rownames(df.overlap), rownames(df.S.NEAR))
only.S.NEAR = setdiff(rownames(df.S.NEAR), rownames(df.overlap))
length(common) # 1647
length(only.OVERLAP) # 0, good, it should be that way.
length(only.S.NEAR) # 24
# For the common ones, get the ones whose gene annotation mismatches.
common.gene.mismatch = common[df.overlap[common, "genes"] != df.S.NEAR[common, "genes"]]
length(common.gene.mismatch) # Only 23 genes 
data.frame(overlap=df.overlap[common.gene.mismatch, "genes"], S.NEAR=df.S.NEAR[common.gene.mismatch, "genes"],
    stringsAsFactors=FALSE)
# Looks good, the S.NEAR have the same gene as OVERLAP plus one more gene next to it.
df.overlap[common.gene.mismatch[1],]
df.S.NEAR[common.gene.mismatch[1],]
dfG[dfG$gene %in% c("Solyc01g005070", "Solyc01g005080"),]
# Note: distance to gene, in (), is to START of gene, even though END of gene is closer.
# It looks good.
# Check an only.S.NEAR gene.
df.S.NEAR[only.S.NEAR[1],]
dfG[dfG$gene == "Solyc01g005560",]
# Looks good.

# Expand the near-distance to 3000 and we should get a lot more, including a bunch with
# multiple genes per marker.
T.df = dfMarkers
S.df = dfG
T.pos = list(id="Hid", start="HampPos1", end="HampPos2")
S.pos = list(id="id", start="start", end="end")
dist = list(method="S.NEAR", closest=0, start.down=3000, end.up=3000)
cols = list(list(col="genes", before="prmSeqL", maxMatch=0, join=TRUE, joinSep=";", format="{+gene}({#T})"))
df = mergeOnMatches(T.df, S.df, T.pos, S.pos, dist, cols)
sum(df$genes != "") # 1815 vs 1671
df.S.NEAR = df[df$genes != "",]
rownames(df.S.NEAR) = paste(df.S.NEAR$HampPos1, df.S.NEAR$HampPos2, sep="_")
common = intersect(rownames(df.overlap), rownames(df.S.NEAR))
only.OVERLAP = setdiff(rownames(df.overlap), rownames(df.S.NEAR))
only.S.NEAR = setdiff(rownames(df.S.NEAR), rownames(df.overlap))
length(common) # 1647
length(only.OVERLAP) # 0, good, it should be that way.
length(only.S.NEAR) # 168
common.gene.mismatch = common[df.overlap[common, "genes"] != df.S.NEAR[common, "genes"]]
length(common.gene.mismatch) # Now 685!
data.frame(name=common.gene.mismatch,
    overlap=df.overlap[common.gene.mismatch, "genes"], S.NEAR=df.S.NEAR[common.gene.mismatch, "genes"],
    stringsAsFactors=FALSE)
# Looks good.  One with 3 genes is 4683900_4684979
df.S.NEAR[only.S.NEAR[1], "genes"]
only.S.NEAR[1]
# 54857_55150 has 3 genes (30/40/50), none of them overlapping.

# Repeat but use closest=1.  We should get only the overlapping genes in the common
# set, and in the unique set we should get only the closest one.
dist = list(method="S.NEAR", closest=1, start.down=3000, end.up=3000)
df = mergeOnMatches(T.df, S.df, T.pos, S.pos, dist, cols)
sum(df$genes != "") # Still 1815
df.S.NEAR = df[df$genes != "",]
rownames(df.S.NEAR) = paste(df.S.NEAR$HampPos1, df.S.NEAR$HampPos2, sep="_")
common = intersect(rownames(df.overlap), rownames(df.S.NEAR))
only.OVERLAP = setdiff(rownames(df.overlap), rownames(df.S.NEAR))
only.S.NEAR = setdiff(rownames(df.S.NEAR), rownames(df.overlap))
length(common) # 1647
length(only.OVERLAP) # 0, good, it should be that way.
length(only.S.NEAR) # 168
common.gene.mismatch = common[df.overlap[common, "genes"] != df.S.NEAR[common, "genes"]]
length(common.gene.mismatch) # Now 0 because if it has an overlap we don't include ANY nearby ones
# We should only see genes 40 and 50 here, because 60 only showed up in the overlap set
df.S.NEAR["4683900_4684979","genes"]
# We should see only nearest one of 30, 40, 50, which is 30 (I looked).
df.S.NEAR["54857_55150","genes"]

# Repeat but use closest=2.  We should get only the overlapping genes in the common
# set, and in the unique set we should get only the closest UPSTREAM *AND* closest
# DOWNSTREAM one.
dist = list(method="S.NEAR", closest=2, start.down=3000, end.up=3000)
df = mergeOnMatches(T.df, S.df, T.pos, S.pos, dist, cols)
sum(df$genes != "") # Still 1815
df.S.NEAR = df[df$genes != "",]
rownames(df.S.NEAR) = paste(df.S.NEAR$HampPos1, df.S.NEAR$HampPos2, sep="_")
common = intersect(rownames(df.overlap), rownames(df.S.NEAR))
only.OVERLAP = setdiff(rownames(df.overlap), rownames(df.S.NEAR))
only.S.NEAR = setdiff(rownames(df.S.NEAR), rownames(df.overlap))
length(common) # 1647
length(only.OVERLAP) # 0, good, it should be that way.
length(only.S.NEAR) # 168
common.gene.mismatch = common[df.overlap[common, "genes"] != df.S.NEAR[common, "genes"]]
length(common.gene.mismatch) # Still 0
# We should only see genes 40 and 50 here, because 60 only showed up in the overlap set
df.S.NEAR["4683900_4684979","genes"]
# We should see only the nearest upstream gene, which is 30, and nearest downstream
# gene, which is 40 (I looked).
df.S.NEAR["54857_55150","genes"]

#############################################
# method = T.NEAR, closest = 0, distance 100bp : amplicon must be within 100 bp of gene (or contain it)
#############################################

# Now reverse it and use T.NEAR, so we are looking for markers that are near genes.
# It should actually generate the same results as with S.NEAR., right?
dist = list(method="T.NEAR", closest=0, start.down=100, end.up=100)
df = mergeOnMatches(T.df, S.df, T.pos, S.pos, dist, cols)
sum(df$genes != "") # 1671
dim(df.overlap) # 1647
df.T.NEAR = df[df$genes != "",]
# Some of the genes overlapping markers should now also be near another marker, so
# the gene list should have an additional gene in it.  Let's use HampPos1_HampPos2
# as row names for df.overlap and df.T.NEAR, which will help us.
rownames(df.overlap) = paste(df.overlap$HampPos1, df.overlap$HampPos2, sep="_")
rownames(df.T.NEAR) = paste(df.T.NEAR$HampPos1, df.T.NEAR$HampPos2, sep="_")
# Get annotated markers common between OVERLAP and T.NEAR, and those specific to each.
common = intersect(rownames(df.overlap), rownames(df.T.NEAR))
only.OVERLAP = setdiff(rownames(df.overlap), rownames(df.T.NEAR))
only.T.NEAR = setdiff(rownames(df.T.NEAR), rownames(df.overlap))
length(common) # 1647
length(only.OVERLAP) # 0, good, it should be that way.
length(only.T.NEAR) # 24
# For the common ones, get the ones whose gene annotation mismatches.
common.gene.mismatch = common[df.overlap[common, "genes"] != df.T.NEAR[common, "genes"]]
length(common.gene.mismatch) # Only 23 genes 
data.frame(overlap=df.overlap[common.gene.mismatch, "genes"], S.NEAR=df.T.NEAR[common.gene.mismatch, "genes"],
    stringsAsFactors=FALSE)
# Looks good, the T.NEAR have the same gene as OVERLAP plus one more gene next to it.
df.overlap[common.gene.mismatch[1],]
df.T.NEAR[common.gene.mismatch[1],]
dfG[dfG$gene %in% c("Solyc01g005070", "Solyc01g005080"),]
# Note: distance to gene, in (), is to START of gene, even though END of gene is closer.
# It looks good.
# Check an only.T.NEAR gene.
df.T.NEAR[only.T.NEAR[1],]
dfG[dfG$gene == "Solyc01g005560",]
# Looks good.

# Expand the near-distance to 3000 and we should get a lot more, including a bunch with
# multiple genes per marker.
T.df = dfMarkers
S.df = dfG
T.pos = list(id="Hid", start="HampPos1", end="HampPos2")
S.pos = list(id="id", start="start", end="end")
dist = list(method="T.NEAR", closest=0, start.down=3000, end.up=3000)
cols = list(list(col="genes", before="prmSeqL", maxMatch=0, join=TRUE, joinSep=";", format="{+gene}({#T})"))
df = mergeOnMatches(T.df, S.df, T.pos, S.pos, dist, cols)
sum(df$genes != "") # 1815 vs 1671
df.T.NEAR = df[df$genes != "",]
rownames(df.T.NEAR) = paste(df.T.NEAR$HampPos1, df.T.NEAR$HampPos2, sep="_")
common = intersect(rownames(df.overlap), rownames(df.T.NEAR))
only.OVERLAP = setdiff(rownames(df.overlap), rownames(df.T.NEAR))
only.T.NEAR = setdiff(rownames(df.T.NEAR), rownames(df.overlap))
length(common) # 1647
length(only.OVERLAP) # 0, good, it should be that way.
length(only.T.NEAR) # 168
common.gene.mismatch = common[df.overlap[common, "genes"] != df.T.NEAR[common, "genes"]]
length(common.gene.mismatch) # Now 685!
data.frame(name=common.gene.mismatch,
    overlap=df.overlap[common.gene.mismatch, "genes"], S.NEAR=df.T.NEAR[common.gene.mismatch, "genes"],
    stringsAsFactors=FALSE)
# Looks good.  One with 3 genes is 4683900_4684979
df.T.NEAR[only.T.NEAR[1], "genes"]
only.T.NEAR[1]
# 54857_55150 has 3 genes (30/40/50), none of them overlapping.

# Repeat but use closest=1.  We should get only the overlapping genes in the common
# set, and in the unique set we should get only the closest marker to each gene,
# which will give a different result than S.NEAR which was closest gene to each
# marker.  With S.NEAR, closest=1 simply reduced the number of genes for each marker,
# but each marker that had genes still had genes.  With T.NEAR, some markers that
# picked up their nearest gene now have that gene picking its closest marker, which
# is a different marker.  So, we should end up with fewer markers with genes.
dist = list(method="T.NEAR", closest=1, start.down=3000, end.up=3000)
Solyc01g005570
df = mergeOnMatches(T.df, S.df, T.pos, S.pos, dist, cols)
sum(df$genes != "") # Now 1691 instead of 1815.  Fewer, as expected.
df.T.NEAR = df[df$genes != "",]
rownames(df.T.NEAR) = paste(df.T.NEAR$HampPos1, df.T.NEAR$HampPos2, sep="_")
common = intersect(rownames(df.overlap), rownames(df.T.NEAR))
only.OVERLAP = setdiff(rownames(df.overlap), rownames(df.T.NEAR))
only.T.NEAR = setdiff(rownames(df.T.NEAR), rownames(df.overlap))
length(common) # still 1647
length(only.OVERLAP) # 0, good, it should be that way.
length(only.T.NEAR) # 44 instead of 168, fewer, as expected.
common.gene.mismatch = common[df.overlap[common, "genes"] != df.T.NEAR[common, "genes"]]
length(common.gene.mismatch) # Now 57 instead of 0.  If a gene overlaps a marker, it gets just that marker, but if not, it gets nearest marker, and so one marker may get several nearest genes even though it overlaps another marker.
# We should only see genes 40 and 50 here, because 60 only showed up in the overlap set
df.T.NEAR["4683900_4684979","genes"]
# Is this marker the closest marker to which genes?  Apparently to 40 and 50.
df.T.NEAR["54857_55150","genes"]

# Repeat but use closest=2.  We should get only the overlapping genes in the common
# set, and in the unique set we should get MORE because now we have both upstream
# and downstream.
dist = list(method="T.NEAR", closest=2, start.down=3000, end.up=3000)
df = mergeOnMatches(T.df, S.df, T.pos, S.pos, dist, cols)
sum(df$genes != "") # 1692 now, we picked up one more, would have thought more than that.
df.T.NEAR = df[df$genes != "",]
rownames(df.T.NEAR) = paste(df.T.NEAR$HampPos1, df.T.NEAR$HampPos2, sep="_")
common = intersect(rownames(df.overlap), rownames(df.T.NEAR))
only.OVERLAP = setdiff(rownames(df.overlap), rownames(df.T.NEAR))
only.T.NEAR = setdiff(rownames(df.T.NEAR), rownames(df.overlap))
length(common) # still 1647
length(only.OVERLAP) # 0, good, it should be that way.
length(only.T.NEAR) # 45 now instead of 4, we gained one.
common.gene.mismatch = common[df.overlap[common, "genes"] != df.T.NEAR[common, "genes"]]
length(common.gene.mismatch) # Now 62 instead of 57, gained 5.
# We should only see genes 40 and 50 here, because 60 only showed up in the overlap set
df.T.NEAR["4683900_4684979","genes"]
# Is this marker the closest marker to which genes?  Apparently to 40 and 50.
df.T.NEAR["54857_55150","genes"]

################################################################################
# End of file.
################################################################################
