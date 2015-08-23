#######################################################################################
# Perform extensive testing of functions in Include_GFFfuncs.R and
# Include_MergeDataUsingPosition.R, using test data in folder test_GFFfuncsAndMergeData.
#######################################################################################

# Include the root R source code file.
R_INCLUDE_DIR = Sys.getenv("R_INCLUDE_DIR")
if (R_INCLUDE_DIR == "")
    stop("Set environment variable R_INCLUDE_DIR or assign that variable to a path here.")
source(paste(R_INCLUDE_DIR, "Include_RootRsourceFile.R", sep=""))

# Source our basic utility functions file.
source(Include_UtilityFunctions)

# Load libraries we need.
addThisLibrary("xlsx")

# Source other files we need.
source(Include_GenomeDb)
source(Include_atGenome)
source(Include_solyGenome)
source(Include_sopeGenome)

# Set working directory.
setwd(paste(genomeDbDir, "Code", sep=PATHSEP))

# Source the files we are testing.
source("Include_GFFfuncs.R")
source("Include_MergeDataUsingPosition.R")

# Read a test file data: a file of markers generated between Soly and Sope on
# the first 14 Mbp of ch01 and ch02.
dfMarkers = read.table("test_GFFfuncsAndMergeData/MarkersOverlapping.test.tsv", header=TRUE, sep="\t", stringsAsFactors=FALSE)
head(dfMarkers)

# Read another test data file: IL introgression positions.
dfILs = read.table("test_GFFfuncsAndMergeData/ILintrogressions.tsv", header=TRUE, sep="\t", stringsAsFactors=FALSE)
dfILs$id = sprintf("SL2.50ch%02d", dfILs$chr)

#######################################################################################
# Testing of Include_GFFfuncs.R.
#######################################################################################

# Read additional test file data: a .gff3 and .gtf gene model file for Soly, truncated
# to first 14 Mbp of ch 01 and ch 02.
dfGff = readFile_GFF3("test_GFFfuncsAndMergeData/test.gff3")
head(dfGff)
dfGtf = readFile_GTF("test_GFFfuncsAndMergeData/test.gtf")
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

writeFile_GFF3_GTF(df1a, "testOutput/gff1.gff3")
writeFile_GFF3_GTF(df2a, "testOutput/gff1.gtf")

df1b = readFile_GFF3_GTF("testOutput/gff1.gff3")
df2b = readFile_GFF3_GTF("testOutput/gff1.gtf")
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

#######################################################################################
# Testing of Include_MergeDataUsingPosition.R.
#######################################################################################

# First test.
t.df = dfMarkers
s.df = dfILs
t.pos = list(id="Hid", start="HampPos1", end="HampPos2")
s.pos = list(id="id", start="start_right", end="end_left")
dist = list(method="OVERLAP")
cols = list(list(col="ILs", before="prmSeqL", maxMatch=0, join="YES", joinSep=";", format="{+IL_segment}({%s})"))
df = mergeOnMatches(t.df, s.df, t.pos, s.pos, dist, cols)
df[1:10,1:11]
df[1500:1510,1:11]
df[1600:1610,1:11]
df[2000:2010,1:11]
tail(df[,1:11], 20)

#############################################
# Check cols arg
#############################################

# Check that we can insert as last column.
cols[[1]]$before = ""
df = mergeOnMatches(t.df, s.df, t.pos, s.pos, dist, cols)
df[1,]

# Check that we can insert just before last column.
cols[[1]]$before = "Pseq2"
df = mergeOnMatches(t.df, s.df, t.pos, s.pos, dist, cols)
df[1,]

# Check that we can insert as first column.
cols[[1]]$before = "NDA"
df = mergeOnMatches(t.df, s.df, t.pos, s.pos, dist, cols)
df[1,]

# Check that we can insert as second column.
cols[[1]]$before = "Hid"
df = mergeOnMatches(t.df, s.df, t.pos, s.pos, dist, cols)
df[1,]

# Restore "before" column.
cols[[1]]$before = "prmSeqL"

#############################################
# Check dist$methods arg
#############################################

dist$method = "OVERLAP"
df = mergeOnMatches(t.df, s.df, t.pos, s.pos, dist, cols)
dim(df)
sum(df$ILs != "") # No surprise, every amplicon is in at least one introgression.

dist$method = "s.TINY" # ILs are not tiny.  Should get empty result = error.
df = mergeOnMatches(t.df, s.df, t.pos, s.pos, dist, cols) # Yup.

dist$method = "t.TINY"
df = mergeOnMatches(t.df, s.df, t.pos, s.pos, dist, cols)
sum(df$ILs != "") # No surprise, all amplicons are INSIDE the introgressions, which are huge.

dist = list(method="s.NEAR", closest=0, start.up=0, end.up=0, start.down=0, end.down=0)
# ILs are huge compared to markers.  Will probably get empty result = error.
df = mergeOnMatches(t.df, s.df, t.pos, s.pos, dist, cols) # Yup.

dist$method = "t.NEAR"
df = mergeOnMatches(t.df, s.df, t.pos, s.pos, dist, cols)
sum(df$ILs != "") # No surprise, all amplicons are INSIDE the introgressions, which are huge.

# We will re-test dist$method with genes.  But now, let's test various output formats.

# Try inserting two different columns of s.df.
cols = list(
    list(col="ILs", before="prmSeqL", maxMatch=0, join="YES", joinSep=";", format="{+IL_segment}({%s})"),
    list(col="ILstart", before="prmSeqL", maxMatch=0, join="YES", joinSep=";", format="{+start_right}")
    )
df = mergeOnMatches(t.df, s.df, t.pos, s.pos, dist, cols)
df[1,]

# Try not joining columns.
dist$method = "OVERLAP"
cols = list(list(col="ILs", before="prmSeqL", maxMatch=0, join="NO", format="{+IL_segment}({%s})"))
df = mergeOnMatches(t.df, s.df, t.pos, s.pos, dist, cols)
df[1,]
unique(df$ILs1)
unique(df$ILs2)
unique(df$ILs3)
unique(df$ILs4)

# Test smaller maxMatch.
cols[[1]]$maxMatch=3
df = mergeOnMatches(t.df, s.df, t.pos, s.pos, dist, cols)
df[1,]

# Join columns with start and end strings.
cols = list(list(col="ILs", before="prmSeqL", maxMatch=0, join="YES", joinSep=";", joinStart="[", joinEnd="]",
    format="{+IL_segment}({%s})"))
df = mergeOnMatches(t.df, s.df, t.pos, s.pos, dist, cols)
df[1,]

#############################################
# Test format options
#############################################

# Test {lb} and {rb}
cols = list(list(col="ILs", before="prmSeqL", maxMatch=0, join="YES", joinSep=";", format="{+IL_segment}({%s})"))
cols[[1]]$format = "{+IL_segment}{lb}{%s}{rb}"
df = mergeOnMatches(t.df, s.df, t.pos, s.pos, dist, cols)
df[100,]
dfILs[dfILs$IL_segment == "IL1-1-3",]
100*(260327-20659)/(4416584-20659)

# Test {#s}
cols[[1]]$format = "{+IL_segment}({#s})"
df = mergeOnMatches(t.df, s.df, t.pos, s.pos, dist, cols)
df[100,]
(260327-20659)

# Test {%t}
cols[[1]]$format = "{+IL_segment}({%t})"
df = mergeOnMatches(t.df, s.df, t.pos, s.pos, dist, cols)
df[100,]
100*(20659-260327)/(261232-260327)

# Test {#t}
cols[[1]]$format = "{+IL_segment}({#t})"
df = mergeOnMatches(t.df, s.df, t.pos, s.pos, dist, cols)
df[100,]
(20659-260327)

# Test {*s*col*val*dgts}
cols[[1]]$format = "{+IL_segment}[{+start_right}]({*s*start_right*1e-3*0}Kbp)"
df = mergeOnMatches(t.df, s.df, t.pos, s.pos, dist, cols)
df[100,]

# Test {*t*col*val*dgts}
cols[[1]]$format = "{+IL_segment}({*t*PampPos1*1e-3*1}Kbp)"
df = mergeOnMatches(t.df, s.df, t.pos, s.pos, dist, cols)
df[100,]

# Test {/s/col/RE/RE.replace}
cols[[1]]$format = "{+IL_segment}({/s/id/SL2.50ch0?/})"
df = mergeOnMatches(t.df, s.df, t.pos, s.pos, dist, cols)
df[100,]

cols[[1]]$format = "{+IL_segment}({/s/id/SL2.50ch0?/chrom})"
df = mergeOnMatches(t.df, s.df, t.pos, s.pos, dist, cols)
df[100,]

#############################################
# Test various bad arguments.
#############################################

df = mergeOnMatches(s.df=s.df, t.pos=t.pos, s.pos=s.pos, dist=dist, cols=cols)
df = mergeOnMatches(t.df=t.df, t.pos=t.pos, s.pos=s.pos, dist=dist, cols=cols)
df = mergeOnMatches(t.df=t.df, s.df=s.df, s.pos=s.pos, dist=dist, cols=cols)
df = mergeOnMatches(t.df=t.df, s.df=s.df, t.pos=t.pos, dist=dist, cols=cols)
df = mergeOnMatches(t.df=t.df, s.df=s.df, t.pos=t.pos, s.pos=s.pos, cols=cols)
df = mergeOnMatches(t.df=t.df, s.df=s.df, t.pos=t.pos, s.pos=s.pos, dist=dist)
df = mergeOnMatches(NULL, s.df, t.pos, s.pos, dist, cols)
df = mergeOnMatches(NA, s.df, t.pos, s.pos, dist, cols)
df = mergeOnMatches("A", s.df, t.pos, s.pos, dist, cols)
df = mergeOnMatches(1:10, s.df, t.pos, s.pos, dist, cols)
df = mergeOnMatches(t.df, NULL, t.pos, s.pos, dist, cols)
df = mergeOnMatches(t.df, NA, t.pos, s.pos, dist, cols)
df = mergeOnMatches(t.df, "A", t.pos, s.pos, dist, cols)
df = mergeOnMatches(t.df, 1:10, t.pos, s.pos, dist, cols)
df = mergeOnMatches(t.df, s.df, NULL, s.pos, dist, cols)
df = mergeOnMatches(t.df, s.df, NA, s.pos, dist, cols)
df = mergeOnMatches(t.df, s.df, list(end=0), s.pos, dist, cols)
df = mergeOnMatches(t.df, s.df, list(start="id"), s.pos, dist, cols)
df = mergeOnMatches(t.df, s.df, list(start="kmer1"), s.pos, dist, cols)
df = mergeOnMatches(t.df, s.df, list(start="Hpct"), s.pos, dist, cols)
df = mergeOnMatches(t.df, s.df, list(end="HampPos1", start="HampPos2"), s.pos, dist, cols)
df = mergeOnMatches(t.df, s.df, list(start="HampPos2", len="kmer1offset"), s.pos, dist, cols)
df = mergeOnMatches(t.df, s.df, t.pos, NULL, dist, cols)
df = mergeOnMatches(t.df, s.df, t.pos, NA, dist, cols)
df = mergeOnMatches(t.df, s.df, t.pos, list(end=0), dist, cols)
df = mergeOnMatches(t.df, s.df, t.pos, list(start="ID"), dist, cols)
df = mergeOnMatches(t.df, s.df, t.pos, list(start="id"), dist, cols)
df = mergeOnMatches(t.df, s.df, t.pos, s.pos, NULL, cols)
df = mergeOnMatches(t.df, s.df, t.pos, s.pos, NA, cols)
df = mergeOnMatches(t.df, s.df, t.pos, s.pos, list(ted=0), cols)
df = mergeOnMatches(t.df, s.df, t.pos, s.pos, list(method=0), cols)
df = mergeOnMatches(t.df, s.df, t.pos, s.pos, list(method="overlap"), cols)
df = mergeOnMatches(t.df, s.df, t.pos, s.pos, list(method="s.NEAR", closest="A"), cols)
df = mergeOnMatches(t.df, s.df, t.pos, s.pos, list(method="s.NEAR", closest=3), cols)
df = mergeOnMatches(t.df, s.df, t.pos, s.pos, list(method="s.NEAR", closest=NA), cols)
df = mergeOnMatches(t.df, s.df, t.pos, s.pos, list(method="s.NEAR", closest=1, start.up="A"), cols)
df = mergeOnMatches(t.df, s.df, t.pos, s.pos, list(method="s.NEAR", closest=1, start.up=7.5), cols)
df = mergeOnMatches(t.df, s.df, t.pos, s.pos, list(method="s.NEAR", closest=1, start.up=NA), cols)
df = mergeOnMatches(t.df, s.df, t.pos, s.pos, list(method="s.NEAR", closest=1, start.down=7.5), cols)
df = mergeOnMatches(t.df, s.df, t.pos, s.pos, list(method="s.NEAR", closest=1, end.up=7.5), cols)
df = mergeOnMatches(t.df, s.df, t.pos, s.pos, list(method="s.NEAR", closest=1, end.down=7.5), cols)
df = mergeOnMatches(t.df, s.df, t.pos, s.pos, dist)
df = mergeOnMatches(t.df, s.df, t.pos, s.pos, dist, NULL)
df = mergeOnMatches(t.df, s.df, t.pos, s.pos, dist, NA)
df = mergeOnMatches(t.df, s.df, t.pos, s.pos, dist, 3)
df = mergeOnMatches(t.df, s.df, t.pos, s.pos, dist, list())
df = mergeOnMatches(t.df, s.df, t.pos, s.pos, dist, list(5))
df = mergeOnMatches(t.df, s.df, t.pos, s.pos, dist, list(list()))
df = mergeOnMatches(t.df, s.df, t.pos, s.pos, dist, list(list(col=NA)))
df = mergeOnMatches(t.df, s.df, t.pos, s.pos, dist, list(list(col="NDA")))
df = mergeOnMatches(t.df, s.df, t.pos, s.pos, dist, list(list(col="X")))
df = mergeOnMatches(t.df, s.df, t.pos, s.pos, dist, list(list(col="X", before="X")))
df = mergeOnMatches(t.df, s.df, t.pos, s.pos, dist, list(list(col="X", before="")))
df = mergeOnMatches(t.df, s.df, t.pos, s.pos, dist, list(list(col="X", before="", format=NA)))
df = mergeOnMatches(t.df, s.df, t.pos, s.pos, dist, list(list(col="X", before="", format="{junk}")))
df = mergeOnMatches(t.df, s.df, t.pos, s.pos, dist, list(list(col="X", before="", format="X", maxMatch="A")))
df = mergeOnMatches(t.df, s.df, t.pos, s.pos, dist, list(list(col="X", before="", format="X", maxMatch=NA)))
df = mergeOnMatches(t.df, s.df, t.pos, s.pos, dist, list(list(col="X", before="", format="X", maxMatch=-1)))
df = mergeOnMatches(t.df, s.df, t.pos, s.pos, dist, list(list(col="X", before="", format="X", joinStart=t.df)))
df = mergeOnMatches(t.df, s.df, t.pos, s.pos, dist, list(list(col="X", before="", format="X", joinEnd=t.df)))

#############################################
# Test using a GFF3 gene file with various dist$methods.
#############################################

# Read .gff3 file of genes in same chromosomes as the sequence used to create the markers in dfMarkers.
dfGff = readFile_GFF3("test_GFFfuncsAndMergeData/test.gff3")
head(dfGff)
dfGff = clean_GFF3(dfGff)
dfG = selectFeatures(dfGff, "gene")
dfG = convertAttrsToCols(dfG, includeAttrs="Alias")
dfG = dfG[, colnames(dfG) %in% c("Alias", "seqname", "start", "end", "strand")]
colnames(dfG) = c("id", "start", "end", "strand", "gene")
dfG = dfG[, c("gene", "id", "start", "end", "strand")]

# Now annotate markers with genes using different methods.

# method = OVERLAP
t.df = dfMarkers
s.df = dfG
t.pos = list(id="Hid", start="HampPos1", end="HampPos2")
s.pos = list(id="id", start="start", end="end")
dist = list(method="OVERLAP")
cols = list(list(col="genes", before="prmSeqL", maxMatch=0, join="YES", joinSep=";", format="{+gene}({#s})"))
df = mergeOnMatches(t.df, s.df, t.pos, s.pos, dist, cols)
sum(df$genes != "")
df[1:10,1:11]
df[1500:1510,1:11]
df[1600:1610,1:11]
tail(df[,1:11], 20)
df[1,]
head(dfG)
# Correct

# Try #t in place of #s
dist = list(method="OVERLAP")
cols = list(list(col="genes", before="prmSeqL", maxMatch=0, join="YES", joinSep=";", format="{+gene}({#t})"))
df = mergeOnMatches(t.df, s.df, t.pos, s.pos, dist, cols)
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

# method = s.TINY : amplicon must contain entire gene
dist = list(method="s.TINY")
df = mergeOnMatches(t.df, s.df, t.pos, s.pos, dist, cols)
sum(df$genes != "") # Not very many (8), sounds good.
df[df$genes != "",1:11]
df[1488,]
dfG[dfG$gene == "Solyc01g009550",]
# Correct

# method = t.TINY : gene must contain entire amplicon
dist = list(method="t.TINY")
df = mergeOnMatches(t.df, s.df, t.pos, s.pos, dist, cols)
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
# method = s.NEAR, closest = 0, distance 100bp : gene must be within 100 bp of amplicon (or contain it)
#############################################

#               start.up: G.start must be no less than M.start-twice max gene start-start distance
#               start.down: G.start must be no more than M.end+100 bp
#               end.up: G.end must be no less than M.start-100 bp
#               end.down: G.end must be no more than M.end+twice max gene end-end distance

dist = list(method="s.NEAR", closest=0, start.down=100, end.up=100)
df = mergeOnMatches(t.df, s.df, t.pos, s.pos, dist, cols)
sum(df$genes != "") # 1671
dim(df.overlap) # 1647
df.s.near = df[df$genes != "",]
# Some of the markers overlapping genes should now also be near another gene, so
# the marker's gene list should have an additional gene in it.  Let's use HampPos1_HampPos2
# as row names for df.overlap and df.s.near, which will help us.
rownames(df.overlap) = paste(df.overlap$HampPos1, df.overlap$HampPos2, sep="_")
rownames(df.s.near) = paste(df.s.near$HampPos1, df.s.near$HampPos2, sep="_")
# Get annotated markers common between OVERLAP and s.NEAR, and those specific to each.
common = intersect(rownames(df.overlap), rownames(df.s.near))
only.overlap = setdiff(rownames(df.overlap), rownames(df.s.near))
only.s.near = setdiff(rownames(df.s.near), rownames(df.overlap))
length(common) # 1647
length(only.overlap) # 0, good, it should be that way.
length(only.s.near) # 24
# For the common ones, get the ones whose gene annotation mismatches.
common.gene.mismatch = common[df.overlap[common, "genes"] != df.s.near[common, "genes"]]
length(common.gene.mismatch) # Only 23 genes 
data.frame(overlap=df.overlap[common.gene.mismatch, "genes"], s.near=df.s.near[common.gene.mismatch, "genes"],
    stringsAsFactors=FALSE)
# Looks good, the s.NEAR have the same gene as OVERLAP plus one more gene next to it.
df.overlap[common.gene.mismatch[1],]
df.s.near[common.gene.mismatch[1],]
dfG[dfG$gene %in% c("Solyc01g005070", "Solyc01g005080"),]
# Note: distance to gene, in (), is to START of gene, even though END of gene is closer.
# It looks good.
# Check an only.s.near gene.
df.s.near[only.s.near[1],]
dfG[dfG$gene == "Solyc01g005560",]
# Looks good.

# Expand the near-distance to 3000 and we should get a lot more, including a bunch with
# multiple genes per marker.
t.df = dfMarkers
s.df = dfG
t.pos = list(id="Hid", start="HampPos1", end="HampPos2")
s.pos = list(id="id", start="start", end="end")
dist = list(method="s.NEAR", closest=0, start.down=3000, end.up=3000)
cols = list(list(col="genes", before="prmSeqL", maxMatch=0, join="YES", joinSep=";", format="{+gene}({#t})"))
df = mergeOnMatches(t.df, s.df, t.pos, s.pos, dist, cols)
sum(df$genes != "") # 1815 vs 1671
df.s.near = df[df$genes != "",]
rownames(df.s.near) = paste(df.s.near$HampPos1, df.s.near$HampPos2, sep="_")
common = intersect(rownames(df.overlap), rownames(df.s.near))
only.overlap = setdiff(rownames(df.overlap), rownames(df.s.near))
only.s.near = setdiff(rownames(df.s.near), rownames(df.overlap))
length(common) # 1647
length(only.overlap) # 0, good, it should be that way.
length(only.s.near) # 168
common.gene.mismatch = common[df.overlap[common, "genes"] != df.s.near[common, "genes"]]
length(common.gene.mismatch) # Now 685!
data.frame(name=common.gene.mismatch,
    overlap=df.overlap[common.gene.mismatch, "genes"], s.near=df.s.near[common.gene.mismatch, "genes"],
    stringsAsFactors=FALSE)
# Looks good.  One with 3 genes is 4683900_4684979
df.s.near[only.s.near[1], "genes"]
only.s.near[1]
# 54857_55150 has 3 genes (30/40/50), none of them overlapping.

# Repeat but use closest=1.  We should get only the overlapping genes in the common
# set, and in the unique set we should get only the closest one.
dist = list(method="s.NEAR", closest=1, start.down=3000, end.up=3000)
df = mergeOnMatches(t.df, s.df, t.pos, s.pos, dist, cols)
sum(df$genes != "") # Still 1815
df.s.near = df[df$genes != "",]
rownames(df.s.near) = paste(df.s.near$HampPos1, df.s.near$HampPos2, sep="_")
common = intersect(rownames(df.overlap), rownames(df.s.near))
only.overlap = setdiff(rownames(df.overlap), rownames(df.s.near))
only.s.near = setdiff(rownames(df.s.near), rownames(df.overlap))
length(common) # 1647
length(only.overlap) # 0, good, it should be that way.
length(only.s.near) # 168
common.gene.mismatch = common[df.overlap[common, "genes"] != df.s.near[common, "genes"]]
length(common.gene.mismatch) # Now 0 because if it has an overlap we don't include ANY nearby ones
# We should only see genes 40 and 50 here, because 60 only showed up in the overlap set
df.s.near["4683900_4684979","genes"]
# We should see only nearest one of 30, 40, 50, which is 30 (I looked).
df.s.near["54857_55150","genes"]

# Repeat but use closest=2.  We should get only the overlapping genes in the common
# set, and in the unique set we should get only the closest UPSTREAM *AND* closest
# DOWNSTREAM one.
dist = list(method="s.NEAR", closest=2, start.down=3000, end.up=3000)
df = mergeOnMatches(t.df, s.df, t.pos, s.pos, dist, cols)
sum(df$genes != "") # Still 1815
df.s.near = df[df$genes != "",]
rownames(df.s.near) = paste(df.s.near$HampPos1, df.s.near$HampPos2, sep="_")
common = intersect(rownames(df.overlap), rownames(df.s.near))
only.overlap = setdiff(rownames(df.overlap), rownames(df.s.near))
only.s.near = setdiff(rownames(df.s.near), rownames(df.overlap))
length(common) # 1647
length(only.overlap) # 0, good, it should be that way.
length(only.s.near) # 168
common.gene.mismatch = common[df.overlap[common, "genes"] != df.s.near[common, "genes"]]
length(common.gene.mismatch) # Still 0
# We should only see genes 40 and 50 here, because 60 only showed up in the overlap set
df.s.near["4683900_4684979","genes"]
# We should see only the nearest upstream gene, which is 30, and nearest downstream
# gene, which is 40 (I looked).
df.s.near["54857_55150","genes"]

#############################################
# method = t.NEAR, closest = 0, distance 100bp : amplicon must be within 100 bp of gene (or contain it)
#############################################

# Now reverse it and use t.NEAR, so we are looking for markers that are near genes.
# It should actually generate the same results as with s.NEAR., right?
dist = list(method="t.NEAR", closest=0, start.down=100, end.up=100)
df = mergeOnMatches(t.df, s.df, t.pos, s.pos, dist, cols)
sum(df$genes != "") # 1671
dim(df.overlap) # 1647
df.t.near = df[df$genes != "",]
# Some of the genes overlapping markers should now also be near another marker, so
# the gene list should have an additional gene in it.  Let's use HampPos1_HampPos2
# as row names for df.overlap and df.t.near, which will help us.
rownames(df.overlap) = paste(df.overlap$HampPos1, df.overlap$HampPos2, sep="_")
rownames(df.t.near) = paste(df.t.near$HampPos1, df.t.near$HampPos2, sep="_")
# Get annotated markers common between OVERLAP and t.NEAR, and those specific to each.
common = intersect(rownames(df.overlap), rownames(df.t.near))
only.overlap = setdiff(rownames(df.overlap), rownames(df.t.near))
only.t.near = setdiff(rownames(df.t.near), rownames(df.overlap))
length(common) # 1647
length(only.overlap) # 0, good, it should be that way.
length(only.t.near) # 24
# For the common ones, get the ones whose gene annotation mismatches.
common.gene.mismatch = common[df.overlap[common, "genes"] != df.t.near[common, "genes"]]
length(common.gene.mismatch) # Only 23 genes 
data.frame(overlap=df.overlap[common.gene.mismatch, "genes"], s.near=df.t.near[common.gene.mismatch, "genes"],
    stringsAsFactors=FALSE)
# Looks good, the t.NEAR have the same gene as OVERLAP plus one more gene next to it.
df.overlap[common.gene.mismatch[1],]
df.t.near[common.gene.mismatch[1],]
dfG[dfG$gene %in% c("Solyc01g005070", "Solyc01g005080"),]
# Note: distance to gene, in (), is to START of gene, even though END of gene is closer.
# It looks good.
# Check an only.t.near gene.
df.t.near[only.t.near[1],]
dfG[dfG$gene == "Solyc01g005560",]
# Looks good.

# Expand the near-distance to 3000 and we should get a lot more, including a bunch with
# multiple genes per marker.
t.df = dfMarkers
s.df = dfG
t.pos = list(id="Hid", start="HampPos1", end="HampPos2")
s.pos = list(id="id", start="start", end="end")
dist = list(method="t.NEAR", closest=0, start.down=3000, end.up=3000)
cols = list(list(col="genes", before="prmSeqL", maxMatch=0, join="YES", joinSep=";", format="{+gene}({#t})"))
df = mergeOnMatches(t.df, s.df, t.pos, s.pos, dist, cols)
sum(df$genes != "") # 1815 vs 1671
df.t.near = df[df$genes != "",]
rownames(df.t.near) = paste(df.t.near$HampPos1, df.t.near$HampPos2, sep="_")
common = intersect(rownames(df.overlap), rownames(df.t.near))
only.overlap = setdiff(rownames(df.overlap), rownames(df.t.near))
only.t.near = setdiff(rownames(df.t.near), rownames(df.overlap))
length(common) # 1647
length(only.overlap) # 0, good, it should be that way.
length(only.t.near) # 168
common.gene.mismatch = common[df.overlap[common, "genes"] != df.t.near[common, "genes"]]
length(common.gene.mismatch) # Now 685!
data.frame(name=common.gene.mismatch,
    overlap=df.overlap[common.gene.mismatch, "genes"], s.near=df.t.near[common.gene.mismatch, "genes"],
    stringsAsFactors=FALSE)
# Looks good.  One with 3 genes is 4683900_4684979
df.t.near[only.t.near[1], "genes"]
only.t.near[1]
# 54857_55150 has 3 genes (30/40/50), none of them overlapping.

# Repeat but use closest=1.  We should get only the overlapping genes in the common
# set, and in the unique set we should get only the closest marker to each gene,
# which will give a different result than s.NEAR which was closest gene to each
# marker.  With s.NEAR closest=1 simply reduced the number of genes for each marker,
# but each marker that had genes still had genes.  With t.NEAR, some markers that
# picked up their nearest gene now have that gene picking its closest marker, which
# is a different marker.  So, we should end up with fewer markers with genes.
dist = list(method="t.NEAR", closest=1, start.down=3000, end.up=3000)
Solyc01g005570
df = mergeOnMatches(t.df, s.df, t.pos, s.pos, dist, cols)
sum(df$genes != "") # Now 1691 instead of 1815.  Fewer, as expected.
df.t.near = df[df$genes != "",]
rownames(df.t.near) = paste(df.t.near$HampPos1, df.t.near$HampPos2, sep="_")
common = intersect(rownames(df.overlap), rownames(df.t.near))
only.overlap = setdiff(rownames(df.overlap), rownames(df.t.near))
only.t.near = setdiff(rownames(df.t.near), rownames(df.overlap))
length(common) # still 1647
length(only.overlap) # 0, good, it should be that way.
length(only.t.near) # 44 instead of 168, fewer, as expected.
common.gene.mismatch = common[df.overlap[common, "genes"] != df.t.near[common, "genes"]]
length(common.gene.mismatch) # Now 57 instead of 0.  If a gene overlaps a marker, it gets just that marker, but if not, it gets nearest marker, and so one marker may get several nearest genes even though it overlaps another marker.
# We should only see genes 40 and 50 here, because 60 only showed up in the overlap set
df.t.near["4683900_4684979","genes"]
# Is this marker the closest marker to which genes?  Apparently to 40 and 50.
df.t.near["54857_55150","genes"]

# Repeat but use closest=2.  We should get only the overlapping genes in the common
# set, and in the unique set we should get MORE because now we have both upstream
# and downstream.
dist = list(method="t.NEAR", closest=2, start.down=3000, end.up=3000)
df = mergeOnMatches(t.df, s.df, t.pos, s.pos, dist, cols)
sum(df$genes != "") # 1692 now, we picked up one more, would have thought more than that.
df.t.near = df[df$genes != "",]
rownames(df.t.near) = paste(df.t.near$HampPos1, df.t.near$HampPos2, sep="_")
common = intersect(rownames(df.overlap), rownames(df.t.near))
only.overlap = setdiff(rownames(df.overlap), rownames(df.t.near))
only.t.near = setdiff(rownames(df.t.near), rownames(df.overlap))
length(common) # still 1647
length(only.overlap) # 0, good, it should be that way.
length(only.t.near) # 45 now instead of 4, we gained one.
common.gene.mismatch = common[df.overlap[common, "genes"] != df.t.near[common, "genes"]]
length(common.gene.mismatch) # Now 62 instead of 57, gained 5.
# We should only see genes 40 and 50 here, because 60 only showed up in the overlap set
df.t.near["4683900_4684979","genes"]
# Is this marker the closest marker to which genes?  Apparently to 40 and 50.
df.t.near["54857_55150","genes"]



# Arguments:
#   t.df: target data frame to which columns will be added and containing positions.
#   s.df: source data frame containing data to be added to t.df and containing positions.
#   t.pos: list giving column names in t.df that define a position.  Members:
#           "start" (required) : start position or main position.
#           "id", "len", "end" (optional) : sequence ID, length, and end position in bp.
#   s.pos: like t.pos, for s.df.
#   dist: list defining how to find position matches between t.df and s.df.  Members:
#           "method" (required) : one of "OVERLAP", "s.TINY", "t.TINY", "t.NEAR", "t.NEAR"
#           "closest" (optional) : for "s/t.NEAR", 0 (all), 1 (nearest 1), or 2 (nearest 2)
#           "start.up", "start.down", "end.up", "end.down" (optional) : for "x.NEAR" (x=s/t):
#                   start.up: x.start must be no less than y.start-start.up
#                   start.down: x.start must be no more than y.end+start.down
#                   end.up: x.end must be no less than y.start-end.up
#                   end.down: x.end must be no more than y.end+end.down
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
