#######################################################################################
# This file contains R definitions and functions for working with .gff3 and .gtf files.
# Author: Ted Toal
# Date: 2015
# Brady Lab, UC Davis
#######################################################################################

#######################################################################################
# Read a .gff3 file.
#
# Arguments:
#   filename: Name of .gff3 file to read.
#
# Returns: data frame containing the data, with columns "seqname", "source", "feature",
#   "start", "end", "score", "strand", "frame", "attributes".
#######################################################################################
readFile_GFF3 = function(filename)
    {
    df = read.table(filename, header=FALSE, sep="\t", quote='"', comment.char="#", stringsAsFactors=FALSE)
    colnames(df) = c("seqname", "source", "feature", "start", "end", "score", "strand", "frame", "attributes")
    return(df)
    }

#######################################################################################
# Read a .gtf file.
#
# Arguments:
#   filename: Name of .gtf file to read.
#
# Returns: data frame containing the data, with columns "seqname", "source", "feature",
#   "start", "end", "score", "strand", "frame", "attributes".
#######################################################################################
readFile_GTF = function(filename)
    {
    df = read.table(filename, header=FALSE, sep="\t", quote="", comment.char="#", stringsAsFactors=FALSE)
    colnames(df) = c("seqname", "source", "feature", "start", "end", "score", "strand", "frame", "attributes")
    return(df)
    }

#######################################################################################
# Write a .gff3 or .gtf file.
#
# Arguments:
#   df: data frame to write, containing either .gff3 or .gtf data.
#   filename: Name of file to write.
#
# Returns: Nothing.
#######################################################################################
writeFile_GFF3_GTF = function(df, filename)
    {
    write.table(df, file=filename, sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)
    }

#######################################################################################
# Clean up .gff3 data.  If data is .gtf, it is more .gff3 after this.
#
# Arguments:
#   df: data frame of .gff3 or .gtf data.
#
# Returns: data frame of cleaned-up .gff3 or .gtf data.
#
# Changes made:
#   1. Attributes column spaces around ";" characters are removed.
#   2. Attributes column multiple consecutive ";" characters are changed to a single ";".
#   3. Attributes column ";" at start or end is removed.
#   4. Attributes column uses an "=" character rather than a space character.
#   5. Attributes column does not use double quotes around text attribute values.
#   6. Feature name "5UTR" is changed to "five_prime_UTR"
#   7. Feature name "3UTR" is changed to "three_prime_UTR"
#
# Certain changes are not made here but need to be made as part of the GTF format:
#   1. Change CDS to include stop codon bases.
#   2. Do not allow start and stop codon features to be split by introns.
#######################################################################################
clean_GFF3 = function(df)
    {
    df$feature[df$feature == "5UTR"] = "five_prime_UTR"
    df$feature[df$feature == "3UTR"] = "three_prime_UTR"
    if (!any(colnames(df) == "attributes"))
        stop("clean_GFF3: data frame has no 'attributes' column")
    df$attributes = gsub(" +;", ";", df$attributes)
    df$attributes = gsub("; +", ";", df$attributes)
    df$attributes = gsub(";;+", ";", df$attributes)
    df$attributes = gsub("^;", "", df$attributes)
    df$attributes = gsub(";$", "", df$attributes)
    df$attributes = gsub('(^|;)([a-zA-Z0-9_]+)[= ]?("?)([^;"]+)("?)', '\\1\\2=\\4', df$attributes)
    rownames(df) = NULL
    return(df)
    }

#######################################################################################
# Clean up .gtf data.  If data is .gff3, it is more .gtf after this.
#
# Arguments:
#   df: data frame of .gtf or .gff3 data.
#
# Returns: data frame of cleaned-up .gtf or .gff3 data.
#
# Changes made:
#   1. Attributes column spaces before ";" characters are removed.
#   2. Attributes column multiple consecutive space characters after ";" characters
#       are changed to a single space.
#   3. Attributes column multiple consecutive ";" characters are changed to a single ";".
#   4. Attributes column ";" at start or end is removed.
#   5. Attributes column single space is added after ";" characters with no following space.
#   6. Attributes column uses a space character rather than "=".
#   7. Attributes column uses double quotes around text attribute values.
#   8. Feature name "five_prime_UTR" is changed to "5UTR"
#   9. Feature name "three_prime_UTR" is changed to "3UTR"
#
# Certain changes are not made here but need to be made as part of the GFF3 format:
#   1. Add frame field values of 0,1,2 for start/stop codons.
#   2. Change CDS to not include stop codon bases.
#   3. Allow start and stop codon features to be split by introns.
#######################################################################################
clean_GTF = function(df)
    {
    df = clean_GFF3(df) # Go the other way first.
    df$feature[df$feature == "five_prime_UTR"] = "5UTR"
    df$feature[df$feature == "three_prime_UTR"] = "3UTR"
    if (!any(colnames(df) == "attributes"))
        stop("clean_GTF: data frame has no 'attributes' column")
    df$attributes = gsub(";", "; ", df$attributes)
    df$attributes = gsub('(^|; )([a-zA-Z0-9_]+)[= ]?"?([^;"]+)"?', '\\1\\2 "\\3"', df$attributes)
    df$attributes = gsub('(^|; )([a-zA-Z0-9_]+) "([0-9.+\\-]+)"', '\\1\\2 \\3', df$attributes)
    rownames(df) = NULL
    return(df)
    }

#######################################################################################
# Return a vector of the unique feature names in a .gff3 or .gtf data frame.
#
# Arguments:
#   df: data frame of .gff3 or .gtf data.
#
# Returns: vector of features.
#######################################################################################
getFeatures = function(df)
    {
    return(unique(df$feature))
    }

#######################################################################################
# Remove .gff3 or .gtf data from a data frame that does not belong to one of a specified
# list of features.
#
# Arguments:
#   df: data frame of .gff3 or .gtf data.
#   selectFeatures: vector of names of features to be retained, others are removed.
#
# Returns: modified data frame.
#######################################################################################
selectFeatures = function(df, selectFeatures)
    {
    df = df[df$feature %in% selectFeatures,]
    rownames(df) = NULL
    return(df)
    }

#######################################################################################
# Convert .gff3 or .gtf attribute data to separate data frame columns.
#
# Arguments:
#   df: data frame of .gff3 or .gtf data, cleaned with clean_GFF3() or clean_GTF().
#   includeAttrs: names of attributes to be made into columns with that name, NULL means
#       to make columns out of every attribute name that is found.
#   excludeAttrs: names of attributes to NOT be made into columns even if in includeAttrs,
#       NULL means to exclude none.
#   removeAttrCol: TRUE to remove the attrCol column from the data frame.  If not removed,
#       the attributes column no longer contains the data moved to separate columns.
#   attrCol: name of attributes column, normally "attributes".
#   newAttrCols: names of new data columns, or NULL to name the columns with the attribute
#       name.  Must be same length and same order as number of new columns to be added.
#   missingAttrVals: values to use in rows lacking the attribute, or NULL to use NA.
#       Must be same length and same order as number of new columns to be added or length 1 (in which
#       case that one value is used for missing attributes in ALL added columns).
#   maxAttrCols: maximum number of attribute data columns ever to be added.
#
# Returns: modified data frame.
#
# Note: The new columns have NA for any row whose attributes did not include that column.
#######################################################################################
convertAttrsToCols = function(df, includeAttrs=NULL, excludeAttrs=NULL, removeAttrCol=TRUE,
    attrCol="attributes", newAttrCols=NULL, missingAttrVals=NULL, maxAttrCols=100)
    {
    if (!any(colnames(df) == attrCol))
        stop("convertAttrsToCols: no such column: ", attrCol)

    N1 = sum(grepl("; ", df[[attrCol]], fixed=TRUE))
    N2 = sum(grepl(";", df[[attrCol]], fixed=TRUE))
    isGTF = (N1 > N2/2)
    if (isGTF)
        df = clean_GFF3(df)

    attrs = strsplit(df[[attrCol]], ";", fixed=TRUE)
    names(attrs) = NULL
    attrNamesVals = sapply(attrs, function(V)
        {
        if (length(V) > 0)
            {
            L = strsplit(V, "=", fixed=TRUE)
            V = sapply(L, "[", 2)
            names(V) = sapply(L, "[", 1)
            }
        return(V)
        }, simplify=FALSE)
    attrNames = unique(names(unlist(attrNamesVals)))
    if (length(attrNames) > maxAttrCols+length(includeAttrs)+length(excludeAttrs))
        stop("convertAttrsToCols: There are ", length(attrNames), " attribute columns to be added, too many")

    if (is.null(includeAttrs))
        includeAttrs = attrNames
    if (!is.null(excludeAttrs))
        includeAttrs = setdiff(includeAttrs, excludeAttrs)
    N = length(includeAttrs)

    if (is.null(newAttrCols))
        newAttrCols = includeAttrs
    else
        {
        if (length(newAttrCols) != N)
            stop("convertAttrsToCols: there are ", N,
                "columns to be included but newAttrCols has only ", length(newAttrCols))
        }
    names(newAttrCols) = includeAttrs

    if (is.null(missingAttrVals))
        missingAttrVals = NA
    if (length(missingAttrVals) == 1)
        missingAttrVals = rep(missingAttrVals, N)
    else
        {
        if (length(missingAttrVals) != N)
            stop("convertAttrsToCols: there are ", N,
                "columns to be included but missingAttrVals has only ", length(missingAttrVals))
        }
    names(missingAttrVals) = includeAttrs

    includeAttrs = intersect(includeAttrs, attrNames)
    newAttrCols = newAttrCols[includeAttrs]
    missingAttrVals = missingAttrVals[includeAttrs]

    rownames(df) = NULL
    for (attr in includeAttrs)
        {
        df[, newAttrCols[attr]] = sapply(attrNamesVals, "[", attr)
        isNA = is.na(df[, newAttrCols[attr]])
        df[isNA, newAttrCols[attr]] = missingAttrVals[attr]
        }

    if (removeAttrCol)
        df = df[, colnames(df) != attrCol, drop=FALSE]
    else
        {
        attrNamesVals = sapply(attrNamesVals, function(V) return(V[!names(V) %in% includeAttrs]), simplify=FALSE)
        df[[attrCol]] = sapply(attrNamesVals, function(V) return(paste(names(V), V, sep="=", collapse=";")))
        if (isGTF)
            df = clean_GTF(df)
        }
    rownames(df) = NULL
    return(df)
    }

#######################################################################################
# Convert separate data frame columns to attributes in the attributes column of a .gff3
# or .gtf data frame.
#
# Arguments:
#   df: data frame of .gff3 or .gtf data with extra columns to merge into attributes.
#   cols: names of df columns to be made into attributes.  The column name is the
#       attribute name.
#   newAttrNames: names to give to the new attributes, or NULL to use 'cols' as names.
#       Must be same length and same order as 'cols'.
#   noAttrValues: values which, when present in a row, indicate that the attribute of
#       that column should not be added at the attributes list for that row.  If NULL,
#       an NA value means to not add an attribute to the row.  Must be same length and
#       same order as 'cols', or length 1, which means that that value applies to all
#       'cols'.
#   attrCol: name of attributes column, normally "attributes".
#   merge: TRUE to merge the data in the "cols" columns into the existing attributes,
#       if any, and FALSE to replace the existing attributes with the new data.
#   remove: TRUE to remove the "cols" columns from the data frame.
#
# Returns: modified data frame in .gff3 format.
#
# Note: When a column in "cols" contains NA, that attribute is not added to the
# attributes column, and if merge is TRUE, the value of the attribute already in the
# attributes column remains unchanged.
#######################################################################################
convertColsToAttrs = function(df, cols, newAttrNames=NULL, noAttrValues=NULL,
    attrCol="attributes", merge=TRUE, remove=TRUE)
    {
    if (is.null(cols) || length(cols) == 0)
        stop("convertColsToAttrs: no column names specified in 'cols' argument")
    N = length(cols)

    if (is.null(newAttrNames))
        newAttrNames = cols
    if (length(newAttrNames) != N)
        stop("convertColsToAttrs: there are ", N, " columns to convert but newAttrNames has only ", N)
    names(newAttrNames) = cols

    if (is.null(noAttrValues))
        noAttrValues = NA
    if (length(noAttrValues) == 1)
        noAttrValues = rep(noAttrValues, N)
    else
        {
        if (length(noAttrValues) != N)
            stop("convertColsToAttrs: there are ", N,
                "columns to be converted but noAttrValues has only ", length(noAttrValues))
        }
    names(noAttrValues) = cols

    if (remove || !any(colnames(df) == attrCol))
        df[[attrCol]] = ""
    else
        df = clean_GFF3(df)

    attrs = strsplit(df[[attrCol]], ";", fixed=TRUE)
    names(attrs) = NULL
    attrNamesVals = sapply(attrs, function(V)
        {
        if (length(V) > 0)
            {
            L = strsplit(V, "=", fixed=TRUE)
            V = sapply(L, "[", 2)
            names(V) = sapply(L, "[", 1)
            }
        return(V)
        }, simplify=FALSE)

    for (col in cols)
        {
        attrName = newAttrNames[col]
        noAttrValue = noAttrValues[col]
        d = df[,col]
        if (is.na(noAttrValue))
            idxs = which(!is.na(d))
        else
            idxs = which(d != noAttrValue)
        attrNamesVals[idxs] = sapply(idxs, function(i)
            {
            V = attrNamesVals[[i]]
            V[attrName] = d[i]
            return(V)
            }, simplify=FALSE)
        }

    df[[attrCol]] = sapply(attrNamesVals, function(V) return(paste(names(V), V, sep="=", collapse=";")))

    if (remove)
        df = df[, !colnames(df) %in% cols, drop=FALSE]
    rownames(df) = NULL
    return(df)
    }

#######################################################################################
#######################################################################################
# Example data from some .gff3 files for reference.
#######################################################################################
#######################################################################################

# Arabidopsis TAIR10 gene model file example:
#Chr3	TAIR10	gene	3026055	3027234	.	+	.	ID=AT3G09860;Note=protein_coding_gene;Name=AT3G09860
#Chr3	TAIR10	mRNA	3026055	3027234	.	+	.	ID=AT3G09860.1;Parent=AT3G09860;Name=AT3G09860.1;Index=1
#Chr3	TAIR10	protein	3026140	3027000	.	+	.	ID=AT3G09860.1-Protein;Name=AT3G09860.1;Derives_from=AT3G09860.1
#Chr3	TAIR10	exon	3026055	3026291	.	+	.	Parent=AT3G09860.1
#Chr3	TAIR10	five_prime_UTR	3026055	3026139	.	+	.	Parent=AT3G09860.1
#Chr3	TAIR10	CDS	3026140	3026291	.	+	0	Parent=AT3G09860.1,AT3G09860.1-Protein;
#Chr3	TAIR10	exon	3026759	3026837	.	+	.	Parent=AT3G09860.1
#Chr3	TAIR10	CDS	3026759	3026837	.	+	1	Parent=AT3G09860.1,AT3G09860.1-Protein;
#Chr3	TAIR10	exon	3026935	3027234	.	+	.	Parent=AT3G09860.1
#Chr3	TAIR10	CDS	3026935	3027000	.	+	0	Parent=AT3G09860.1,AT3G09860.1-Protein;
#Chr3	TAIR10	three_prime_UTR	3027001	3027234	.	+	.	Parent=AT3G09860.1

# Heinz ITAG2.3 gene model file example:
#SL2.40ch07	ITAG_eugene	gene	3314726	3318115	.	-	.	Alias=Solyc07g008440;ID=gene:Solyc07g008440.2;Name=Solyc07g008440.2;from_BOGAS=1;length=3390
#SL2.40ch07	ITAG_eugene	mRNA	3314726	3318115	.	-	.	ID=mRNA:Solyc07g008440.2.1;Name=Solyc07g008440.2.1;Note=Purine permease family protein (AHRD V1 **-* D7MCV2_ARALY)%3B contains Interpro domain(s)  IPR004853  Protein of unknown function DUF250 ;Ontology_term=GO:0005345;Parent=gene:Solyc07g008440.2;from_BOGAS=1;interpro2go_term=GO:0016020;length=3390;nb_exon=2;sifter_term=GO:0005345
#SL2.40ch07	ITAG_eugene	exon	3314726	3315117	.	-	.	ID=exon:Solyc07g008440.2.1.2;Parent=mRNA:Solyc07g008440.2.1;from_BOGAS=1
#SL2.40ch07	ITAG_eugene	three_prime_UTR	3314726	3314766	.	-	.	ID=three_prime_UTR:Solyc07g008440.2.1.0;Parent=mRNA:Solyc07g008440.2.1;from_BOGAS=1
#SL2.40ch07	ITAG_eugene	CDS	3314767	3315117	.	-	0	ID=CDS:Solyc07g008440.2.1.2;Parent=mRNA:Solyc07g008440.2.1;from_BOGAS=1
#SL2.40ch07	ITAG_eugene	intron	3315118	3317308	.	-	.	ID=intron:Solyc07g008440.2.1.2;Parent=mRNA:Solyc07g008440.2.1;from_BOGAS=1
#SL2.40ch07	ITAG_eugene	exon	3317309	3318115	.	-	0	ID=exon:Solyc07g008440.2.1.1;Parent=mRNA:Solyc07g008440.2.1;from_BOGAS=1
#SL2.40ch07	ITAG_eugene	CDS	3317309	3318061	.	-	0	ID=CDS:Solyc07g008440.2.1.1;Parent=mRNA:Solyc07g008440.2.1;from_BOGAS=1
#SL2.40ch07	ITAG_eugene	five_prime_UTR	3318062	3318115	.	-	.	ID=five_prime_UTR:Solyc07g008440.2.1.0;Parent=mRNA:Solyc07g008440.2.1;from_BOGAS=1

# Heinz ITAG2.4 gene model file example:
#SL2.50ch07	ITAG_eugene	gene	3314726	3318115	.	-	.	ID=gene:Solyc07g008440.2;Name=Solyc07g008440.2;Alias=Solyc07g008440;from_BOGAS=1;length=3390
#SL2.50ch07	ITAG_eugene	mRNA	3314726	3318115	.	-	.	ID=mRNA:Solyc07g008440.2.1;Name=Solyc07g008440.2.1;Parent=gene:Solyc07g008440.2;Note=Purine permease family protein (AHRD V1 **-* D7MCV2_ARALY)%3B contains Interpro domain(s)  IPR004853  Protein of unknown function DUF250 ;Ontology_term=GO:0005345;from_BOGAS=1;interpro2go_term=GO:0016020;length=3390;nb_exon=2;sifter_term=GO:0005345
#SL2.50ch07	ITAG_eugene	exon	3314726	3315117	.	-	.	ID=exon:Solyc07g008440.2.1.2;Parent=mRNA:Solyc07g008440.2.1;from_BOGAS=1
#SL2.50ch07	ITAG_eugene	three_prime_UTR	3314726	3314766	.	-	.	ID=three_prime_UTR:Solyc07g008440.2.1.0;Parent=mRNA:Solyc07g008440.2.1;from_BOGAS=1
#SL2.50ch07	ITAG_eugene	CDS	3314767	3315117	.	-	0	ID=CDS:Solyc07g008440.2.1.2;Parent=mRNA:Solyc07g008440.2.1;from_BOGAS=1
#SL2.50ch07	ITAG_eugene	intron	3315118	3317308	.	-	.	ID=intron:Solyc07g008440.2.1.2;Parent=mRNA:Solyc07g008440.2.1;from_BOGAS=1
#SL2.50ch07	ITAG_eugene	exon	3317309	3318115	.	-	0	ID=exon:Solyc07g008440.2.1.1;Parent=mRNA:Solyc07g008440.2.1;from_BOGAS=1
#SL2.50ch07	ITAG_eugene	CDS	3317309	3318061	.	-	0	ID=CDS:Solyc07g008440.2.1.1;Parent=mRNA:Solyc07g008440.2.1;from_BOGAS=1
#SL2.50ch07	ITAG_eugene	five_prime_UTR	3318062	3318115	.	-	.	ID=five_prime_UTR:Solyc07g008440.2.1.0;Parent=mRNA:Solyc07g008440.2.1;from_BOGAS=1

# Penn V2 scaffolds gene model file example:
#scaffold680.1	AUGUSTUS	transcript	14931	21707	0.03	-	.	Name=scaffold680.1.g168.t1;ID=462805;Note=%20unknown%20:%20%20no%20hits%20%20(original%20descriptio..
#scaffold680.1	AUGUSTUS	transcription_end_site	14931	14931	.	-	.	Parent=462805;ID=462806
#scaffold680.1	AUGUSTUS	exon	14931	15169	.	-	.	Parent=462805;ID=462807
#scaffold680.1	AUGUSTUS	stop_codon	14978	14980	.	-	0	Parent=462805;ID=462808
#scaffold680.1	AUGUSTUS	intron	15170	16937	0.44	-	.	Parent=462805;ID=462809
#scaffold680.1	AUGUSTUS	intron	17202	17359	0.59	-	.	Parent=462805;ID=462810
#scaffold680.1	AUGUSTUS	intron	17613	18187	0.31	-	.	Parent=462805;ID=462811
#scaffold680.1	AUGUSTUS	intron	18276	18364	0.37	-	.	Parent=462805;ID=462812
#scaffold680.1	AUGUSTUS	CDS	14978	15169	1	-	0	Parent=462805;ID=462813
#scaffold680.1	AUGUSTUS	CDS	16938	17201	0.3	-	0	Parent=462805;ID=462813
#scaffold680.1	AUGUSTUS	CDS	17360	17612	0.63	-	1	Parent=462805;ID=462813
#scaffold680.1	AUGUSTUS	CDS	18188	18275	0.32	-	2	Parent=462805;ID=462813
#scaffold680.1	AUGUSTUS	CDS	18365	18842	0.5	-	0	Parent=462805;ID=462813
#scaffold680.1	AUGUSTUS	exon	16938	17201	.	-	.	Parent=462805;ID=462814
#scaffold680.1	AUGUSTUS	exon	17360	17612	.	-	.	Parent=462805;ID=462815
#scaffold680.1	AUGUSTUS	exon	18188	18275	.	-	.	Parent=462805;ID=462816
#scaffold680.1	AUGUSTUS	exon	18365	18848	.	-	.	Parent=462805;ID=462817
#scaffold680.1	AUGUSTUS	start_codon	18840	18842	.	-	0	Parent=462805;ID=462818
#scaffold680.1	AUGUSTUS	exon	21564	21707	.	-	.	Parent=462805;ID=462819
#scaffold680.1	AUGUSTUS	transcription_start_site	21707	21707	.	-	.	Parent=462805;ID=462820

# Penn V2.0 genome gene model file example:
#Spenn-ch00	AUGUSTUS	gene	5691	9072	0.84	-	.	ID=Sopen00g001010;Name=Sopen00g001010;
#Spenn-ch00	AUGUSTUS	mRNA	5691	9072	0.84	-	.	ID=Sopen00g001010.1;Name=Sopen00g001010.1;Parent=Sopen00g001010;Note=hypothetical protein;
#Spenn-ch00	AUGUSTUS	exon	5691	5899	.	-	.	ID=exon:Sopen00g001010.1.1;Parent=Sopen00g001010.1;
#Spenn-ch00	AUGUSTUS	CDS	5723	5899	1	-	0	ID=cds:Sopen00g001010.1.1;Parent=Sopen00g001010.1;
#Spenn-ch00	AUGUSTUS	intron	5900	8915	1	-	.	ID=intron:Sopen00g001010.1.1;Parent=Sopen00g001010.1;
#Spenn-ch00	AUGUSTUS	CDS	8916	9002	0.99	-	0	ID=cds:Sopen00g001010.1.2;Parent=Sopen00g001010.1;
#Spenn-ch00	AUGUSTUS	exon	8916	9072	.	-	.	ID=exon:Sopen00g001010.1.2;Parent=Sopen00g001010.1;
#Spenn-ch00	AUGUSTUS	gene	11246	12192	0.17	-	.	ID=Sopen00g001020;Name=Sopen00g001020;
#Spenn-ch00	AUGUSTUS	mRNA	11246	12192	0.17	-	.	ID=Sopen00g001020.1;Name=Sopen00g001020.1;Parent=Sopen00g001020;Note=hypothetical protein;
#Spenn-ch00	AUGUSTUS	exon	11246	11553	.	-	.	ID=exon:Sopen00g001020.1.1;Parent=Sopen00g001020.1;
#Spenn-ch00	AUGUSTUS	CDS	11373	11553	0.88	-	1	ID=cds:Sopen00g001020.1.1;Parent=Sopen00g001020.1;
#Spenn-ch00	AUGUSTUS	intron	11554	12021	0.88	-	.	ID=intron:Sopen00g001020.1.1;Parent=Sopen00g001020.1;
#Spenn-ch00	AUGUSTUS	CDS	12022	12134	0.89	-	0	ID=cds:Sopen00g001020.1.2;Parent=Sopen00g001020.1;
#Spenn-ch00	AUGUSTUS	exon	12022	12192	.	-	.	ID=exon:Sopen00g001020.1.2;Parent=Sopen00g001020.1;

#######################################################################################
#######################################################################################
# GTF file format.
# From http://mblab.wustl.edu/GTF22.html
#######################################################################################
#######################################################################################
#
# GTF2.2: A Gene Annotation Format (Revised Ensembl GTF)
#
# GTF Field Definitions Examples Scripts and Resources Introduction
#
# GTF stands for Gene transfer format. It borrows from GFF, but has additional
# structure that warrants a separate definition and format name. Structure is as
# GFF, so the fields are:
#   <seqname> <source> <feature> <start> <end> <score> <strand> <frame> [attributes] [comments]
#
# Here is a simple example with 3 translated exons. Order of rows is not important.
#
# 381 Twinscan  CDS          380   401   .   +   0  gene_id "001"; transcript_id "001.1";
# 381 Twinscan  CDS          501   650   .   +   2  gene_id "001"; transcript_id "001.1";
# 381 Twinscan  CDS          700   707   .   +   2  gene_id "001"; transcript_id "001.1";
# 381 Twinscan  start_codon  380   382   .   +   0  gene_id "001"; transcript_id "001.1";
# 381 Twinscan  stop_codon   708   710   .   +   0  gene_id "001"; transcript_id "001.1";
# The whitespace in this example is provided only for readability. In GTF, fields
# must be separated by a single TAB and no white space.
#
# GTF Field Definitions
#
# <seqname> The name of the sequence. Commonly, this is the chromosome ID or
# contig ID. Note that the coordinates used must be unique within each sequence
# name in all GTFs for an annotation set.
#
# <source> The source column should be a unique label indicating where the
# annotations came from --- typically the name of either a prediction program or
# a public database.
#
# <feature> The following feature types are required: "CDS", "start_codon",
# "stop_codon".  The features "5UTR", "3UTR", "inter", "inter_CNS", "intron_CNS"
# and "exon" are optional.  All other features will be ignored.  The types must
# have the correct capitalization shown here.
#
# CDS represents the coding sequence starting with the first translated codon and
# proceeding to the last translated codon. Unlike Genbank annotation, the stop
# codon is not included in the CDS for the terminal exon. The optional feature
# "5UTR" represents regions from the transcription start site or beginning of the
# known 5' UTR to the base before the start codon of the transcript. If this
# region is interrupted by introns then each exon or partial exon is annotated as
# a separate 5UTR feature. Similarly, "3UTR" represents regions after the stop
# codon and before the polyadenylation site or end of the known 3' untranslated
# region. Note that the UTR features can only be used to annotate portions of mRNA
# genes, not non-coding RNA genes.
#
# The feature "exon" more generically describes any transcribed exon. Therefore,
# exon boundaries will be the transcription start site, splice donor, splice
# acceptor and poly-adenylation site. The start or stop codon will not necessarily
# lie on an exon boundary.
#
# The "start_codon" feature is up to 3bp long in total and is included in the
# coordinates for the "CDS" features. The "stop_codon" feature similarly is up to
# 3bp long and is excluded from the coordinates for the "3UTR" features, if used.
#
# The "start_codon" and "stop_codon" features are not required to be atomic; they
# may be interrupted by valid splice sites. A split start or stop codon appears as
# two distinct features. All "start_codon" and "stop_codon" features must have a
# 0,1,2 in the <frame> field indicating which part of the codon is represented by
# this feature. Contiguous start and stop codons will always have frame 0.
#
# The "inter" feature describes an intergenic region, one which is by almost all
# accounts not transcribed. The "inter_CNS" feature describes an intergenic
# conserved noncoding sequence region. All of these should have an empty
# transcript_id attribute, since they are not transcribed and do not belong to any
# transcript. The "intron_CNS" feature describes a conserved noncoding sequence
# region within an intron of a transcript, and should have a transcript_id
# associated with it.
#
# <start> <end> Integer start and end coordinates of the feature relative to the
# beginning of the sequence named in <seqname>.  <start> must be less than or
# equal to <end>. Sequence numbering starts at 1. Values of <start> and <end> that
# extend outside the reference sequence are technically acceptable, but they are
# discouraged.
#
# <score> The score field indicates a degree of confidence in the feature's
# existence and coordinates. The value of this field has no global scale but may
# have relative significance when the <source> field indicates the prediction
# program used to create this annotation. It may be a floating point number or
# integer, and not necessary and may be replaced with a dot.
#
# <frame> 0 indicates that the feature begins with a whole codon at the 5' most
# base. 1 means that there is one extra base (the third base of a codon) before
# the first whole codon and 2 means that there are two extra bases (the second and
# third bases of the codon) before the first codon. Note that for reverse strand
# features, the 5' most base is the <end> coordinate. Here are the details excised
# from the GFF spec. Important: Note comment on reverse strand.
#
# '0' indicates that the specified region is in frame, i.e. that its first base
# corresponds to the first base of a codon. '1' indicates that there is one extra
# base, i.e. that the second base of the region corresponds to the first base of a
# codon, and '2' means that the third base of the region is the first base of a
# codon. If the strand is '-', then the first base of the region is value of
# <end>, because the corresponding coding region will run from <end> to <start> on
# the reverse strand. Frame is calculated as (3 - ((length-frame) mod 3)) mod 3.
# (length-frame) is the length of the previous feature starting at the first whole
# codon (and thus the frame subtracted out). (length-frame) mod 3 is the number of
# bases on the 3' end beyond the last whole codon of the previous feature.
# 3-((length-frame) mod 3) is the number of bases left in the codon after removing
# those that are represented at the 3' end of the feature. (3-((length-frame) mod
# 3)) mod 3 changes a 3 to a 0, since three bases makes a whole codon, and 1 and 2
# are left unchanged.
#
# [attributes] All nine features have the same two mandatory attributes at the end
# of the record:
#   gene_id value;
# A globally unique identifier for the genomic locus of the transcript. If empty,
# no gene is associated with this feature.
#   transcript_id value;
# A globally unique identifier for the predicted transcript. If empty, no transcript
# is associated with this feature. These attributes are designed for handling multiple
# transcripts from the same genomic region. Any other attributes or comments must
# appear after these two and will be ignored. Attributes must end in a semicolon
# which must then be separated from the start of any subsequent attribute by
# exactly one space character (NOT a tab character).
#
# Textual attributes should be surrounded by doublequotes.
#
# These attributes are required even for non-mRNA transcribed regions such as
# "inter" and "inter_CNS" features.
#
# [comments] Comments begin with a hash ('#') and continue to the end of the line.
# Nothing beyond a hash will be parsed. These may occur anywhere in the file,
# including at the end of a feature line.
#
# Examples
#
# Here is an example of a gene on the negative strand including UTR regions.
# Larger coordinates are 5' of smaller coordinates. Thus, the start codon is 3 bp
# with largest coordinates among all those bp that fall within the CDS regions.
# Note that the stop codon lies between the 3UTR and the CDS
#
# 140    Twinscan    inter         5141     8522     .    -    .    gene_id ""; transcript_id "";
# 140    Twinscan    inter_CNS     8523     9711     .    -    .    gene_id ""; transcript_id "";
# 140    Twinscan    inter         9712     13182    .    -    .    gene_id ""; transcript_id "";
# 140    Twinscan    3UTR          65149    65487    .    -    .    gene_id "140.000"; transcript_id "140.000.1";
# 140    Twinscan    3UTR          66823    66992    .    -    .    gene_id "140.000"; transcript_id "140.000.1";
# 140    Twinscan    stop_codon    66993    66995    .    -    0    gene_id "140.000"; transcript_id "140.000.1";
# 140    Twinscan    CDS           66996    66999    .    -    1    gene_id "140.000"; transcript_id "140.000.1";
# 140    Twinscan    intron_CNS    70103    70151    .    -    .    gene_id "140.000"; transcript_id "140.000.1";
# 140    Twinscan    CDS           70207    70294    .    -    2    gene_id "140.000"; transcript_id "140.000.1";
# 140    Twinscan    CDS           71696    71807    .    -    0    gene_id "140.000"; transcript_id "140.000.1";
# 140    Twinscan    start_codon   71805    71806    .    -    0    gene_id "140.000"; transcript_id "140.000.1";
# 140    Twinscan    start_codon   73222    73222    .    -    2    gene_id "140.000"; transcript_id "140.000.1";
# 140    Twinscan    CDS           73222    73222    .    -    0    gene_id "140.000"; transcript_id "140.000.1";
# 140    Twinscan    5UTR          73223    73504    .    -    .    gene_id "140.000"; transcript_id "140.000.1";
#
# Note the frames of the coding exons. For example:
#
# The first CDS (from 71807 to 71696) always has frame zero. Frame of the 1st CDS
# =0, length =112.  (3-((length - frame) mod 3)) mod 3  = 2, the frame of the 2nd
# CDS. Frame of the 2nd CDS=2, length=88. (3-((length - frame) mod 3)) mod 3  = 1,
# the frame of the terminal CDS. Alternatively, the frame of terminal CDS can be
# calculated without the rest of the gene. Length of the terminal CDS=4. length
# mod 3 =1, the frame of the terminal CDS. Note the split start codon. The second
# start codon region has a frame of 2, since it is the second base, and has an
# accompanying CDS feature, since CDS always includes the start codon.
#
# Here is an example in which the "exon" feature is used. It is a 5 exon gene with
# 3 translated exons.
#
# 381 Twinscan  exon         150   200   .   +   .  gene_id "381.000"; transcript_id "381.000.1";
# 381 Twinscan  exon         300   401   .   +   .  gene_id "381.000"; transcript_id "381.000.1";
# 381 Twinscan  CDS          380   401   .   +   0  gene_id "381.000"; transcript_id "381.000.1";
# 381 Twinscan  exon         501   650   .   +   .  gene_id "381.000"; transcript_id "381.000.1";
# 381 Twinscan  CDS          501   650   .   +   2  gene_id "381.000"; transcript_id "381.000.1";
# 381 Twinscan  exon         700   800   .   +   .  gene_id "381.000"; transcript_id "381.000.1";
# 381 Twinscan  CDS          700   707   .   +   2  gene_id "381.000"; transcript_id "381.000.1";
# 381 Twinscan  exon         900  1000   .   +   .  gene_id "381.000"; transcript_id "381.000.1";
# 381 Twinscan  start_codon  380   382   .   +   0  gene_id "381.000"; transcript_id "381.000.1";
# 381 Twinscan  stop_codon   708   710   .   +   0  gene_id "381.000"; transcript_id "381.000.1";
#
# Scripts and Resources
#
# Several Perl scripts have been written for checking, parsing, correcting, and
# comparing GTF-formatted annotations. Most of the important ones are included the
# Eval package, which comes equipped with a GTF parsing Perl package GTF.pm.
#
#  Eval Software Eval Documentation The Eval documentation contains a complete
#  code-level documentation of GTF.pm, suitable for able Perl programmers to
#  create and parse GTF files.
#
# The script validate_gtf.pl included in the Eval package is particularly useful
# for checking that your GTF annotation is consistent and well-formed. Here are
# some more useful links: GFF Specification at Sanger Brent Lab Homepage Top 

#######################################################################################
# End of file.
#######################################################################################
