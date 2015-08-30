#######################################################################################
# See usage below for description.
# Author: Ted Toal
# Date: 2015
# Brady Lab, UC Davis
#######################################################################################

# Enclose everything in braces so stop statements will work correctly.
{

# Pathname separator.
PATHSEP = ifelse(grepl("/", Sys.getenv("HOME")), "/", "\\")

# For some reason this is needed when using RScript.
library("methods")

################################################################################
# Process program arguments.
################################################################################

# Get arguments.
testing = FALSE
#testing = TRUE # For testing only.
{
if (!testing)
    args = commandArgs(TRUE)
# For testing only:
else
    {
    args = "annotate.template"
    }
}

Nexpected = 1
if (length(args) != 1)
    {
    usage = c(
        "Read a parameter file that specifies parameters for reading one or two data files",
        "containing columns of data, and parameters for converting the files to other formats",
        "and modifying their column data.  See annotate.template file for sample parameter file.",
        "",
        "Usage: Rscript annotateFile.R <paramFile>",
        "",
        "Arguments:",
        "   <paramFile> : Path of parameter file, see example parameter file annotate.template."
        )
    for (S in usage)
        cat(S, "\n", sep="")
    stop("Try again with correct number of arguments")
    }

cat("annotateFile.R arguments:\n")
paramFile = args[1]
cat("  paramFile: ", paramFile, "\n")
if (!file.exists(paramFile))
    stop("File doesn't exist: ", paramFile)


# Get directory where this file resides.
XSEP = ifelse(PATHSEP == "\\", "\\\\", PATHSEP)
RE = paste("^.*--file=(([^", XSEP, "]*", XSEP, ")*)[^", XSEP, "]+$", sep="")
allArgs = commandArgs()
thisDir = sub(RE, "\\1", allArgs[grepl("--file=", allArgs)])

# Source the necessary include files from the same directory containing this file.
source(paste(thisDir, "Include_GFFfuncs.R", sep=""))
source(paste(thisDir, "Include_MergeDataUsingPosition.R", sep=""))

################################################################################
# Functions.
################################################################################

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
# Read parameter file and extract relevant lines.
################################################################################

# Read the parameter file lines.
S = readLines(paramFile)

# Remove comments.
hasQuote = grepl('^[^#]*"', S)
S[!hasQuote] = sub("#.*$", "", S[!hasQuote])
S[hasQuote] = sub('^([^#]*"[^"]*".*)#.*$', "\\1", S[hasQuote])

# Trim whitespace from end of line, but keep whitespace at start of line.
S = sub("[ \t]+$", "", S)

# Remove blank lines.
S = S[!grepl("^[ \t]*$", S)]

# Keep a copy of the lines as they appear at this point, for error messages.
Orig = S

# Convert := to =.
S = sub(":=", "=", S)

# Remove whitespace around =.
S = sub("[ \t]+=[ \t]+", "=", S)

# Remove quote marks.
S = sub('="(.*)"', "=\\1", S)
S = sub("='(.*)'", "=\\1", S)

################################################################################
# Convert the parameter lines into a structured list.
################################################################################

# Lines ending in ":" are section header lines that mark the start of sections.

# Each section will be a sublist of a main list L, named with the section name.
# Sections can be nested, and this creates sub-sublists.  

# Since the same section name can appear multiple times ("column" section), when
# this happens, insert two additional sublists between these multiple-instance
# lists and their parent, as this example shows:
#   First occurrence of section "column" for new column "gene":
#       L
#         L[["column"]]:
#           .indent = 4
#           .outer = "mergeCols"
#           name = "gene"
#           etc.
#   After second occurrence of section "column" for new column "position":
#       L
#         L[["column"]]:
#           .indent = 4
#           .outer = "mergeCols"
#           .sublists
#             .sublists[[1]]:
#               name = "gene"
#               etc.
#             .sublists[[2]]:
#               name = "position"
#               etc.

# Lines not ending in ":" are parameter lines consisting of a parameter name and
# value, and these are placed in the sublist for the section they are in, with
# the parameter name becoming the list element name, and value becoming the
# element value.
#
# The amount of indent on each line determines whether that line belongs in
# the most recent preceding section; if the indent is more than that preceding
# section line, it belongs in that section, else it belongs in the section
# before that IF the indent is more than that one, etc., and if there is no
# unindented section header, a line that doesn't belong in any of the sections
# becomes an element in the main list L.  To track indent level, each section's
# sublist has an element named .indent that contains the length of the indent
# string for that section header, and another member .outer that is the name
# of the section that contains it, or ".outermost" if it is an outermost section.
# Tabs are converted into four spaces for the purposes of counting indent amount.

isSectionHeader = grepl(":$", S)
indents = sub("[^ \t].*$", "", S)
indents = gsub("\t", "    ", indents)
indentAmt = nchar(indents)
S = sub("^[ \t]*", "", S)

# Loop for each line and insert the line into list L.  Initially, all sublists
# appear as a separate sublist inside L, and are not nested as they should be.
# We will nest them properly after this.
L = list()
curSection = "outermost"
L[[curSection]]$.indent = -1
L[[curSection]]$.outer = "outermost"
for (i in 1:length(S))
    {
    # First determine which section is the parent section that contains this line,
    # and make that section be "curSection".
    while (indentAmt[i] <= L[[curSection]]$.indent)
        curSection = L[[curSection]]$.outer
    #cat("Parent section: ", curSection, "\n")

    # If this line is a section header line, start a new list in L, named with the
    # section name.  If this section name already exists within L, use two nested
    # lists as shown in example above.
    if (isSectionHeader[i])
        {
        newSection = sub("^[ \t]*([^:]+):$", "\\1", S[i])
        if (grepl("[^A-Za-z0-9_.]", newSection))
            stop("Invalid characters in section name '", newSection, "' in line: ", Orig[i])

        # Basic: add a new section sublist.
        if (!any(names(L) == newSection))
            {
            L2 = list()
            L2$.indent = indentAmt[i]
            L2$.outer = curSection
            L[[newSection]] = L2
            }
        # Repeated section name: if indexed sublist not yet created, create it.
        else if (is.null(L[[newSection]]$.sublists))
            {
            L2 = L[[newSection]]
            .indent = L2$.indent
            .outer = L2$.outer
            L2$.indent = NULL
            L2$.outer = NULL
            L[[newSection]] = list(.indent=.indent, .outer=.outer, .sublists=list(L2, list()))
            }
        # Repeated section name: if indexed sublist has been created, insert a new empty sublist.
        else
            L[[newSection]]$.sublists = c(L[[newSection]]$.sublists, list())

        # Record new section name as current section.
        curSection = newSection
        #cat("New section: ", curSection, "\n")
        }

    # If this line is a parameter line, add a new element to the list for the current
    # section, with the parameter name as the element name and parameter value as the value.
    else
        {
        pos = regexpr("=", S[i])
        if (length(pos) == 0)
            stop("Expected to find := on line but did not: ", Orig[i])
        if (pos[1] == 1)
            stop("Missing parameter name on line: ", Orig[i])
        paramName = substring(S[i], 1, pos[1]-1)
        paramVal = substring(S[i], pos[1]+1)
        if (grepl("[^A-Za-z0-9_.]", paramName))
            stop("Invalid characters in parameter name '", paramName, "' in line: ", Orig[i])

        # If the current section is an indexed section (e.g. the "column" example
        # above), add the new parameter to the last sublist of the current section
        # .sublists member.  Otherwise, just add it to the current section sublist.
        if (is.null(L[[curSection]]$.sublists))
            L[[curSection]][[paramName]] = paramVal
        else
            L[[curSection]]$.sublists[[length(L[[curSection]]$.sublists)]][[paramName]] = paramVal
        }
    }
# str(L) # Is it good?  Especially, L$column

# Now loop through each list member of L and move it into its proper containing
# parent if it isn't at the outermost level.  For this to work, we must do it
# from the innermost lists to the outermost ones, so that each parent can always
# be found at the outer level of L because it hasn't been moved into a nested
# sublist yet.  So, work from largest .indent to smallest.  Get rid of the
# .indent and .outer members, and also, collapse .sublists member up into its
# parent, so for example we have "L$mergeCols$column[[i]]" rather than
# "L$mergeCols$column$.sublists[[i]]".
L$outermost = NULL
indents = sapply(L, function(L) { if (!is.list(L)) return(NA); return(L$.indent) })
outer = sapply(L, function(L) L$.outer)
ord = order(indents, decreasing=TRUE)
indents = indents[ord]
outer = outer[ord]
subsections = names(indents)
for (subsection in subsections)
    {
    outer = L[[subsection]]$.outer
    L[[subsection]]$.indent = NULL
    L[[subsection]]$.outer = NULL
    L2 = L[[subsection]]$.sublists
    if (!is.null(L2))
        {
        L[[subsection]]$.sublists = NULL
        L[[subsection]] = L2
        } 
    if (outer != "outermost")
        {
        L[[outer]][[subsection]] = L[[subsection]]
        L[[subsection]] = NULL
        }
    }
# str(L) # Now how does it look?

################################################################################
# Verify that we have the parameters we expect and need.
################################################################################

# Define the expected parameters by defining a list with the same structure as
# the expected list.  Member names in this "expected" list must be the expected
# member name OR, for a list member, that name with ".match" appended.  If a
# list member is missing, an empty list is added as that member.  When ".match"
# is appended, it means that list is ignored if inputFileS path is empty (no
# match/merge operation is to be done).  In that case the list is deleted.
# The value of each non-list member indicates the action to take to evaluate the
# parameter value:
#       "rmv.if.empty" : remove the parameter if it is an empty string value 
#       "error.if.empty" : it is an error if the parameter is an empty string or missing
#       "rmv.parent.if.empty" : remove enclosing list if parameter is empty string
#       "(A,B,C,...)" : the parameter value must equal either A, B, C, ...
#           (ignoring case) or, if it is missing, set it to A.
#       "NA.if.missing" : set the parameter to NA if it is missing, else leave as-is.
#       any other string : set the parameter to that string if the parameter is an
#           empty string or missing.

expected = list(
    inputFileT = list(
        path = "error.if.empty",
        type = "(tsv,csv,gff3,gtf)",
        columns = "rmv.if.empty",
        column1 = "rmv.if.empty",
        quote = "",
        comment = "",
        keepColumns = "rmv.if.empty"
        ),
    inputFileS.match = list(
        path = "rmv.parent.if.empty",
        type = "(tsv,csv,gff3,gtf)",
        columns = "rmv.if.empty",
        column1 = "rmv.if.empty",
        quote = "",
        comment = "",
        keepColumns = "rmv.if.empty"
        ),
    attrExtractT = list(
        attributesColumn = "rmv.parent.if.empty",
        keepFeatures = "rmv.if.empty",
        extractAttrs = "rmv.if.empty",
        excludeAttrs = "rmv.if.empty",
        newAttrColumns = "rmv.if.empty",
        missingAttrValues = "NA.if.missing",
        removeAttrColumns = "(NO,YES)"
        ),
    attrExtractS.match = list(
        attributesColumn = "rmv.parent.if.empty",
        keepFeatures = "rmv.if.empty",
        extractAttrs = "rmv.if.empty",
        excludeAttrs = "rmv.if.empty",
        newAttrColumns = "rmv.if.empty",
        missingAttrValues = "NA.if.missing",
        removeAttrColumns = "(NO,YES)"
        ),
    positionT.match = list(
        start = "error.if.empty",
        end = "rmv.if.empty",
        len = "rmv.if.empty",
        id = "rmv.if.empty"
        ),
    positionS.match = list(
        start = "error.if.empty",
        end = "rmv.if.empty",
        len = "rmv.if.empty",
        id = "rmv.if.empty"
        ),
    match.match = list(
        method = "(OVERLAP,S.TINY,T.TINY,S.NEAR,T.NEAR)",
        start.up = "rmv.if.empty",
        start.down = "rmv.if.empty",
        end.up = "rmv.if.empty",
        end.down = "rmv.if.empty",
        closest = "rmv.if.empty"
        ),
    mergeCols.match = list(
        column = list(
            list(
                name = "rmv.parent.if.empty",
                before = "rmv.if.empty",
                maxMatch = "rmv.if.empty",
                join = "(YES,NO)",
                joinPfx = "",
                joinSfx = "",
                joinSep = ",",
                format = "error.if.empty"
                )
            )
        ),
    attrCreation = list(
        columnsForAttrs = "rmv.parent.if.empty",
        newAttrNames = "",
        noAttrValues = "NA.if.missing",
        attributesColumn = "attributes",
        merge = "(YES,NO)",
        remove = "(YES,NO)"
        ),
    outputFile = list(
        path = "error.if.empty",
        type = "(tsv,csv,gff3,gtf)",
        columns = "rmv.if.empty",
        header = "rmv.if.empty",
        column1 = "rmv.if.empty",
        quote = "(NO,YES)"
        )
    )

# Use the "expected" structure to verify the parameter data.
#
# Define a recursive function for testing each level of parameter data in L.
# Remove ".match" from the names of sublists after handling it properly.
# Parameters or their parent sections are deleted as requested by "expected" list.
#
# Arguments:
#   expected.sublist: a sublist of the "expected" list, or the list itself.
#   L.sublist: the sublist of L at the same level as expected.sublist.
#   fullName: names of sections within which L.sublist is embedded, "" if none.
#   doMatchMerge: global parameter, TRUE if match/merge is being done.
#
# Returns: the modified L.sublist.

doMatchMerge = (!is.null(L$inputFileS) && inputFileS != "")
checkParams = function(expected.sublist, L.sublist, fullName)
    {
    #cat("fullName:", fullName, "\n")
    # If the first element of expected.sublist is a list and if expected.sublist
    # has names, the names are section names.
    if (is.list(expected.sublist[[1]]) && !is.null(names(expected.sublist)))
        {
        # Loop for each section name in expected.sublist
        for (section in names(expected.sublist))
            {
            #cat("section:", section, " in ", fullName, "\n")
            # Strip ".match" from the section name.
            removeIfNoMatchMerge = grepl("\\.match$", section)
            sectionName = sub("\\.match$", "", section)
            # Add the sectionName to L.sublist as an empty list if not already there.
            if (is.null(L.sublist[[sectionName]]))
                {
                L.sublist[[sectionName]] = list()
                #cat("Remove section", sectionName, " of ", fullName, "\n")
                }
            # If ".match" and not doing match/merge, remove the entire section.
            if (removeIfNoMatchMerge && !doMatchMerge)
                {
                L.sublist[[sectionName]] = NULL
                #cat("Remove section", sectionName, " of ", fullName, "\n")
                }
            # Process the section's members.
            else
                {
                newFullName = ifelse(fullName == "", sectionName, paste(fullName, "$", sectionName, sep=""))
                L.sublist[[sectionName]] = checkParams(expected.sublist[[section]], L.sublist[[sectionName]], newFullName)
                if (is.null(L.sublist[[sectionName]]))
                    {
                    #cat("Remove section", sectionName, " of ", fullName, "\n")
                    }
                }
            }
        }
    # Else if the first element of expected.sublist is a list then expected.sublist
    # must not have names, so expected.sublist is an unnamed section of parameters
    # (mergeCols.column is our only case like this).
    else if (is.list(expected.sublist[[1]]))
        {
        # Make sure L.sublist is unnamed.
        if (!is.null(names(L.sublist)))
            stop("Unexpected parameter names found within section '", fullName, "'")
        # There is only one element in the unnamed sublist of expected.sublist, but
        # but there can be any number of elements in L.sublist.
        # Loop for each member of L.sublist, backwards so if an element is removed,
        # it doesn't shift the index of the remaining elements in the loop.
        for (i in length(L.sublist):1)
            {
            #cat("section:", i, " in ", fullName, "\n")
            # Process the section's members.
            newFullName = paste(fullName, "[", i, "]", sep="")
            newSublist = checkParams(expected.sublist[[1]], L.sublist[[i]], newFullName)
            L.sublist[[i]] = newSublist
            if (is.null(newSublist))
                {
                #cat("Remove section", i, " of ", fullName, "\n")
                }
            }
        }
    # Else expected.sublist names are parameter names.
    else
        {
        # Loop for each parameter name in expected.sublist
        for (param in names(expected.sublist))
            {
            handling = expected.sublist[[param]]
            #cat("param:", param, " in ", fullName, " handling:", handling, "\n")
            if (handling == "rmv.if.empty")
                {
                if (!is.null(L.sublist[[param]]) && L.sublist[[param]] == "")
                    {
                    L.sublist[[param]] = NULL
                    #cat(" Parameter removed\n")
                    }
                }
            else if (handling == "error.if.empty")
                {
                if (is.null(L.sublist[[param]]) || L.sublist[[param]] == "")
                    stop("Expected a value for parameter '", param, "' in section '", fullName, "'")
                }
            else if (handling == "rmv.parent.if.empty")
                {
                if (!is.null(L.sublist[[param]]) && L.sublist[[param]] == "")
                    {
                    #cat(" Removing parent\n")
                    return(NULL)
                    }
                }
            else if (handling == "NA.if.missing")
                {
                if (is.null(L.sublist[[param]]))
                    {
                    L.sublist[[param]] = NA
                    #cat(" Parameter set to NA\n")
                    }
                }
            else if (grepl("^\\(", handling))
                {
                handling = gsub("[()]", "", handling)
                handling = unlist(strsplit(handling, ",", fixed=TRUE))
                if (is.null(L.sublist[[param]]))
                    L.sublist[[param]] = handling[1]
                else if (!toupper(L.sublist[[param]]) %in% toupper(handling))
                    stop("Expected parameter '", param, "' in section '", fullName, "' to be one of ",
                        expected.sublist[[param]])
                }
            else
                {
                if (is.null(L.sublist[[param]]) || L.sublist[[param]] == "")
                    L.sublist[[param]] = handling
                }
            }
        }
    return(L.sublist)
    }

# Now call it to check the parameters.
L = checkParams(expected, L, "", !is.null(L$inputFileS))

################################################################################
################################################################################
# Ready.  Do it.
################################################################################
################################################################################

################################################################################
# Function to split a comma-separated string apart into a vector of separate
# strings.  If S is NULL, NULL is returned.  If any substring is NAstring, it
# is changed to NA.
################################################################################

splitCommaList = function(S, NAstring=NULL)
    {
    if (is.null(S))
        return(NULL)
    S = unlist(strsplit(S, ",", fixed=TRUE))
    if (!is.null(NAstring))
        S[S == NAstring] = NA
    return(S)
    }

################################################################################
# Function to read an input file.
################################################################################

readInputFile = function(sectionName, params)
    {
    if (!file.exists(params$path))
        stop(sectionName, ": ", params$path, " does not exist")

    # gff3
    if (params$type == "gff3")
        {
        df = readFile_GFF3(params$path)
        df = clean_GFF3(df)
        }
    # gtf
    else if (params$type == "gtf")
        {
        df = readFile_GTF(params$path)
        df = clean_GTF(df)
        }
    # tsv or csv
    else
        {
        # Set up args to call read function.
        args = list(file=params$path, stringsAsFactors=FALSE)
        if (params$type == "tsv")
        args$header = is.null(params$columns)
        if (!args$header)
            args$col.names = params$columns
        if (params$type == "tsv")
            args$sep = "\t"
        args$quote = params$quote
        args$comment.char = params$comment
        # Args are set up, now read the file.
        if (params$type == "tsv")
            df = do.call(read.table, args)
        else
            df = do.call(read.csv, args)
        if (nrow(df) == 0)
            stop(sectionName, " had no data")
        if (ncol(df) == 0)
            stop(sectionName, " had no columns")
        # columns
        if (!is.null(params$columns))
            {
            columns = params$columns
            if (length(columns) != ncol(df))
                stop("Expected ", sectionName, " to have ", length(columns), " columns but it had ",
                    ncol(df), " which were ", paste(colnames(df), collapse=","))
            colnames(df) = columns
            }
        # column1
        if (!is.null(params$column1))
            df[[params$column1]] = rownames(df)
        rownames(df) = NULL
        }

    # keepColumns
    if (!is.null(params$keepColumns))
        {
        columns = params$keepColumns
        missing = columns[!columns %in% colnames(df)]
        if (length(missing) > 0)
            stop(sectionName, " does not have these 'keepColumns': ", paste(missing, collapse=","))
        df = df[,columns]
        }
    return(df)
    }

################################################################################
# Read inputFileT.
################################################################################

T.df = readInputFile("inputFileT", L$inputFileT)

################################################################################
# Read inputFileS if it was specified.
################################################################################

S.df = NULL
if (doMatchMerge)
    S.df = readInputFile("inputFileS", L$inputFileS)

################################################################################
# Extract attributes from inputFileT.
################################################################################

if (!is.null(L$attrExtractT))
    {
    keepFeatures = splitCommaList(L$attrExtractT$keepFeatures)
    extractAttrs = splitCommaList(L$attrExtractT$extractAttrs)
    excludeAttrs = splitCommaList(L$attrExtractT$excludeAttrs)
    newAttrColumns = splitCommaList(L$attrExtractT$newAttrColumns)
    missingAttrValues = splitCommaList(L$attrExtractT$missingAttrValues, NAstring="R.NA")

    if (!is.null(L$attrExtractT$keepFeatures))
        T.df = selectFeatures(T.df, keepFeatures)

    T.df = convertAttrsToCols(T.df, extractAttrs, excludeAttrs,
        L$attrExtractT$removeAttrColumns,
        L$attrExtractT$attributesColumn, newAttrColumns, missingAttrValues)
    }

################################################################################
# Extract attributes from inputFileS.
################################################################################

if (doMatchMerge && !is.null(L$attrExtractS))
    {
    keepFeatures = splitCommaList(L$attrExtractS$keepFeatures)
    extractAttrs = splitCommaList(L$attrExtractS$extractAttrs)
    excludeAttrs = splitCommaList(L$attrExtractS$excludeAttrs)
    newAttrColumns = splitCommaList(L$attrExtractS$newAttrColumns)
    missingAttrValues = splitCommaList(L$attrExtractS$missingAttrValues, NAstring="R.NA")

    if (!is.null(L$attrExtractS$keepFeatures))
        S.df = selectFeatures(S.df, keepFeatures)

    S.df = convertAttrsToCols(S.df, extractAttrs, excludeAttrs,
        L$attrExtractS$removeAttrColumns,
        L$attrExtractS$attributesColumn, newAttrColumns, missingAttrValues)
    }

################################################################################
# Match and merge.
################################################################################

if (doMatchMerge && !is.null(L$attrExtractS))
    {
    T.df = mergeOnMatches(T.df, S.df, L$positionT, L$positionS, L$match, L$mergeCols)
    }

################################################################################
# Create attributes in output file.
################################################################################

if (!is.null(L$attrCreation))
    {
    columnsForAttrs = splitCommaList(L$attrCreation$columnsForAttrs)
    newAttrNames = splitCommaList(L$attrCreation$newAttrNames)
    noAttrValues = splitCommaList(L$attrCreation$noAttrValues, NAstring="R.NA")

    T.df = convertColsToAttrs(T.df, columnsForAttrs, newAttrNames, noAttrValues,
        L$attrCreation$attributesColumn, L$attrCreation$merge, L$attrCreation$remove)
    }

################################################################################
# Create output file.
################################################################################

# always
    {
    # Reduce output to only the requested columns.
    columns = splitCommaList(L$outputFile$columns)
    if (!is.null(columns))
        {
        missing = columns[!columns %in% colnames(T.df)]
        if (any(missing))
            stop("No such columns in output data frame: ", paste(missing, collapse=","), "\n")
        T.df = T.df[, columns]
        }

    # Eliminate any possible row names.
    rownames(T.df) = NULL

    # Write the data to the output file.

    # tsv or csv
    if (L$outputFile$type == "tsv" || L$outputFile$type == "csv")
        {
        col.names = L$outputFile$header
        row.names = FALSE

        if (!is.null(L$outputFile$column1))
            {
            row.names = TRUE
            if (col.names)
                col.names = NA
            column1 = L$outputFile$column1
            if (column1 != "#")
                {
                if (!any(colnames(T.df) == column1))
                    stop("No such column in output data frame: ", column1)
                if (any(duplicated(T.df[[column1]])))
                    stop("Duplicate values are not allowed in output data frame column: ", column1)
                }
            }

        # tsv
        if (L$outputFile$type == "tsv")
            {
            write.table(T.df, L$outputFile$path, quote=L$outputFile$quote,
                col.names=col.names, row.names=row.names, sep="\t")
            }
        # csv
        else
            {
            write.csv(T.df, L$outputFile$path, quote=L$outputFile$quote,
                col.names=col.names, row.names=row.names)
            }
        }
    # gff3 or gtf
    else
        {
        writeFile_GFF3_GTF(T.df, L$outputFile$path)
        }
    }

################################################################################
# Finished.
################################################################################

cat("Finished conversion and/or merge operation.\n")
}

################################################################################
# End of file.
################################################################################
