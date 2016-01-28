################################################################################
# See usage below for description.
# Author: Ted Toal
# Date: 2015
# Brady Lab, UC Davis
################################################################################

# Enclose everything in braces so stop statements will work correctly.
{

# Pathname separator.
#PATHSEP = ifelse(grepl("/", Sys.getenv("HOME")), "/", "\\")
PATHSEP = "/"

# Get directory where this file resides.
XSEP = ifelse(PATHSEP == "\\", "\\\\", PATHSEP)
RE = paste("^.*--file=(([^", XSEP, "]*", XSEP, ")*)[^", XSEP, "]+$", sep="")
args = commandArgs(FALSE)
thisDir = sub(RE, "\\1", args[grepl("--file=", args)])
#thisDir = "~/Documents/UCDavis/BradyLab/IGGPIPE/IGGPIPE/code/R/" # For testing only.

# Source the necessary include files from the same directory containing this file.
source(paste(thisDir, "Include_Common.R", sep=""))
source(paste(thisDir, "Include_GFFfuncs.R", sep=""))
source(paste(thisDir, "Include_MergeDataUsingPosition.R", sep=""))

################################################################################
# Process program arguments.
################################################################################

# Get arguments.
testing = 0
#testing = 1 # For testing only.
#testing = 2 # For testing only.
#testing = 3 # For testing only.
#testing = 4 # For testing only.
#testing = 5 # For testing only.
#testing = 6 # For testing only.
#testing = 7 # For testing only.
{
if (testing == 0)
    args = commandArgs(TRUE)
else
    {
    # For testing only:
    if (testing == 1)
        args = "annotate.template"
    else if (testing == 2)
        args = "annotate/test_gff3_to_tsv.gene_positions"
    else if (testing == 3)
        args = "annotate/test_tsv_to_gff3.markers"
    else if (testing == 4)
        args = "annotate/HPintrogressions_to_gff3"
    else if (testing == 5)
        args = "annotate/test_addILs_column.markers"
    else if (testing == 6)
        args = "annotate/test_add_genesColumn.markers"
    else if (testing == 7)
        args = "annotate/test_add_isInNearColumn.markers"
    else
        stop("Unknown value for 'testing'")
    }
}

# Get the arguments specific to this program.
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

# In some cases a section name can be a number ("mergeCols" subsections).  Here
# we treat that no differently than named sections (the section name is a number).

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
if (any(grepl(" ", indents)) && any(grepl("\t", indents)))
    cat("Warning: parameter file lines are indented with a mixture of spaces and tab characters.\n",
        "This may cause confusion over how much each line is indented.\n")
indents = gsub("\t", "    ", indents)
indentAmt = nchar(indents)
S = sub("^[ \t]*", "", S)

# Define a recursive function to process the parameter lines of a section and turn
# them into a list.  The list element names are the section names or the parameter
# names, and the list element values are sublists (for section names) or parameter
# values (for parameter names).
#
# Arguments:
#   Sidx: index into parameter line vector S of next line to process
#   indent: amount of indent of the section that contains lines starting at S[Sidx].
#       If a line is indented by this amount or less, it is a line that does not
#       belong within this section, so return from this function.
#
# Returns:
#   A list with members "Sidx" (index of line whose indent caused return from
#   this function) and "L" (the sublist containing the contents of this section).
collectSectionIntoList = function(Sidx, indent)
    {
    # Make empty list for the section.
    L = list()

    # Loop for each line until end of lines or until a line is found which is not
    # a member of this section.
    # Is the next line a member of this section?
    while (Sidx <= length(S) && indentAmt[Sidx] > indent)
        {
        # If next line is a section line, process a sub-section and add its sublist
        # to L, named with the section name.
        if (isSectionHeader[Sidx])
            {
            newSectionName = sub("^[ \t]*([^:]+):$", "\\1", S[Sidx])
            if (grepl("[^A-Za-z0-9_.]", newSectionName))
                stop("Invalid characters in section name '", newSectionName, "' in line: ", Orig[Sidx])

            # Process a subsection and add its sublist to the current section's list.
            L2 = collectSectionIntoList(Sidx+1, indentAmt[Sidx])
            Sidx = L2$Sidx
            L[[newSectionName]] = L2$L
            }
        # Else this line is a parameter line, so add a new element to L, named with
        # the parameter name and with value equal to the parameter value.
        else
            {
            pos = regexpr("=", S[Sidx])
            if (length(pos) == 0)
                stop("Expected to find := on line but did not: ", Orig[Sidx])
            if (pos[1] == 1)
                stop("Missing parameter name on line: ", Orig[Sidx])
            paramName = substring(S[Sidx], 1, pos[1]-1)
            paramVal = substring(S[Sidx], pos[1]+1)
            if (grepl("[^A-Za-z0-9_.]", paramName))
                stop("Invalid characters in parameter name '", paramName, "' in line: ", Orig[Sidx])

            L[[paramName]] = paramVal

            # Next S line.
            Sidx = Sidx + 1
            }
        }
    # Return the current value of Sidx and the section list L that we made.
    return(list(Sidx=Sidx, L=L))
    }

# Get the list at the root level, starting with the first line, and using an
# indent of -1 so all lines with no indent (indent 0) will be part of this
# section.
L = collectSectionIntoList(1, -1)$L

# str(L) # Is it good?

################################################################################
# Verify that we have the parameters we expect and need.
################################################################################

# Define the expected parameters by defining a list with the same structure as
# the expected input data list L.  Member names in this "expected" list must be
# the expected member name, either a section name containing a sublist or a
# parameter name containing a parameter value.  If it is a sublist member, that
# sublist name may optionally have ".match" appended.  If an expected sublist
# member is missing in the input list L, an empty sublist of that name will be
# added to it.  When ".match" is appended to a sublist name in the "expected"
# list, it means that the expected sublist in L is ignored if inputFileS path
# is empty (no match/merge operation is to be done).  In that case the sublist
# is deleted from L.  If any sublist (nested or not) has no names (e.g.
# expected$mergeCols.match), this means there may be any number of instances
# of that sublist at that position in L, and their names don't matter and don't
# even have to exist. For parameter members of "expected", the member value
# indicates the action to take to evaluate the parameter value in L:
#       "rmv.if.empty" : remove the parameter if it is an empty string value 
#       "error.if.empty" : it is an error if the parameter is an empty string or missing
#       "rmv.parent.if.empty" : remove enclosing list if parameter is empty string or missing
#       "NA.if.missing" : set the parameter to NA if it is missing, else leave as-is.
#       "logical.YES" or "logical.NO" : parameter must be YES/NO/TRUE/FALSE/1/0
#           (ignoring case) and it is set to R TRUE or FALSE, or if missing is set
#           to the value after "logical.".
#       "(A,B,C,...)" : the parameter value must equal either A, B, C, ...
#           (ignoring case) or, if it is missing, set it to A.
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
        keepColumns = "rmv.if.empty",
        addColNames = "rmv.if.empty",
        addColValues = "NA.if.missing"
        ),
    inputFileS.match = list(
        path = "rmv.parent.if.empty",
        type = "(tsv,csv,gff3,gtf)",
        columns = "rmv.if.empty",
        column1 = "rmv.if.empty",
        quote = "",
        comment = "",
        keepColumns = "rmv.if.empty",
        addColNames = "rmv.if.empty",
        addColValues = "NA.if.missing"
        ),
    attrExtractT = list(
        attributesColumn = "rmv.parent.if.empty",
        keepFeatures = "rmv.if.empty",
        extractAttrs = "rmv.if.empty",
        excludeAttrs = "rmv.if.empty",
        newAttrColumns = "rmv.if.empty",
        missingAttrValues = "NA.if.missing",
        removeAttrColumns = "logical.NO"
        ),
    attrExtractS.match = list(
        attributesColumn = "rmv.parent.if.empty",
        keepFeatures = "rmv.if.empty",
        extractAttrs = "rmv.if.empty",
        excludeAttrs = "rmv.if.empty",
        newAttrColumns = "rmv.if.empty",
        missingAttrValues = "NA.if.missing",
        removeAttrColumns = "logical.NO"
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
        closest = "rmv.if.empty",
        start.up = "rmv.if.empty",
        start.down = "rmv.if.empty",
        end.up = "rmv.if.empty",
        end.down = "rmv.if.empty"
        ),
    mergeCols.match = list(
        list(
            col = "rmv.parent.if.empty",
            before = "rmv.if.empty",
            maxMatch = "rmv.if.empty",
            join = "logical.YES",
            joinPfx = "",
            joinSfx = "",
            joinSep = ",",
            format = "error.if.empty"
            )
        ),
    attrCreation = list(
        columnsForAttrs = "rmv.parent.if.empty",
        newAttrNames = "rmv.if.empty",
        noAttrValues = "NA.if.missing",
        attributesColumn = "attributes",
        merge = "logical.YES",
        remove = "logical.YES"
        ),
    outputFile = list(
        path = "error.if.empty",
        type = "(tsv,csv,gff3,gtf)",
        outColNames = "rmv.if.empty",
        newColNames = "rmv.if.empty",
        header = "logical.YES",
        column1 = "rmv.if.empty",
        quote = "logical.NO"
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
doMatchMerge = (!is.null(L$inputFileS) && L$inputFileS != "")

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
    # Else if the first element of expected.sublist is a list (without names),
    # then L.sublist can have any number of members like this (mergeCols.column
    # is our only case like this) and their names are optional.
    else if (is.list(expected.sublist[[1]]))
        {
        # There should only be one element in expected.sublist.
        if (length(expected.sublist) != 1)
            stop("Programming error, 'expected' sublist without names should be the only member of its parent")

        # There is only one element in the unnamed sublist of expected.sublist, but
        # but there can be any number of elements in L.sublist.
        # Loop for each member of L.sublist, backwards so if an element is removed,
        # it doesn't shift the index of the remaining elements in the loop.
        for (i in length(L.sublist):1)
            {
            #cat("section:", i, " in ", fullName, "\n")

            # Get name or number for L sublist.
            if (is.null(names(L.sublist)))
                sectionName = i
            else
                sectionName = names(L.sublist)[i]

            # Get fullName for the sublist.
            newFullName = paste(fullName, "[", sectionName, "]", sep="")

            # Process the section's members.
            newSublist = checkParams(expected.sublist[[1]], L.sublist[[sectionName]], newFullName)
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
                    stop("Expected a value for parameter '", param, "' in section '", fullName, "'", call.=FALSE)
                }
            else if (handling == "rmv.parent.if.empty")
                {
                if (is.null(L.sublist[[param]]) || L.sublist[[param]] == "")
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
            else if (handling == "logical.YES" || handling == "logical.NO")
                {
                val = L.sublist[[param]]
                if (is.null(val) || val == "")
                    val = sub("logical\\.", "", handling)
                val = toupper(val)
                if (val %in% c("YES", "TRUE", "1"))
                    L.sublist[[param]] = TRUE
                else if (val %in% c("NO", "FALSE", "0"))
                    L.sublist[[param]] = FALSE
                else
                    stop("Expected a logical value (YES or NO) for parameter ", param)
                }
            else if (grepl("^\\(", handling))
                {
                handling = gsub("[()]", "", handling)
                handling = unlist(strsplit(handling, ",", fixed=TRUE))
                if (is.null(L.sublist[[param]]))
                    L.sublist[[param]] = handling[1]
                else if (!toupper(L.sublist[[param]]) %in% toupper(handling))
                    stop("Expected parameter '", param, "' in section '", fullName, "' to be one of ",
                        expected.sublist[[param]], call.=FALSE)
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
L = checkParams(expected, L, "")

################################################################################
################################################################################
# Define helper functions.
################################################################################
################################################################################

#######################################################################################
# Function to parse outputFile$addColValues format strings and generate a vector of
# data formatted accordingly.
#
# Arguments:
#   format: a outputFile$addColValues format string.
#   T.df: the T.df data frame.
#
# Returns: vector of data formatted as per "format", and of length equal to nrow(T.df).
#######################################################################################
formatNewColData = function(format, T.df)
    {
    chk.col = function(col)
        {
        if (!any(col == colnames(T.df)))
            stop("'", col, "' is not one of the T data frame columns: ", paste(colnames(T.df), collapse=":"), call.=FALSE)
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

    # Split on "{" or "}".
    V = unlist(strsplit(format, "[{}]"))
    if (length(V) %% 2 == 1)
        V = c(V, "")
    if (length(V) == 0) stop("formatData: software error, V is empty")
    V = matrix(V, ncol=2, byrow=TRUE, dimnames=list(NULL, c("verbatim", "format")))

    # Generate formatted strings from the format codes.
    fmtstrs = rep("", nrow(T.df))
    numUniqueID = 0 # 0 means don't add unique add to end.
    for (i in 1:nrow(V))
        {
        code = V[i,2]

        # Add the verbatim text.
        fmtstrs = paste(fmtstrs, V[i, "verbatim"], sep="")

        # Get first character of the code after "{"
        firstChar = substring(code, 1, 1)

        # Process the code.

        # {com} : a verbatim comma
        if (code == "com")
            fmtstrs = paste(fmtstrs, ",", sep="")

        # {lb} : a verbatim left brace
        else if (code == "lb")
            fmtstrs = paste(fmtstrs, "{", sep="")

        # {rb} : a verbatim right brace
        else if (code == "rb")
            fmtstrs = paste(fmtstrs, "}", sep="")

        # {+col} : the value in T.df column "col"
        else if (firstChar == "+")
            {
            col = sub("^\\+(.*)$", "\\1", code)
            chk.col(col)
            fmtstrs = paste(fmtstrs, T.df[[col]], sep="")
            }

        # {/col/RE/RE.replace} : regular expression search/replace of T.df column "col"
        else if (firstChar == "/")
            {
            s = unlist(strsplit(code, "/", fixed=TRUE))
            if (length(s) == 3)
                s[4] = ""
            if (length(s) != 4)
                error("expected {/col/RE/RE.replace} but got: ", code)
            col = s[2]
            chk.col(col)
            RE = s[3]
            RE.replace = s[4]
            S = sub.catch(RE, RE.replace, T.df[[col]])
            fmtstrs = paste(fmtstrs, S, sep="")
            }

        # {#n}
        else if (firstChar == "#")
            {
            numUniqueID = as.integer(substring(code, 2))
            if (numUniqueID <= 0)
                stop("Invalid outputFile$addColValues {#n} code, must be > 0: '", code, "'")
            }

        # else if not empty, it is an invalid code.
        else if (code != "")
            {
            stop("Invalid outputFile$addColValues format code: '", code, "'")
            }
        }

    # If requested to do so, add a unique ID string.
    if (numUniqueID > 0)
        {
        idxs = unlist(split(1:length(fmtstrs), fmtstrs))
        IDs = unlist(as.list(tapply(1:length(fmtstrs), fmtstrs, function(V) 1:length(V))))
        if (numUniqueID > 0)
            {
            Nadd0s = numUniqueID - nchar(IDs)
            Nadd0s[Nadd0s < 0] = 0
            maxAdd0s = max(Nadd0s)
            if (maxAdd0s > 0)
                {
                zeroes = paste(rep("0", maxAdd0s), collapse="")
                IDs = paste(substring(zeroes, 1, Nadd0s), IDs, sep="")
                }
            }
        fmtstrs[idxs] = paste(fmtstrs[idxs], IDs, sep="")
        }

    # Return result.
    return(fmtstrs)
    }

################################################################################
# Function to split a comma-separated string apart into a vector of separate
# strings.  If S is NULL, NULL is returned.  If any substring is NAstring, it
# is changed to NA.
################################################################################

splitCommaList = function(S, NAstring=NULL)
    {
    if (is.null(S))
        return(NULL)
    if (!is.na(S) && S != "")
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
        stop(sectionName, ": ", params$path, " does not exist", call.=FALSE)

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
            args$col.names = splitCommaList(params$columns)
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
            stop(sectionName, " had no data", call.=FALSE)
        if (ncol(df) == 0)
            stop(sectionName, " had no columns", call.=FALSE)
        # columns
        if (!is.null(params$columns))
            {
            columns = splitCommaList(params$columns)
            if (length(columns) != ncol(df))
                stop("Expected ", sectionName, " to have ", length(columns), " columns but it had ",
                    ncol(df), " which were ", paste(colnames(df), collapse=":"), call.=FALSE)
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
        keepColumns = splitCommaList(params$keepColumns)
        missing = keepColumns[!keepColumns %in% colnames(df)]
        if (length(missing) > 0)
            stop(sectionName, " does not have these 'keepColumns':\n",
                paste(missing, collapse=":"), "\nIt has:\n",
                paste(colnames(df), collapse=":"), call.=FALSE)
        df = df[,keepColumns]
        }

    # addColNames, addColValues.
    addColNames = splitCommaList(params$addColNames)
    addColValues = splitCommaList(params$addColValues, NAstring="R.NA")
    if (!is.null(addColNames))
        {
        if (length(addColValues) == 1)
            addColValues = rep(addColValues, length(addColNames))
        else if (length(addColValues) != length(addColNames))
            stop("Expected ", length(addColNames), " column values in addColValues, section '",
                sectionName, "' but found ", length(addColValues))
        names(addColValues) = addColNames
        for (col in addColNames)
            df[[col]] = formatNewColData(addColValues[col], df)
        }

    return(df)
    }

################################################################################
################################################################################
# Ready.  Do it.
################################################################################
################################################################################

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

if (doMatchMerge)
    {
    # mergeOnMatches = function(T.df, S.df, T.pos, S.pos, match, mergeCols)
    T.pos = L$positionT
    S.pos = L$positionS
    match = L$match
    mergeCols = L$mergeCols
    T.df = mergeOnMatches(T.df, S.df, T.pos, S.pos, match, mergeCols)
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
    outColNames = splitCommaList(L$outputFile$outColNames)
    if (!is.null(outColNames))
        {
        missing = outColNames[!outColNames %in% colnames(T.df)]
        if (length(missing) > 0)
            stop("No such columns in output data frame: ", paste(missing, collapse=":"))
        T.df = T.df[, outColNames]
        }

    # Rename columns if requested.
    newColNames = splitCommaList(L$outputFile$newColNames)
    if (!is.null(newColNames))
        {
        if (length(newColNames) != ncol(T.df))
            stop("Expected ", ncol(T.df), " column names in newColNames but found ", length(newColNames))
        colnames(T.df) = newColNames
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
            write.table.winSafe(T.df, L$outputFile$path, quote=L$outputFile$quote,
                col.names=col.names, row.names=row.names, sep="\t")
            }
        # csv
        else
            {
            write.csv.winSafe(T.df, L$outputFile$path, quote=L$outputFile$quote,
                col.names=col.names, row.names=row.names)
            }
        }
    # gff3 or gtf
    else
        {
        if (ncol(T.df) != 9)
            stop("GFF3 and GTF files must have 9 columns, this has ", ncol(T.df))
        writeFile_GFF3_GTF(T.df, L$outputFile$path)
        }
    }

cat("Output file ", L$outputFile$path, " written\n")

################################################################################
# Finished.
################################################################################

cat("Finished conversion and/or merge operation.\n")
}

################################################################################
# End of file.
################################################################################
