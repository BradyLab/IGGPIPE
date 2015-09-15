###############################################################################
# General-purpose utility functions by Ted Toal.
#
# To make use of these utility functions, do:
# use lib "/path-to-directory-containing-the-file";
# use TedsLibrary ':all';
###############################################################################

package TedsLibrary;

use strict;
use warnings;

require Exporter;
our @ISA = qw(Exporter);

###############################################################################
# EXPORT all module global variables and all module functions when ":all"
# is specified.
###############################################################################
our %EXPORT_TAGS = ( all => [ qw(
    %nt_to_aa %aaInfo1 %aaInfo3 %aaInfo
    getComputer readTabFile readTabFileNoHeader readTabFileWithHeader writeTabFile
    getTime log10 log2 signif
    createLogFile closeLogFile printLog showFileAndLogOutput
    eolChomp commify makeNumberPretty interpolateTildes findFirstMatch indexOf
    ucArray trArray sumArray countTrueArray
    R_table R_unique R_intersect R_union R_setdiff
    complement reverseComplement entropy KLdivergence
    ntToAA AA1_to_AA3 AA1_ToNames AA3_to_AA1 AA3_ToNames
    AA1_GetPolarity AA1_GetCharge AA1_GetHydropathy PAM10
    averageInWindow
    findORF lenORF countCodons countCodons_Array checkGoodCodingSeq
    readGenBankFile_genomeSeq readGenBankFile_CDSseq readGenBankFile_countCodons
    readFastaSeq
    ) ] );
our @EXPORT_OK = ( @{ $EXPORT_TAGS{'all'} } );
our @EXPORT = qw();
our $VERSION = '0.01';

# Preloaded methods go here.

###############################################################################
# Return the name of the computer we are running on, obtained from the OS
# hostname command.
###############################################################################
sub getComputer {
    my $computer = `hostname`;
    $computer = eolChomp($computer);
    $computer =~ s/\..*//;
    return($computer);
}

###############################################################################
# Read a tab-separated file containing a header line.
# Arguments:
#   $fileName: name of file to read.
#   $keyName: (optional argument) name of column heading to be used as the
#       key to the hash that stores the data.  If this argument is "" or is
#       not specified, the data is stored in an array rather than a hash.
#   $commentChar: (optional argument) character indicating a comment line.
#       If this argument is "" or is not specified, comment lines are
#       assumed not to exist.
# Returns: list ($data, $columns):
#   $data: reference to a hash ($keyName specified) or an array ($keyName
#       not specified) whose elements are references to hashes whose keys
#       are the column names and whose values are the column values read
#       from the file.
#   $columns: reference to array of column names.
###############################################################################
sub readTabFile {
    my $fileName = shift;
    my $keyName = "";
    if (scalar(@_) > 0) { $keyName = shift; }
    my $commentChar = "";
    if (scalar(@_) > 0) { $commentChar = shift; }
    open(FILE_DATA, "<" . $fileName) or die("Can't open file $fileName: $!");
    # Get header line and extract column names.
    my $hdr = <FILE_DATA>;
    $hdr = eolChomp($hdr);
    if ($commentChar ne "") {
        while ($hdr =~ m/^$commentChar/) { $hdr = <FILE_DATA>; $hdr = eolChomp($hdr); }
    }
    $hdr =~ s/[\x0A\x0D]+//g;
    $hdr =~ s/\"//g;
    my @columnNames = split("\t", $hdr, -1);
    my %columnIdx;
    @columnIdx{@columnNames} = (0..$#columnNames); # Hash to convert column name to column index.
    if ($keyName ne "" && !defined($columnIdx{$keyName})) {
        die("readTabFile: $keyName does not exist as a column heading in file $fileName");
    }
    my @a;
    my %h;
    while (defined(my $line = <FILE_DATA>)) {
        $line = eolChomp($line);
        if ($commentChar ne "") {
            while ($line =~ m/^$commentChar/) { $line = <FILE_DATA>; $line = eolChomp($line); }
        }
        $line =~ s/[\x0A\x0D]+//g;
        my @elements = split("\t", $line, -1);
        my %h2;
        foreach my $col (@columnNames) {
            my $s = $elements[$columnIdx{$col}];
            $s =~ s/^"//;
            $s =~ s/"$//;
            $h2{$col} = $s;
        }
        if ($keyName eq "") { push(@a, \%h2); }
        elsif (defined($h{$h2{$keyName}})) { die("Key column value $h2{$keyName} not unique in file $fileName"); }
        else { $h{$h2{$keyName}} = \%h2; }
    }
    close(FILE_DATA);
    if ($keyName eq "") { return(\@a, \@columnNames); }
    return(\%h, \@columnNames);
}

###############################################################################
# Read a tab-separated file that does NOT contain a header line.
# Arguments:
#   $fileName: name of file to read.
#   $commentChar: (optional argument) character indicating a comment line.
#       If this argument is "" or is not specified, comment lines are
#       assumed not to exist.
# Returns: reference to an array whose elements are references (one per
# line or row of file data) to arrays of the values read for each row.
###############################################################################
sub readTabFileNoHeader {
    my $fileName = shift;
    my $commentChar = "";
    if (scalar(@_) > 0) { $commentChar = shift; }
    open(FILE_DATA, "<" . $fileName) or die("Can't open file $fileName: $!");
    my @a;
    while (defined(my $line = <FILE_DATA>)) {
        $line = eolChomp($line);
        if ($commentChar ne "") {
            while ($line =~ m/^$commentChar/) { $line = <FILE_DATA>; $line = eolChomp($line); }
        }
        $line =~ s/[\x0A\x0D]+//g;
        my @elements = split("\t", $line, -1);
        push(@a, \@elements);
    }
    close(FILE_DATA);
    return(\@a);
}

###############################################################################
# Read a tab-separated file that contains a header line.  Similar to readTabFile
# except that the return value is a hash of references to arrays.
# Arguments:
#   $fileName: name of file to read.
#   $commentChar: (optional argument) character indicating a comment line.
#       If this argument is "" or is not specified, comment lines are
#       assumed not to exist.
# Returns: list ($data, $columns):
#   $data: reference to a hash whose keys are the column names and whose values
#       are references to arrays of values in those columns.
#   $columns: reference to array of column names.
###############################################################################
sub readTabFileWithHeader {
    my $fileName = shift;
    my $commentChar = "";
    if (scalar(@_) > 0) { $commentChar = shift; }
    open(FILE_DATA, "<" . $fileName) or die("Can't open file $fileName: $!");
    # Get header line and extract column names.
    my $row = 0;
    my $hdr = <FILE_DATA>;
    $hdr = eolChomp($hdr);
    $row++;
    if ($commentChar ne "") {
        while ($hdr =~ m/^$commentChar/) { $hdr = <FILE_DATA>; $hdr = eolChomp($hdr); $row++; }
    }
    $hdr =~ s/[\x0A\x0D]+//g;
    $hdr =~ s/\"//g;
    my @columnNames = split("\t", $hdr, -1);
    my $numColumnNames = @columnNames;
    my %h;
    while (defined(my $line = <FILE_DATA>)) {
        $line = eolChomp($line);
        $row++;
        if ($commentChar ne "") {
            while ($line =~ m/^$commentChar/) { $line = <FILE_DATA>; $line = eolChomp($line); $row++; }
        }
        $line =~ s/[\x0A\x0D]+//g;
        my @elements = split("\t", $line, -1);
        if (scalar(@elements) != $numColumnNames) { die("Number of data items on row $row not $numColumnNames"); }
        for (my $i = 0; $i < $numColumnNames; $i++) { push(@{$h{$columnNames[$i]}}, $elements[$i]); }
    }
    close(FILE_DATA);
    return(\%h, \@columnNames);
}

###############################################################################
# Output values to a tab-separated file, with a header line included.
# Three different ways of passing the data are supported:
# 1. An array of references to hashes, with each hash having the SAME keys,
#    which are the column headings of the tab-separated file.  The hash
#    values are the data written to the file.
# 2. A hash of references to arrays.  The hash keys are the column headings
#    of the tab-separated file, and each referenced array must be the same
#    length.  The arrays contain the data written to the file.
# 3. A hash of references to hashes.  The referenced hashes must all have
#    the same keys and these are the column headings, as in 1 above.  In
#    this case, there is no array to determine order of the rows that are
#    written, so the row order is determined by sorting the main hash keys.
# Arguments:
#   $fileName: name of file to generate.
#   $dataRef: a reference to an array of references to hashes OR a reference
#       to a hash of references to arrays OR a reference to a hash of
#       references to hashes.
#   $descRef: reference to a hash of descriptions, using the same or a subset
#       of the keys in the hash or hashes referenced by $dataRef.  The keys
#       of this hash determine which columns are output.  The values of the
#       keys determine the column order.  See below for more.
# The hash keys obtained from the %$descRef hash are sorted by the values in
# that hash and written to the output file header line as column names.
# The %$descRef values might be descriptions used in a different form of output
# of the hash, to explain the column meanings.  If those descriptions start
# with numbers or letters (1. 2. 3. or a. b. c. for example), then that
# provides a way to sort the keys = column headings in the output file into
# a desired order.  If no descriptions are needed, the %$descRef values can
# simply be numbers or letters (sort uses cmp, not <=>, however).
###############################################################################
sub writeTabFile {
    my $fileName = shift; 
    my $dataRef = shift;
    my $descRef = shift;
    open(FILE_DATA, ">" . $fileName) or die("Can't create file $fileName: $!");
    my @colKeys = keys(%$descRef);
    @colKeys = sort {$$descRef{$a} cmp $$descRef{$b}} @colKeys;
    print(FILE_DATA join("\t", @colKeys), "\n");
    # 1. Array of references to hashes
    if (ref($dataRef) ne "HASH") {
        foreach my $href (@$dataRef) {
            my %h = %$href;
            print(FILE_DATA join("\t", @h{@colKeys}), "\n");
        }
    # 2. Hash of references to arrays
    } elsif (defined($$dataRef{$colKeys[0]})) {
        foreach my $colKey (@colKeys) {
            if (!exists($$dataRef{$colKey})) { die("Column $colKey doesn't exist in \$dataRef"); }
        }
        my $len = scalar(@{$$dataRef{$colKeys[0]}});
        for (my $i = 0; $i < $len; $i++) {
            my $line = "";
            foreach my $colKey (@colKeys) {
                $line = $line . "\t" . $$dataRef{$colKey}->[$i];
            }
            $line = substr($line, 1); # Remove first tab.
            print(FILE_DATA $line, "\n");
        }
    # 3. Hash of references to hashes
    } else {
        my @rowKeys = sort { $a cmp $b } keys(%$dataRef);
        foreach my $rowKey (@rowKeys) {
            # There should be a one-line notation for this but I can't find it.
            # print(FILE_DATA join("\t", @{$$dataRef{$rowKey}}->{@colKeys}), "\n");  DOESN'T WORK
            my $line = "";
            foreach my $colKey (@colKeys) {
                $line = $line . "\t" . $$dataRef{$rowKey}->{$colKey};
            }
            $line = substr($line, 1); # Remove first tab.
            print(FILE_DATA $line, "\n");
        }
    }
    close(FILE_DATA);
}

###############################################################################
# Return a time string giving current date/time.
###############################################################################
sub getTime {
    my ($sec, $min, $hour, $mday, $mon, $year, $wday, $yday, $isdst) = localtime(time);
    my @month_abbr = qw(Jan Feb Mar Apr May Jun Jul Aug Sep Oct Nov Dec);
    my $time = sprintf "%4d-%s-%02d %02d:%02d:%02d", $year+1900, $month_abbr[$mon], $mday, $hour, $min, $sec;
    return($time);
}

###############################################################################
# Return the log base 10 of the argument.
###############################################################################
sub log10 {
    my $n = shift;
    return(log($n)/log(10));
}

###############################################################################
# Return the log base 2 of the argument.
###############################################################################
sub log2 {
    my $n = shift;
    return(log($n)/log(2));
}

###############################################################################
# Return a string representation of $v, rounded to $n significant digits.
# Arguments:
#   $v: number to convert to string.
#   $n: number of significant digits desired.
# Returns: string representation of $v in base 10.
###############################################################################
sub signif {
    my $v = shift;
    my $n = shift;
    return(sprintf("%.*g", $n, $v));
}

###############################################################################
# Create a log file.
# Arguments:
#   $logFileName: name of the log file.
###############################################################################
our $logFileName;
our $logFileHandle;
sub createLogFile {
    $logFileName = shift;
    open($logFileHandle, ">" . "$logFileName") or die("Can't create file $logFileName: $!");
}

###############################################################################
# Close the log file created by createLogFile.
###############################################################################
sub closeLogFile {
    close($logFileHandle);
}

###############################################################################
# Print arguments to standard output and to log file.
###############################################################################
sub printLog {
    my @args = @_;
    print("@args");
    print($logFileHandle "@args");
    $logFileHandle->flush;
    *STDOUT->flush;
}

###############################################################################
# Read the argument file as text and call printLog on each line, then
# delete the file.
###############################################################################
sub showFileAndLogOutput {
    my $filename = shift;
    open(DATAFILE, "<" . $filename) or die("Can't open '$filename': $!");
    my $line;
    while (defined($line = <DATAFILE>)) {
        printLog("$line"); # $line has \n
    }
    close(DATAFILE);
    unlink($filename);
}

###############################################################################
# Remove end-of-line characters from a string.
# Arguments:
#   $s: string whose end-of-line characters are to be removed.
# Returns: modified string.
# Note: this is like chomp(), but chomp() removes whatever $/ is set to, which
# is often "\n" and therefore it fails to remove "\r\n" off of lines that have
# it.  This function removes all \r and \n characters from the end of the line.
# Remember to assign the return value (unlike chomp, which doesn't require that).
###############################################################################
sub eolChomp {
    my $s = shift;
    my $i = 1;
    my $len = length($s);
    while ($i <= $len) {
        my $ch = substr($s, -$i, 1);
        if (($ch ne "\r") && ($ch ne "\n")) { last; }
        $i++;
    }
    return(substr($s, 0, $len + 1 - $i));
}

###############################################################################
# Commify a number provided as an argument.  Perl Cookbook, 2.17, p. 64.
# Arguments:
#   $N: number to be commified.
# Returns: commified number string.
###############################################################################
sub commify {
    my $N = shift;
    my $text = reverse $N;
    $text =~ s/(\d\d\d)(?=\d)(?!\d*\.)/$1,/g;
    return(scalar reverse $text);
}

###############################################################################
# Make and return a string containing a description followed by a value.
# Arguments:
#   $descr: the description
#   $val: the value
#   $percent100: value representing 100%.
#   $pctReplacement: replace % at end of description with this.
# Returns: pretty number string.
#
# If $descr ends with %, replace the % with $pctReplacement.  If $val is purely an
# integer, commify it.  If $descr ended with %, append in parentheses after $val its
# percent of $percent100.  Return the resulting description/value string.
###############################################################################
sub makeNumberPretty {
    my $descr = shift;
    my $val = shift;
    my $percent100 = shift;
    my $pctReplacement = shift;
    my $endedPct = ($descr =~ s/%$/$pctReplacement/);
    my $valStr = $val;
    if ($val =~ m/^[0-9]+$/) {
        $valStr = commify($val);
        if ($endedPct) {
            my $pct = 100 * $val / $percent100;
            $pct = sprintf("%.1f", $pct);
            $valStr = $valStr . " (" . $pct . "%)";
         }
    }
    return($descr . " => " . $valStr);
}

###############################################################################
# Apply interpolation to ~{names} in a string.
# Arguments:
#   $s: string to be tilde-interpolated, containing ~{name}'s.
#   %ithref: reference to hash containing keys equal to the ~names and
#       values being the replacement text for ~{name}.
# Returns: tilde-interpolated string.
###############################################################################
sub interpolateTildes {
    my $s = shift;
    my $ithref = shift;
    my %ith = %$ithref;
    foreach my $key (keys(%ith)) {
        $s =~ s/~{$key}/$ith{$key}/g;
    }
    return($s);
}

###############################################################################
# Search an array of strings for a particular string and return its index.
# Arguments:
#   $s: string to search for.
#   $arrRef: reference to array of strings to search.
# Returns: index in @arr of first occurrence of $s, or -1 if not found.
###############################################################################
sub findFirstMatch {
    my $s = shift;
    my $arrRef = shift;
    for (my $i = 0; $i < scalar(@$arrRef); $i++) {
        if ($s eq $$arrRef[$i]) { return($i); }
    }
    return(-1);
}

###############################################################################
# Search an array of numbers for a particular value and return its index.
# Arguments:
#   $N: value to search for.
#   $arrRef: reference to array of numbers to search.
# Returns: index in @arr of first occurrence of $N, or -1 if not found.
###############################################################################
sub indexOf {
    my $N = shift;
    my $arrRef = shift;
    for (my $i = 0; $i < scalar(@$arrRef); $i++) {
        if ($N == $$arrRef[$i]) { return($i); }
    }
    return(-1);
}

###############################################################################
# Apply the uc() function to an array of strings.
# Arguments:
#   $arrRef: reference to array of strings to which to apply uc().
# Returns: Nothing.  The strings are modified in-place in the array.
###############################################################################
sub ucArray {
    my $arrRef = shift;
    for (my $i = 0; $i < scalar(@$arrRef); $i++) { $$arrRef[$i] = uc($$arrRef[$i]); }
}

###############################################################################
# Apply the tr operator to an array of strings.
# Arguments:
#   $arrRef: reference to array of strings to which to apply tr.
#   $searchList: the list of symbols to search for.
#   $replaceList: the corresponding list of symbols to replace.
# Returns: Nothing.  The strings are modified in-place in the array.
###############################################################################
sub trArray {
    my $arrRef = shift;
    my $searchList = shift;
    my $replaceList = shift;
    for (my $i = 0; $i < scalar(@$arrRef); $i++) { eval("\$\$arrRef[$i] =~ tr/$searchList/$replaceList/"); }
}

###############################################################################
# Return the sum of an array of values.
# Arguments:
#   $arrRef: reference to array of values to sum.
# Returns: sum of values in the array.
###############################################################################
sub sumArray {
    my $arrRef = shift;
    my $sum = 0;
    foreach my $val (@$arrRef) { $sum += $val; }
    return($sum);
}

###############################################################################
# Return the number of TRUE values in an array.  In Perl, there is no boolean
# value per se; non-zero non-blank non-undefined is treated as TRUE.
# Arguments:
#   $arrRef: reference to array of boolean values to count.
# Returns: number of TRUE values in the array.
###############################################################################
sub countTrueArray {
    my $arrRef = shift;
    my $count = 0;
    foreach my $val (@$arrRef) { if ($val) { $count++; } }
    return($count);
}

###############################################################################
# Count the number of occurrences of each value in an array.  This emulates the
# R "table" function.
# Arguments:
#   $arrRef: reference to array of values to count.
# Returns: reference to hash whose keys are the array values and whose values
# are the number of times those array values occurred.
###############################################################################
sub R_table {
    my $arrRef = shift;
    my %h;
    foreach my $val (@$arrRef) { $h{$val}++; }
    return(\%h);
}

###############################################################################
# Return an array of the unique values found in one or more supplied arrays.
# This emulates the R "unique" function.
# Arguments:
#   $arrRef, ...: one or more references to arrays of values to find unique ones.
# Returns: reference to array of values unique across all supplied arrays.
###############################################################################
sub R_unique {
    my %h;
    while (scalar(@_) > 0) {
        my $arrRef = shift;
        foreach my $val (@$arrRef) { $h{$val} = 1; }
    }
    return(\keys(%h));
}

###############################################################################
# Return an array of the intersection of the values found in two supplied arrays.
# This emulates the R "intersect" function.
# Arguments:
#   $arrRef1, $arrRef2: references to array of values to find intersection.
# Returns: reference to array of values found in both @$arrRef1 and @$arrRef2.
###############################################################################
sub R_intersect {
    my $arrRef1 = shift;
    my $arrRef2 = shift;
    my %h;
    foreach my $val (@$arrRef1) { $h{$val} = 1; }
    my @a;
    foreach my $val (@$arrRef2) { if (exists($h{$val})) { push(@a, $val); } }
    return(\@a);
}

###############################################################################
# Return an array of the union of the values found in two supplied arrays.
# This emulates the R "union" function.
# Arguments:
#   $arrRef1, $arrRef2: references to array of values to find union.
# Returns: reference to array of values found in @$arrRef1 and/or @$arrRef2.
###############################################################################
sub R_union {
    my $arrRef1 = shift;
    my $arrRef2 = shift;
    my %h;
    foreach my $val (@$arrRef1) { $h{$val} = 1; }
    foreach my $val (@$arrRef2) { $h{$val} = 1; }
    return(\keys(%h));
}

###############################################################################
# Return an array of the set difference of the values found in two supplied arrays.
# This emulates the R "setdiff" function.
# Arguments:
#   $arrRef1, $arrRef2: references to array of values to find set difference.
# Returns: reference to array of values found in @$arrRef1 but not in @$arrRef2.
###############################################################################
sub R_setdiff {
    my $arrRef1 = shift;
    my $arrRef2 = shift;
    my %h;
    foreach my $val (@$arrRef1) { $h{$val} = 1; }
    foreach my $val (@$arrRef2) { delete($h{$val}); }
    return(\keys(%h));
}

###############################################################################
# Compute complement of a DNA or RNA sequence.
# Arguments:
#   $seqRef: sequence to complement, upper and/or lower case letters.
#   $rna: optional argument, if present with 1 value it signals that $seq
#       is RNA rather than DNA.  By default DNA or RNA are identified
#       automatically by searching $seq for T or U.  If it contains both,
#       it is treated as DNA although the complement of U is taken as A.
# Returns: complemented sequence.  Lower and upper case are preserved.
#
# IUPAC symbols are handled as follows:
#
# IUPAC Base(s)     Complement
# ----- -------     ----------
#   A   Adenine     T
#   T   Thymine     A
#   U   Uracil      A
#   C   Cytosine    G
#   G   Guanine     C
#   R   A or G      Y
#   Y   C or T/U    R
#   S   G or C      W
#   W   A or T/U    S
#   K   G or T/U    M
#   M   A or C      K
#   B   C, G or T/U V
#   V   A, C or G   B
#   D   A, G or T/U H
#   H   A, C or T/U D
#   N   any base    N
#   .   gap         .
#   -   gap         -
#   All other symbols are left unchanged.
###############################################################################
sub complement {
    my $seq = shift;
    my $isRNA = ($seq =~ m/U/) || ((scalar(@_) > 0) && ($_[0] == 1));
    if (!$isRNA) {
        $seq =~ tr/ATUCGRYSWKMBVDHNatucgryswkmbvdhn/TAAGCYRWSMKVBHDNtaagcyrwsmkvbhdn/;
    } else {
        $seq =~ tr/ATUCGRYSWKMBVDHNatucgryswkmbvdhn/UAAGCYRWSMKVBHDNuaagcyrwsmkvbhdn/;
    }
    return($seq);
}

###############################################################################
# Compute reverse complement of a DNA or RNA sequence.
# Arguments:
#   $seq: sequence to reverse complement, upper and/or lower case letters.
#   $rna: optional argument, if present with 1 value it signals that $seq
#       is RNA rather than DNA.  By default DNA or RNA are identified
#       automatically by searching $seq for T or U.  If it contains both,
#       it is treated as DNA although the complement of U is taken as A.
# Returns: reverse-complemented sequence.  Lower and upper case are preserved.
###############################################################################
sub reverseComplement {
    my $seq = shift;
    my $isRNA = ($seq =~ m/U/) || ((scalar(@_) > 0) && ($_[0] == 1));
    $seq = reverse($seq);
    $seq = complement($seq, $isRNA);
    return($seq);
}

###############################################################################
# Compute the entropy of a symbol string.
# Arguments:
#   $seq: symbol string whose entropy is to be computed.  Upper and lower
#       case symbols are considered distinct.
# Returns: computed entropy, in bits.
###############################################################################
sub entropy {
    my $seq = shift;
    my @s = split("", $seq);
    my %counts;
    my $sum = 0;
    foreach my $ch (@s) { $counts{$ch}++; $sum++; }
    my $H = 0;
    foreach my $ch (keys(%counts)) { my $p = $counts{$ch}/$sum; $H += -$p*log2($p); }
    return($H);
}

###############################################################################
# Compute the Kullback–Leibler divergence of a pair of discrete distributions.
# Arguments:
#   $P: reference to a hash representing a discrete distribution.  The
#       values of the hash must all be non-negative, and when divided by
#       the sum of the values, represent probabilities.
#   $Q: reference to a second hash representing a discrete distribution.
#       Again, the values of the hash must all be non-negative, and when
#       divided by the sum of the values, represent probabilities.  Also,
#       a $Q hash value is only allowed to be 0 when the corresponding $P
#       value is also 0, but this function ignores all entries where this
#       is violated.
# Returns: a list ($D, $error):
#   $D: computed K-L divergence, in bits.
#   $error: string describing error, or empty string if no error.
#
# The $P and $Q hashes represent distributions over the same probability
# sample space, and must use the same keys, except when one does not use
# a key because its value for that key is 0.  Note that the K-L divergence
# is undefined when there is any sample space element with probability 0 in
# $Q but non-zero in $P.  We return an error when that happens, but ignore
# the fact and compute $D anyway.
#
# The K-L divergence is NOT a distance, because it is non-symmetric; the
# divergence of $P to $Q is different than from $Q to $P, generally.
###############################################################################
sub KLdivergence {
    my $P = shift;
    my $Q = shift;

    # Check values, sum them, and compute part of $D.
    my $sumP = 0;
    my $sumQ = 0;
    my $D = 0;
    my $numBad = 0;
    foreach my $key (keys(%$P)) {
        my $p = $$P{$key};
        my $q = defined($$Q{$key}) ? $$Q{$key} : 0;
        if ($p < 0 || $q < 0) { return(0, "Negative element(s)"); }
        if ($p > 0) {
            if ($q == 0) { $numBad++; }
            else {
                $sumP += $p;
                $sumQ += $q;
                $D += $p * log2($p/$q);
            }
        }
    }
    if ($sumP == 0) { return(0, "\$P sums to 0"); }
    my $error = "";
    if ($numBad > 0) { $error = "\$Q element is 0 but \$P is not: $numBad times"; }

    # Final adjustment is done.  Do the math.  This is necessary because we chose to
    # compute $D in the loop above rather than waiting until we had $sumP and $sumQ and
    # then looping a second time and computing $D by dividing $p and $q by those sums.
    $D = $D / $sumP + log2($sumQ / $sumP);
    return($D, $error);
}

###############################################################################
# Genetic code (from Wikipedia):
#                                       2nd base
#                   U               C               A               G
#   1st base
#       U       UUU:Phe/F       UCU:Ser/S       UAU:Tyr/Y       UGU:Cys/C
#       U       UUC:Phe/F       UCC:Ser/S       UAC:Tyr/Y       UGC:Cys/C
#       U       UUA:Leu/L       UCA:Ser/S       UAA:Ochre Stop  UGA:Opal Stop
#       U       UUG:Leu/L       UCG:Ser/S       UAG:Amber Stop  UGG:Trp/W    
#       C       CUU:Leu/L       CCU:Pro/P       CAU:His/H       CGU:Arg/R
#       C       CUC:Leu/L       CCC:Pro/P       CAC:His/H       CGC:Arg/R
#       C       CUA:Leu/L       CCA:Pro/P       CAA:Gln/Q       CGA:Arg/R
#       C       CUG:Leu/L       CCG:Pro/P       CAG:Gln/Q       CGG:Arg/R
#       A       AUU:Ile/I       ACU:Thr/T       AAU:Asn/N       AGU:Ser/S
#       A       AUC:Ile/I       ACC:Thr/T       AAC:Asn/N       AGC:Ser/S
#       A       AUA:Ile/I       ACA:Thr/T       AAA:Lys/K       AGA:Arg/R
#       A       AUG:Met/M*      ACG:Thr/T       AAG:Lys/K       AGG:Arg/R
#       G       GUU:Val/V       GCU:Ala/A       GAU:Asp/D       GGU:Gly/G
#       G       GUC:Val/V       GCC:Ala/A       GAC:Asp/D       GGC:Gly/G
#       G       GUA:Val/V       GCA:Ala/A       GAA:Glu/E       GGA:Gly/G
#       G       GUG:Val/V       GCG:Ala/A       GAG:Glu/E       GGG:Gly/G
# *AUG both codes for methionine and serves as an initiation site.
#
# Amino acids (from Wikipedia):
#
# Amino Acid    3   1   polarity    charge (pH 7.4)     hydropathy index
# ------------- -----   --------    ---------------     ----------------
# Alanine       Ala A   nonpolar    neutral                 1.8
# Arginine      Arg R   polar       positive                −4.5
# Asparagine    Asn N   polar       neutral                 −3.5
# Aspartic acid Asp D   polar       negative                −3.5
# Cysteine      Cys C   polar       neutral                 2.5
# Glutamic acid Glu E   polar       negative                −3.5
# Glutamine     Gln Q   polar       neutral                 −3.5
# Glycine       Gly G   nonpolar    neutral                 −0.4
# Histidine     His H   polar       pos(10%), neut(90%)     -3.2
# Isoleucine    Ile I   nonpolar    neutral                 4.5
# Leucine       Leu L   nonpolar    neutral                 3.8
# Lysine        Lys K   polar       positive                −3.9
# Methionine    Met M   nonpolar    neutral                 1.9
# Phenylalanine Phe F   nonpolar    neutral                 2.8
# Proline       Pro P   nonpolar    neutral                 −1.6
# Serine        Ser S   polar       neutral                 −0.8
# Threonine     Thr T   polar       neutral                 −0.7
# Tryptophan    Trp W   nonpolar    neutral                 −0.9
# Tyrosine      Tyr Y   polar       neutral                 −1.3
# Valine        Val V   nonpolar    neutral                 4.2
###############################################################################

###############################################################################
# Export a hash whose keys are RNA nucleotide triplets (in upper case) and
# whose values are the corresponding single-letter amino acid code (in
# upper case) translated using the standard genetic code.
###############################################################################
our %nt_to_aa = (
    UUU => "F",     UCU => "S",     UAU => "Y",     UGU => "C",
    UUC => "F",     UCC => "S",     UAC => "Y",     UGC => "C",
    UUA => "L",     UCA => "S",     UAA => "*",     UGA => "*",
    UUG => "L",     UCG => "S",     UAG => "*",     UGG => "W",    
    CUU => "L",     CCU => "P",     CAU => "H",     CGU => "R",
    CUC => "L",     CCC => "P",     CAC => "H",     CGC => "R",
    CUA => "L",     CCA => "P",     CAA => "Q",     CGA => "R",
    CUG => "L",     CCG => "P",     CAG => "Q",     CGG => "R",
    AUU => "I",     ACU => "T",     AAU => "N",     AGU => "S",
    AUC => "I",     ACC => "T",     AAC => "N",     AGC => "S",
    AUA => "I",     ACA => "T",     AAA => "K",     AGA => "R",
    AUG => "M",     ACG => "T",     AAG => "K",     AGG => "R",
    GUU => "V",     GCU => "A",     GAU => "D",     GGU => "G",
    GUC => "V",     GCC => "A",     GAC => "D",     GGC => "G",
    GUA => "V",     GCA => "A",     GAA => "E",     GGA => "G",
    GUG => "V",     GCG => "A",     GAG => "E",     GGG => "G"
    );

###############################################################################
# Export hashes whose keys are amino acid full names, 3-letter names, or
# 1-letter names, and whose values are references to hashes containing
# amino acid information (all sub-hashes have the same keys).
###############################################################################
our %aaInfo1 = (
A=>{name=>"alanine",       name3=>"ALA", name1=>"A", polarity=>"nonpolar", charge=>"neutral",  hydropathy=> 1.8},
R=>{name=>"arginine",      name3=>"ARG", name1=>"R", polarity=>"polar",    charge=>"positive", hydropathy=> -4.5},
N=>{name=>"asparagine",    name3=>"ASN", name1=>"N", polarity=>"polar",    charge=>"neutral",  hydropathy=> -3.5},
D=>{name=>"aspartic acid", name3=>"ASP", name1=>"D", polarity=>"polar",    charge=>"negative", hydropathy=> -3.5},
C=>{name=>"cysteine",      name3=>"CYS", name1=>"C", polarity=>"polar",    charge=>"neutral",  hydropathy=> 2.5},
E=>{name=>"glutamic acid", name3=>"GLU", name1=>"E", polarity=>"polar",    charge=>"negative", hydropathy=> -3.5},
Q=>{name=>"glutamine",     name3=>"GLN", name1=>"Q", polarity=>"polar",    charge=>"neutral",  hydropathy=> -3.5},
G=>{name=>"glycine",       name3=>"GLY", name1=>"G", polarity=>"nonpolar", charge=>"neutral",  hydropathy=> -0.4},
H=>{name=>"histidine",     name3=>"HIS", name1=>"H", polarity=>"polar",    charge=>"p10n90",   hydropathy=> -3.2},
I=>{name=>"isoleucine",    name3=>"ILE", name1=>"I", polarity=>"nonpolar", charge=>"neutral",  hydropathy=> 4.5},
L=>{name=>"leucine",       name3=>"LEU", name1=>"L", polarity=>"nonpolar", charge=>"neutral",  hydropathy=> 3.8},
K=>{name=>"lysine",        name3=>"LYS", name1=>"K", polarity=>"polar",    charge=>"positive", hydropathy=> -3.9},
M=>{name=>"methionine",    name3=>"MET", name1=>"M", polarity=>"nonpolar", charge=>"neutral",  hydropathy=> 1.9},
F=>{name=>"phenylalanine", name3=>"PHE", name1=>"F", polarity=>"nonpolar", charge=>"neutral",  hydropathy=> 2.8},
P=>{name=>"proline",       name3=>"PRO", name1=>"P", polarity=>"nonpolar", charge=>"neutral",  hydropathy=> -1.6},
S=>{name=>"serine",        name3=>"SER", name1=>"S", polarity=>"polar",    charge=>"neutral",  hydropathy=> -0.8},
T=>{name=>"threonine",     name3=>"THR", name1=>"T", polarity=>"polar",    charge=>"neutral",  hydropathy=> -0.7},
W=>{name=>"tryptophan",    name3=>"TRP", name1=>"W", polarity=>"nonpolar", charge=>"neutral",  hydropathy=> -0.9},
Y=>{name=>"tyrosine",      name3=>"TYR", name1=>"Y", polarity=>"polar",    charge=>"neutral",  hydropathy=> -1.3},
V=>{name=>"valine",        name3=>"VAL", name1=>"V", polarity=>"nonpolar", charge=>"neutral",  hydropathy=> 4.2},
"*"=>{name=>"STOP",        name3=>"STP", name1=>"*", polarity=>"N/A",      charge=>"N/A",      hydropathy=> "N/A"}
    );

our %aaInfo3 = (
    ALA => $aaInfo1{A}, ARG => $aaInfo1{R}, ASN => $aaInfo1{N}, ASP => $aaInfo1{D},
    CYS => $aaInfo1{C}, GLU => $aaInfo1{E}, GLN => $aaInfo1{Q}, GLY => $aaInfo1{G},
    HIS => $aaInfo1{H}, ILE => $aaInfo1{I}, LEU => $aaInfo1{L}, LYS => $aaInfo1{K},
    MET => $aaInfo1{M}, PHE => $aaInfo1{F}, PRO => $aaInfo1{P}, SER => $aaInfo1{S},
    THR => $aaInfo1{T}, TRP => $aaInfo1{W}, TYR => $aaInfo1{Y}, VAL => $aaInfo1{V},
    STP => $aaInfo1{"*"}
    );

our %aaInfo = (
    alanine => $aaInfo1{A},         arginine => $aaInfo1{R},        asparagine => $aaInfo1{N}, 
    aspartic => $aaInfo1{D},        cysteine => $aaInfo1{C},        glutamic => $aaInfo1{E}, 
    glutamine => $aaInfo1{Q},       glycine => $aaInfo1{G},         histidine => $aaInfo1{H}, 
    isoleucine => $aaInfo1{I},      leucine => $aaInfo1{L},         lysine => $aaInfo1{K}, 
    methionine => $aaInfo1{M},      phenylalanine => $aaInfo1{F},   proline => $aaInfo1{P}, 
    serine => $aaInfo1{S},          threonine => $aaInfo1{T},       tryptophan => $aaInfo1{W}, 
    tyrosine => $aaInfo1{Y},        valine => $aaInfo1{V},          STOP => $aaInfo1{"*"}
    );

###############################################################################
# Export two-level hash containing the PAM10 amino acid substitution matrix.
# This matrix was downloaded from the NCBI website, and included the following
# comments:
#
#   This matrix was produced by "pam" Version 1.0.6 [28-Jul-93]
#   PAM 10 substitution matrix, scale = ln(2)/2 = 0.346574
#   Expected score = -8.27, Entropy = 3.43 bits
#   Lowest score = -23, Highest score = 13
#
# Here we exclude the special codes B, Z, and X.
###############################################################################
our %PAM10 = (
  A => {A=>  7, R=>-10, N=> -7, D=> -6, C=>-10, Q=> -7, E=> -5, G=> -4, H=>-11, I=> -8,
        L=> -9, K=>-10, M=> -8, F=>-12, P=> -4, S=> -3, T=> -3, W=>-20, Y=>-11, V=> -5},
  R => {A=>-10, R=>  9, N=> -9, D=>-17, C=>-11, Q=> -4, E=>-15, G=>-13, H=> -4, I=> -8,
        L=>-12, K=> -2, M=> -7, F=>-12, P=> -7, S=> -6, T=>-10, W=> -5, Y=>-14, V=>-11},
  N => {A=> -7, R=> -9, N=>  9, D=> -1, C=>-17, Q=> -7, E=> -5, G=> -6, H=> -2, I=> -8,
        L=>-10, K=> -4, M=>-15, F=>-12, P=> -9, S=> -2, T=> -5, W=>-11, Y=> -7, V=>-12},
  D => {A=> -6, R=>-17, N=> -1, D=>  8, C=>-21, Q=> -6, E=>  0, G=> -6, H=> -7, I=>-11,
        L=>-19, K=> -8, M=>-17, F=>-21, P=>-12, S=> -7, T=> -8, W=>-21, Y=>-17, V=>-11},
  C => {A=>-10, R=>-11, N=>-17, D=>-21, C=> 10, Q=>-20, E=>-20, G=>-13, H=>-10, I=> -9,
        L=>-21, K=>-20, M=>-20, F=>-19, P=>-11, S=> -6, T=>-11, W=>-22, Y=> -7, V=> -9},
  Q => {A=> -7, R=> -4, N=> -7, D=> -6, C=>-20, Q=>  9, E=> -1, G=>-10, H=> -2, I=>-11,
        L=> -8, K=> -6, M=> -7, F=>-19, P=> -6, S=> -8, T=> -9, W=>-19, Y=>-18, V=>-10},
  E => {A=> -5, R=>-15, N=> -5, D=>  0, C=>-20, Q=> -1, E=>  8, G=> -7, H=> -9, I=> -8,
        L=>-13, K=> -7, M=>-10, F=>-20, P=> -9, S=> -7, T=> -9, W=>-23, Y=>-11, V=>-10},
  G => {A=> -4, R=>-13, N=> -6, D=> -6, C=>-13, Q=>-10, E=> -7, G=>  7, H=>-13, I=>-17,
        L=>-14, K=>-10, M=>-12, F=>-12, P=>-10, S=> -4, T=>-10, W=>-21, Y=>-20, V=> -9},
  H => {A=>-11, R=> -4, N=> -2, D=> -7, C=>-10, Q=> -2, E=> -9, G=>-13, H=> 10, I=>-13,
        L=> -9, K=>-10, M=>-17, F=> -9, P=> -7, S=> -9, T=>-11, W=>-10, Y=> -6, V=> -9},
  I => {A=> -8, R=> -8, N=> -8, D=>-11, C=> -9, Q=>-11, E=> -8, G=>-17, H=>-13, I=>  9,
        L=> -4, K=> -9, M=> -3, F=> -5, P=>-12, S=>-10, T=> -5, W=>-20, Y=> -9, V=> -1},
  L => {A=> -9, R=>-12, N=>-10, D=>-19, C=>-21, Q=> -8, E=>-13, G=>-14, H=> -9, I=> -4,
        L=>  7, K=>-11, M=> -2, F=> -5, P=>-10, S=>-12, T=>-10, W=> -9, Y=>-10, V=> -5},
  K => {A=>-10, R=> -2, N=> -4, D=> -8, C=>-20, Q=> -6, E=> -7, G=>-10, H=>-10, I=> -9,
        L=>-11, K=>  7, M=> -4, F=>-20, P=>-10, S=> -7, T=> -6, W=>-18, Y=>-12, V=>-13},
  M => {A=> -8, R=> -7, N=>-15, D=>-17, C=>-20, Q=> -7, E=>-10, G=>-12, H=>-17, I=> -3,
        L=> -2, K=> -4, M=> 12, F=> -7, P=>-11, S=> -8, T=> -7, W=>-19, Y=>-17, V=> -4},
  F => {A=>-12, R=>-12, N=>-12, D=>-21, C=>-19, Q=>-19, E=>-20, G=>-12, H=> -9, I=> -5,
        L=> -5, K=>-20, M=> -7, F=>  9, P=>-13, S=> -9, T=>-12, W=> -7, Y=> -1, V=>-12},
  P => {A=> -4, R=> -7, N=> -9, D=>-12, C=>-11, Q=> -6, E=> -9, G=>-10, H=> -7, I=>-12,
        L=>-10, K=>-10, M=>-11, F=>-13, P=>  8, S=> -4, T=> -7, W=>-20, Y=>-20, V=> -9},
  S => {A=> -3, R=> -6, N=> -2, D=> -7, C=> -6, Q=> -8, E=> -7, G=> -4, H=> -9, I=>-10,
        L=>-12, K=> -7, M=> -8, F=> -9, P=> -4, S=>  7, T=> -2, W=> -8, Y=>-10, V=>-10},
  T => {A=> -3, R=>-10, N=> -5, D=> -8, C=>-11, Q=> -9, E=> -9, G=>-10, H=>-11, I=> -5,
        L=>-10, K=> -6, M=> -7, F=>-12, P=> -7, S=> -2, T=>  8, W=>-19, Y=> -9, V=> -6},
  W => {A=>-20, R=> -5, N=>-11, D=>-21, C=>-22, Q=>-19, E=>-23, G=>-21, H=>-10, I=>-20,
        L=> -9, K=>-18, M=>-19, F=> -7, P=>-20, S=> -8, T=>-19, W=> 13, Y=> -8, V=>-22},
  Y => {A=>-11, R=>-14, N=> -7, D=>-17, C=> -7, Q=>-18, E=>-11, G=>-20, H=> -6, I=> -9,
        L=>-10, K=>-12, M=>-17, F=> -1, P=>-20, S=>-10, T=> -9, W=> -8, Y=> 10, V=>-10},
  V => {A=> -5, R=>-11, N=>-12, D=>-11, C=> -9, Q=>-10, E=>-10, G=> -9, H=> -9, I=> -1,
        L=> -5, K=>-13, M=> -4, F=>-12, P=> -9, S=>-10, T=> -6, W=>-22, Y=>-10, V=>  8}
    );

###############################################################################
# Translate a nucleotide sequence from DNA or RNA to amino acids using the
# standard genetic code.
# Arguments:
#   $seq: sequence to translate, may contain T or U (or both, they are treated identically).
#       Upper/lower case are also treated identically.  
# Returns: translated sequence consisting of single-letter amino acid codes in upper
# case.  STOP codons are translated to *, and non-ATUCG letters are translated to X.
# If $seq length is not a multiple of 3, one or two extra letters at the end are IGNORED.
###############################################################################
sub ntToAA {
    my $seq = shift;
    $seq = uc($seq);
    $seq =~ tr/T/U/;
    my $numCodons = int(length($seq) / 3);
    my $AAseq = "";
    for (my $i = 0; $i < $numCodons; $i++) {
        my $codon = substr($seq, $i*3, 3);
        if (defined($nt_to_aa{$codon})) { $AAseq .= $nt_to_aa{$codon}; }
        else { $AAseq .= "X"; }
    }
    return($AAseq);
}

###############################################################################
# Translate a string of 1-letter amino acid codes into an array of 3-letter
# amino acid codes.
# Arguments:
#   $seq1: amino acid letter sequence to translate.  Upper/lower case are identical.  
# Returns: reference to array of 3-letter amino acid codes, or if $seq1 is
# one letter, a string is returned.
###############################################################################
sub AA1_to_AA3 {
    my $seq1 = shift;
    $seq1 = uc($seq1);
    my @s = split("", $seq1);
    if (scalar(@s) == 1) { return($aaInfo1{$seq1}{name3}); }
    my @names;
    foreach my $AA1 (@s) { push(@names, $aaInfo1{$AA1}{name3}); }
    return(\@names);
}

###############################################################################
# Just like AA1_to_AA3 except result is full AA names.
###############################################################################
sub AA1_ToNames {
    my $seq1 = shift;
    $seq1 = uc($seq1);
    my @s = split("", $seq1);
    if (scalar(@s) == 1) { return($aaInfo1{$seq1}{name}); }
    my @names;
    foreach my $AA1 (@s) { push(@names, $aaInfo1{$AA1}{name}); }
    return(\@names);
}

###############################################################################
# Translate an array of 3-letter amino acid codes into a string of 1-letter
# amino acid codes.
# Arguments:
#   $seq: reference to array of 3-letter amino acid names to translate.
#       Upper/lower case are identical.  
# Returns: string of corresponding 1-letter amino acid codes.
###############################################################################
sub AA3_to_AA1 {
    my $seq3 = shift;
    my $AA1seq = "";
    foreach my $AA3 (@$seq3) { $AA1seq .= $aaInfo3{uc($AA3)}{name1}; }
    return($AA1seq);
}

###############################################################################
# Just like AA1_ToNames except argument is reference to array of 3-letter
# amino acid codes.
###############################################################################
sub AA3_ToNames {
    my $seq3 = shift;
    my @names;
    foreach my $AA3 (@$seq3) { push(@names, $aaInfo3{uc($AA3)}{name}); }
    return(\@names);
}

###############################################################################
# Get the polarity of the amino acids whose 1-letter codes are given by the
# argument sequence.
# Arguments:
#   $seq1: amino acid letter sequence.  Upper/lower case are identical.  
# Returns: reference to array of amino acid polarities, or if $seq1 is
# one letter, a single polarity value is returned.  The polarity values
# are either "polar" or "nonpolar".
###############################################################################
sub AA1_GetPolarity {
    my $seq1 = shift;
    $seq1 = uc($seq1);
    my @s = split("", $seq1);
    if (scalar(@s) == 1) { return($aaInfo1{$seq1}{polarity}); }
    my @polarities;
    foreach my $AA1 (@s) { push(@polarities, $aaInfo1{$AA1}{polarity}); }
    return(\@polarities);
}

###############################################################################
# Get the charge of the amino acids whose 1-letter codes are given by the
# argument sequence.
# Arguments:
#   $seq1: amino acid letter sequence.  Upper/lower case are identical.  
# Returns: reference to array of amino acid charges at pH 7.4, or if $seq1
# is one letter, a single charge value is returned.  The charge values are
# either "positive", "negative", "neutral", or, for histidine, "p10n90"
# indicating 10% positive and 90% negative.
###############################################################################
sub AA1_GetCharge {
    my $seq1 = shift;
    $seq1 = uc($seq1);
    my @s = split("", $seq1);
    if (scalar(@s) == 1) { return($aaInfo1{$seq1}{charge}); }
    my @charges;
    foreach my $AA1 (@s) { push(@charges, $aaInfo1{$AA1}{charge}); }
    return(\@charges);
}

###############################################################################
# Get the hydropathy (aka hydrophobicity) index of the amino acids whose
# 1-letter codes are given by the argument sequence.
# Arguments:
#   $seq1: amino acid letter sequence.  Upper/lower case are identical.  
# Returns: reference to array of amino acid hydropathy indices, or if $seq1
# is one letter, a single hydropathy index is returned.  The hydropathy
# indices are floating point numbers.
#
# Note: See Kyte J, Doolittle RF (May 1982). "A simple method for
#       displaying the hydropathic character of a protein".
# Journal of Molecular Biology 157 (1): 105–32.
# doi:10.1016/0022-2836(82)90515-0. PMID 7108955
###############################################################################
sub AA1_GetHydropathy {
    my $seq1 = shift;
    $seq1 = uc($seq1);
    my @s = split("", $seq1);
    if (scalar(@s) == 1) { return($aaInfo1{$seq1}{hydropathy}); }
    my @hydropathies;
    foreach my $AA1 (@s) { push(@hydropathies, $aaInfo1{$AA1}{hydropathy}); }
    return(\@hydropathies);
}

###############################################################################
# Compute averages of an array of values across sliding windows.
# Arguments:
#   $values: reference to array of values to average across sliding windows
#       of indexes into the array.
#   $winSize: window size, i.e. number of values in $values array to average
#       together to produce a single output value.  Must be <= size of the
#       $values array.
#   $winDelta: amount to move the window between each averaged output value.
#       A value of 1 means the number of output values will be N-$winSize+1,
#       where N is the size of the $values array.  The number of output
#       values will be int((N-$winSize+$winDelta)/$winDelta).  Must be >= 1.
# Returns: reference to array of averages of $values over sliding windows.
# The first array entry is the average of $values[0]..@
###############################################################################
sub averageInWindow {
    my $values = shift;
    my $winSize = shift;
    my $winDelta = shift;
    my $N = @$values;
    if ($winSize > $N) { die("averageInWindow requires $winSize <= # of values"); }
    $winDelta = int($winDelta);
    if ($winDelta < 1) { die("averageInWindow requires $winDelta >= 1"); }
    my @averages;

    # Compute average for the first window.
    my $sum = 0;
    for (my $j = 0; $j < $winSize; $j++) { $sum += $values->[$j]; }
    push(@averages, $sum/$winSize);

    # For subsequent windows, repeating the above is efficient if 2*$winDelta >= $winSize,
    # but otherwise, it is more efficient to simply update the previous window sum by
    # removing the first $winDelta values and adding $winDelta new values to the end.
    if (2*$winDelta >= $winSize) {
        for (my $i = $winDelta; $i + $winSize <= $N; $i += $winDelta) {
            $sum = 0;
            for (my $j = 0; $j < $winSize; $j++) { $sum += $values->[$i+$j]; }
            push(@averages, $sum/$winSize);
        }
    } else {
        for (my $i = $winDelta; $i + $winSize <= $N; $i += $winDelta) {
            for (my $j = 1; $j <= $winDelta; $j++) {
                $sum += $values->[$i+$winSize-$j] - $values->[$i-$j];
            }
            push(@averages, $sum/$winSize);
        }
    }
    return(\@averages);
}

###############################################################################
# Search a nucleotide sequence from DNA or RNA to find the ORFs in it,
# within a specified reading frame.
# Arguments:
#   $seq: sequence to search, may contain T or U (or both, they are treated
#       identically).  Upper/lower case are also treated identically.
#   $frame: 0, 1, or 2, indicating which of the three reading frames is to
#       be searched.
#   $START: (optional) starting codon.  If not specified, "AUG" is used.
#       May contain T or U. 
# Returns: reference to array of ORF starting locations, using 0-based
# indexes into $seq to identify the beginning of an ORF (the AUG codon).
# It will be true that every returned array element modulo 3 == $frame.
###############################################################################
sub findORF {
    my $seq = shift;
    my $frame = shift;
    my $START = (scalar(@_) > 0) ? shift : "AUG";
    if ($frame > 0) { $seq = substr($seq, $frame); }
    $seq = uc($seq);
    $seq =~ tr/T/U/;
    $START = uc($START);
    $START =~ tr/T/U/;
    my $numCodons = int(length($seq) / 3);
    my @ORFs = ();
    for (my $i = 0; $i < $numCodons; $i++) {
        my $codon = substr($seq, $i*3, 3);
        if ($codon eq $START) { push(@ORFs, $i*3+$frame); }
    }
    return(\@ORFs);
}

###############################################################################
# Search a nucleotide sequence that starts with an ORF (e.g. with AUG) to
# locate the ORF stop codon.  Only the reading frame of the starting codon
# is searched.
# Arguments:
#   $seq: sequence to search, may contain T or U (or both, they are treated
#       identically).  Upper/lower case are also treated identically.
# Returns: length of ORF in codons, including the starting codon and ending
# stop codon, or 0 if stop codon not found.  A return value of 1 does not
# occur.
###############################################################################
sub lenORF {
    my $seq = shift;
    $seq = uc($seq);
    $seq =~ tr/T/U/;
    my $numCodons = int(length($seq) / 3);
    for (my $i = 1; $i < $numCodons; $i++) {
        my $codon = substr($seq, $i*3, 3);
        if (defined($nt_to_aa{$codon}) && $nt_to_aa{$codon} eq "*") { return($i+1); }
    }
    return(0);
}

###############################################################################
# Count codons and first-codons in a nucleotide coding sequence.
# Arguments:
#   $seq: sequence to be counted, may contain T or U (or both, they are
#       treated identically).  Upper/lower case are also treated identically.
#   $codonCountsRef: reference to hash of codon counters, key is codon and
#       value is count of number of times that codon was seen.
#   $firstCodonCountsRef: reference to another hash of codon counters, key
#       is FIRST codon in the sequence and value is count of times it was seen.
# Returns: Total number of codons found in the sequence.  Hashes %$codonCountsRef
# and %$firstCodonCountsRef have updated counts, preserving previous counts
# and adding to them.
###############################################################################
sub countCodons {
    my $seq = shift;
    my $codonCountsRef = shift;
    my $firstCodonCountsRef = shift;
    $seq = uc($seq);
    $seq =~ tr/T/U/;
    my $numCodons = int(length($seq) / 3);
    if ($numCodons > 0) {
        my $codon = substr($seq, 0, 3);
        $$firstCodonCountsRef{$codon}++;
    }       
    for (my $i = 0; $i < $numCodons; $i++) {
        my $codon = substr($seq, $i*3, 3);
        $$codonCountsRef{$codon}++;
    }
    return($numCodons);
}

###############################################################################
# Count codons and first-codons in an array of nucleotide coding sequences,
# using countCodons() above for each sequence in the array.
# Arguments:
#   $arrRef: reference to array of sequences to be counted.
#   $codonCountsRef: reference to hash of codon counters, key is codon and
#       value is count of number of times that codon was seen.
#   $firstCodonCountsRef: reference to another hash of codon counters, key
#       is FIRST codon in the sequence and value is count of times it was seen.
# Returns: Total number of codons found in the sequence.  Hashes %$codonCountsRef
# and %$firstCodonCountsRef have updated counts, preserving previous counts
# and adding to them.
###############################################################################
sub countCodons_Array {
    my $arrRef = shift;
    my $codonCountsRef = shift;
    my $firstCodonCountsRef = shift;
    my $numCodons = 0;
    foreach my $seq (@$arrRef) { $numCodons += countCodons($seq, $codonCountsRef, $firstCodonCountsRef); }
    return($numCodons);
}

###############################################################################
# Check a nucleotide coding sequence to see if it appears to be a good
# coding sequence.
# Arguments:
#   $seq: sequence to be checked, may contain T or U (or both, they are
#       treated identically).  Upper/lower case are also treated identically.
#   $STARTsRef: (optional) reference to array of permissible start codons.
#       If unspecified, "AUG" is used.  May contain T or U.
# Returns: empty string if no errors found, else a string describing the
# first problem noted.
###############################################################################
sub checkGoodCodingSeq {
    my $seq = shift;
    my $STARTsRef = (scalar(@_) > 0) ? shift : ["AUG"];
    $seq = uc($seq);
    $seq =~ tr/T/U/;
    my @STARTs = @$STARTsRef;
    ucArray(\@STARTs);
    trArray(\@STARTs, "T", "U");
    my $len = length($seq);
    if ($len == 0) { return("Sequence has 0 length"); }
    if (($len % 3) != 0) { return("Length not multiple of 3"); }
    my $firstCodon = substr($seq, 0, 3);
    if (findFirstMatch($firstCodon, \@STARTs) == -1) { return("First codon not a start codon: $firstCodon"); }
    my $lastCodon = substr($seq, $len-3);
    if (findFirstMatch($lastCodon, ["UGA","UAG","UAA"]) == -1) { return("Last codon not a stop codon: $lastCodon"); }
    return("");
}

###############################################################################
# Read the genome sequence from a GenBank file.
# Arguments:
#   $fileName: name of GenBank file to read.
# Returns: a list ($genome, $error):
#   $genome: A single long string that is the genome sequence, empty if error.
#   $error: An error string if genome could not be read, else empty string.
###############################################################################
sub readGenBankFile_genomeSeq {
    my $fileName = shift;
    if (!open(IN, "<$fileName")) { return("", "Error opening $fileName for reading"); }
    my $line;
    my $lineNum = 0;
    my $foundSeq = 0;
    while (!$foundSeq and my $line = <IN>) { $foundSeq = ($line =~ m/^ORIGIN/); }
    if (!$foundSeq) { return("", "Genome sequence not found in $fileName"); }
    my $genome = "";
    while (my $line = <IN>) {
        $lineNum++;
        $line = eolChomp($line);
        last if ($line !~ m/^ *\d/);
        my $dnaSeq = uc($line);
        $dnaSeq =~ s/[ 0-9]//g;
        $genome .= $dnaSeq;
    }
    close(IN);
    if (length($genome) == 0) { return("", "Genome sequence from $fileName is empty"); }
    return($genome, ""); 
}

###############################################################################
# Read the CDS sequence positions from a GenBank file, extract the sequences
# from the genome sequence, and return an array of CDS sequences.
# Arguments:
#   $fileName: name of GenBank file to read.
#   $genome: genome from the GenBank file.
#   $STARTsRef: (optional) reference to array of permissible start codons.
#       If unspecified, "AUG" is used.  May contain T or U.
# Returns: a list ($cdsRef, $errorRef):
#   $cdsRef: reference to an array of CDS sequences, empty array if unable to read.
#   $errorsRef: reference to an array of error messages, empty array if no errors.
###############################################################################
sub readGenBankFile_CDSseq {
    my $fileName = shift;
    my $genome = shift;
    my $STARTsRef = (scalar(@_) > 0) ? shift : ["AUG"];
    if (!open(IN, "<$fileName")) { return([], ["Error opening $fileName for reading"]); }
    my @cds;
    my @errors;
    my $line;
    my $lineNum = 0;
    while (my $line = <IN>) {
        $lineNum++;
        $line = eolChomp($line);
        # Ignore all but CDS lines.
        if ($line =~ m/^ +CDS +/) {
            # Parse CDS line to get start/end position and complement flag.
            #     CDS             5234..5530
            #     CDS             complement(5683..6459)
            #     CDS             complement(join(4225844..4226284,4227001..4227165))
            my ($compl, $start, $end) = ($line =~ m/^ +CDS +(complement)*\(*(\d+)\.\.(\d+)\)*$/);
            $compl = "" if (!defined($compl)); # Don't error if "complement" not found
            # Test for problems with the data.
            if ($line =~ m/^ +CDS +.*join\(/) {
                push(@errors, sprintf("Line %6d: joins are not supported", $lineNum));
            }
            elsif (!defined($start) or !defined($end) or $start < 0 or $end < $start) {
                push(@errors, sprintf("Line %6d: Bad start/end position:\n  %s", $lineNum, $line));
            }
            elsif ($end > length($genome)) {
                push(@errors, sprintf("Line %6d: End position beyond genome:\n  %s", $lineNum, $line));
            }
            else {
                # Retrieve the sequence and check for problems with it.
                my $seq = substr($genome, $start-1, $end-$start+1);
                if ($compl) { $seq = reverseComplement($seq); }
                my $seqError = checkGoodCodingSeq($seq, $STARTsRef);
                if ($seqError ne "") {
                    push(@errors, sprintf("Line %6d: Not a good coding sequence: %s", $lineNum, $seqError));
                }
                else { push(@cds, $seq); }
            }
        }
        last if ($line =~ m/^ORIGIN/); # Start of genome sequence.
    }
    close(IN);
    return(\@cds, \@errors);
}

###############################################################################
# Read CDS sequences from a GenBank file and count their codon distribution.
# Arguments:
#   $fileName: name of GenBank file to read.
#   $codonCountsRef: reference to hash of codon counters, key is codon and
#       value is count of number of times that codon was seen.
#   $firstCodonCountsRef: reference to another hash of codon counters, key
#       is FIRST codon in the sequence and value is count of times it was seen.
#   $STARTsRef: (optional) reference to array of permissible start codons.
#       If unspecified, "AUG" is used.  May contain T or U.
# Returns: a list ($genomeLen, $codonCount, $errorsRef):
#   $genomeLen: length of the genome that was read, in nucleotide bases.
#   $codonCount: Total number of codons found in the sequence.
#   $errorsRef: reference to an array of error messages, empty array if no errors.
# On return, hashes %$codonCountsRef and %$firstCodonCountsRef have updated
# counts, preserving previous counts and adding to them.
###############################################################################
sub readGenBankFile_countCodons {
    my $fileName = shift;
    my $codonCountsRef = shift;
    my $firstCodonCountsRef = shift;
    my $STARTsRef = (scalar(@_) > 0) ? shift : ["AUG"];

    # Read the genome from the GenBank file.
    my ($genome, $error) = readGenBankFile_genomeSeq($fileName);
    if ($error ne "") { return(0, 0, [$error]); }
    my $genomeLen = length($genome);
    
    # Read the CDS sequences from the GenBank file.
    my ($cdsRef, $errorsRef) = readGenBankFile_CDSseq($fileName, $genome, $STARTsRef);
    if (scalar(@$cdsRef) == 0) {
        push(@$errorsRef, "No CDS's found in $fileName");
        return($genomeLen, 0, $errorsRef);
    }
    
    # Accumulate codon stats for the CDS's.
    my $codonCount = countCodons_Array($cdsRef, $codonCountsRef, $firstCodonCountsRef);
    if ($codonCount == 0) { push(@$errorsRef, "No codons found in file $fileName"); }
    return($genomeLen, $codonCount, $errorsRef);
}

###############################################################################
# Read a sequence from a FASTA file.  This can also be used to read quality
# scores from a FASTA-like file.
# Arguments:
#	$file: handle of open file from which to read.
#   $nextLine: the next line from the file, or undef if none.
#   $qual: if present, indicates the file contains quality scores (optional arg).
# Returns: list ($seqName, $seqArgs, $seq, $nextLine):
#	$seqName: Name of the sequence (identifier just after ">") or undef if
#       end of file reached.
#   $seqArgs: Other arguments following sequence identifier, which is taken to end at a
#       space, tab, or | character.
#   $seq: REFERENCE to the sequence itself.  For a quality score file, this is a
#       reference to the concatenation of all quality scores.
#   $nextLine: The next line from the file, or undef if none.
#   $seqLen: length of $seq.
###############################################################################
sub readFastaSeq {
    my $file = shift;
    my $nextLine = shift;
    my $qual = "";
    if (scalar(@_) > 0) { $qual = shift; }
    if (!defined($nextLine)) { $nextLine = <$file>; }
    if (!defined($nextLine)) { return(undef, undef, undef, undef); }
    $nextLine = eolChomp($nextLine);
    if ($nextLine !~ m/^>/) { die("Expected > as next line in FASTA, got $nextLine"); }
    my ($seqName, $seqArgs) = ($nextLine =~ m/^>([^ \t|]*)(.*)$/);
    my $seq = "";
    while (defined($nextLine = <$file>)) {
        $nextLine = eolChomp($nextLine);
        last if ($nextLine =~ m/^>/);
        if ($qual) { $seq .= " "; }
        $seq .= $nextLine;
    }
	return($seqName, $seqArgs, \$seq, $nextLine, length($seq));
}

###############################################################################
# Required 1 at end of a module.
###############################################################################
1;
__END__
# Module documentation.

=head1 NAME

TedsLibrary - Ted Toal's general-purpose Perl functions.

=head1 SYNOPSIS

  use TedsLibrary ':all';

=head1 DESCRIPTION

Refer to the source code comments in TedsLibrary.pm for details.

=head2 EXPORT

None by default.



=head1 SEE ALSO

Source code in file TedsLibrary.pm.

=head1 AUTHOR

Ted Toal, E<lt>twtoal@ucdavis.edu<gt>

=head1 COPYRIGHT AND LICENSE

Copyright (C) 2011 by Ted Toal

This library is free software; you can redistribute it and/or modify
it under the same terms as Perl itself, either Perl version 5.10.0 or,
at your option, any later version of Perl 5 you may have available.


=cut
