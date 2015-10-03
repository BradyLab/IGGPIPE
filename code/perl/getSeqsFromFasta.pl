#!/usr/bin/perl 
###########################################################################
# Read a FASTA file and extract specified sequences from it, writing them
# to separate files.  This also supports FASTA quality files, which have
# quality numbers instead of bases in them.
###########################################################################
use strict;
use warnings;
use File::Basename;
use lib dirname($0); # Look for TedsLibrary.pm in same directory as this .pl file.
use TedsLibrary ':all';
use Getopt::Long;

# Process options.
Getopt::Long::Configure('bundling');
Getopt::Long::Configure('bundling_override');
my %options;
GetOptions(\%options,
    'h', 'help', 'c|contains', 'a|anywhere', 'e|end', 'i|idfile=s', 'o|output=s',
    'p|prefix=s', 's|suffix=s', 'r|revComp', 'R=s', 'O|omit', 'n|name=s',
    't=s', 'U=s', 'u=s', 'q|qual', 'd', 'delete=s', 'l|lineSize=i', 'v|verbose');

# Compact help text.
my $S0 = basename($0);
my @compactHelpText = (
	"Extract sequences from FASTA or quality file into FASTA output files.",
	"Each extracted sequence is written to an individual output file unless -o",
	"Usage: perl $S0 [options] <inputFastaFile> [<seq>...]",
	" <inputFastaFile>        : FASTA file from which to read sequences.",
	" <seq>                   : Sequences to extract, in format:",
	"    [!]<seqId>[:<seqStart>..<seqEnd>]   (Optional ! means to reverse-complement.)",
	" [options]:",
	"  --help                 : Show more verbose usage info.",
	"  -c or --contains       : Match sequence IDs CONTAINING <seqId>.",
	"                           (Default is to match <seqId> exactly).",
	"  -a or --anywhere       : <seqId>'s may match ANY PART of sequence ID line, ignore case.",
	"  -e or --end            : Error if sequence shorter than requested length.",
	"  -i=FIL or --idfile=FIL : Text file FIL is a list of <seq>'s to extract.",
	"  -o=FIL or --output=FIL : Output all gene sequences to file FIL.",
	"                           (Sequences in output file are NOT in same order as input!)",
	"  -p=PFX or --prefix=PFX : The output FASTA file(s) are named PFX+seqStr+SFX",
	"  -s=SFX or --suffix=SFX : The output FASTA file(s) are named PFX+seqStr+SFX",
	"  -r or --revcomp        : Reverse complement ALL output sequences even if no !",
	"  -R=PTN                 : Reverse complement if regular expression matches ID line.",
	"  -O or --omit           : Do not write the sequence ID lines to the output file(s).",
	"  -n=NAME or --name=NAME : Use NAME for %ISID% (see -U), for making an output ID line.",
	"  -t=PTN                 : Regular expression to convert the input sequence ID into",
	"                           %ISID%, for making an output ID line, w/parenthesized subexpr.",
    "  -U=template            : Template string for generating the output sequence names.",
    "                           Replace %ISID%, %ISARGS%, %SEQ%, %REVCOMP%, %SUBRGN% (see -hh).",
	"  -u=PTN                 : Search and replace regular expression to convert the template",
	"                           string specified by -U into an output ID line, with internal /",
	"  -q or --qual           : Sequence is quality values rather than base values.",
	"  -d                     : Delete all '*' characters from the sequence.",
	"  --delete=ch            : Delete all ch characters from the sequence.",
	"  -l=# or --lineSize=#   : Number of sequence characters per line in output, 0=all 1 line.",
	"  -v or -verbose         : Verbose debug output.",
	"Example: perl $S0 ../ITAG2.3_cdna.fa Solyc04g064820:2500..3000 -o=G64820.fa",
	""
    );

# Long help text.
my @longHelpText = (
	"Extract sequences from FASTA sequence or quality file into one or more FASTA output files.",
	" Each sequence with specified <seqId> and optional subsequence position is",
	" extracted from <inputFastFile> and written to an individual output file unless",
	" -o is used to output to one file.",
	"Usage: perl $S0 [options] <inputFastaFile> [[!]<seq>...]",
	" <inputFastaFile>        : FASTA file from which to read sequences.",
	" <seq>                   : zero or more sequence IDs to extract from <inputFastaFile>.",
	"    Each <seq> is of the form <seqId>[:<seqStart>..<seqEnd>], i.e. a sequence",
	"    name (<seqId>) optionally followed by a subsequence (:<seqStart>..<seqEnd>).",
	"    If none specified, ALL sequences are extracted to output file(s) unless -i.",
	"    If '!' precedes a <seq>, the sequence is reverse-complemented (surround !<seq>",
	"    with ' marks so bash shell doesn't get upset by the !).",
	"    If a sub-sequence is specified, this means to extract a subsequence from within",
	"    the sequence.  E.g. 'Chr7:356..9876' means to extract bases 356 through 9876 of",
	"    the sequence whose name is Chr7 (first base is numbered 1).",
	" [options]:",
	"  -h                     : Show less verbose usage info.",
	"  --help                 : Show this usage info.",
	"  -c or --contains       : Specified <seqId> matches sequence IDs CONTAINING <seqId>.",
	"    (<seqId>s minus a leading ! are actually treated as Perl regular expressions.)",
    "    (Searches use the i modifier for case-insensitive searches.)",
	"    (If unspecified, each sequence ID matching <seqId> EXACTLY is extracted.)",
	"  -a or --anywhere       : Specified <seqId>'s may match any part of the sequence ID",
	"    or the text following that ID on the FASTA ID line, and case is ignored.",
	"  -e or --end            : If a sequence is shorter than the requested length, halt",
	"    with an error.  Default is to extract to the end of the sequence with no error.",
	"  -i=FIL or --idfile=FIL : Text file FIL is a list of <seq>'s to extract,",
	"    with optional ! in front of a name to indicate reverse complement wanted.",
	"    Each <seq> is of the form '<seqId>[:<seqStart>..<seqEnd>]' ([] = optional).",
	"  -o=FIL or --output=FIL : Output all gene sequences to file FIL.  If this",
	"    option is not specified, each gene sequence is output to its own file named",
	"    PFX+seqStr+SFX where seqStr is <seqId>, or if <seqId> is followed by",
	"    a subsequence specified, seqStr is <seqId>_<start>_<end> where <start> and",
	"    <end> are the subsequence start/end numbers (see -p and -s options below).",
	"    Note that the sequences in the output file are generally NOT in the same order",
	"    as the <seq>'s in the request list!!!",
	"  -p=PFX or --prefix=PFX : The output FASTA file(s) are named PFX+seqStr+SFX",
	"    (If unspecified, PFX is empty.  You might choose to make it a dir, e.g. 'out/'.)",
	"  -s=SFX or --suffix=SFX : The output FASTA file(s) are named PFX+seqStr+SFX",
	"    (If unspecified, SFX is '.fasta'.)",
	"  -r or --revcomp        : Reverse complement all output sequences.",
	"  -R=PTN                 : PTN is a Perl regular expression pattern that is matched",
	"    against each FASTA sequence ID line, and if it matches, that sequence will be",
	"    reverse complemented before being written to the output file.  The pattern must",
	"    not include outer // delimiters.",
	"  -O or --omit           : Do not write the sequence ID lines to the output file(s)",
	"    Note that if the -o option below is used, the sequences in the output file are",
	"    generally NOT in the same order as the <seq>'s in the request list!!!",
	"  -n=NAME or --name=NAME : Use NAME for %ISID% (see -U), for making an output ID line.",
	"  -t=PTN                 : Perform pattern search/replace against the input sequence",
    "    ID.  PTN is a Perl regular expression pattern that is matched against the ID.  It",
	"    must include a parenthesized subexpression whose value will be used for %ISID%",
	"    for replacement within the -T template for generating the output sequence ID lines.",
	"    The pattern must not include the outer // delimiters.",
    "  -U=template            : Template string for generating the output sequence names.",
    "    Within this template string make these replacements:",
    "      %ISID% with the input FASTA sequence ID or -n argument.",
    "      %ISARGS% with input FASTA sequence arguments after ID.",
    "      %SEQ% with <seq> argument used to select sequence.",
    "      %REVCOMP% with ' revcomp:yesORno' (rev. compl. state)",
    "      %SUBRGN% with ' sub_region:start#-end#' if <seq> has :<seqStart>..<seqEnd> in it,",
    "        else this is an empty string.",
    "    After replacements, apply the -u search/replace pattern.",
    "    Default is '%ISID%%REVCOMP%%SUBRGN%%ISARGS%'",
	"  -u=PTN                 : Search and replace regular expression to convert the template",
	"    string specified by -U into an output ID line.  This is a Perl substitution operator",
	"    regular expression that searches and replaces characters in the template string.  The",
	"    pattern must not include the outer // delimiters but it must include an inner /",
	"    delimiter, e.g. 'revcomp:/COMP:'   Default is an empty string (use template as-is).",
	"  -q or --qual           : The sequence is quality values rather than base values, so",
	"    the -r and -R options only reverse the number order and do not complement, and if",
	"    a subsequence was specified, the indexes for the subsequence are indexes into the",
	"    array of quality scores.",
	"  -d                     : Delete all '*' characters from the sequence",
	"  -sdelete=ch            : Delete all ch characters from the sequence",
	"  -l=# or --lineSize=#   : Number of sequence characters per sequence line in output,",
	"    or, for -q quality files, number of quality values per line in output.  Use 0 to",
	"    put all sequence characters or values on one line.",
	"  -v or -verbose         : Verbose debug output.",
	"Example 1: perl $S0 ../ITAG2.3_cdna.fa '!Solyc04g064820' -s=.fa",
	"Example 2: perl $S0 ../ITAG2.3_cdna.fa Solyc04g064820:2500..3000 -s=.fa",
	"Example 3: perl $S0 -nM82 Pseudo_M82.fasta -r -s=.rev.fa",
	"Example 4: perl $S0 -t 'Solyc..g(.{6}...[a-z]+)' -l80 primers.fasta -r -o=myData.fa", # Shows user two apostrophes, which he must type.
    "Example 5: perl $S0 ../ITAG2.3_cdna.fa Solyc04g064820 -s=.fa -U '%ISID%' -u 'Soly(.*)...\$/\$1'",
	""
    );

# Get script arguments, display help.
my $n = @ARGV;
if (exists($options{help})) { foreach my $s (@longHelpText) { print("$s\n"); } }
elsif (exists($options{h}) || $n < 1) { foreach my $s (@compactHelpText) { print("$s\n"); } }
if ($n < 1) { die("Insufficient arguments"); }
my ($inputFastaFile, @seqs) = @ARGV;
 
# Number of sequence characters on a sequence line in FASTA file.
my $numSeqValsPerLine = (exists($options{l}) ? $options{l} : 60);

# Verbose flag.
my $verbose = exists($options{v});

# Is this a sequence quality scores file?
my $qual = exists($options{q});

# Template string.
my $template = (exists($options{U}) ? $options{U} : "%ISID%%REVCOMP%%SUBRGN%%ISARGS%");

# Template string substitution strings.
my $templateSearch;
my $templateReplace;
if (exists($options{u})) {
    $templateSearch = $options{u};
    $templateReplace = $options{u};
    $templateSearch =~ s|/.*$||;
    $templateReplace =~ s|^.*/||;
    $templateReplace = '"' . $templateReplace . '"';
}

# Read -i file if any specified, and append sequence Ids to seqs list.
if (exists($options{i}) && $options{i}) {
    my $seqFileName = $options{i};
    if ($verbose) { print("Reading input file $seqFileName\n"); }
    open(my $seqFile, "<", $seqFileName) or die("Can't open -i file '$seqFileName': $!");
    my $line;
    while (defined($line = <$seqFile>)) {
        $line = eolChomp($line);
        push(@seqs, $line);
    }
    close($seqFile);
}

# Parse sequence Ids to get sequence name, start, end, and reverse complement flags.
# Use the reqSeqIds as keys to hash %seqs, whose values are themselves hashes with keys
# "seq", "start", "end", and "revcomp".  The values of those four keys are arrays for
# that reqSeqId, containing:
#   seq: the full <seq> string specified by the user.
#   start: the start position that was specified, 0 if none
#   end: the end position that was specified, 0 if none
#   revcomp: the reverse complement flag
# Note that the same seqid might be specified more than once!  We use start=end=0
# to mean the entire sequence is to be output.
my %seqs;
if ($verbose) { print("Parsing sequence IDs\n"); }
foreach my $seq (@seqs) {
    my ($rev, $reqSeqId, $subSeq) = ($seq =~ m/^(!?)([^:]+)(:[0-9]+..[0-9]+)?$/);
    #if (!defined($reqSeqId)) { print("seq=$seq rev=$rev reqSeqId=$reqSeqId subSeq=$subSeq\n"); }
    #print("rev=$rev reqSeqId=$reqSeqId subSeq=$subSeq\n");
    my ($start, $mid, $end) = (0, 0, 0); # Don't care about $mid but it is captured.
    if (defined($subSeq)) {
        ($start, $mid, $end) = ($subSeq =~ m/^:([0-9]+)(\.\.|\-|:)([0-9]+)$/); # Can we allow #-# or #:# also?
        if ($start < 1) { die("Start of subsequence cannot be < 1"); }
        if ($end < $start) { die("End of subsequence cannot be < start"); }
    }
    push(@{$seqs{$reqSeqId}->{seq}}, $seq);
    push(@{$seqs{$reqSeqId}->{revcomp}}, ($rev eq "!"));
    push(@{$seqs{$reqSeqId}->{start}}, $start);
    push(@{$seqs{$reqSeqId}->{end}}, $end);
}
my @reqSeqIds = keys(%seqs);
my %hSeqFound = map { $_ => 0 } @reqSeqIds;
my $numSeqs = @reqSeqIds;
my $extractAllSeqs = ($numSeqs == 0);

# Open FASTA file.
open(my $inputFasta, "<", $inputFastaFile) or die("Can't open '$inputFastaFile': $!");

# If all sequences go to a single output file, create that file.
my $outFilename = "";
if (exists($options{o})) { 
    $outFilename = $options{o};
    open(OUTFILE, ">" . $outFilename) or die("Can't create '$outFilename': $!");
}

# Loop reading sequences and outputting any in the requested sequence name list to an
# output file.
my ($seqId, $seqArgs, $seqRef, $seqNext, $seqLen);
$seqNext = <$inputFasta>;
my $numOutSeqs = 0;
if ($verbose) { print("Starting read of FASTA file $inputFastaFile\n"); }
while (1) {
    ($seqId, $seqArgs, $seqRef, $seqNext, $seqLen) = readFastaSeq($inputFasta, $seqNext, $qual);
    if (!defined($seqId)) { last; }
    if ($verbose) { print("Processing sequence ID $seqId\n"); }

    # Set $outputSeq true if we decide the sequence is to be output.
    my $outputSeq = 0;

    # Recall from above that the %seqs hash contains the sequences requested by the user.
    # See comments above for description of the sub-keys of %seqs.  Below we will assign
    # to %subseqs the %seqs member corresponding to $seqId when $outputSeq is true.  (It
    # is plural %subseqs instead of %subseq only because the "seq", "start", "end", and
    # "revcomp" key values are arrays).
    my %subseqs;

    # Output the sequence (and/or sub-sequences of it) to a file if all sequences are to
    # be extracted, or if an exact match was requested and the sequence name is in the
    # hash of sequence names, or if no exact match and the sequence name matches any of
    # the specified names.
    if ($verbose) { print("Checking to see which requested sequences are in this sequence ID\n"); }
    if ($extractAllSeqs) {
        $outputSeq = 1;
        $subseqs{revcomp}->[0] = ""; # Don't reverse complement it.
        $subseqs{start}->[0] = 0; # Output entire sequence.
        $subseqs{end}->[0] = 0; # Output entire sequence.
    } else {
        if (!exists($options{c})) {
            $outputSeq = exists($hSeqFound{$seqId});
            if ($outputSeq) { %subseqs = %{$seqs{$seqId}}; }
            $hSeqFound{$seqId} = 1;
        } elsif (exists($options{a})) {
            my $fullIDline = $seqId . $seqArgs;
            foreach my $reqSeqId (@reqSeqIds) {
                if ($fullIDline =~ m/$reqSeqId/i) {
                    $outputSeq = 1;
                    %subseqs = %{$seqs{$reqSeqId}};
                    $hSeqFound{$reqSeqId} = 1;
                    last;
                }
            }
        } else {
            foreach my $reqSeqId (@reqSeqIds) {
                if ($seqId =~ m/$reqSeqId/i) {
                    $outputSeq = 1;
                    %subseqs = %{$seqs{$reqSeqId}};
                    $hSeqFound{$reqSeqId} = 1;
                    last;
                }
            }
        }
    }

    # If the sequence (and/or sub-sequences of it) is to be output, do it.
    if (!$outputSeq) {
        if ($verbose) { print("No sequences to output\n"); }
    } else {
        if ($verbose) { print("Outputting requested sequence IDs\n"); }
        my $originalSeqLine = $seqId . $seqArgs;

        # Get global reverse complement flag, for all subsequences of this sequence.
        my $revCompAll = (exists($options{r})) || (exists($options{R}) && $originalSeqLine =~ m/$options{R}/);

        # Get the default value for %ISID% replacement in the template.
        my $ISID = $seqId;

        # Change it according to -t and -n options.
        if (exists($options{t})) { $ISID =~ s/$options{t}/$1/; }
        if (exists($options{n})) { $ISID = $options{n}; }

        # Replace %ISID% and %ISARGS% in template.
        my $templatePartiallySubstituted = $template;
        $templatePartiallySubstituted =~ s/%ISID%/$ISID/;
        $templatePartiallySubstituted =~ s/%ISARGS%/$seqArgs/;

        # Determine whether or not %SEQ5, %REVCOMP%, and %SUBRGN% occur in the
        # template, so we can skip the search/replace if it isn't necessary.  We
        # may be doing this on many thousands of strings, it can be slow.
        my $Replace_SEQ = ($templatePartiallySubstituted =~ m/%SEQ%/);
        my $Replace_REVCOMP = ($templatePartiallySubstituted =~ m/%REVCOMP%/);
        my $Replace_SUBRGN = ($templatePartiallySubstituted =~ m/%SUBRGN%/);

        # Delete * (e.g. STOP at end), or whatever user wants deleted.
        if (exists($options{d})) { $$seqRef =~ s/\*//g; }
        if (exists($options{delete})) { $$seqRef =~ s/$options{delete}//g; }

        # Loop for each subsequence to be output.
        my $numSubseqs = @{$subseqs{revcomp}};
        if ($verbose) { print("Outputting $numSubseqs subsequences\n"); }
        for (my $i = 0; $i < $numSubseqs; $i++) {
            if ($verbose) { print("Outputting subsequence $i\n"); }
            # Get the default value for %SUBRGN% replacement in the template.
            my $SUBRGN = "";
            # Get the individual output filename.
            my $individualOutFilename = $seqId;
            # Get start and end positions and length.
            my ($start, $end) = ($subseqs{start}->[$i], $subseqs{end}->[$i]);
            my $sublen = $end - $start + 1;
            # Get reverse complement flag for this subsequence.
            my $doRev = ($subseqs{revcomp}->[$i]) || $revCompAll;
            # Get the default value for %REVCOMP% replacement in the template.
            my $REVCOMP = " revcomp:no";
            if ($doRev) { $REVCOMP = " revcomp:yes"; }
            # Sequence to be output.  Don't initialize to $$seqRef, which might be very long.
            my $thisSequence;

            # If not a quality score file, do basic processing.
            if (!$qual) {
                # Extract sub-sequence if requested.
                if ($start > 0) {
                    my $len = $seqLen;
                    if ($len < $end) {
                        if (exists($options{e})) { die("Sequence $seqId length $len is shorter than requested end position $end"); }
                        $end = $len;
                        my $sublen = $end - $start + 1;
                    }
                    $thisSequence = substr($$seqRef, $start-1, $sublen);
                    # Include sub-range in sequence args.
                    $SUBRGN = " sub_region:$start-$end";
                    $individualOutFilename = $individualOutFilename . "_" . $start . "_" . $end;
                } else { $thisSequence = $$seqRef; }

                # Reverse complement the sequence if unconditional -r or if the specified -r
                # pattern matches or if the reversal [^] flag was specified for the sequence.
                if ($doRev) { $thisSequence = reverseComplement($thisSequence); }
                # Perform replacements in the template string.
                my $outseqIDline;
                if (!exists($options{O})) {
                    $outseqIDline = $templatePartiallySubstituted;
                    if ($Replace_SEQ) { $outseqIDline =~ s/%SEQ%/$subseqs{seq}->[$i]/; }
                    if ($Replace_REVCOMP) { $outseqIDline =~ s/%REVCOMP%/$REVCOMP/; }
                    if ($Replace_SUBRGN) { $outseqIDline =~ s/%SUBRGN%/$SUBRGN/; }
                    # See http://stackoverflow.com/questions/392643/how-to-use-a-variable-in-the-replacement-side-of-the-perl-substitution-operator
                    if (exists($options{u})) {
                        $outseqIDline =~ s/$templateSearch/$templateReplace/ee;
                    }
                }
                # If every sequence is written to its own file, create that file.
                my $outFilename1 = "";
                if ($outFilename eq "") {
                    $outFilename1 = $individualOutFilename;
                    if (exists($options{p})) { $outFilename1 = $options{p} . $outFilename1; }
                    if (exists($options{s})) { $outFilename1 .= $options{s};
                    } else { $outFilename1 .= ".fasta"; }
                    open(OUTFILE, ">" . $outFilename1) or die("Can't create '$outFilename1': $!");
                }
                # Output the sequence.
                if (!exists($options{O})) { print(OUTFILE ">" . $outseqIDline . "\n"); }
                if ($numSeqValsPerLine < 1) {
                    print(OUTFILE $thisSequence, "\n");
                } else {
                    for (my $i = 0; $i < length($thisSequence); $i += $numSeqValsPerLine) {
                        print(OUTFILE substr($thisSequence, $i, $numSeqValsPerLine), "\n");
                    }
                }
                $numOutSeqs++;
                # Close the file if individual files per sequence.
                if ($outFilename eq "") {
                    close(OUTFILE);
                    print("Created file $outFilename1 containing sequence for $seqId\n");
                }

            # If this is a quality score file, subsequent processing is different.
            } else {
                # Split the sequence string into an array of quality score strings.
                my @arrSeq = split(" ", $$seqRef);
                # Extract sub-sequence if requested.
                if ($start > 0) {
                    my $len = scalar(@arrSeq);
                    if ($len < $end) {
                        if (exists($options{e})) { die("Sequence $seqId length $len is shorter than requested end position $end"); }
                        $end = $len;
                    }
                    @arrSeq = @arrSeq[$start-1 .. $end-1];
                    # Include sub-range in sequence args.
                    $SUBRGN = " sub_region:$start-$end";
                    $individualOutFilename = $individualOutFilename . "_" . $start . "_" . $end;
                }
                # Reverse the scores if unconditional -r or if the specified -r
                # pattern matches or if the reversal [^] flag was specified for the sequence.
                if ($doRev) { @arrSeq = reverse(@arrSeq); }
                # Perform replacements in the template string.
                my $outseqIDline;
                if (!exists($options{O})) {
                    $outseqIDline = $templatePartiallySubstituted;
                    if ($Replace_SEQ) { $outseqIDline =~ s/%SEQ%/$subseqs{seq}->[$i]/; }
                    if ($Replace_REVCOMP) { $outseqIDline =~ s/%REVCOMP%/$REVCOMP/; }
                    if ($Replace_SUBRGN) { $outseqIDline =~ s/%SUBRGN%/$SUBRGN/; }
                    if (!exists($options{u})) {
                        $outseqIDline =~ s/$templateSearch/$templateReplace/;
                    }
                }
                # If every sequence is written to its own file, create that file.
                my $outFilename1 = "";
                if ($outFilename eq "") {
                    $outFilename1 = $individualOutFilename;
                    if (exists($options{p})) { $outFilename1 = $options{p} . $outFilename1; }
                    if (exists($options{s})) { $outFilename1 .= $options{s};
                    } else { $outFilename1 .= ".fasta"; }
                    open(OUTFILE, ">" . $outFilename1) or die("Can't create '$outFilename1': $!");
                }
                # Output the sequence.
                if (!exists($options{O})) { print(OUTFILE ">" . $outseqIDline . "\n"); }
                if ($numSeqValsPerLine < 1) {
                    print(OUTFILE join(" ", @arrSeq), "\n");
                } else {
                    for (my $i = 0; $i < scalar(@arrSeq); $i += $numSeqValsPerLine) {
                        my $j = $i+$numSeqValsPerLine-1;
                        if ($j >= scalar(@arrSeq)) { $j = scalar(@arrSeq)-1; }
                        $thisSequence = join(" ", @arrSeq[$i..$j]);
                        print(OUTFILE $thisSequence, "\n");
                    }
                }
                $numOutSeqs++;
                # Close the file if individual files per sequence.
                if ($outFilename eq "") {
                    close(OUTFILE);
                    print("Created file $outFilename1 containing sequence for $seqId\n");
                }
            }
        }
    }
}
close($inputFasta);
if ($outFilename ne "") {
    close(OUTFILE);
    print("Created file $outFilename containing $numOutSeqs sequences\n");
}

if ($verbose) { print("Checking for sequence IDs not found.\n"); }
foreach my $seqId (keys(%hSeqFound)) {
    if (!$hSeqFound{$seqId}) {
        print("***** NOT FOUND: $seqId\n");
    }
}
