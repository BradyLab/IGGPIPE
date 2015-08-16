#!/usr/bin/perl 
###########################################################################
# Read a FASTA file, extract its sequence IDs and their lengths, and write
# the data to a file.  Output format is controlled by arguments, with the
# default format:
# id	len
# SL2.40ch00	16439357
# SL2.40ch01	68062687
###########################################################################
use strict;
use warnings;
use File::Basename;
use lib dirname($0); # Look for TedsLibrary.pm in same directory as this .pl file.
use TedsLibrary ':all';
use Getopt::Long;

# Hash with keys whose names are #word# words specified for the header argument to be
# replaced by the specified hash values.
my %headerReplaceHash = ("t" => "\t", "#" => "#");

# Hash with keys whose names are #word# words specified for the data argument to be
# replaced by actual sequence values.
my %dataReplaceHash = (id => "", len => "", end => "", args => "", seq => "", "/" => "/",
    "t" => "\t", "#" => "#");

# Default options.
my $defaultHeader = "id#t#len";
my $defaultData = "#id##t##len#";

# Process options.
Getopt::Long::Configure('bundling');
Getopt::Long::Configure('bundling_override');
my %options;
GetOptions(\%options, 'h|help', 'header=s', 'data=s');

# Get script arguments.
my $n = @ARGV;
if (exists($options{h}) || $n < 2) {
	print("Extract sequence IDs from FASTA file into text file.\n");
	print("Usage: perl $0 [options] <inputFastaFile> <outputTextFile>\n");
	print(" <inputFastaFile> : FASTA file from which to read sequences.\n");
	print(" <outputTextFile> : Text file to which to write sequence IDs/lengths.\n");
	print(" [options]:\n");
	print("  -h or --help: Show this usage info.\n");
	print("  --header=HDR: Header text to write to output file.  If empty string, no\n");
	print("    header is written.  If unspecified, HDR='$defaultHeader' is used.  #t# is\n");
	print("    replaced by tab and ### by one #.\n");
	print("  --data=DAT  : Data format of each text line.  If unspecified, the default\n");
	print("    is DAT='$defaultData'.  #id# is replaced by the sequence id, #len# by its\n");
	print("    length, #end# by its length minus 1, #args# by any arguments following\n");
	print("    sequence id, #seq# by the sequence itself, #/# by /, #t# by tab, and ### by\n");
	print("    one #.  Also, /regexp/ is replaced by the first matching substring when\n");
	print("    regexp is matched against the sequence ID.\n");
	print("Example 1: perl $0 ITAG2.3_genomic.fa ITAG2.3_genomic.idlens\n");
	print("Example 2: perl $0 --header \"id\tlen\tchr\" --data \"#id#\t#len#\t/SL2.40ch([0-9]{2})/\" ITAG2.3_genomic.fa ITAG2.3_genomic.idlens\n");
	die;
}
my ($inputFastaFile, $outputTextFile) = @ARGV;
my $header = $defaultHeader;
if (exists($options{header})) { $header = $options{header}; }
my $data = $defaultData;
if (exists($options{data})) { $data = $options{data}; }

###########################################################################
# Replace /regexp/ sequences within a string by applying the regexp to a
# specified string and using the first matched subexpression as the
# replacement value.
# Arguments:
#	$s: string in which replacements are to be made.
#   $matchStr: the string to apply each regexp to.
# Returns: substituted string.
###########################################################################
sub replaceRegexps {
    my $s = shift;
    my $matchStr = shift;
    while (1) {
        my ($regexp) = ($s =~ m@/([^/]*)/@);
        last if (!defined($regexp));
        my ($replacement) = ($matchStr =~ m/$regexp/);
        if (!defined($replacement)) { $replacement = ""; }
        $s =~ s@/[^/]*/@$replacement@g;
    }
    return($s);
}

###########################################################################
# Replace #word# sequences within a string by values specified by a hash,
# where the hash key is "word" and the hash value is the replacement value.
# Arguments:
#	$s: string in which replacements are to be made.
#   $href: reference to hash containing replacements.
# Returns: substituted string.
###########################################################################
sub replaceWords {
    my $s = shift;
    my $href = shift;
    my %h = %$href;
    my @words = keys(%h);
    foreach my $word (@words) {
        my $replacement = $h{$word};
        $s =~ s/#$word#/$replacement/g;
    }
    return($s);
}

# Open FASTA file.
open(my $inputFasta, "<", $inputFastaFile) or die("Can't open '$inputFastaFile': $!");

# Create output file.
open(OUTFILE, ">" . $outputTextFile) or die("Can't create '$outputTextFile': $!");
$header = replaceWords($header, \%headerReplaceHash);
print(OUTFILE "$header\n");

# Loop reading sequences and outputting one text line for each one.
my ($seqName, $seqArgs, $seq, $seqNext);
$seqNext = <$inputFasta>;
while (1) {
    ($seqName, $seqArgs, $seq, $seqNext) = readFastaSeq($inputFasta, $seqNext);
    if (!defined($seqName)) { last; }
    my $s = replaceRegexps($data, $seqName);
    $dataReplaceHash{id} = $seqName;
    $dataReplaceHash{len} = length($seq);
    $dataReplaceHash{end} = length($seq)-1;
    $dataReplaceHash{args} = $seqArgs;
    $dataReplaceHash{seq} = $seq;
    $s = replaceWords($s, \%dataReplaceHash);
    print(OUTFILE "$s\n");
}
close($inputFasta);
close(OUTFILE);
