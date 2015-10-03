#!/usr/bin/perl -w 
###########################################################################
# See usage below.
###########################################################################
use strict;
use warnings;

# Get arguments.
my $n = @ARGV;
if ($n < 3) {
	print("Read N lists of unique k-mers obtained from N genomes, intersect them,\n");
	print("intersection size, and write intersecting k-mers to an intersection file.\n");
	print("The k-mer lists are assumed to be in alphabetical order.  If Jellyfish is\n");
	print("used to obtain the k-mers, it seems to output k-mers in alphabetical order\n");
	print("if the number of them is not too big, but when it gets large, it stops doing\n");
	print("it alphabetically, so sort them afterwards to make sure they are alphabetical.\n");
	print("Usage: perl $0 <outFile> <inFile1> <inFile2> ... <inFileN>\n");
	print(" <inFile1> : Text file containing genome #1 unique k-mers, one per line.\n");
	print(" <inFile2> ... <inFileN> : Same, for genomes #2..N.\n");
	print(" <outFile> : Output file.\n");
	die;
}
my ($outFile, @inFiles) = @ARGV;
my $inFile;
my $N = @inFiles;
my $Nm1 = $N-1;

# Open the k-mer files.
my @INFILES;
my $INFILE;
foreach my $i (0..$Nm1) {
    $inFile = $inFiles[$i];
    open($INFILES[$i], "<" . $inFile) or die("Can't open '$inFile': $!");
}

# Create output file.
open(ISECTFILE, ">" . $outFile) or die("Can't create '$outFile': $!");

# Read through the k-mer files consecutively, maintaining a matched position
# between the N files by comparing the k-mers and reading from the file which is
# alphabetically lower.  Keep track of number of k-mers seen in each input file.
my @kmers; # Current k-mer in each of N input files.
my $kmer;
my @Nkmers; # Number of k-mers read from each of N input files.
my $Nkmer;
foreach my $i (0..$Nm1) { $Nkmers[$i] = 0; }
my $readAllKmers = 1; # This is 1 if a k-mer needs to be read from each of the N files.
# Now loop until end of one file is reached.  For each k-mer common to all N
# files, write the k-mer to the output file and count the number of such k-mers.
my $nIsect = 0;
my $endOfAfile = 0;
while (!$endOfAfile) {
    # If readAllKmers, read through all files and get the first k-mer from each
    # one.  If the end of any file is reached, exit loop.
    if ($readAllKmers) {
        foreach my $i (0..$Nm1) {
            $INFILE = $INFILES[$i];
            $kmer = <$INFILE>;
            if (!defined($kmer)) {
                $endOfAfile = 1;
            } else {
                $kmer = eolChomp($kmer);
                $kmers[$i] = $kmer;
                $Nkmers[$i]++;
            }
        }
        last if ($endOfAfile);
        $readAllKmers = 0;
    }

    # Find largest of the k-mers, alphabetically.
    my $maxKmer = $kmers[0];
    foreach my $i (1..$Nm1) { if ($maxKmer lt $kmers[$i]) { $maxKmer = $kmers[$i]; } }

    # Start with allEqual set to 1, and then for each k-mer that is less than
    # the largest k-mer, set allEqual to 0 and read k-mers from file until one
    # is equal to or greater than the largest k-mer or end of file is reached.
    # If end of file, that is the end of the outer loop, set endOfAfile to 1.
    my $allEqual = 1; # Assume they are all equal to maxKmer.
    foreach my $i (0..$Nm1) {
        if ($kmers[$i] ne $maxKmer) {
            $allEqual = 0;
            $INFILE = $INFILES[$i];
            while ($kmers[$i] lt $maxKmer) {
                if (!defined($kmer = <$INFILE>)) {
                    $endOfAfile = 1;
                    last;
                } else {
                    $kmer = eolChomp($kmer);
                    $kmers[$i] = $kmer;
                    $Nkmers[$i]++;
                }
            }
        }
    }

    # If all k-mers were equal (allEqual is 1), that k-mer is in the intersection
    # set so write it out.  Then, set readAllKmers to 1 to re-read k-mers from
    # every file.
    if ($allEqual) {
        $nIsect++;
        print(ISECTFILE "$maxKmer\n");
        $readAllKmers = 1;
    }   
}

# Now loop until the end of ALL FILES is reached, simply reading and counting k-mers.
foreach my $i (0..$Nm1) {
    $INFILE = $INFILES[$i];
    while (defined($kmer = <$INFILE>)) {
        $Nkmers[$i]++;
    }
}

# Close the input files.
foreach my $i (0..$Nm1) { close($INFILES[$i]); }

# Close the output file.
close(ISECTFILE);

# Report stats.
foreach my $i (1..$N) { print(" Genome $i k-mers: $Nkmers[$i-1]\n"); }
print(" Intersection k-mers: $nIsect\n");
exit;
