/*
    Find k-mer positions in FASTA sequences.  See usage() below.

    This program uses a simple scheme for storing k-mers, requiring 2 bits of
    memory per POSSIBLE k-mer, for a total of 2^(2k+1) bits or 2^(2k-2) bytes
    of memory required.  For k=14, 64 MB are required.
    
    K-mers are represented (temporarily, not for storing in 2-bit array) in
    binary using 2 bits per base, with the following encoding so that
    ones-complementing a base gives its complementary base:
        A - 00
        T - 11
        C - 01
        G - 10

    An option is available to intersect two k-mer lists and find positions of the
    resulting intersection k-mers.

    Algorithm:
        1. Read the k-mer file and initialize a table.  The table has 2^2k
            entries, one for each POSSIBLE k-mer, and each entry is 2 bits.
            The values are as follows:
                0 - this is not a k-mer from the input k-mer text file (else it
                    is one).
                1 - k-mer not yet seen in input FASTA sequence
                2 - k-mer seen once in input FASTA sequence
                3 - k-mer seen more than once in input FASTA sequence
        2. If k-mer intersection option is specified, read the second k-mer file
            and intersect it with the k-mer table initialized above.  Change any
            k-mer whose value is 1 to a value of 2 to indicate that the k-mer
            appears in both files.  When finished, reset all k-mers whose value
            is 1 to a value of 0, and all k-mers whose value is 2 to a value of 1.
        3. Read input FASTA sequence file, one sequence at a time, run through
            the sequence and generate each k-mer and reverse-complement k-mer,
            look them up in the k-mer table, and if found adjust the table entry
            and write the k-mer and its position to the k-mer position output file.
        4. When finished, list any k-mers from the text file that were not seen
            in the input sequences or were seen more than once.

        A k-mer that is its own reverse complement is not counted twice when seen.
*/
#include <iostream>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <math.h>
#include <fstream>
#include <vector>
#include "stdhdr.h"
#include "Strings.h"
#include "InputOutput.h"

#define NUM_BASES 4
#define BITS_PER_BASE 2
#define BASE_BIT_MASK 3
#define VAL_A   0
#define VAL_T   3
#define VAL_C   1
#define VAL_G   2

#define TBL_BITS_PER_KMER   2   // In k-mer table, # bits for each k-mer
#define TBL_KMER_BIT_MASK   3   // Mask for those bits = 2^TBL_BITS_PER_KMER-1
#define TBL_KMERS_PER_BYTE  4   // Number of k-mers in each table byte = 8/TBL_BITS_PER_KMER
#define TBL_KMER_IDX_SHIFT  2   // Bit shift count for a kmer to get k-mer table index = log2 TBL_KMERS_PER_BYTE
#define TBL_KMER_IDX_MASK   3   // Mask to get index into a table byte of the k-mer = 2^TBL_KMER_IDX_SHIFT-1

#define TBL_UNKNOWN_KMER    0
#define TBL_KMER_UNSEEN     1
#define TBL_KMER_SEEN_ONCE  2
#define TBL_KMER_SEEN_MANY  3

#define BITS_PER_BYTE 8

#define MIN_K   1
#define MAX_K   16
typedef ubyte4 kmer; // 32 bits holds 16-mer

// Define types to hold a set of sequence contigs, each consisting of a string of ATCG's
// encoded using 2 bits per base, with any non-ATCG character splitting the sequence into
// multiple contigs.  (I.e. N's break up the FASTA sequence into contigs).
typedef std::vector<ubyte> seqcontigbytes;
struct seqcontig
    {
    unsigned startPos;
    unsigned contigLen;
    seqcontigbytes* contigBytes;
    };
typedef std::vector<struct seqcontig> seqcontigs;
// Function to discard a dynamically allocated seqcontigs struct.
static void discardSeqContigs(seqcontigs* seqContigs)
    {
    for (unsigned i = 0; i < seqContigs->size(); i++)
        {
        delete (*seqContigs)[i].contigBytes;
        (*seqContigs)[i].contigBytes =  NULL;
        }
    delete seqContigs;
    }

// Show program usage information and exit program.
static void usage(void)
    {
    static const char* const usage[] =
        {
        "NAME",
        "       findMers - locate k-mers in FASTA sequences",
        "",
        "SYNOPSIS",
        "       findMers [OPTIONS] <FASTA_FILE> [<KMER_FILE> [<KMER_POS_FILE>]]",
        "",
        "DESCRIPTION",
        "       Find the position(s) within FASTA file sequences of a simple text-file",
        "       list of k-mers.  The intent is that the k-mers in the list are known",
        "       to be unique, but this will be double-checked and reported if not true,",
        "       and even if not unique, all positions of each k-mer will be found.",
        "",
        "       The k-mers are required to be DNA base k-mers, i.e. alphabet ATCG.",
        "       Also, if the k-mers are unique, none should be the reverse complement",
        "       of another.  There may be additional data on each line of the k-mer",
        "       file (or the -i option's k-mer file), separated from the k-mer by at",
        "       least one space or tab character, and that data is ignored/discarded.",
        "",
        "       A k-mer position output file is written, with each line being a k-mer,",
        "       a tab, a sequence name (from the input FASTA sequence file), a tab, a",
        "       1-based k-mer start position, a tab, and a '+' or '-' indicating whether",
        "       the k-mer was found on the + or - strand.  The starting position is the",
        "       5' position on the strand where the k-mer is found.  The first line is a",
        "       header line.",
        "",
        "       OPTIONS",
        "           -v verbosity",
        "               -v0 (default) for quiet mode",
        "               -v1 for basic information output,",
        "               -v2 for that and output indicating # Mb processed",
        "               -v3 for that and verbose operational output",
        "",
        "           -i <INTERSECT_KMER_FILE>",
        "               Intersect k-mers in this file with those in <KMER_FILE> and only",
        "               report position(s) of the intersect k-mers",
        "",
        "           -n <UNSEEN_KMERS_FILE>",
        "               Output to this file a list of k-mers from input file that were not",
        "               seen in FASTA seqs",
        "",
        "           -m <MULTIPLE_KMERS_FILE>",
        "               Output to this file a list of k-mers from input file that were seen",
        "               more than once in FASTA seqs",
        "",
        "           -f <CONTIG_FILE>",
        "               Output to this file a list of contig positions and lengths for the",
        "               FASTA seqs",
        "",
        "       FASTA_FILE",
        "               Input FASTA sequences file name",
        "",
        "       KMER_FILE",
        "               Optional input k-mer text file name (size of k determined from first k-mer)",
        "               This option might not be specified if the reason for running the program is",
        "               to produce the <CONTIG_FILE> output file.",
        "",
        "       KMER_POS_FILE",
        "               Optional name of text file to receive tab-separated k-mer position",
        "               information.  If not specified, no such file is produced.",
        "",
        "AUTHOR",
        "       Written by Ted Toal.",
        NULL
        };
    for (const char *const *p = usage; *p != NULL; p++)
        cout << *p << "\n";
    exit(EXIT_SUCCESS);
    }

// Convert character c to a VAL_ constant in b if it is ATCG, and return true.
// Otherwise return false.  The character is converted to upper case.
static bool cvtBase(char c, unsigned& b)
    {
    c = toupper(c);
    if (c == 'A')
        b = VAL_A;
    else if (c == 'T')
        b = VAL_T;
    else if (c == 'C')
        b = VAL_C;
    else if (c == 'G')
        b = VAL_G;
    else
        return(false);
    return(true);
    }

// Convert k-mer string s into a binary k-mer in Kmer and return true, or return
// false if any k-mer character is not ATCG.  s length must be k.  If false is
// returned, badChar contains the non-ATCG character that was encountered.
static bool cvtKmerStrToBinary(unsigned k, const char* s, kmer& Kmer, char& badChar)
    {
    Kmer = 0;
    unsigned b;
    for (unsigned i = 0; i < k; i++)
        {
        if (!cvtBase(*s++, b))
            return(false);
        Kmer = (Kmer << BITS_PER_BASE) | b;
        }
    return(*s == 0);
    }

// Convert binary k-mer Kmer into a string k-mer in buf.  k-mer length is k,
// buf length must be at least k+1.
static void cvtKmerBinaryToStr(unsigned k, kmer Kmer, char* buf)
    {
    unsigned b;
    char c;
    for (unsigned i = 0; i < k; i++)
        {
        b = Kmer & BASE_BIT_MASK;
        Kmer >>= BITS_PER_BASE;
        if (b == VAL_A)
            c = 'A';
        else if (b == VAL_T)
            c = 'T';
        else if (b == VAL_C)
            c = 'C';
        else
            c = 'G';
        buf[k-i-1] = c;
        }
    buf[k] = 0;
    }

// Read k-mer file and initialize k-mer table.  Exit with fatal error if unable
// to read the file.  If successful, the value of k is returned (k-mer size) in
// argument k, the number of k-mers in the file in argument numKmers, the size
// of the k-mer table in bytes in argument tblSizeBytes, and the return value
// is a pointer to the k-mer table.
static ubyte* readKmerFile(const char* KMERfile, unsigned& k, unsigned& numKmers,
    unsigned& tblSizeBytes, unsigned verbosity)
    {
    kmer tblIdx;
    unsigned tblShift;

    // Read first k-mer to get k.
    InputFile kfile;
    if (!kfile.Open(KMERfile))
        {
        cout << "Can't open file " << KMERfile << "\n";
        exit(EXIT_FAILURE);
        }
    char line[MAX_K*100];
    kfile.getline(line, sizeof(line));
    if (!kfile.good())
        {
        cout << "No k-mers in file " << KMERfile << "\n";
        exit(EXIT_FAILURE);
        }
    char* token = strtok(line, " \t"); // Discard anything after first space/tab on the line.
    k = strlen(token);
    if (k < MIN_K || k > MAX_K)
        {
        cout << "First k-mer in file " << KMERfile << " gives out-of-range size k=" << k << "\n";
        exit(EXIT_FAILURE);
        }

    // Create the table.
    tblSizeBytes = (unsigned) ((ubyte8)pow(NUM_BASES, k) * TBL_BITS_PER_KMER / BITS_PER_BYTE);
    ubyte* kmerTbl = new ubyte[tblSizeBytes];
    std::fill_n(kmerTbl, tblSizeBytes, 0); // Note TBL_UNKNOWN_KMER = 0

    // Read the k-mers and initialize the table.
    numKmers = 0;
    while (kfile.good()) // Until end of file.
        {
        // First convert k-mer string in 'token' into binary k-mer in Kmer.
        kmer Kmer;
        char badChar;
        if (!cvtKmerStrToBinary(k, token, Kmer, badChar))
            {
            cout << "Unexpected k-mer character '" << badChar << "' in file " << KMERfile << "\n";
            exit(EXIT_FAILURE);
            }
        numKmers++;
        // Now initialize the table entry.
        tblIdx = Kmer >> TBL_KMER_IDX_SHIFT;
        tblShift = (Kmer & TBL_KMER_IDX_MASK)*TBL_BITS_PER_KMER;
        kmerTbl[tblIdx] |= TBL_KMER_UNSEEN << tblShift;
        // Read next k-mer.
        kfile.getline(line, sizeof(line));
        // Discard anything after first space/tab on the line.
        token = strtok(line, " \t");
        }
    if (!kfile.eof())
        {
        cout << "Unknown error reading file " << KMERfile << "\n";
        exit(EXIT_FAILURE);
        }
    // Close the k-mer file.
    kfile.Close();
    if (verbosity >= 1)
        cout << "Number of " << k << "-mers in file " << KMERfile << " is " << numKmers << "\n";
    return(kmerTbl);
    }

// Like readKmerFile(), but intersect the k-mers that are read from the file with
// the k-mers already in the table.  Set k-mers not in this file to 0 in the table.
// Exit with fatal error if unable to read the file.  If successful, the number of
// k-mers read from the file is returned in argument numKmers2, and the return value
// is the number of k-mers remaining in the table.
static unsigned intersectKmerFile(const char* KMER2file, unsigned k, ubyte* kmerTbl,
    unsigned tblSizeBytes, unsigned verbosity, unsigned& numKmers2)
    {
    kmer tblIdx;
    unsigned tblShift;

    // Read first k-mer to get k and make sure it matches argument k value.
    InputFile kfile2;
    if (!kfile2.Open(KMER2file))
        {
        cout << "Can't open file " << KMER2file << "\n";
        exit(EXIT_FAILURE);
        }
    char line[MAX_K*100];
    kfile2.getline(line, sizeof(line));
    if (!kfile2.good())
        {
        cout << "No k-mers in file " << KMER2file << "\n";
        exit(EXIT_FAILURE);
        }
    char* token = strtok(line, " \t"); // Discard anything after first space/tab on the line.
    unsigned k2 = strlen(token);
    if (k2 != k)
        {
        cout << "First k-mer in file " << KMER2file << " does not have k=" << k << " bases\n";
        exit(EXIT_FAILURE);
        }

    // Read the k-mers and change value of each that is 1 in the table to 2.
    numKmers2 = 0;
    unsigned numIsectKmers = 0;
    while (kfile2.good()) // Until end of file.
        {
        // First convert k-mer string in 'token' into binary k-mer in Kmer.
        kmer Kmer;
        char badChar;
        if (!cvtKmerStrToBinary(k, token, Kmer, badChar))
            {
            cout << "Unexpected k-mer character '" << badChar << "' in file " << KMER2file << "\n";
            exit(EXIT_FAILURE);
            }
        numKmers2++;
        // Test the table entry to see if it is 1 = TBL_KMER_UNSEEN.
        tblIdx = Kmer >> TBL_KMER_IDX_SHIFT;
        tblShift = (Kmer & TBL_KMER_IDX_MASK)*TBL_BITS_PER_KMER;
        unsigned v = (kmerTbl[tblIdx] >> tblShift) & TBL_KMER_BIT_MASK;
        if (v == TBL_KMER_UNSEEN)
            {
            // The k-mer is in the table, change the table entry to 2 = TBL_KMER_SEEN_ONCE.
            kmerTbl[tblIdx] = (kmerTbl[tblIdx] & ~(TBL_KMER_BIT_MASK << tblShift)) |
                (TBL_KMER_SEEN_ONCE << tblShift);
            ++numIsectKmers;
            }
        // Read next k-mer.
        kfile2.getline(line, sizeof(line));
        // Discard anything after first space/tab on the line.
        token = strtok(line, " \t");
        }
    if (!kfile2.eof())
        {
        cout << "Unknown error reading file " << KMER2file << "\n";
        exit(EXIT_FAILURE);
        }
    // Close the k-mer file.
    kfile2.Close();
    if (verbosity >= 1)
        {
        cout << "Number of " << k << "-mers in file " << KMER2file << " is " << numKmers2 << "\n";
        cout << "Number of " << k << "-mers in intersection is " << numIsectKmers << "\n";
        }

    // Go back through the table and change k-mers with value 1 = TBL_KMER_UNSEEN to value
    // 0 = TBL_UNKNOWN_KMER and with value 2 = TBL_KMER_SEEN_ONCE to 1 = TBL_KMER_UNSEEN.
    kmer tblSizeKmers = tblSizeBytes * TBL_KMERS_PER_BYTE;
    tblIdx = 0;
    tblShift = 0;
    for (kmer Kmer = 0; Kmer < tblSizeKmers; Kmer++)
        {
        unsigned v = (kmerTbl[tblIdx] >> tblShift) & TBL_KMER_BIT_MASK;
        if (v == TBL_KMER_UNSEEN || v == TBL_KMER_SEEN_ONCE)
            kmerTbl[tblIdx] = (kmerTbl[tblIdx] & ~(TBL_KMER_BIT_MASK << tblShift)) |
                ((v-1) << tblShift);

        tblShift += TBL_BITS_PER_KMER;
        if (tblShift == BITS_PER_BYTE)
            {
            tblShift = 0;
            tblIdx++;
            }
        }

    return(numIsectKmers);
    }

// Read one sequence from FASTA sequence file and store its contigs in a vector.
// Arguments:
//      FASTAfile: name of the FASTA file being read.
//      ffile: open file handle for the fASTA file.
//      line: pointer to buffer containing first line of sequence to read (ID line).
//      lineSize: size of line buffer.
//      totalNumBases: counter of total number of bases read from FASTA file.
//      seqID: place to return the sequence ID for the sequence that is read.
// Returns: pointer to a seqcontigs vector containing the contigs making up the
// sequence.  Use discardSeqContigs() to discard it when finished.
// On return, line contains the sequence line of the next sequence in the file,
// if any, and totalNumBases has been updated.
static seqcontigs* readFASTAseq(const char* FASTAfile, InputFile& ffile,
    char* line, size_t lineSize, ubyte8& totalNumBases, S& seqID, unsigned verbosity)
    {
    // Make sure line contains a sequence ID line.
    if (line[0] != '>')
        {
        cout << "Expected FASTA sequence ID line in " << FASTAfile << " but got:\n" << line << "\n";
        exit(EXIT_FAILURE);
        }
    // Extract sequence ID.
    for (int i = 1; line[i] != 0; i++)
        {
        if (line[i] == ' ' || line[i] == '\t' || line[i] == '|')
            break;
        seqID += line[i];
        }
    if (seqID.length() == 0)
        {
        cout << "Bad FASTA sequence ID in " << FASTAfile << ":\n" << line << "\n";
        exit(EXIT_FAILURE);
        }
    if (verbosity >= 3)
        cout << "Seq ID " << seqID << "\n";

    // Read the entire sequence and convert the bases to 2 bits.  If we
    // encounter non-ATCG characters (not just N but anything, assume there
    // may be IUPAC base codes too, which we ignore), split the sequence
    // into a new contig at that point.  Accumulate the sequence bits in
    // byte vectors, one vector per sequence contig.  Also track the contig
    // starting position and length.
    seqcontigs* seqContigs = new seqcontigs;
    unsigned curContigNum = 0;
    unsigned curContigByteNum = 0;
    unsigned curContigBitNum = 0;
    unsigned seqPos = 0;
    ubyte8 numBases = 0;
    ffile.getline(line, lineSize); // Read first sequence line.
    if (!ffile.good() && line[0] != 0)
        ffile.clear(); // Long line will fill buffer before \n and set error flag.
    while (ffile.good() && line[0] != '>')
        {
        // Process each FASTA sequence character and add it to the current
        // sequence contig if it is ATCG.
        for (int i = 0; line[i] != 0; i++)
            {
            unsigned b;
            numBases++; // Count number of bases.
            totalNumBases++; // Count total number of bases.
            if (verbosity >= 2)
                if (totalNumBases % 10000000 == 0)
                    cout << "  " << (totalNumBases/1000000) << " Mb\n";
            seqPos++; // Advance sequence position, first position is 1.
            if (!cvtBase(line[i], b))
                {
                if (curContigNum < seqContigs->size()) // If currently have a contig's memory allocated...
                    curContigNum++; // Advance to a new contig.
                }
            else
                {
                // If currently don't have a contig's vector allocated, do it.
                if (curContigNum == seqContigs->size())
                    {
                    seqcontig SF;
                    SF.startPos = seqPos;
                    SF.contigLen = 0;
                    SF.contigBytes = new seqcontigbytes;
                    seqContigs->push_back(SF);
                    curContigByteNum = 0;
                    }
                // If don't have the current contig byte added to the vector, do it.
                if (curContigByteNum == (*seqContigs)[curContigNum].contigBytes->size())
                    {
                    (*seqContigs)[curContigNum].contigBytes->push_back(0);
                    curContigBitNum = 0;
                    }
                // Insert the new base into the current contig byte at correct bit position.
                (*(*seqContigs)[curContigNum].contigBytes)[curContigByteNum] |= (ubyte)(b << curContigBitNum);
                (*seqContigs)[curContigNum].contigLen++;
                curContigBitNum += BITS_PER_BASE;
                if (curContigBitNum == BITS_PER_BYTE) // End of byte reached.
                    {
                    curContigByteNum++; // Defer byte allocation until we have something to put in it.
                    curContigBitNum = 0;
                    }
                }
            }
        // Read next sequence line (or ID line of following sequence).
        ffile.getline(line, lineSize);
        if (!ffile.good() && line[0] != 0)
            ffile.clear(); // Long line will fill buffer before \n and set error flag.
        }
    if (verbosity >= 3)
        cout << "Finished reading seq, #contigs = " << seqContigs->size() << " #bases = " << numBases << "\n";
    return(seqContigs);
    }

// Given a list of sequence contigs from a FASTA file sequence, for each base
// of each contig, form binary k-mers and their reverse complements, look them
// up in a k-mer table, and if found modify the table to indicate so and write
// out the k-mer to a k-mer position output file.
// Arguments:
//      k: k-mer size.
//      kmerTbl: pointer to k-mer table.
//      seqID: name of the FASTA sequence being processed.
//      seqContigs: vector of sequence contigs to process.
//      KMERPOSfile: name of output file to which to write k-mers found in the k-mer table.
//      kmerposfile: open file handle for KMERPOSfile.
// On return, the table has been modified to reflect the k-mers that were found
// in the FASTA sequence.
static void processSeqContigs(unsigned k, ubyte *kmerTbl, S& seqID, seqcontigs* seqContigs,
    const char* KMERPOSfile, OutputFile& kmerposfile, unsigned verbosity)
    {
    unsigned numKmersFound = 0;
    kmer kmerMask = ((kmer)1 << (k*BITS_PER_BASE)) - 1; // 1-bits in lowest k bit-pairs.
    unsigned bitShiftUpperKmer = ((k-1)*BITS_PER_BASE); // Left shift to get bits into highest k-mer base.
    for (unsigned i = 0; i < seqContigs->size(); i++)
        {
        unsigned startPos = (*seqContigs)[i].startPos;
        unsigned contigLen = (*seqContigs)[i].contigLen;
        unsigned curContigByteNum = 0;
        unsigned curContigBitNum = 0;
        ubyte curContigByte;
        kmer Kmer = 0; // This will hold a travelling k-mer, travelling along the sequence contig.
        kmer revKmer = 0; // Reverse complement of Kmer.
        // The following loop is where this program will be most of the time.
        // It needs to be very optimized, and I'm afraid it isn't good enough.
        for (unsigned j = 0; j < contigLen; j++)
            {
            // If at first base of a contig byte, get the byte.
            if (curContigBitNum == 0)
                curContigByte = (*(*seqContigs)[i].contigBytes)[curContigByteNum];
            // Get next base from the sequence contig.
            kmer b = (curContigByte >> curContigBitNum);
            kmer bRev = (~b) & BASE_BIT_MASK;
            b &= BASE_BIT_MASK;
            // Shift the new base into the two travelling k-mers.
            Kmer = ((Kmer << BITS_PER_BASE) | b) & kmerMask;
            revKmer = (revKmer >> BITS_PER_BASE) | (bRev << bitShiftUpperKmer);
            // Ignore the first k-1 k-mers as we haven't accumulated enough bases yet.
            if (j >= k-1) // At j == k-1 we have k bases.
                {
                // Compute k-mer table index and shift amount for Kmer.
                kmer tblIdx = Kmer >> TBL_KMER_IDX_SHIFT;
                unsigned tblShift = (Kmer & TBL_KMER_IDX_MASK)*TBL_BITS_PER_KMER;
                // Look for Kmer in the table, and if not there, try revKmer,
                // since the table contains canonical k-mers using the lexically
                // smaller of a k-mer and its reverse complement.
                unsigned v = (kmerTbl[tblIdx] >> tblShift) & TBL_KMER_BIT_MASK;
                if (v != TBL_UNKNOWN_KMER)
                    {
                    numKmersFound++;
                    if (v != TBL_KMER_SEEN_MANY)
                        {
                        // Increment k-mer counter in table entry.
                        v++;
                        kmerTbl[tblIdx] = (kmerTbl[tblIdx] & ~(TBL_KMER_BIT_MASK << tblShift)) |
                            (v << tblShift);
                        }
                    if (KMERPOSfile != NULL)
                        {
                        char buf[MAX_K+1];
                        cvtKmerBinaryToStr(k, Kmer, buf);
                        unsigned pos = (startPos+j-k+1);
                        kmerposfile << buf << "\t" << seqID << "\t" << pos << "\t+\t" << (i+1) << "\t" << (j+1) << "\n";
                        }
                    }
                else
                    {
                    tblIdx = revKmer >> TBL_KMER_IDX_SHIFT;
                    tblShift = (revKmer & TBL_KMER_IDX_MASK)*TBL_BITS_PER_KMER;
                    v = (kmerTbl[tblIdx] >> tblShift) & TBL_KMER_BIT_MASK;
                    if (v != TBL_UNKNOWN_KMER)
                        {
                        numKmersFound++;
                        if (v != TBL_KMER_SEEN_MANY)
                            {
                            // Increment k-mer counter in table entry.
                            v++;
                            kmerTbl[tblIdx] = (kmerTbl[tblIdx] & ~(TBL_KMER_BIT_MASK << tblShift)) |
                                (v << tblShift);
                            }
                        if (KMERPOSfile != NULL)
                            {
                            char buf[MAX_K+1];
                            cvtKmerBinaryToStr(k, revKmer, buf);
                            unsigned pos = (startPos+j);
                            kmerposfile << buf << "\t" << seqID << "\t" << pos << "\t-\t" << (i+1) << "\t" << (j+1) << "\n";
                            }
                        }
                    }
                }
            // Advance to next contig base, and to next contig byte when needed.
            curContigBitNum += BITS_PER_BASE;
            if (curContigBitNum == BITS_PER_BYTE) // End of byte reached.
                {
                curContigByteNum++;
                curContigBitNum = 0;
                }
            }
        }
    if (verbosity >= 3)
        cout << " # k-mers found in seq: " << numKmersFound << "\n";
    }

// Write the position and length of a vector of contigs to contig output file.
// Arguments:
//      seqID: name of the FASTA sequence containing the contigs.
//      seqContigs: vector of sequence contigs to output.
//      CONTIGfile: name of contig file.
//      contigfile: open file handle for CONTIGfile.
static void outputContigInfo(S& seqID, seqcontigs* seqContigs, const char* CONTIGfile,
    OutputFile& contigfile)
    {
    unsigned contigNum = 0;
    for (seqcontigs::iterator it = seqContigs->begin() ; it != seqContigs->end(); ++it)
        {
        contigNum++;
        contigfile << seqID << "\t" << it->startPos << "\t" << contigNum << "\t" << it->contigLen << "\n";
        }
    }

// Read input FASTA sequence file, one sequence at a time, run through the
// sequence and generate each k-mer and reverse-complement k-mer, look them up
// in the k-mer table, and if found adjust the table entry and write the k-mer
// and its position to the k-mer position output file.
// If kmerTbl is NULL, no per-k-mer processing is done, but a contig file is
// still produced if CONTIGfile is not NULL.
// Arguments:
//      k: k-mer size, 0 if kmerTbl is NULL.
//      kmerTbl: pointer to k-mer table, or NULL to not search for k-mers.
//      FASTAfile: name of FASTA file to read.
//      CONTIGfile: name of contig file to create, or NULL to not create it.
//      KMERPOSfile: name of k-mer position file to create, or NULL to not create it.
// On return, the table has been modified to reflect the k-mers that were found
// in the FASTA sequences.
static void processFASTAfile(unsigned k, ubyte *kmerTbl, const char* FASTAfile,
    const char* CONTIGfile, const char* KMERPOSfile, unsigned verbosity)
    {
    // Open FASTA file and read first line = ID line of first sequence.
    InputFile ffile;
    if (!ffile.Open(FASTAfile))
        {
        cout << "Can't open file " << FASTAfile << "\n";
        exit(EXIT_FAILURE);
        }
    char line[1000];
    ffile.getline(line, sizeof(line));
    if (!ffile.good())
        {
        cout << "No sequences in file " << FASTAfile << "\n";
        exit(EXIT_FAILURE);
        }

    // If contig info is to be written to a file, create the file and write header line.
    OutputFile contigfile;
    if (CONTIGfile != NULL)
        {
        if (!contigfile.Create(CONTIGfile))
            {
            cout << "Can't create file " << CONTIGfile << "\n";
            exit(EXIT_FAILURE);
            }
        contigfile << "seqID\tpos\tcontig\tlen\n";
        }

    // Create k-mer position output file and write header line.
    OutputFile kmerposfile;
    if (KMERPOSfile != NULL)
        {
        if (!kmerposfile.Create(KMERPOSfile))
            {
            cout << "Can't create file " << KMERPOSfile << "\n";
            exit(EXIT_FAILURE);
            }
        kmerposfile << "kmer\tseqID\tpos\tstrand\tcontig\tcontigPos\n";
        }

    // Loop processing input FASTA sequences until there are no more.
    ubyte8 totalNumBases = 0;
    unsigned numSeqs = 0;
    while (ffile.good()) // Until end of file.
        {
        // Read the sequence into a vector of sequence contigs.
        S seqID;
        seqcontigs* seqContigs = readFASTAseq(FASTAfile, ffile, line, sizeof(line),
            totalNumBases, seqID, verbosity);
        numSeqs++;

        // Write contig info to file if requested.
        if (CONTIGfile != NULL)
            outputContigInfo(seqID, seqContigs, CONTIGfile, contigfile);

        // Search for the k-mers of each contig in the k-mer table and modify
        // the table and output the k-mers when they are found, if kmerTbl isn't
        // NULL.
        if (kmerTbl != NULL)
            processSeqContigs(k, kmerTbl, seqID, seqContigs, KMERPOSfile, kmerposfile, verbosity);

        // Discard the contigs that were read.
        discardSeqContigs(seqContigs);
        }
    if (!ffile.eof())
        {
        cout << "Unknown error reading file " << FASTAfile << "\n";
        exit(EXIT_FAILURE);
        }

    // Close the FASTA file.
    ffile.Close();
    // Close the k-mer position output file.
    if (KMERPOSfile != NULL)
        kmerposfile.Close();
    // Close the contig info output file.
    if (CONTIGfile != NULL)
        contigfile.Close();

    // Report stats.
    if (verbosity >= 1)
        {
        cout << "Number of FASTA sequences processed: " << numSeqs << "\n";
        cout << "Total number of bases processed: " << totalNumBases << "\n";
        if (CONTIGfile != NULL)
            cout << "FASTA sequence contig info written to file " << CONTIGfile << "\n";
        if (KMERPOSfile != NULL)
            cout << "k-mer positions written to file " << KMERPOSfile << "\n";
        }
    }

// List any k-mers not seen or seen more than once.
static void writeUnseenAndMultipleKmersAndSummary(unsigned k, ubyte* kmerTbl,
    unsigned tblSizeBytes, unsigned verbosity, const char* UNSEEN_KMERSfile, const char* MULTIPLE_KMERSfile)
    {
    kmer tblSizeKmers = tblSizeBytes * TBL_KMERS_PER_BYTE;
    kmer tblIdx = 0;
    unsigned tblShift = 0;
    unsigned numInputKmers = 0;
    unsigned numNotFound = 0;
    unsigned numFoundOnce = 0;
    unsigned numFoundMore = 0;

    // Create unseen k-mer output file.
    OutputFile unseen_kmers_file;
    if (UNSEEN_KMERSfile != NULL)
        {
        if (!unseen_kmers_file.Create(UNSEEN_KMERSfile))
            {
            cout << "Can't create file " << UNSEEN_KMERSfile << "\n";
            exit(EXIT_FAILURE);
            }
        }

    // Create multiple k-mer output file.
    OutputFile multiple_kmers_file;
    if (MULTIPLE_KMERSfile != NULL)
        {
        if (!multiple_kmers_file.Create(MULTIPLE_KMERSfile))
            {
            cout << "Can't create file " << MULTIPLE_KMERSfile << "\n";
            exit(EXIT_FAILURE);
            }
        }

    // Search for unseen and multiple k-mers.
    for (kmer Kmer = 0; Kmer < tblSizeKmers; Kmer++)
        {
        char buf[MAX_K+1];
        // Get the table entry.
        unsigned v = (kmerTbl[tblIdx] >> tblShift) & TBL_KMER_BIT_MASK;
        if (v != TBL_UNKNOWN_KMER)
            numInputKmers++;
        if (v == TBL_KMER_UNSEEN)
            numNotFound++;
        else if (v == TBL_KMER_SEEN_ONCE)
            numFoundOnce++;
        else if (v == TBL_KMER_SEEN_MANY)
            numFoundMore++;
        if (UNSEEN_KMERSfile != NULL && v == TBL_KMER_UNSEEN)
            {
            cvtKmerBinaryToStr(k, Kmer, buf);
            unseen_kmers_file << buf << "\n";
            }
        else if (MULTIPLE_KMERSfile != NULL && v == TBL_KMER_SEEN_MANY)
            {
            cvtKmerBinaryToStr(k, Kmer, buf);
            multiple_kmers_file << buf << "\n";
            }
        // Next k-mer.
        tblShift += TBL_BITS_PER_KMER;
        if (tblShift == BITS_PER_BYTE)
            {
            tblShift = 0;
            tblIdx++;
            }
        }

    // Close the unseen k-mer output file.
    if (UNSEEN_KMERSfile != NULL)
        unseen_kmers_file.Close();

    // Close the multiple k-mer output file.
    if (MULTIPLE_KMERSfile != NULL)
        multiple_kmers_file.Close();

    // Final summary.
    if (verbosity >= 1)
        {
        cout << "Number of input list " << k << "-mers:                                    " << numInputKmers << "\n";
        cout << "Number of input list " << k << "-mers NOT FOUND in FASTA seqs:            " << numNotFound << "\n";
        cout << "Number of input list " << k << "-mers found ONCE in FASTA seqs:           " << numFoundOnce << "\n";
        cout << "Number of input list " << k << "-mers found MORE THAN ONCE in FASTA seqs: " << numFoundMore << "\n";
        }
    }

// Return a pointer to argv[a1][a2+1], or if that is 0, to argv[a1+1] if it exists.  Modify a1 and a2 to advance
// them to the end of the string.  If string cannot be extracted, call usage() and don't return.
static const char* getOptionValueString(int argc, char* const argv[], int& a1, int& a2)
    {
    const char* p;
    a2++;
    if (argv[a1][a2] != 0)
        {
        p = &argv[a1][a2];
        a2 += strlen(p) - 1;
        }
    else
        {
        a1++;
        if (argc <= a1)
            usage();
        p = argv[a1];
        a2 = strlen(p) - 1;
        }
    return(p);
    }

// Main function for findMers program.
int main(int argc, char* const argv[])
    {
    char* q;

    // Option variables.
    unsigned verbosity = 0;                 // -v
    const char* UNSEEN_KMERSfile = NULL;    // -n
    const char* MULTIPLE_KMERSfile = NULL;  // -m
    const char* CONTIGfile = NULL;          // -f
    const char* KMER2file = NULL;           // -i

    // Options.
    int a1; // argv first index.
    int a2; // argv second index.
    for (a1 = 1; a1 < argc && argv[a1][0] == '-'; a1++)
        {
        for (a2 = 1; argv[a1][0] == '-' && argv[a1][a2] != 0; a2++)
            switch (argv[a1][a2])
                {
                case 'v':
                    verbosity = (unsigned) strtoul(getOptionValueString(argc, argv, a1, a2), &q, 0);
                    if (*q != 0)
                        usage();
                    break;
                case 'n':
                    UNSEEN_KMERSfile = getOptionValueString(argc, argv, a1, a2);
                    break;
                case 'm':
                    MULTIPLE_KMERSfile = getOptionValueString(argc, argv, a1, a2);
                    break;
                case 'f':
                    CONTIGfile = getOptionValueString(argc, argv, a1, a2);
                    break;
                case 'i':
                    KMER2file = getOptionValueString(argc, argv, a1, a2);
                    break;
                }
        }

    // Usage.
    if (argc < a1+1)
        usage();

    // Arguments.
    const char* FASTAfile = argv[a1+0];
    const char* KMERfile = NULL;
    const char* KMERPOSfile = NULL;
    if (argc > a1+1)
        KMERfile = argv[a1+1];
    if (argc > a1+2)
        KMERPOSfile = argv[a1+2];

    /*
        1. If KMERfile was specified, read the k-mer file and initialize k-mer table.
    */
    unsigned k = 0;
    unsigned numKmers = 0;
    unsigned tblSizeBytes = 0;
    ubyte* kmerTbl = NULL;
    if (KMERfile != NULL)
        kmerTbl = readKmerFile(KMERfile, k, numKmers, tblSizeBytes, verbosity);

    /*
        2. If k-mer intersection option, read and intersect second k-mer file.
    */
    unsigned numKmers2 = 0;
    unsigned numIsectKmers = 0;
    if (KMERfile != NULL && KMER2file != NULL)
        numIsectKmers = intersectKmerFile(KMER2file, k, kmerTbl, tblSizeBytes, verbosity, numKmers2);

    /*
        3. Read input FASTA sequence file and look up k-mers in table.
    */
    processFASTAfile(k, kmerTbl, FASTAfile, CONTIGfile, KMERPOSfile, verbosity);

    /*
        4. If a k-mer table was produced, display summary info and write out k-mers not seen or seen more than once.
    */
    if (kmerTbl != NULL)
        writeUnseenAndMultipleKmersAndSummary(k, kmerTbl, tblSizeBytes, verbosity, UNSEEN_KMERSfile, MULTIPLE_KMERSfile);

    return 0;
    }
