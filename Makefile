# This is a Makefile for running the IGGPIPE pipeline on genome files to produce
# candidate IGG markers.  Before running this, edit file allParameters.template.
# Run this Makefile with command 'make PARAMS=<allParameters filename>' to get
# further usage instructions, e.g. 'make PARAMS=allParameters.myFile'.

# Delete created files if error occurs creating them.
.DELETE_ON_ERROR:

# Run all commands of a target build using a single shell invocation.
.ONESHELL:

# We will use secondary expansion.  This affects only prerequisities, and then,
# only if the prerequisite includes $$.  If it does, the first (normal) expansion
# changes $$ to a single $.  Then, the second expansion expands that $, which
# would normally be followed by something that included $* (the "stem" match of
# a static pattern rule.  Thus, $$(VAR_NAME_PREFIX$$*) would allow a nested
# variable reference where the variable name ends in a value equal to the stem
# of the static pattern rule.  We use the genome number as the stem value here.
.SECONDEXPANSION:

# Use this shell.
SHELL=/bin/sh

# Some useful variables.
EMPTY :=
SPACE := $(EMPTY) $(EMPTY)
INDENT := $(EMPTY)    $(EMPTY)

# Cancel C-compiler implicit rules so it doesn't try to compile findMers.cpp.
% : %.cpp
%.o : %.cpp
% : %.o

################################################################################
# Include path parameters definition file.
################################################################################
include allPathParameters.ours

################################################################################
# Include main parameter definition file (there might not be one).
################################################################################
include $(PARAMS)

################################################################################
# Default make target.
# If variable PARAMS is not defined, show basic usage info, else make target ALL.
################################################################################
ifeq ($(PARAMS),)
all: basic_usage
else
all: ALL
endif

################################################################################
# 'basic_usage' target, show basic usage info.
################################################################################
ifeq ($(PARAMS),)
basic_usage:
	@echo
	@echo "For basic usage info, use this make command:"
	@echo
	@echo "$(INDENT)make usage"
	@echo
	@echo "Summarizing the usage info, you must specify a parameter file using PARAMS in"
	@echo "the make command, and either specify target ALL to run the entire pipeline or",
	@echo "specify a pipeline step name.  For example:"
	@echo 
	@echo "$(INDENT)make PARAMS=allParameters.test ALL"
	@echo 
endif

################################################################################
# If user invokes 'make usage', show usage info.
################################################################################
usage:
	more help.txt

################################################################################
# If user invokes 'make clean', give him some instructions on how to clean.
################################################################################
clean:
	@echo "Use 'make PARAMS=filename CLEAN=1 ALL' to remove files made with parameter file 'filename'"
	@echo "Use 'make PARAMS=filename CLEAN_OUT_DIR=1 ALL' to remove all output directory files."

################################################################################
# Set variables for selecting genome.
################################################################################
# If variable GENOME was not specified on 'make' command line, define it as 'ALL'
GENOME ?= ALL

# If variable GENOME is not 'ALL' or a valid genome number, exit with error.
ifneq ($(GENOME),ALL)
ifeq ($(GENOME_$(GENOME)),)
$(error GENOME must be either a valid genome number (1..$(N_GENOMES)) or 'ALL'))
endif
endif

################################################################################
# Set variables for cleaning operations.
################################################################################
# If variable CLEAN was defined as anything at all, set it to "CLEAN=1" so it can
# be easily used with recursive invocations of make.
ifneq ($(CLEAN),)
override CLEAN := CLEAN=1
endif

################################################################################
# Set variables for command timing.
################################################################################
# If variable TIME_CMDS is YES, set TIME and REDIR to time commands and redirect
# stderr to stdout.  If not YES, leave these empty.
ifeq ($(TIME_CMDS),YES)
TIME := (time
REDIR := ) 2>&1
endif

################################################################################
# General pipeline targets.
################################################################################

# Target for running or cleaning entire pipeline or cleaning entire output directory.
ifeq ($(CLEAN_OUT_DIR),1)
ALL:
	@echo
	@echo "Removing all files from the output directory."
	@# Try to remove output directory first.  If trashing, entire thing will be in trash.
	-$(CMD_DELETE_WHEN_CLEANING) $(DIR_IGGPIPE_OUT)
	-$(CMD_DELETE_WHEN_CLEANING) $(DIR_GENOME_OUT_DATA)/*
	-$(CMD_DELETE_WHEN_CLEANING) $(DIR_KMERS)/*
	-$(CMD_DELETE_WHEN_CLEANING) $(DIR_PRIMER_DATA)/*
	-$(CMD_DELETE_WHEN_CLEANING) $(DIR_IGGPIPE_OUT)/*
	-$(CMD_DELETE_WHEN_CLEANING) $(DIR_IGGPIPE_OUT)
	@echo "Removed all files from the output directory."
	@echo
else ifeq ($(CLEAN),1)
ALL: all_getSeqInfo all_getKmers all_kmerStats all_sortKmers all_getContigFile \
    kmerIsect all_getGenomicPos all_mergeKmers getCommonUniqueKmers findLCRs \
    findIndelGroups all_getDNAseqs findPrimers all_ePCRtesting removeBadMarkers \
    plotMarkers
	@echo
	@echo "ALL files are cleaned"
	@echo
else
ALL: all_getSeqInfo all_getKmers all_kmerStats all_sortKmers all_getContigFile \
    kmerIsect all_getGenomicPos all_mergeKmers getCommonUniqueKmers findLCRs \
    findIndelGroups all_getDNAseqs findPrimers all_ePCRtesting removeBadMarkers \
    plotMarkers
	@echo
	@echo "ALL files are up to date"
	@echo
endif

# Target for getting or cleaning sequence info for all genomes.
all_getSeqInfo:
	@$(MAKE) PARAMS=$(PARAMS) $(CLEAN) getSeqInfo GENOME=ALL

# Target for getting or cleaning contig files for all genomes.
all_getContigFile:
	@$(MAKE) PARAMS=$(PARAMS) $(CLEAN) getContigFile GENOME=ALL

# Target for getting or cleaning unique k-mers from all genomes.
all_getKmers:
	@$(MAKE) PARAMS=$(PARAMS) $(CLEAN) getKmers GENOME=ALL

# Target for getting or cleaning unique k-mer statistics for all genomes.
all_kmerStats:
	@$(MAKE) PARAMS=$(PARAMS) $(CLEAN) kmerStats GENOME=ALL

# Target for getting or cleaning sorted unique k-mers for all genomes.
all_sortKmers:
	@$(MAKE) PARAMS=$(PARAMS) $(CLEAN) sortKmers GENOME=ALL

# Target for getting or cleaning unique k-mer genomic position for all genomes.
all_getGenomicPos:
	@$(MAKE) PARAMS=$(PARAMS) $(CLEAN) getGenomicPos GENOME=ALL

# Target for merging common unique k-mer genomic positions files for all genomes.
all_mergeKmers:
	@$(MAKE) PARAMS=$(PARAMS) $(CLEAN) mergeKmers GENOME=ALL

# Target for getting or cleaning DNA for making primers for all genomes.
all_getDNAseqs:
	@$(MAKE) PARAMS=$(PARAMS) $(CLEAN) getDNAseqs GENOME=ALL

# Target for testing primers using ePCR on all genomes.
all_ePCRtesting:
	@$(MAKE) PARAMS=$(PARAMS) $(CLEAN) ePCRtesting GENOME=ALL

################################################################################
# Create output directories.
################################################################################

# Target for making the main output directory.
$(DIR_IGGPIPE_OUT):
	@echo
	@echo "Creating directory $(DIR_IGGPIPE_OUT)"
	mkdir -p $(DIR_IGGPIPE_OUT)

# Target for making the genome data output directory.
$(DIR_GENOME_OUT_DATA):
	@echo
	@echo "Creating directory $(DIR_GENOME_OUT_DATA)"
	mkdir -p $(DIR_GENOME_OUT_DATA)

# Target for making the k-mer output directory.
$(DIR_KMERS):
	@echo
	@echo "Creating directory $(DIR_KMERS)"
	mkdir -p $(DIR_KMERS)

# Target for making the primer data directory.
$(DIR_PRIMER_DATA):
	@echo
	@echo "Creating directory $(DIR_PRIMER_DATA)"
	mkdir -p $(DIR_PRIMER_DATA)

################################################################################
# getSeqInfo: Get sequence ID and length information.  Argument: GENOME
################################################################################

# A list of all FASTA files to be analyzed.
FASTA_FILES := $(foreach G,$(GENOME_NUMBERS),$(PATH_GENOME_FASTA_$(G)))

# A list of all .idlen files to be produced.
IDLEN_FILES := $(foreach G,$(GENOME_NUMBERS),$(PFX_GENOME_DATA_FILE)$(G).idlens)

# Target .idlen file(s) for genome GENOME.
ifeq ($(GENOME),ALL)
TARGET_IDLEN := $(IDLEN_FILES)
else
TARGET_IDLEN := $(PFX_GENOME_DATA_FILE)$(GENOME).idlens
endif

# Phony target to make or clean TARGET_IDLEN file(s).
ifeq ($(CLEAN),)
getSeqInfo: $(TARGET_IDLEN)
	@echo
	@echo "getSeqInfo file(s) for genome(s) $(GENOME) are up to date."
else
getSeqInfo:
	@$(CMD_DELETE_WHEN_CLEANING) $(TARGET_IDLEN)
	@echo "getSeqInfo output file(s) for genome(s) $(GENOME) removed."
endif

# IDLEN_FILES is multiple targets, one per genome.
# Here, % is a genome number.

$(IDLEN_FILES) : $(PFX_GENOME_DATA_FILE)%.idlens : $$(PATH_GENOME_FASTA_$$*) | \
        $(DIR_GENOME_OUT_DATA) $(PATH_PERL) $(PATH_EXTRACT_SEQ_IDS)
	@echo
	@echo "*** getSeqInfo PARAMS=$(PARAMS) $(CLEAN) GENOME=$* ***"
	@echo "Extracting sequence IDs and lengths from $< into $@"
	$(TIME) $(PATH_PERL) $(PATH_EXTRACT_SEQ_IDS) $< $@ $(REDIR)
	@echo "Finished."

################################################################################
# getContigFile: Get sequence contig positions and lengths.  Argument: GENOME
################################################################################

# A list of all FASTA files already defined above: FASTA_FILES

# A list of all .contigs files to be produced.
CONTIG_FILES := $(foreach G,$(GENOME_NUMBERS),$(PFX_GENOME_DATA_FILE)$(G).contigs)

# Target .contigs file(s) for genome GENOME.
ifeq ($(GENOME),ALL)
TARGET_CONTIG := $(CONTIG_FILES)
else
TARGET_CONTIG := $(PFX_GENOME_DATA_FILE)$(GENOME).contigs
endif

# Phony target to make or clean TARGET_CONTIG file(s).
ifeq ($(CLEAN),)
getContigFile: $(TARGET_CONTIG)
	@echo
	@echo "getContigFile file(s) for genome(s) $(GENOME) are up to date."
else
getContigFile:
	@$(CMD_DELETE_WHEN_CLEANING) $(TARGET_CONTIG)
	@echo "getContigFile output file(s) for genome(s) $(GENOME) removed."
endif

# CONTIG_FILES is multiple targets, one per genome.
# Here, % is a genome number.

$(CONTIG_FILES) : $(PFX_GENOME_DATA_FILE)%.contigs : $$(PATH_GENOME_FASTA_$$*) | \
        $(DIR_GENOME_OUT_DATA) $(PATH_FINDMERS)
	@echo
	@echo "*** getContigFile PARAMS=$(PARAMS) $(CLEAN) GENOME=$* ***"
	@echo "Extracting contig positions and lengths from $< into $@"
	$(TIME) $(PATH_FINDMERS) $(ARGS_FINDMER) -f $@ $< $(REDIR)
	@echo "Finished."

################################################################################
# getKmers: Get unique k-mers.  Argument: GENOME
################################################################################

# A list of all FASTA files already defined above: FASTA_FILES

# A list of all .kmers files to be produced.
KMERS_FILES := $(foreach G,$(GENOME_NUMBERS),$(PFX_KMERS_DATA_FILE)$(G).kmers)

# Target .kmers file(s) for genome GENOME.
ifeq ($(GENOME),ALL)
TARGET_KMERS := $(KMERS_FILES)
else
TARGET_KMERS := $(PFX_KMERS_DATA_FILE)$(GENOME).kmers
endif

# Phony target to make or clean TARGET_KMERS file(s).
ifeq ($(CLEAN),)
getKmers: $(TARGET_KMERS)
	@echo
	@echo "getKmers file(s) for genome(s) $(GENOME) are up to date."
else
getKmers:
	@$(CMD_DELETE_WHEN_CLEANING) $(TARGET_KMERS)
	@echo "getKmers output file(s) for genome(s) $(GENOME) removed."
endif

# KMERS_FILES is multiple targets, one per genome.
# Here, % is a genome number.

$(KMERS_FILES) : $(PFX_KMERS_DATA_FILE)%.kmers : $$(PATH_GENOME_FASTA_$$*) | \
        $(DIR_KMERS) $(PATH_JELLYFISH)
	@echo
	@echo "*** getKmers PARAMS=$(PARAMS) $(CLEAN) GENOME=$* ***"
	@echo "Extracting unique $(K)-mers from $< into $(PFX_KMERS_DATA_FILE)$*.kmers_*"
	$(TIME) $(PATH_JELLYFISH) count -C -m $(K) -s $(JELLYFISH_HASH_SIZE) -U 1 \
	    -o $(PFX_KMERS_DATA_FILE)$*.kmers $< $(REDIR)
	@echo "Finished."

################################################################################
# kmerStats: Get unique k-mer statistics.  Argument: GENOME
################################################################################

# A list of all .kmers files already defined above: KMERS_FILES

# A list of all .stats files to be produced.
KMERS_STATS_FILES := $(foreach G,$(GENOME_NUMBERS),$(PFX_KMERS_DATA_FILE)$(G).stats)

# Target .stats file(s) for genome GENOME.
ifeq ($(GENOME),ALL)
TARGET_STATS := $(KMERS_STATS_FILES)
else
TARGET_STATS := $(PFX_KMERS_DATA_FILE)$(GENOME).stats
endif

# Phony target to make or clean TARGET_STATS file(s).
ifeq ($(CLEAN),)
kmerStats: $(TARGET_STATS)
	@echo
	@echo "kmerStats file(s) for genome(s) $(GENOME) are up to date."
else
kmerStats:
	@$(CMD_DELETE_WHEN_CLEANING) $(TARGET_STATS)
	@echo "kmerStats output file(s) for genome(s) $(GENOME) removed."
endif

# KMERS_STATS_FILES is multiple targets, one per genome.
# Here, % is a genome number.

$(KMERS_STATS_FILES) : $(PFX_KMERS_DATA_FILE)%.stats : $(PFX_KMERS_DATA_FILE)%.kmers | \
        $(DIR_KMERS)
	@echo
	@echo "*** kmerStats PARAMS=$(PARAMS) $(CLEAN) GENOME=$* ***"
	@echo "Getting statistics for $(K)-mers from $< into $@"
	$(TIME) $(PATH_JELLYFISH) stats -v $< >$@ $(REDIR)
	@echo "Finished."
	@echo
	@echo "Statistics are:"
	@cat $@
	@echo

################################################################################
# sortKmers: Get sorted unique k-mers.  Argument: GENOME
################################################################################

# A list of all .kmers files already defined above: KMERS_FILES

# A list of all .sorted files to be produced.
SORTED_KMERS_FILES := $(foreach G,$(GENOME_NUMBERS),$(PFX_KMERS_DATA_FILE)$(G).sorted)

# Target .sorted file(s) for genome GENOME.
ifeq ($(GENOME),ALL)
TARGET_SORTED_KMERS := $(SORTED_KMERS_FILES)
else
TARGET_SORTED_KMERS := $(PFX_KMERS_DATA_FILE)$(GENOME).sorted
endif

# Phony target to make or clean TARGET_SORTED_KMERS file(s).
ifeq ($(CLEAN),)
sortKmers: $(TARGET_SORTED_KMERS)
	@echo
	@echo "sortKmers file(s) for genome(s) $(GENOME) are up to date."
else
sortKmers:
	@$(CMD_DELETE_WHEN_CLEANING) $(TARGET_SORTED_KMERS)
	@echo "sortKmers output file(s) for genome(s) $(GENOME) removed."
endif

# SORTED_KMERS_FILES is multiple targets, one per genome.
# Here, % is a genome number.

$(SORTED_KMERS_FILES) : $(PFX_KMERS_DATA_FILE)%.sorted : $(PFX_KMERS_DATA_FILE)%.kmers | \
        $(DIR_KMERS)
	@echo
	@echo "*** sortKmers PARAMS=$(PARAMS) $(CLEAN) GENOME=$* ***"
	@echo "Sorting text-based $(K)-mers from $< into $@"
	@echo
	@echo "Extracting text-based $(K)-mers from binary file $< into $@.tmp"
	$(TIME) $(PATH_JELLYFISH) dump -c -o $@.tmp $< $(REDIR)
	@echo "Sorting text-based $(K)-mers from $@.tmp into $@.tsv."
	$(TIME) sort $@.tmp >$@.tsv $(REDIR)
	@echo "Removing $@.tmp"
	$(TIME) rm $@.tmp $(REDIR)
	@echo "Extracting field 1 ($(K)-mer) from $@.tsv to $@"
	$(TIME) cut -f 1 -d " " $@.tsv >$@ $(REDIR)
	@rm $@.tsv
	@echo "Finished."

################################################################################
# kmerIsect: Get intersection of unique k-mers in all genomes.
################################################################################

# A list of all sorted k-mer files already defined above: SORTED_KMERS_FILES

# Phony target to make or clean PATH_ISECT_KMERS file(s).
ifeq ($(CLEAN),)
kmerIsect: $(PATH_ISECT_KMERS)
	@echo
	@echo "kmerIsect files are up to date."
else
kmerIsect:
	@$(CMD_DELETE_WHEN_CLEANING) $(PATH_ISECT_KMERS)
	@echo "kmerIsect output file(s) removed."
endif

# PATH_ISECT_KMERS target.

$(PATH_ISECT_KMERS) : $(SORTED_KMERS_FILES) | $(DIR_KMERS) $(PATH_PERL) $(PATH_KMER_ISECT)
	@echo
	@echo "*** kmerIsect PARAMS=$(PARAMS) $(CLEAN) ***"
	@echo "Intersecting unique $(K)-mers from .sorted files into $@"
	$(TIME) $(PATH_PERL) $(PATH_KMER_ISECT) $@ $^ $(REDIR)
	@echo "Finished."

################################################################################
# getGenomicPos: Get genomic position of common unique k-mers.  Argument: GENOME
################################################################################

# A list of all FASTA files already defined above: FASTA_FILES

# A list of all .isect files to be produced.
ISECT_KMER_FILES := $(foreach G,$(GENOME_NUMBERS),$(PFX_KMERS_DATA_FILE)$(G).isect)

# Target .isect file(s) for genome GENOME.
ifeq ($(GENOME),ALL)
TARGET_ISECT := $(ISECT_KMER_FILES)
else
TARGET_ISECT := $(PFX_KMERS_DATA_FILE)$(GENOME).isect
endif

# Phony target to make or clean TARGET_ISECT file(s).
ifeq ($(CLEAN),)
getGenomicPos: $(TARGET_ISECT)
	@echo
	@echo "getGenomicPos file(s) for genome(s) $(GENOME) are up to date."
else
getGenomicPos:
	@$(CMD_DELETE_WHEN_CLEANING) $(TARGET_ISECT)
	@echo "getGenomicPos output file(s) for genome(s) $(GENOME) removed."
endif

# ISECT_KMER_FILES is multiple targets, one per genome.
# Here, % is a genome number.

$(ISECT_KMER_FILES) : $(PFX_KMERS_DATA_FILE)%.isect : $(PATH_ISECT_KMERS) $$(PATH_GENOME_FASTA_$$*) | \
        $(DIR_KMERS) $(PATH_FINDMERS)
	@echo
	@echo "*** getGenomicPos PARAMS=$(PARAMS) $(CLEAN) GENOME=$* ***"
	@echo "Finding genomic positions of common unique $(K)-mers in $< into $@"
	@echo
	@echo "Adding genomic positions to common unique $(K)-mers from $< into $@.tmp"
	$(TIME) $(PATH_FINDMERS) $(ARGS_FINDMER) -v2 $(PATH_GENOME_FASTA_$*) $< $@.tmp $(REDIR)
	@echo "Sorting by $(K)-mer from $@.tmp into $@, removing header line"
	$(TIME) tail -n +2 $@.tmp | sort >$@ $(REDIR)
	@echo "Removing $@.tmp"
	$(TIME) rm $@.tmp $(REDIR)
	@echo "Finished."

################################################################################
# mergeKmers: Merge files of common unique k-mers having genomic positions into
# a single file containing the data for all genomes.  Argument: GENOME
################################################################################

# We want to join corresponding lines of multiple genome files, with the k-mer
# in the first column being the join key (the files are sorted by k-mer, so we
# are actually just joining the lines together, except that the k-mer appears
# only one time as the first column).  The files to be joined have suffix
# ".isect", and the joined files have suffix ".merge".
# This is problematic because the 'join' command can only join two files, but we
# have as many files to be joined as we have genomes.  We will join one genome
# at a time to the previously joined file, starting by joining the second genome
# k-mer positions file to the reference genome k-mer positions file.  When this
# is called with GENOME=1 (reference genome), we will just COPY the reference
# genome positions file (_1.isect) to a new file (_1.merge), which
# will be the join file for the second genome.  The second join joins the
# _1.merge file and the _2.isect file to form the _2.merge
# file, etc.

# A list of all .merge files to be produced.
MERGE_KMER_FILES := $(foreach G,$(GENOME_NUMBERS),$(PFX_KMERS_DATA_FILE)$(G).merge)

# Target .merge file(s) for genome GENOME.
ifeq ($(GENOME),ALL)
TARGET_MERGE := $(MERGE_KMER_FILES)
else
TARGET_MERGE := $(PFX_KMERS_DATA_FILE)$(GENOME).merge
endif

# Phony target to make or clean TARGET_MERGE file(s).
ifeq ($(CLEAN),)
mergeKmers: $(TARGET_MERGE)
	@echo
	@echo "mergeKmers files for genome(s) $(GENOME) are up to date."
else
mergeKmers:
	@$(CMD_DELETE_WHEN_CLEANING) $(TARGET_MERGE)
	@echo "mergeKmers output file(s) for genome(s) $(GENOME) removed."
endif

# Define the target file recipe for genome 1 separately, it is different.
TARGET_MERGE_1 := $(word 1,$(MERGE_KMER_FILES))

$(TARGET_MERGE_1) : $(PFX_KMERS_DATA_FILE)1.isect
	@echo
	@echo "*** mergeKmers PARAMS=$(PARAMS) $(CLEAN) GENOME=$* ***"
	@echo "Copy common unique $(K)-mers for genome 1 from $< to $@"
	$(TIME) cp $< $@ $(REDIR)
	@echo "Finished."

# Define the target files for all other genomes, they are done the same way.
TARGET_MERGE_OTHERS := $(filter-out $(TARGET_MERGE_1),$(MERGE_KMER_FILES))

# Define variable G_PREV to be equal to GENOME_NUMBERS with an extra "0" prepended
# to the beginning, so that $(word $(GENOME),$(G_PREV)) will give the number of the
# genome PRECEDING genome $(GENOME).
G_PREV := 0 $(GENOME_NUMBERS)

# MERGE_KMER_FILES is multiple targets, one per genome, but we EXCLUDE GENOME 1.
# Here, % is a genome number.

$(TARGET_MERGE_OTHERS) : $(PFX_KMERS_DATA_FILE)%.merge : $(PFX_KMERS_DATA_FILE)%.isect \
        $(PFX_KMERS_DATA_FILE)$$(word %,$$(G_PREV)).merge | $(PATH_RSCRIPT)
	@echo
	@echo "*** mergeKmers PARAMS=$(PARAMS) $(CLEAN) GENOME=$* ***"
	@echo "Merge common unique $(K)-mers for genomes $(word $*,$(G_PREV)) and $* to $@"
	$(TIME) join -t '	' $(PFX_KMERS_DATA_FILE)$(word $*,$(G_PREV)).merge $< >$@ $(REDIR)
	@echo "Finished."

################################################################################
# getCommonUniqueKmers: Sort merged kmers file to obtain sorted common unique
# k-mers file.
################################################################################

# The last file in the chain of merged files from the merge operation preceding this.
# This file is the dependent file we use as input to findLCRs (after sorting it).
UNSORTED_COMMON_UNIQUE_KMERS := $(PFX_KMERS_DATA_FILE)$(N_GENOMES).merge

# Phony target to make or clean PATH_COMMON_UNIQUE_KMERS file.
ifeq ($(CLEAN),)
getCommonUniqueKmers: $(PATH_COMMON_UNIQUE_KMERS)
	@echo
	@echo "Common unique k-mers file is up to date."
else
getCommonUniqueKmers:
	@$(CMD_DELETE_WHEN_CLEANING) $(PATH_COMMON_UNIQUE_KMERS)
	@echo "Common unique k-mers file removed."
endif

# Target for the sorted common unique k-mers file.
$(PATH_COMMON_UNIQUE_KMERS) : $(UNSORTED_COMMON_UNIQUE_KMERS)
	@echo "*** findLCRs PARAMS=$(PARAMS) $(CLEAN) ***"
	@echo
	@echo "Sort merged common unique $(K)-mers by reference genome position from $< into $(PATH_COMMON_UNIQUE_KMERS)"
	$(TIME) sort -k 2,2 -k 5,5n -k 3,3n $< >$(PATH_COMMON_UNIQUE_KMERS) $(REDIR)
	@echo "Finished."

################################################################################
# findLCRs: Find locally conserved regions using common unique k-mers with
# genomic positions that have been merged into a single file sorted by reference
# genome position.
################################################################################

# Phony target to make or clean PATH_LCR_FILE and PATH_BAD_KMERS_FILE files.
ifeq ($(CLEAN),)
findLCRs: $(PATH_LCR_FILE) $(PATH_BAD_KMERS_FILE)
	@echo
	@echo "findLCRs files are up to date."
else
findLCRs:
	@$(CMD_DELETE_WHEN_CLEANING) $(PATH_LCR_FILE) $(PATH_BAD_KMERS_FILE)
	@echo "findLCRs output file(s) removed."
endif

# PATH_LCR_FILE and PATH_BAD_KMERS_FILE targets are built at the same time.
# Define target names with % in them, so we can create a pattern target.  The %
# tells make that all the target files are made at one time by the recipe.  Here
# the only thing in common between the names is "$(SFX_LCR_FILE).tsv", and
# $(SFX_LCR_FILE) might be empty, but % at least matches .tsv.
PTN_LCR_FILE := $(DIR_MAIN_OUTPUT)/$(PFX_LCR_FILE)%
PTN_BAD_KMERS_FILE := $(DIR_MAIN_OUTPUT)/$(PFX_BAD_KMERS_FILE)%

# Use the patterns in a pattern target.  Dependent file is the SORTED version of
# the common unique k-mers file.
$(PTN_LCR_FILE) $(PTN_BAD_KMERS_FILE) : $(PATH_COMMON_UNIQUE_KMERS) | \
        $(DIR_IGGPIPE_OUT) $(PATH_RSCRIPT) $(PATH_FIND_LCRS)
	@echo "*** findLCRs PARAMS=$(PARAMS) $(CLEAN) ***"
	@echo
	@echo "Find locally conserved regions using common unique $(K)-mers from $< into $@"
	$(TIME) $(PATH_RSCRIPT) $(PATH_FIND_LCRS) $(WD) $(PATH_COMMON_UNIQUE_KMERS) \
	    $(GENOME_LETTERS_SQUISHED) $(KMIN) $(LMIN) $(DMIN) $(DMAX) \
	    $(PATH_LCR_FILE) \
	    $(PATH_BAD_KMERS_FILE) \
	    $(INVESTIGATE_FINDLCRS) $(REDIR)
	@echo "Finished."

################################################################################
# findIndelGroups: Analyze locally conserved regions to find insertions/deletions.
################################################################################

# A list of all .idlens files is already defined above: IDLEN_FILES

# Phony target to make or clean PATH_OVERLAPPING_INDEL_GROUPS_FILE and PATH_NONOVERLAPPING_INDEL_GROUPS_FILE files.
ifeq ($(CLEAN),)
findIndelGroups: $(PATH_OVERLAPPING_INDEL_GROUPS_FILE) $(PATH_NONOVERLAPPING_INDEL_GROUPS_FILE)
	@echo
	@echo "findIndelGroups files are up to date."
else
findIndelGroups:
	@$(CMD_DELETE_WHEN_CLEANING) $(PATH_OVERLAPPING_INDEL_GROUPS_FILE) $(PATH_NONOVERLAPPING_INDEL_GROUPS_FILE)
	@echo "findIndelGroups output file(s) removed."
endif

# PATH_OVERLAPPING_INDEL_GROUPS_FILE and PATH_NONOVERLAPPING_INDEL_GROUPS_FILE targets are built at the same time.
# Define target names with % in them, so we can create a pattern target.  The %
# tells make that all the target files are made at one time by the recipe.  Here
# the only thing in common between the names is "$(SFX_INDEL_GROUPS_FILE).tsv", and
# $(SFX_INDEL_GROUPS_FILE) might be empty, but % at least matches .tsv.
PTN_OVERLAPPING_INDELS_FILE := $(DIR_MAIN_OUTPUT)/$(PFX_OVERLAPPING_INDEL_GROUPS_FILE)%
PTN_NONOVERLAPPING_INDELS_FILE := $(DIR_MAIN_OUTPUT)/$(PFX_NONOVERLAPPING_INDEL_GROUPS_FILE)%

# Use the patterns in a pattern target.
$(PTN_OVERLAPPING_INDELS_FILE) $(PTN_NONOVERLAPPING_INDELS_FILE) : $(PATH_LCR_FILE) $(IDLEN_FILES) | \
        $(DIR_IGGPIPE_OUT) $(PATH_RSCRIPT) $(PATH_FIND_INDEL_GROUPS)
	@echo
	@echo "*** findIndelGroups PARAMS=$(PARAMS) $(CLEAN) ***"
	@echo "Find Indel Groups using locally conserved regions in $< and write them to two output files"
	$(TIME) $(PATH_RSCRIPT) $(PATH_FIND_INDEL_GROUPS) $(WD) \
	    $(PATH_LCR_FILE) \
	    $(AMIN) $(AMAX) $(ADMIN) $(ADMAX) $(NDAMIN) $(MINFLANK) $(OVERLAP_REMOVAL) \
	    $(PATH_OVERLAPPING_INDEL_GROUPS_FILE) $(PATH_NONOVERLAPPING_INDEL_GROUPS_FILE) \
	    $(GENOME_LETTERS_SQUISHED) $(INVESTIGATE_ANALYZELCRS) \
	    $(IDLEN_FILES) $(REDIR)
	@echo "Finished."

################################################################################
# getDNAseqs: Extract DNA sequence around Indel Groups, for making primers.
################################################################################

# A list of all FASTA files already defined above: FASTA_FILES

# Likewise, a list of all .contigs files is already defined above: CONTIG_FILES

# A list of all .dnaseqs files to be produced.
DNA_SEQ_FILES := $(foreach G,$(GENOME_NUMBERS),$(PFX_DNA_SEQS_PATH)_$(G).dnaseqs)

# Target .dnaseqs file(s) for genome GENOME.
ifeq ($(GENOME),ALL)
TARGET_DNASEQ := $(DNA_SEQ_FILES)
else
TARGET_DNASEQ := $(PFX_DNA_SEQS_PATH)_$(GENOME).dnaseqs
endif

# Phony target to make or clean TARGET_DNASEQ file(s).
ifeq ($(CLEAN),)
getDNAseqs: $(TARGET_DNASEQ)
	@echo
	@echo "getDNAseqs file(s) for genome(s) $(GENOME) are up to date."
else
getDNAseqs:
	@$(CMD_DELETE_WHEN_CLEANING) $(TARGET_DNASEQ) $(DIR_GENOME_OUT_DATA)/extract*.txt $(DIR_GENOME_OUT_DATA)/seqs*.txt
	@echo "getDNAseqs output file(s) for genome(s) $(GENOME) removed."
endif

# DNA_SEQ_FILES is multiple targets, one per genome.
# Here, % is a genome number.

$(DNA_SEQ_FILES) : $(PFX_DNA_SEQS_PATH)_%.dnaseqs : $$(PATH_GENOME_FASTA_$$*) \
        $(PFX_GENOME_DATA_FILE)$$*.contigs $(PATH_OVERLAPPING_INDEL_GROUPS_FILE) | \
        $(DIR_IGGPIPE_OUT) $(DIR_GENOME_OUT_DATA) \
        $(PATH_RSCRIPT) $(PATH_GET_DNA_SEQS) $(PATH_PERL) $(PATH_GET_SEQS_FASTA)
	@echo
	@echo "*** getDNAseqs PARAMS=$(PARAMS) $(CLEAN) GENOME=$* ***"
	@echo "Extract DNA sequence around Indel Groups and write to $@"
	$(TIME) $(PATH_RSCRIPT) $(PATH_GET_DNA_SEQS) $(WD) \
	    $(PATH_OVERLAPPING_INDEL_GROUPS_FILE) $* \
	    $@ \
		$(DIR_GENOME_OUT_DATA) $(EXTENSION_LEN) \
		$(PATH_PERL) $(PATH_GET_SEQS_FASTA) \
		$(PATH_GENOME_FASTA_$*) \
		$(PFX_GENOME_DATA_FILE)$*.contigs $(INVESTIGATE_GETDNASEQS) $(REDIR)
	@echo "Finished."

################################################################################
# findPrimers: Run primer3 to search for primers around Indel Groups.
################################################################################

# A list of all .dnaseqs files is already defined above: DNA_SEQ_FILES

# Phony target to make or clean PATH_NONVALIDATED_MARKER_FILE file.
ifeq ($(CLEAN),)
findPrimers: $(PATH_NONVALIDATED_MARKER_FILE)
	@echo
	@echo "findPrimers files are up to date."
else
findPrimers:
	@$(CMD_DELETE_WHEN_CLEANING) $(PATH_NONVALIDATED_MARKER_FILE) $(PATH_PRIMER3_IN) $(PATH_PRIMER3_OUT)
	@echo "findPrimers output file(s) removed."
endif

# PATH_NONVALIDATED_MARKER_FILE target.

$(PATH_NONVALIDATED_MARKER_FILE) : $(PATH_OVERLAPPING_INDEL_GROUPS_FILE) $(DNA_SEQ_FILES) | $(PATH_RSCRIPT) \
        $(DIR_PRIMER_DATA) $(PATH_FIND_PRIMERS) $(PATH_PRIMER3CORE) $(PATH_PRIMER3_SETTINGS)
	@echo
	@echo "*** findPrimers PARAMS=$(PARAMS) $(CLEAN) ***"
	@echo "Find primers around Indel Groups in $< and write to $@"
	$(TIME) $(PATH_RSCRIPT) $(PATH_FIND_PRIMERS) $(WD) \
		$(PATH_OVERLAPPING_INDEL_GROUPS_FILE) \
		$(PFX_DNA_SEQS_PATH) \
		$(PATH_NONVALIDATED_MARKER_FILE) \
		$(PATH_PRIMER3CORE) $(PATH_PRIMER3_SETTINGS) $(PATH_PRIMER3CONFIG) \
		$(PATH_PRIMER3_IN) $(PATH_PRIMER3_OUT) \
		$(INVESTIGATE_FINDPRIMERS) $(REDIR)
	@echo "Finished."

################################################################################
# ePCRtesting: Run e-PCR to test all primer pairs for proper amplicon size.
################################################################################

# A list of all FASTA files already defined above: FASTA_FILES

# A list of all .bad.tsv files to be produced.
BAD_MARKER_FILES := $(foreach G,$(GENOME_NUMBERS),$(PFX_BAD_MARKER_ERROR_PATH)_$(G).bad.tsv)

# Target .bad.tsv file(s) for genome GENOME.
ifeq ($(GENOME),ALL)
TARGET_BAD_MARKER := $(BAD_MARKER_FILES)
else
TARGET_BAD_MARKER := $(PFX_BAD_MARKER_ERROR_PATH)_$(word $(GENOME),$(GENOME_NUMBERS)).bad.tsv
endif

# Phony target to make or clean TARGET_BAD_MARKER file(s).
ifeq ($(CLEAN),)
ePCRtesting: $(TARGET_BAD_MARKER)
	@echo
	@echo "ePCRtesting file(s) for genome(s) $(GENOME) are up to date."
else
ePCRtesting:
	@$(CMD_DELETE_WHEN_CLEANING) $(TARGET_BAD_MARKER) $(DIR_GENOME_OUT_DATA)/*.epcr.in $(DIR_GENOME_OUT_DATA)/*.epcr.out
	@echo "ePCRtesting output file(s) for genome(s) $(GENOME) removed."
endif

# BAD_MARKER_FILES is multiple targets, one per genome.
# Here, % is a genome number.

$(BAD_MARKER_FILES) : $(PFX_BAD_MARKER_ERROR_PATH)_%.bad.tsv : $$(PATH_GENOME_FASTA_$$*) \
        $(PATH_NONVALIDATED_MARKER_FILE) | $(DIR_IGGPIPE_OUT) $(DIR_GENOME_OUT_DATA) $(DIR_PRIMER_DATA) \
        $(PATH_RSCRIPT) $(PATH_EPCR_TESTING) $(PATH_EPCR)
	@echo
	@echo "*** ePCRtesting PARAMS=$(PARAMS) $(CLEAN) GENOME=$* ***"
	@echo "Use e-PCR to test marker primer pairs and write errors to $@"
	$(TIME) $(PATH_RSCRIPT) $(PATH_EPCR_TESTING) $(WD) \
	    $(PATH_NONVALIDATED_MARKER_FILE) $* $@ \
	    $(DIR_GENOME_OUT_DATA) $(PATH_EPCR) \
	    $(EPCR_MAX_DEV) $(EPCR_WORD_SIZE) $(EPCR_MAX_MISMATCH) $(EPCR_MAX_GAPS) \
		$(PATH_GENOME_FASTA_$*) $(INVESTIGATE_EPCRTESTING) $(REDIR)
	@echo "Finished."

################################################################################
# removeBadMarkers: Read bad marker files that failed e-PCR, then read full
# marker file, remove the bad markers from it, and write new files of good ones,
# one set overlapping, the other not.
################################################################################

# A list of all .bad.tsv files is already defined above: BAD_MARKER_FILES

# Phony target to make or clean PATH_OVERLAPPING_MARKERS_FILE and PATH_NONOVERLAPPING_MARKERS_FILE.
ifeq ($(CLEAN),)
removeBadMarkers: $(PATH_OVERLAPPING_MARKERS_FILE) $(PATH_NONOVERLAPPING_MARKERS_FILE)
	@echo
	@echo "removeBadMarkers files are up to date."
else
removeBadMarkers:
	@$(CMD_DELETE_WHEN_CLEANING) $(PATH_OVERLAPPING_MARKERS_FILE) $(PATH_NONOVERLAPPING_MARKERS_FILE)
	@echo "removeBadMarkers output file(s) removed."
endif

# PATH_OVERLAPPING_MARKERS_FILE and PATH_NONOVERLAPPING_MARKERS_FILE targets are built at the same time.
# Define target names with % in them, so we can create a pattern target.  The %
# tells make that all the target files are made at one time by the recipe.  Here
# the only thing in common between the names is "$(SFX_GOOD_MARKER_FILE).tsv", and
# $(SFX_GOOD_MARKER_FILE) might be empty, but % at least matches .tsv.
PTN_OVERLAPPING_MARKERS_FILE := $(DIR_MAIN_OUTPUT)/$(PFX_OVERLAPPING_MARKERS_FILE)%
PTN_NONOVERLAPPING_MARKERS_FILE := $(DIR_MAIN_OUTPUT)/$(PFX_NONOVERLAPPING_MARKERS_FILE)%

# Use the patterns in a pattern target.
$(PTN_OVERLAPPING_MARKERS_FILE) $(PTN_NONOVERLAPPING_MARKERS_FILE) : $(PATH_NONVALIDATED_MARKER_FILE) $(BAD_MARKER_FILES) | \
        $(PATH_RSCRIPT) $(PATH_RMV_BAD_MARKERS)
	@echo
	@echo "*** removeBadMarkers PARAMS=$(PARAMS) $(CLEAN) ***"
	@echo "Remove markers identified by e-PCR as bad from $< and write good ones to two output files"
	$(TIME) $(PATH_RSCRIPT) $(PATH_RMV_BAD_MARKERS) $(WD) $(OVERLAP_REMOVAL) \
	    $(PATH_NONVALIDATED_MARKER_FILE) \
	    $(PFX_BAD_MARKER_ERROR_PATH) \
		$(PATH_OVERLAPPING_MARKERS_FILE) \
		$(PATH_NONOVERLAPPING_MARKERS_FILE) \
		$(INVESTIGATE_RMVBADMARKERS) $(REDIR)
	@echo "Finished."

################################################################################
# plotMarkers: Make density plots of final "good" candidate IGG markers, both
# overlapping and non-overlapping.
################################################################################

# A list of all .idlens files is already defined above: IDLEN_FILES

# The full name of the counts plot output file.
MARKER_COUNTS_FILE := $(PFX_MARKER_COUNTS_PATH).plot.pdf

# A list of all .png files to be produced.
MARKER_DENSITY_FILES := $(foreach G,$(GENOME_LETTERS),$(PFX_MARKER_DENSITY_PATH)_$(G).plot.png)

# Phony target to make or clean MARKER_COUNTS_FILE and MARKER_DENSITY_FILES files.
ifeq ($(CLEAN),)
plotMarkers: $(MARKER_COUNTS_FILE) $(MARKER_DENSITY_FILES)
	@echo
	@echo "plotMarkers files are up to date."
else
plotMarkers:
	@$(CMD_DELETE_WHEN_CLEANING) $(MARKER_COUNTS_FILE) $(MARKER_DENSITY_FILES)
	@echo "plotMarkers output file(s) removed."
endif

# MARKER_COUNTS_FILE and MARKER_DENSITY_FILES targets are built at the same time.
# Define target names with % in them, so we can create a pattern target.  The %
# tells make that all the target files are made at one time by the recipe.  Here
# the only thing in common between the names is ".plot.", which is what % matches.
PTN_COUNTS_FILE := $(PFX_MARKER_COUNTS_PATH)%pdf
PTN_DENSITY_FILES := $(foreach G,$(GENOME_LETTERS),$(PFX_MARKER_DENSITY_PATH)_$(G)%png)

# Use the patterns in a pattern target.
$(PTN_COUNTS_FILE) $(PTN_DENSITY_FILES) : $(PATH_OVERLAPPING_MARKERS_FILE) $(PATH_NONOVERLAPPING_MARKERS_FILE) \
        $(IDLEN_FILES) | $(PATH_RSCRIPT) $(PATH_PLOT_MARKERS)
	@echo
	@echo "*** plotMarkers PARAMS=$(PARAMS) $(CLEAN) ***"
	@echo "Make density plots of 'good' candidate IGG markers to file $@"
	$(TIME) $(PATH_RSCRIPT) $(PATH_PLOT_MARKERS) $(WD) $(PLOT_NDAMIN) $(PLOT_ALPHA) \
	    $(PATH_OVERLAPPING_MARKERS_FILE) \
	    $(PATH_NONOVERLAPPING_MARKERS_FILE) \
		$(PFX_MARKER_COUNTS_PATH) \
		$(PFX_MARKER_DENSITY_PATH) \
		$(IDLEN_FILES) $(REDIR)
	@echo "Finished."

################################################################################
# InDels: Read input file and perform alignments, then search them for InDels.
################################################################################

# A list of all FASTA files already defined above: FASTA_FILES

# Do this to avoid case annoyance.
Indels: InDels

# Do this to avoid case annoyance.
indels: InDels

# Phony target to make or clean PATH_INDELS_OUTPUT_FILE file.
# If variable PARAMS is not defined, show basic usage info, else make InDels
# output file target.
ifeq ($(PARAMS),)
InDels:
	@echo
	@echo "You must specify a PARAMS file:"
	@echo
	@echo "$(INDENT)make PARAMS=<allParametersFile> InDels"
	@echo 
else ifeq ($(CLEAN),)
InDels: $(PATH_INDELS_OUTPUT_FILE)
	@echo
	@echo "InDels files are up to date."
else
InDels:
	@$(CMD_DELETE_WHEN_CLEANING) $(PATH_INDELS_OUTPUT_FILE)
	@echo "InDels output file(s) removed."
endif

# PATH_INDELS_OUTPUT_FILE target.

$(PATH_INDELS_OUTPUT_FILE) : $(PATH_INDELS_INPUT_FILE) $(FASTA_FILES) | $(DIR_GENOME_OUT_DATA) \
        $(PATH_RSCRIPT) $(PATH_ALIGN_AND_GET_INDELS) $(PATH_PERL) $(PATH_GET_SEQS_FASTA) $(PATH_ALIGNER)
	@echo
	@echo "*** InDels PARAMS=$(PARAMS) $(CLEAN) ***"
	@echo "Align sequences of $< and find InDels and write them to $@"
	$(TIME) $(PATH_RSCRIPT) $(PATH_ALIGN_AND_GET_INDELS) $(WD) \
	    $(PATH_INDELS_INPUT_FILE) $(PATH_INDELS_OUTPUT_FILE) \
	    $(DIR_GENOME_OUT_DATA) \
		$(PATH_PERL) $(PATH_GET_SEQS_FASTA) \
		$(PATH_ALIGNER) $(INVESTIGATE_ALIGN_AND_GET_INDELS) \
		$(FASTA_FILES) $(REDIR)
	@echo "Finished."

################################################################################
# plotInDels: Make plots of information about InDels found within InDel groups.
################################################################################

# Do this to avoid case annoyance.
plotIndels: plotInDels

# Phony target to make or clean PATH_INDELS_PLOT_FILE file.
# If variable PARAMS is not defined, show basic usage info, else make plotInDels
# output file target.
ifeq ($(PARAMS),)
plotInDels:
	@echo
	@echo "You must specify a PARAMS file:"
	@echo
	@echo "$(INDENT)make PARAMS=<allParametersFile> plotInDels"
	@echo 
else ifeq ($(CLEAN),)
plotInDels: $(PATH_INDELS_PLOT_FILE)
	@echo
	@echo "plotInDels file is up to date."
else
plotInDels:
	@$(CMD_DELETE_WHEN_CLEANING) $(PATH_INDELS_PLOT_FILE)
	@echo "plotInDels output file(s) removed."
endif

# PATH_INDELS_PLOT_FILE target.

$(PATH_INDELS_PLOT_FILE) : $(PATH_INDELS_OUTPUT_FILE) | $(PATH_RSCRIPT) $(PATH_PLOT_INDELS)
	@echo
	@echo "*** plotInDels PARAMS=$(PARAMS) $(CLEAN) ***"
	@echo "Make plots of InDel information to file $@"
	$(TIME) $(PATH_RSCRIPT) $(PATH_PLOT_INDELS) $(WD) $(PATH_INDELS_OUTPUT_FILE) \
	    $(PATH_INDELS_PLOT_FILE) $(REDIR)
	@echo "Finished."

################################################################################
# End of file.
################################################################################
