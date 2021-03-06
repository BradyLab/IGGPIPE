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
# a static pattern rule).  Thus, $$(VAR_NAME_PREFIX$$*) would allow a nested
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
# Include VERSION.txt.
################################################################################
include VERSION.txt

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
# Target version: show version number from VERSION.txt file.
################################################################################
version:
	@echo "IGGPIPE Version $(VERSION)"

################################################################################
# target basic_usage: show basic usage info.
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
# Target usage: show usage info.
################################################################################
usage:
	less help.txt
help:
	less help.txt

################################################################################
# Target clean: give user some instructions on how to clean away files.
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
$(error GENOME must be either a valid genome number (1..$(N_GENOMES)) or 'ALL')
endif
endif

################################################################################
# Set variables for cleaning operations.
################################################################################
# If variable CLEAN was defined as anything at all, set it to "1".
ifneq ($(CLEAN),)
override CLEAN := 1
endif

################################################################################
# Set variables for command timing.
################################################################################
# If variable TIME_CMDS is YES, set TIME and REDIR to time commands and redirect
# stderr to stdout.  If not YES, leave these empty.
ifeq ($(TIME_CMDS),YES)
TIME := ($(CMD_TIME)
REDIR := ) 2>&1
endif

################################################################################
# ALL target.
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
	@echo
	@echo "Removed all files from the output directory."
	@echo
else ifeq ($(CLEAN),1)
ALL: getSeqInfo getKmers kmersToText kmerStats getContigFile getGenomicPosIsect mergeKmers \
	sortCommonUniqueKmers findLCRs findIndelGroups getDNAseqsForPrimers findPrimers \
	ePCRtesting removeBadMarkers plotMarkers
	@echo
	@echo "ALL files are cleaned"
	@echo
else
ALL: getSeqInfo getKmers kmersToText kmerStats getContigFile getGenomicPosIsect mergeKmers \
	sortCommonUniqueKmers findLCRs findIndelGroups getDNAseqsForPrimers findPrimers \
	ePCRtesting removeBadMarkers plotMarkers
	@echo
	@echo "ALL files are up to date"
	@echo
endif

################################################################################
# Targets to create output directories.
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
# Target getSeqInfo: Get sequence ID and length information.  Argument: GENOME
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
.PHONY: getSeqInfo
ifeq ($(CLEAN),)
getSeqInfo: $(TARGET_IDLEN)
	@echo
	@echo "getSeqInfo file(s) for genome(s) $(GENOME) are up to date."
	@echo
else
getSeqInfo:
	@$(CMD_DELETE_WHEN_CLEANING) $(TARGET_IDLEN)
	@echo "getSeqInfo output file(s) for genome(s) $(GENOME) removed."
	@echo
endif

# IDLEN_FILES is multiple targets, one per genome.
# Here, % is a genome number.

$(IDLEN_FILES) : $(PFX_GENOME_DATA_FILE)%.idlens : $$(PATH_GENOME_FASTA_$$*) | \
        $(DIR_GENOME_OUT_DATA) $(PATH_EXTRACT_SEQ_IDS)
	@echo
	@echo "*** getSeqInfo PARAMS=$(PARAMS) GENOME=$* ***"
	@echo "Extracting sequence IDs and lengths from $< into $@"
	$(TIME) \
		$(CMD_PERL) $(PATH_EXTRACT_SEQ_IDS) $< $@ $(REDIR) ### getSeqInfo_$*
	@echo "Finished."

################################################################################
# Target getContigFile: Get sequence contig positions and lengths.
# Argument: GENOME
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
.PHONY: getContigFile
ifeq ($(CLEAN),)
getContigFile: $(TARGET_CONTIG)
	@echo
	@echo "getContigFile file(s) for genome(s) $(GENOME) are up to date."
	@echo
else
getContigFile:
	@$(CMD_DELETE_WHEN_CLEANING) $(TARGET_CONTIG)
	@echo "getContigFile output file(s) for genome(s) $(GENOME) removed."
	@echo
endif

# CONTIG_FILES is multiple targets, one per genome.
# Here, % is a genome number.

$(CONTIG_FILES) : $(PFX_GENOME_DATA_FILE)%.contigs : $$(PATH_GENOME_FASTA_$$*) | \
        $(DIR_GENOME_OUT_DATA) $(PATH_FINDMERS)
	@echo
	@echo "*** getContigFile PARAMS=$(PARAMS) GENOME=$* ***"
	@echo "Extracting contig positions and lengths from $< into $@"
	$(TIME) \
		$(PATH_FINDMERS) -f $@ $< $(REDIR) ### getContigFile_$*
	@echo "Finished."

################################################################################
# Target getKmers: Get unique k-mers.  Argument: GENOME
################################################################################

# A list of all FASTA files already defined above: FASTA_FILES

# A list of all binary .kmers files to be produced.
KMERS_BINARY_FILES := $(foreach G,$(GENOME_NUMBERS),$(PFX_KMERS_DATA_FILE)$(G).kmers)

# Target binary .kmers file(s) for genome GENOME.
ifeq ($(GENOME),ALL)
TARGET_BINARY_KMERS := $(KMERS_BINARY_FILES)
else
TARGET_BINARY_KMERS := $(PFX_KMERS_DATA_FILE)$(GENOME).kmers
endif

# Target .kmers.txt file(s) for genome GENOME.
ifeq ($(GENOME),ALL)
TARGET_TEXT_KMERS := $(KMERS_TEXT_FILES)
else
TARGET_TEXT_KMERS := $(PFX_KMERS_DATA_FILE)$(GENOME).kmers.txt
endif

# Phony target to make or clean TARGET_BINARY_KMERS file(s).
.PHONY: getKmers
ifeq ($(CLEAN),)
getKmers: $(TARGET_TEXT_KMERS)
	@echo
	@echo "getKmers file(s) for genome(s) $(GENOME) are up to date."
	@echo
else
getKmers:
	@$(CMD_DELETE_WHEN_CLEANING) $(TARGET_BINARY_KMERS)
	@echo "getKmers output file(s) for genome(s) $(GENOME) removed."
	@echo
endif

# KMERS_BINARY_FILES is multiple targets, one per genome.
# Here, % is a genome number.

$(KMERS_BINARY_FILES) : $(PFX_KMERS_DATA_FILE)%.kmers : $$(PATH_GENOME_FASTA_$$*) | $(DIR_KMERS)
	@echo
	@echo "*** getKmers PARAMS=$(PARAMS) GENOME=$* ***"
	@echo "Extracting unique $(K)-mers from $< into binary file $@"
	$(TIME) \
		$(CMD_JELLYFISH) count -C -m $(K) -s $(JELLYFISH_HASH_SIZE) -U 1 -o $@ $< $(REDIR) ### getKmers_$*
	@echo "Finished."

################################################################################
# Target kmerStats: Get unique k-mer statistics.  Argument: GENOME
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
.PHONY: kmerStats
ifeq ($(CLEAN),)
kmerStats: $(TARGET_STATS)
	@echo
	@echo "kmerStats file(s) for genome(s) $(GENOME) are up to date."
	@echo
else
kmerStats:
	@$(CMD_DELETE_WHEN_CLEANING) $(TARGET_STATS)
	@echo "kmerStats output file(s) for genome(s) $(GENOME) removed."
	@echo
endif

# KMERS_STATS_FILES is multiple targets, one per genome.
# Here, % is a genome number.

$(KMERS_STATS_FILES) : $(PFX_KMERS_DATA_FILE)%.stats : $(PFX_KMERS_DATA_FILE)%.kmers | \
        $(DIR_KMERS)
	@echo
	@echo "*** kmerStats PARAMS=$(PARAMS) GENOME=$* ***"
	@echo "Getting statistics for $(K)-mers from $< into $@"
	$(TIME) \
		$(CMD_JELLYFISH) stats -v $< >$@ $(REDIR) ### jellyfish_$*
	@echo "Finished."
	@echo
	@echo "Statistics are:"
	@cat $@
	@echo

################################################################################
# Target kmersToText: Convert binary k-mers to text k-mers.  Argument: GENOME
################################################################################

# A list of all binary .kmers files already defined above: KMERS_BINARY_FILES

# A list of all .kmers.txt files to be produced.
KMERS_TEXT_FILES := $(foreach G,$(GENOME_NUMBERS),$(PFX_KMERS_DATA_FILE)$(G).kmers.txt)

# Target .kmers.txt file(s) for genome GENOME.
ifeq ($(GENOME),ALL)
TARGET_TEXT_KMERS := $(KMERS_TEXT_FILES)
else
TARGET_TEXT_KMERS := $(PFX_KMERS_DATA_FILE)$(GENOME).kmers.txt
endif

# Phony target to make or clean TARGET_TEXT_KMERS file(s).
.PHONY: kmersToText
ifeq ($(CLEAN),)
kmersToText: $(TARGET_TEXT_KMERS)
	@echo
	@echo "kmersToText file(s) for genome(s) $(GENOME) are up to date."
	@echo
else
kmersToText:
	@$(CMD_DELETE_WHEN_CLEANING) $(TARGET_TEXT_KMERS)
	@echo "kmersToText output file(s) for genome(s) $(GENOME) removed."
	@echo
endif

# KMERS_TEXT_FILES is multiple targets, one per genome.
# Here, % is a genome number.

$(KMERS_TEXT_FILES) : $(PFX_KMERS_DATA_FILE)%.kmers.txt : $(PFX_KMERS_DATA_FILE)$$*.kmers | $(DIR_KMERS)
	@echo
	@echo "*** kmersToText PARAMS=$(PARAMS) GENOME=$* ***"
	@echo "Converting $(K)-mers from binary file $< into text format in file $@"
	$(TIME) \
		$(CMD_JELLYFISH) dump -c -o $@ $< $(REDIR) ### jellydump_$*
	@echo "Finished."

################################################################################
# Target getGenomicPosIsect: Intersect unique k-mers to get common unique
# k-mers, and add genomic positions to them.  Argument: GENOME
################################################################################

# A list of all FASTA files already defined above: FASTA_FILES

# A list of all .kmers.txt files preceded by "-i "
KMERS_TEXT_FILES_i := $(foreach G,$(GENOME_NUMBERS),-i $(PFX_KMERS_DATA_FILE)$(G).kmers.txt)

# A list of all .isect files to be produced.
ISECT_KMER_FILES := $(foreach G,$(GENOME_NUMBERS),$(PFX_KMERS_DATA_FILE)$(G).isect)

# Target .isect file(s) for genome GENOME.
ifeq ($(GENOME),ALL)
TARGET_ISECT := $(ISECT_KMER_FILES)
else
TARGET_ISECT := $(PFX_KMERS_DATA_FILE)$(GENOME).isect
endif

# Phony target to make or clean TARGET_ISECT file(s).
.PHONY: getGenomicPosIsect
ifeq ($(CLEAN),)
getGenomicPosIsect: $(TARGET_ISECT)
	@echo
	@echo "getGenomicPosIsect file(s) for genome(s) $(GENOME) are up to date."
	@echo
else
getGenomicPosIsect:
	@$(CMD_DELETE_WHEN_CLEANING) $(TARGET_ISECT)
	@echo "getGenomicPosIsect output file(s) for genome(s) $(GENOME) removed."
	@echo
endif

# ISECT_KMER_FILES is multiple targets, one per genome.
# Here, % is a genome number.

$(ISECT_KMER_FILES) : $(PFX_KMERS_DATA_FILE)%.isect : $(KMERS_TEXT_FILES) $$(PATH_GENOME_FASTA_$$*) | \
        $(DIR_KMERS) $(PATH_FINDMERS)
	@echo
	@echo "*** getGenomicPosIsect PARAMS=$(PARAMS) GENOME=$* ***"
	@echo "Intersect unique k-mers in $< and get genomic positions of common unique $(K)-mers into $@"
	@echo
	@echo "Adding genomic positions to common unique $(K)-mers from $< into $@.pos"
	$(TIME) \
		$(PATH_FINDMERS) -v2 $(KMERS_TEXT_FILES_i) $(PATH_GENOME_FASTA_$*) $(K) $@.pos $(REDIR) ### getGenomicPosIsect_$*
	@echo "Sorting by $(K)-mer from $@.pos into $@, removing header line"
	$(TIME) \
		tail -n +2 $@.pos | sort >$@ $(REDIR) ### sort_kmer_$*
	@echo "Removing $@.pos"
	$(TIME) \
		rm $@.pos $(REDIR) ### rm_$*
	@echo "Finished."

################################################################################
# Target mergeKmers: Merge files of common unique k-mers having genomic
# positions into a single file containing the data for all genomes.
# Argument: GENOME
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
.PHONY: mergeKmers
ifeq ($(CLEAN),)
mergeKmers: $(TARGET_MERGE)
	@echo
	@echo "mergeKmers files for genome(s) $(GENOME) are up to date."
	@echo
else
mergeKmers:
	@$(CMD_DELETE_WHEN_CLEANING) $(TARGET_MERGE)
	@echo "mergeKmers output file(s) for genome(s) $(GENOME) removed."
	@echo
endif

# Define the target file recipe for genome 1 separately, it is different.
TARGET_MERGE_1 := $(word 1,$(MERGE_KMER_FILES))

$(TARGET_MERGE_1) : $(PFX_KMERS_DATA_FILE)1.isect
	@echo
	@echo "*** mergeKmers PARAMS=$(PARAMS) ***"
	@echo "Copy common unique $(K)-mers for genome 1 from $< to $@"
	$(TIME) \
		cp $< $@ $(REDIR) ### cp
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
        $(PFX_KMERS_DATA_FILE)$$(word %,$$(G_PREV)).merge
	@echo
	@echo "*** mergeKmers PARAMS=$(PARAMS) GENOME=$* ***"
	@echo "Merge common unique $(K)-mers for genomes $(word $*,$(G_PREV)) and $* to $@"
	$(TIME) \
		join -t '	' $(PFX_KMERS_DATA_FILE)$(word $*,$(G_PREV)).merge $< >$@ $(REDIR) ### join_$*
	@echo "Finished."

################################################################################
# Target sortCommonUniqueKmers: Sort merged kmers file to obtain sorted common
# unique k-mers file.
################################################################################

# The last file in the chain of merged files from the merge operation preceding this.
# This file is the dependent file we use as input to findLCRs (after sorting it).
UNSORTED_COMMON_UNIQUE_KMERS := $(PFX_KMERS_DATA_FILE)$(N_GENOMES).merge

# Phony target to make or clean PATH_COMMON_UNIQUE_KMERS file.
.PHONY: sortCommonUniqueKmers
ifeq ($(CLEAN),)
sortCommonUniqueKmers: $(PATH_COMMON_UNIQUE_KMERS)
	@echo
	@echo "Common unique k-mers file is up to date."
	@echo
else
sortCommonUniqueKmers:
	@$(CMD_DELETE_WHEN_CLEANING) $(PATH_COMMON_UNIQUE_KMERS)
	@echo "Common unique k-mers file removed."
	@echo
endif

# Target for the sorted common unique k-mers file.
$(PATH_COMMON_UNIQUE_KMERS) : $(UNSORTED_COMMON_UNIQUE_KMERS)
	@echo "*** sortCommonUniqueKmers PARAMS=$(PARAMS) ***"
	@echo
	@echo "Sort merged common unique $(K)-mers by reference genome position from $< into $(PATH_COMMON_UNIQUE_KMERS)"
	$(TIME) \
		sort -k 2,2 -k 5,5n -k 3,3n $< >$(PATH_COMMON_UNIQUE_KMERS) $(REDIR) ### sort_merged
	@echo "Finished."

################################################################################
# Target findLCRs: Find locally conserved regions using common unique k-mers
# with genomic positions that have been merged into a single file sorted by
# reference genome position.
################################################################################

# Phony target to make or clean PATH_LCR_FILE and PATH_BAD_KMERS_FILE files.
.PHONY: findLCRs
ifeq ($(CLEAN),)
findLCRs: $(PATH_LCR_FILE) $(PATH_BAD_KMERS_FILE)
	@echo
	@echo "findLCRs files are up to date."
	@echo
else
findLCRs:
	@$(CMD_DELETE_WHEN_CLEANING) $(PATH_LCR_FILE) $(PATH_BAD_KMERS_FILE)
	@echo "findLCRs output file(s) removed."
	@echo
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
        $(DIR_IGGPIPE_OUT) $(PATH_FIND_LCRS)
	@echo "*** findLCRs PARAMS=$(PARAMS) ***"
	@echo
	@echo "Find locally conserved regions using common unique $(K)-mers from $< into $@"
	$(TIME) \
		$(CMD_RSCRIPT) $(PATH_FIND_LCRS) $(WD) $(PATH_COMMON_UNIQUE_KMERS) \
	    $(GENOME_LETTERS_SQUISHED) $(KMIN) $(LMIN) $(DMIN) $(DMAX) \
	    $(PATH_LCR_FILE) \
	    $(PATH_BAD_KMERS_FILE) \
	    $(INVESTIGATE_FINDLCRS) $(REDIR) ### findLCRs
	@echo "Finished."

################################################################################
# Target findIndelGroups: Analyze locally conserved regions to find
# insertions/deletions.
################################################################################

# A list of all .idlens files is already defined above: IDLEN_FILES

# Phony target to make or clean PATH_OVERLAPPING_INDEL_GROUPS_FILE and PATH_NONOVERLAPPING_INDEL_GROUPS_FILE files.
.PHONY: findIndelGroups
ifeq ($(CLEAN),)
findIndelGroups: $(PATH_OVERLAPPING_INDEL_GROUPS_FILE) $(PATH_NONOVERLAPPING_INDEL_GROUPS_FILE)
	@echo
	@echo "findIndelGroups files are up to date."
	@echo
else
findIndelGroups:
	@$(CMD_DELETE_WHEN_CLEANING) $(PATH_OVERLAPPING_INDEL_GROUPS_FILE) $(PATH_NONOVERLAPPING_INDEL_GROUPS_FILE)
	@echo "findIndelGroups output file(s) removed."
	@echo
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
        $(DIR_IGGPIPE_OUT) $(PATH_FIND_INDEL_GROUPS)
	@echo
	@echo "*** findIndelGroups PARAMS=$(PARAMS) ***"
	@echo "Find Indel Groups using locally conserved regions in $< and write them to two output files"
	$(TIME) \
		$(CMD_RSCRIPT) $(PATH_FIND_INDEL_GROUPS) $(WD) \
	    $(PATH_LCR_FILE) \
	    $(AMIN) $(AMAX) $(ADMIN) $(ADMAX) $(NDAMIN) $(MINFLANK) $(OVERLAP_REMOVAL) \
	    $(PATH_OVERLAPPING_INDEL_GROUPS_FILE) $(PATH_NONOVERLAPPING_INDEL_GROUPS_FILE) \
	    $(GENOME_LETTERS_SQUISHED) $(INVESTIGATE_ANALYZELCRS) \
	    $(IDLEN_FILES) $(REDIR) ### findIndelGroups
	@echo "Finished."

################################################################################
# Target getDNAseqsForPrimers: Extract DNA sequence around Indel Groups, for
# making primers.
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
.PHONY: getDNAseqsForPrimers
ifeq ($(CLEAN),)
getDNAseqsForPrimers: $(TARGET_DNASEQ)
	@echo
	@echo "getDNAseqsForPrimers file(s) for genome(s) $(GENOME) are up to date."
	@echo
else
getDNAseqsForPrimers:
	@$(CMD_DELETE_WHEN_CLEANING) $(TARGET_DNASEQ) $(DIR_GENOME_OUT_DATA)/extract*.txt $(DIR_GENOME_OUT_DATA)/seqs*.txt
	@echo "getDNAseqsForPrimers output file(s) for genome(s) $(GENOME) removed."
	@echo
endif

# DNA_SEQ_FILES is multiple targets, one per genome.
# Here, % is a genome number.

$(DNA_SEQ_FILES) : $(PFX_DNA_SEQS_PATH)_%.dnaseqs : $$(PATH_GENOME_FASTA_$$*) \
        $(PFX_GENOME_DATA_FILE)$$*.contigs $(PATH_OVERLAPPING_INDEL_GROUPS_FILE) | \
        $(DIR_IGGPIPE_OUT) $(DIR_GENOME_OUT_DATA) \
        $(PATH_GET_DNA_SEQS_FOR_PRIMERS) $(PATH_GET_SEQS_FASTA)
	@echo
	@echo "*** getDNAseqsForPrimers PARAMS=$(PARAMS) GENOME=$* ***"
	@echo "Extract DNA sequence around Indel Groups and write to $@"
	$(TIME) \
		$(CMD_RSCRIPT) $(PATH_GET_DNA_SEQS_FOR_PRIMERS) $(WD) \
	    $(PATH_OVERLAPPING_INDEL_GROUPS_FILE) $* \
	    $@ \
		$(DIR_GENOME_OUT_DATA) $(EXTENSION_LEN) \
		$(CMD_PERL) $(PATH_GET_SEQS_FASTA) \
		$(PATH_GENOME_FASTA_$*) \
		$(PFX_GENOME_DATA_FILE)$*.contigs $(INVESTIGATE_GETDNASEQS) $(REDIR) ### getDNAseqsForPrimers_$*
	@echo "Finished."

################################################################################
# Target findPrimers: Run primer3 to search for primers around Indel Groups.
################################################################################

# A list of all .dnaseqs files is already defined above: DNA_SEQ_FILES

# Phony target to make or clean PATH_NONVALIDATED_MARKER_FILE file.
.PHONY: findPrimers
ifeq ($(CLEAN),)
findPrimers: $(PATH_NONVALIDATED_MARKER_FILE)
	@echo
	@echo "findPrimers files are up to date."
	@echo
else
findPrimers:
	@$(CMD_DELETE_WHEN_CLEANING) $(PATH_NONVALIDATED_MARKER_FILE) $(PATH_PRIMER3_IN) $(PATH_PRIMER3_OUT)
	@echo "findPrimers output file(s) removed."
	@echo
endif

# PATH_NONVALIDATED_MARKER_FILE target.

$(PATH_NONVALIDATED_MARKER_FILE) : $(PATH_OVERLAPPING_INDEL_GROUPS_FILE) $(DNA_SEQ_FILES) \
        $(PATH_PRIMER3_SETTINGS) | $(DIR_PRIMER_DATA) $(PATH_FIND_PRIMERS)
	@echo
	@echo "*** findPrimers PARAMS=$(PARAMS) ***"
	@echo "Find primers around Indel Groups in $< and write to $@"
	$(TIME) \
		$(CMD_RSCRIPT) $(PATH_FIND_PRIMERS) $(WD) \
		$(PATH_OVERLAPPING_INDEL_GROUPS_FILE) \
		$(PFX_DNA_SEQS_PATH) \
		$(PATH_NONVALIDATED_MARKER_FILE) \
		$(CMD_PRIMER3CORE) $(PATH_PRIMER3_SETTINGS) $(DIR_PRIMER3CONFIG) \
		$(PATH_PRIMER3_IN) $(PATH_PRIMER3_OUT) \
		$(INVESTIGATE_FINDPRIMERS) $(REDIR) ### findPrimers
	@echo "Finished."

################################################################################
# Target ePCRtesting: Run e-PCR to test all primer pairs for proper amplicon
# size.
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
.PHONY: ePCRtesting
ifeq ($(CLEAN),)
ePCRtesting: $(TARGET_BAD_MARKER)
	@echo
	@echo "ePCRtesting file(s) for genome(s) $(GENOME) are up to date."
	@echo
else
ePCRtesting:
	@$(CMD_DELETE_WHEN_CLEANING) $(TARGET_BAD_MARKER) $(DIR_GENOME_OUT_DATA)/*.epcr.in $(DIR_GENOME_OUT_DATA)/*.epcr.out
	@echo "ePCRtesting output file(s) for genome(s) $(GENOME) removed."
	@echo
endif

# BAD_MARKER_FILES is multiple targets, one per genome.
# Here, % is a genome number.

$(BAD_MARKER_FILES) : $(PFX_BAD_MARKER_ERROR_PATH)_%.bad.tsv : $$(PATH_GENOME_FASTA_$$*) \
        $(PATH_NONVALIDATED_MARKER_FILE) | $(DIR_IGGPIPE_OUT) $(DIR_GENOME_OUT_DATA) $(DIR_PRIMER_DATA) \
        $(PATH_EPCR_TESTING)
	@echo
	@echo "*** ePCRtesting PARAMS=$(PARAMS) GENOME=$* ***"
	@echo "Use e-PCR to test marker primer pairs and write errors to $@"
	$(TIME) \
		$(CMD_RSCRIPT) $(PATH_EPCR_TESTING) $(WD) \
	    $(PATH_NONVALIDATED_MARKER_FILE) $* $@ \
	    $(DIR_GENOME_OUT_DATA) $(CMD_EPCR) \
	    $(EPCR_MAX_DEV) $(EPCR_WORD_SIZE) $(EPCR_MAX_MISMATCH) $(EPCR_MAX_GAPS) \
		$(PATH_GENOME_FASTA_$*) $(INVESTIGATE_EPCRTESTING) $(REDIR) ### ePCRtesting_$*
	@echo "Finished."

################################################################################
# Target removeBadMarkers: Read bad marker files that failed e-PCR, then read
# full marker file, remove the bad markers from it, and write new files of good
# ones, one set overlapping, the other non-overlapping.  ID numbers are added to
# the markers.
################################################################################

# A list of all .bad.tsv files is already defined above: BAD_MARKER_FILES

# Phony target to make or clean PATH_OVERLAPPING_MARKERS_FILE and PATH_NONOVERLAPPING_MARKERS_FILE.
.PHONY: removeBadMarkers
ifeq ($(CLEAN),)
removeBadMarkers: $(PATH_OVERLAPPING_MARKERS_FILE) $(PATH_NONOVERLAPPING_MARKERS_FILE)
	@echo
	@echo "removeBadMarkers files are up to date."
	@echo
else
removeBadMarkers:
	@$(CMD_DELETE_WHEN_CLEANING) $(PATH_OVERLAPPING_MARKERS_FILE) $(PATH_NONOVERLAPPING_MARKERS_FILE)
	@echo "removeBadMarkers output file(s) removed."
	@echo
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
        $(PATH_RMV_BAD_MARKERS)
	@echo
	@echo "*** removeBadMarkers PARAMS=$(PARAMS) ***"
	@echo "Remove markers identified by e-PCR as bad from $< and write good ones to two output files"
	$(TIME) \
		$(CMD_RSCRIPT) $(PATH_RMV_BAD_MARKERS) $(WD) \
	    $(OVERLAP_REMOVAL) $(ID_PREFIX) \
	    $(PATH_NONVALIDATED_MARKER_FILE) \
	    $(PFX_BAD_MARKER_ERROR_PATH) \
		$(PATH_OVERLAPPING_MARKERS_FILE) \
		$(PATH_NONOVERLAPPING_MARKERS_FILE) \
		$(INVESTIGATE_RMVBADMARKERS) $(REDIR) ### removeBadMarkers
	@echo "Finished."

################################################################################
# Target plotMarkers: Make density plots of final "good" candidate IGG markers,
# both overlapping and non-overlapping.
################################################################################

# A list of all .idlens files is already defined above: IDLEN_FILES

# The full name of the counts plot output file.
MARKER_COUNTS_FILE := $(PFX_MARKER_COUNTS_PATH).plot.pdf

# A list of all .png files to be produced.
MARKER_DENSITY_FILES := $(foreach G,$(GENOME_LETTERS),$(PFX_MARKER_DENSITY_PATH)_$(G).plot.png)

# Phony target to make or clean MARKER_COUNTS_FILE and MARKER_DENSITY_FILES files.
.PHONY: plotMarkers
ifeq ($(CLEAN),)
plotMarkers: $(MARKER_COUNTS_FILE) $(MARKER_DENSITY_FILES)
	@echo
	@echo "plotMarkers files are up to date."
	@echo
else
plotMarkers:
	@$(CMD_DELETE_WHEN_CLEANING) $(MARKER_COUNTS_FILE) $(MARKER_DENSITY_FILES)
	@echo "plotMarkers output file(s) removed."
	@echo
endif

# MARKER_COUNTS_FILE and MARKER_DENSITY_FILES target.

$(MARKER_COUNTS_FILE) $(MARKER_DENSITY_FILES) : $(PATH_OVERLAPPING_MARKERS_FILE) $(PATH_NONOVERLAPPING_MARKERS_FILE) \
        $(IDLEN_FILES) | $(PATH_PLOT_MARKERS)
	@echo
	@echo "*** plotMarkers PARAMS=$(PARAMS) ***"
	@echo "Make density plots of 'good' candidate IGG markers to output files"
	$(TIME) \
		$(CMD_RSCRIPT) $(PATH_PLOT_MARKERS) $(WD) $(PLOT_NDAMIN) $(PLOT_ALPHA) \
	    $(PATH_OVERLAPPING_MARKERS_FILE) \
	    $(PATH_NONOVERLAPPING_MARKERS_FILE) \
		$(PFX_MARKER_COUNTS_PATH) \
		$(PFX_MARKER_DENSITY_PATH) \
		$(IDLEN_FILES) $(REDIR) ### plotMarkers
	@echo "Finished."

################################################################################
# Target LCRsToLCBs: Read LCRs data file and convert to LCBs.
################################################################################

# Phony target to make or clean PATH_LCB_FILE file.
# If variable PARAMS is not defined, show basic usage info, else make
# LCRsToLCBs output file target.
.PHONY: LCRsToLCBs
ifeq ($(PARAMS),)
LCRsToLCBs:
	@echo
	@echo "You must specify a PARAMS file:"
	@echo
	@echo "$(INDENT)make PARAMS=<allParametersFile> LCRsToLCBs"
	@echo
else ifeq ($(CLEAN),)
LCRsToLCBs: $(PATH_LCB_FILE)
	@echo
	@echo "LCBs data file is up to date."
	@echo
else
LCRsToLCBs:
	@$(CMD_DELETE_WHEN_CLEANING) $(PATH_LCB_FILE)
	@echo "LCBs data file removed."
	@echo
endif

# PATH_LCB_FILE target.

$(PATH_LCB_FILE) : $(PATH_LCR_FILE) | $(DIR_GENOME_OUT_DATA)
	@echo
	@echo "*** LCRsToLCBs PARAMS=$(PARAMS) ***"
	@echo "Convert LCRs in $< to LCBs and write to $@"
	$(TIME) \
		$(CMD_RSCRIPT) $(PATH_LCRS_TO_LCBS) $(WD) $(PATH_LCR_FILE) $(PATH_LCB_FILE) $(REDIR) ### LCRsToLCBs
	@echo "Finished."

################################################################################
# Target getDNAseqsForIndelsSNPs: Read input file and extract DNA sequences from
# each genome, then write combined data to output file.
################################################################################

# A list of all FASTA files already defined above: FASTA_FILES

# Phony target to make or clean PATH_INDELS_SNPS_SEQS_FILE file.
# If variable PARAMS is not defined, show basic usage info, else make
# getDNAseqsForIndelsSNPs output file target.
.PHONY: getDNAseqsForIndelsSNPs
ifeq ($(PARAMS),)
getDNAseqsForIndelsSNPs:
	@echo
	@echo "You must specify a PARAMS file:"
	@echo
	@echo "$(INDENT)make PARAMS=<allParametersFile> getDNAseqsForIndelsSNPs"
	@echo
else ifeq ($(CLEAN),)
getDNAseqsForIndelsSNPs: $(PATH_INDELS_SNPS_SEQS_FILE)
	@echo
	@echo "Indels and SNPs sequence file is up to date."
	@echo
else
getDNAseqsForIndelsSNPs:
	@$(CMD_DELETE_WHEN_CLEANING) $(PATH_INDELS_SNPS_SEQS_FILE)
	@echo "Indels and SNPs sequence file removed."
	@echo
endif

# PATH_INDELS_SNPS_SEQS_FILE target.

$(PATH_INDELS_SNPS_SEQS_FILE) : $(PATH_INDELS_SNPS_INPUT_FILE) $(FASTA_FILES) | \
        $(DIR_GENOME_OUT_DATA) $(PATH_GET_DNA_SEQS_FOR_INDELS_SNPS) $(PATH_GET_SEQS_FASTA)
	@echo
	@echo "*** getDNAseqsForIndelsSNPs PARAMS=$(PARAMS) ***"
	@echo "Get DNA sequences for $< and write merged data to $@"
	$(TIME) \
		$(CMD_RSCRIPT) $(PATH_GET_DNA_SEQS_FOR_INDELS_SNPS) $(WD) $(K) \
	    $(PATH_INDELS_SNPS_INPUT_FILE) $(PATH_INDELS_SNPS_SEQS_FILE) \
	    $(DIR_GENOME_OUT_DATA) \
		$(CMD_PERL) $(PATH_GET_SEQS_FASTA) \
		$(INVESTIGATE_GET_DNA_SEQS_FOR_INDELS_SNPS) \
		$(FASTA_FILES) $(REDIR) ### getDNAseqsForIndelsSNPs
	@echo "Finished."

################################################################################
# Target IndelsSNPs: Read input file and perform alignments, then search them
# for Indels and SNPs.
################################################################################

# Do this to avoid case annoyance.
indelsSNPs: IndelsSNPs

# Phony target to make or clean PATH_INDELS_OUTPUT_FILE and PATH_SNPS_OUTPUT_FILE files.
# If variable PARAMS is not defined, show basic usage info, else make IndelsSNPs
# output file target.
.PHONY: IndelsSNPs
ifeq ($(PARAMS),)
IndelsSNPs:
	@echo
	@echo "You must specify a PARAMS file:"
	@echo
	@echo "$(INDENT)make PARAMS=<allParametersFile> IndelsSNPs"
	@echo
else ifeq ($(CLEAN),)
IndelsSNPs: $(PATH_INDELS_OUTPUT_FILE) $(PATH_SNPS_OUTPUT_FILE)
	@echo
	@echo "Indels and SNPs files are up to date."
	@echo
else
IndelsSNPs:
	@$(CMD_DELETE_WHEN_CLEANING) $(PATH_INDELS_OUTPUT_FILE) $(PATH_SNPS_OUTPUT_FILE)
	@echo "Indels and SNPs output file(s) removed."
	@echo
endif

# PATH_INDELS_OUTPUT_FILE and PATH_SNPS_OUTPUT_FILE target.

$(PATH_INDELS_OUTPUT_FILE) $(PATH_SNPS_OUTPUT_FILE) : $(PATH_INDELS_SNPS_SEQS_FILE) | \
        $(PATH_ALIGN_AND_GET_INDELS_SNPS)
	@echo
	@echo "*** IndelsSNPs PARAMS=$(PARAMS) ***"
	@echo "Align sequences of $< and find Indels and SNPs and write them to output files"
	$(TIME) \
		$(CMD_RSCRIPT) $(PATH_ALIGN_AND_GET_INDELS_SNPS) $(WD) \
	    $(PATH_INDELS_SNPS_SEQS_FILE) $(PATH_INDELS_OUTPUT_FILE) $(PATH_SNPS_OUTPUT_FILE) \
		$(CMD_ALIGNER) $(MAX_INDELS_PER_KBP) $(MAX_SNPS_PER_KBP) $(SCRAMBLE_SEQUENCE) \
		$(INVESTIGATE_ALIGN_AND_GET_INDELS_SNPS) $(REDIR) ### IndelsSNPs
	@echo "Finished."

################################################################################
# Target plotIndels: Make plots of information about Indels found within Indel
# groups.
################################################################################

# Do this to avoid case annoyance.
plotindels: plotIndels

# Phony target to make or clean PATH_INDELS_PLOT_FILE file.
# If variable PARAMS is not defined, show basic usage info, else make plotIndels
# output file target.
.PHONY: plotIndels
ifeq ($(PARAMS),)
plotIndels:
	@echo
	@echo "You must specify a PARAMS file:"
	@echo
	@echo "$(INDENT)make PARAMS=<allParametersFile> plotIndels"
	@echo
else ifeq ($(CLEAN),)
plotIndels: $(PATH_INDELS_PLOT_FILE)
	@echo
	@echo "plotIndels file is up to date."
	@echo
else
plotIndels:
	@$(CMD_DELETE_WHEN_CLEANING) $(PATH_INDELS_PLOT_FILE)
	@echo "plotIndels output file(s) removed."
	@echo
endif

# PATH_INDELS_PLOT_FILE target.

$(PATH_INDELS_PLOT_FILE) : $(PATH_INDELS_OUTPUT_FILE) | $(PATH_PLOT_INDELS)
	@echo
	@echo "*** plotIndels PARAMS=$(PARAMS) ***"
	@echo "Make plots of Indel information to file $@"
	$(TIME) \
		$(CMD_RSCRIPT) $(PATH_PLOT_INDELS) $(WD) $(PATH_INDELS_OUTPUT_FILE) \
	    $(PATH_INDELS_PLOT_FILE) $(REDIR) ### plotIndels
	@echo "Finished."

################################################################################
# Target documents: Convert .asciidoc documents to .html files using asciidoc
# application (developer only).
#
# Target README: Update README document
# Target INSTALL: Update INSTALL document
# Target RUN: Update RUN document
#
# Requirements:
#   1. "python" is installed on your system (https://www.python.org/downloads/)
#   2. "asciidoc" tools are installed (http://asciidoc.org/INSTALL.html)
# Make a symbolic link in a directory on your path that is named asciidoc,
# pointing to the "asciidoc.py" file in the asciidoc install directory.
# For example, if asciidoc-8.6.9 were downloaded into ~/src, and if ~/bin
# is on your path, this would make a symbolic link:
#   ln -s ~/src/asciidoc-8.6.9/asciidoc.py ~/bin/asciidoc
# Also make sure the file "asciidoc.py" is executable:
#   chmod +x ~/src/asciidoc-8.6.9/asciidoc.py
################################################################################

documents: README INSTALL RUN

README: README.html

README.html: README.asciidoc
	asciidoc -a VERSION=$(VERSION) -b html -o README.html README.asciidoc

INSTALL: INSTALL.html

INSTALL.html: INSTALL.asciidoc
	asciidoc -a VERSION=$(VERSION) -b html -o INSTALL.html INSTALL.asciidoc

RUN: RUN.html

RUN.html: RUN.asciidoc
	asciidoc -a VERSION=$(VERSION) -b html -o RUN.html RUN.asciidoc

################################################################################
# End of file.
################################################################################
