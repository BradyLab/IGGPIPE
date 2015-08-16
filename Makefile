# This is a Makefile for running the SCARF pipeline on genome files to produce
# candidate SCAR markers.  Before running this, edit file allParameters.template.
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
# Include parameter definition file (there might not be one).
################################################################################
include $(PARAMS)

################################################################################
# Default make target.
# If variable PARAMS is not defined, show basic usage info, else set target ALL.
################################################################################
ifeq ($(PARAMS),)
all:
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
else
all: ALL
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

# Define variables GENOME_x := #, where x = a genome letter, # is genome number.
# For example, if GENOME_LETTERS are "H P" we have GENOME_H := 1, GENOME_P := 2
$(foreach G,$(GENOME_NUMBERS),$(eval GENOME_$(word $(G),$(GENOME_LETTERS)) := $(G)))

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
ifneq ($(CLEAN_OUT_DIR),1)
ALL: all_getSeqInfo all_getKmers all_kmerStats all_sortKmers all_getContigFile \
    kmerIsect all_getGenomicPos all_splitKmers findLCRs findINDELs \
    all_getDNAseqs findPrimers all_ePCRtesting removeBadMarkers plotMarkers
	@echo
ifeq ($(CLEAN),)
	@echo "ALL files are up to date"
else
	@echo "ALL files are cleaned"
endif
else
ALL:
	@echo
	@echo "Removing all files from the output directory."
	@# Try to remove output directory first.  If trashing, entire thing will be in trash.
	-$(CMD_DELETE_WHEN_CLEANING) $(DIR_SCARF_OUT)
	-$(CMD_DELETE_WHEN_CLEANING) $(DIR_GENOME_OUT_DATA)/*
	-$(CMD_DELETE_WHEN_CLEANING) $(DIR_SPLIT_KMERS)/*
	-$(CMD_DELETE_WHEN_CLEANING) $(DIR_KMERS)/*
	-$(CMD_DELETE_WHEN_CLEANING) $(DIR_PRIMER_DATA)/*
	-$(CMD_DELETE_WHEN_CLEANING) $(DIR_SCARF_OUT)/*
	-$(CMD_DELETE_WHEN_CLEANING) $(DIR_SCARF_OUT)
	@echo "Removed all files from the output directory."
endif
	@echo

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

# Target for getting or cleaning split of common unique k-mer genomic positions file for all genomes.
all_splitKmers:
	@$(MAKE) PARAMS=$(PARAMS) $(CLEAN) splitKmers GENOME=ALL

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
$(DIR_SCARF_OUT):
	@echo
	@echo "Creating directory $(DIR_SCARF_OUT)"
	mkdir -p $(DIR_SCARF_OUT)

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

# Target for making the split k-mer data output directory.
$(DIR_SPLIT_KMERS):
	@echo
	@echo "Creating directory $(DIR_SPLIT_KMERS)"
	mkdir -p $(DIR_SPLIT_KMERS)

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
	@echo "getSeqInfo files for genome(s) $(GENOME) are up to date."
else
getSeqInfo:
	@$(CMD_DELETE_WHEN_CLEANING) $(TARGET_IDLEN)
	@echo "getSeqInfo output files removed."
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
	@echo "getContigFile files for genome(s) $(GENOME) are up to date."
else
getContigFile:
	@$(CMD_DELETE_WHEN_CLEANING) $(TARGET_CONTIG)
	@echo "getContigFile output files removed."
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

# A list of all .kmers_0 files to be produced.
KMERS_0_FILES := $(foreach G,$(GENOME_NUMBERS),$(PFX_KMERS_DATA_FILE)$(G).kmers_0)

# Target .kmers_0 and .kmers_* file(s) for genome GENOME.
ifeq ($(GENOME),ALL)
TARGET_KMERS := $(KMERS_0_FILES)
KMERS_ALL_FILES := $(foreach G,$(GENOME_NUMBERS),$(PFX_KMERS_DATA_FILE)$(G).kmers_*)
else
TARGET_KMERS := $(PFX_KMERS_DATA_FILE)$(GENOME).kmers_0
KMERS_ALL_FILES := $(PFX_KMERS_DATA_FILE)$(GENOME).kmers_*
endif

# Phony target to make or clean TARGET_KMERS file(s).
ifeq ($(CLEAN),)
getKmers: $(TARGET_KMERS)
	@echo
	@echo "getKmers files for genome(s) $(GENOME) are up to date."
else
getKmers:
	@$(CMD_DELETE_WHEN_CLEANING) $(KMERS_ALL_FILES)
	@echo "getKmers output files removed."
endif

# KMERS_0_FILES is multiple targets, one per genome.
# This is a little troublesome because the actual files end with _0 (and perhaps
# _1, etc.), while the ARGUMENT to jellyfish does not include _0.
# Here, % is a genome number.

$(KMERS_0_FILES) : $(PFX_KMERS_DATA_FILE)%.kmers_0 : $$(PATH_GENOME_FASTA_$$*) | \
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

# A list of all .kmers_0 files already defined above: KMERS_0_FILES

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
	@echo "kmerStats files for genome(s) $(GENOME) are up to date."
else
kmerStats:
	@$(CMD_DELETE_WHEN_CLEANING) $(TARGET_STATS)
	@echo "kmerStats output files removed."
endif

# KMERS_STATS_FILES is multiple targets, one per genome.
# Here, % is a genome number.

$(KMERS_STATS_FILES) : $(PFX_KMERS_DATA_FILE)%.stats : $(PFX_KMERS_DATA_FILE)%.kmers_0 | \
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

# A list of all .kmers_0 files already defined above: KMERS_0_FILES

# A list of all .isect.sorted files to be produced.
SORTED_KMERS_FILES := $(foreach G,$(GENOME_NUMBERS),$(PFX_KMERS_DATA_FILE)$(G).isect.sorted)

# Target .isect.sorted file(s) for genome GENOME.
ifeq ($(GENOME),ALL)
TARGET_SORTED_KMERS := $(SORTED_KMERS_FILES)
else
TARGET_SORTED_KMERS := $(PFX_KMERS_DATA_FILE)$(GENOME).isect.sorted
endif

# Phony target to make or clean TARGET_SORTED_KMERS file(s).
ifeq ($(CLEAN),)
sortKmers: $(TARGET_SORTED_KMERS)
	@echo
	@echo "sortKmers files for genome(s) $(GENOME) are up to date."
else
sortKmers:
	@$(CMD_DELETE_WHEN_CLEANING) $(TARGET_SORTED_KMERS)
	@echo "sortKmers output files removed."
endif

# SORTED_KMERS_FILES is multiple targets, one per genome.
# Here, % is a genome number.

$(SORTED_KMERS_FILES) : $(PFX_KMERS_DATA_FILE)%.isect.sorted : $(PFX_KMERS_DATA_FILE)%.kmers_0 | \
        $(DIR_KMERS)
	@echo
	@echo "*** sortKmers PARAMS=$(PARAMS) $(CLEAN) GENOME=$* ***"
	@echo "Sorting text-based $(K)-mers from $< into $@"
	@echo
	@echo "Extracting text-based $(K)-mers from binary file $< into $@.tmp"
	$(TIME) $(PATH_JELLYFISH) dump -c -o $@.tmp $< $(REDIR)
	@echo "Sorting text-based $(K)-mers from $@.tmp into $@.tsv.  This can take a long time."
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
	@echo "kmerIsect output files removed."
endif

# .isect target.

$(PATH_ISECT_KMERS) : $(SORTED_KMERS_FILES) | $(DIR_KMERS) $(PATH_PERL) $(PATH_KMER_ISECT)
	@echo
	@echo "*** kmerIsect PARAMS=$(PARAMS) $(CLEAN) ***"
	@echo "Intersecting unique $(K)-mers from .isect.sorted files into $@"
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
	@echo "getGenomicPos files for genome(s) $(GENOME) are up to date."
else
getGenomicPos:
	@$(CMD_DELETE_WHEN_CLEANING) $(TARGET_ISECT)
	@echo "getGenomicPos output files removed."
endif

# ISECT_KMER_FILES is multiple targets, one per genome.
# Here, % is a genome number.

$(ISECT_KMER_FILES) : $(PFX_KMERS_DATA_FILE)%.isect : $(PATH_ISECT_KMERS) $$(PATH_GENOME_FASTA_$$*) | \
        $(DIR_KMERS) $(PATH_FINDMERS)
	@echo
	@echo "*** getGenomicPos PARAMS=$(PARAMS) $(CLEAN) GENOME=$* ***"
	@echo "Finding genomic positions of common unique $(K)-mers in $< into $@"
	$(TIME) $(PATH_FINDMERS) $(ARGS_FINDMER) -v2 $(PATH_GENOME_FASTA_$*) $< $@ $(REDIR)
	@echo "Finished."

################################################################################
# splitKmers: Get genomic position of common unique k-mers.  Argument: GENOME
################################################################################

# A list of all FASTA files already defined above: FASTA_FILES

# This is problematic because the output file names include the sequence IDs,
# which we don't know here.  So, what we will do is make an empty target file
# with suffix <genome>.isect.split using 'touch' to record the time we finished
# creating the actual split k-mer output files.

SPLIT_KMER_EMPTY_TARGET := $(foreach G,$(GENOME_NUMBERS),$(PFX_SPLIT_UNIQ_KMERS_FILE)$(G).isect.split)
SPLIT_KMER_ALL_TARGETS := $(foreach G,$(GENOME_NUMBERS),$(PFX_SPLIT_UNIQ_KMERS_FILE)$(G)_*.isect.split)

# Target .isect.split file(s) for genome GENOME.
ifeq ($(GENOME),ALL)
TARGET_SPLIT := $(SPLIT_KMER_EMPTY_TARGET)
else
TARGET_SPLIT := $(PFX_SPLIT_UNIQ_KMERS_FILE)$(GENOME).isect.split
endif

# Phony target to make or clean TARGET_SPLIT file(s).
ifeq ($(CLEAN),)
splitKmers: $(TARGET_SPLIT)
	@echo
	@echo "splitKmers files for genome(s) $(GENOME) are up to date."
else
splitKmers:
	@$(CMD_DELETE_WHEN_CLEANING) $(TARGET_SPLIT) $(SPLIT_KMER_ALL_TARGETS)
	@echo "splitKmers output files removed."
endif

# These variables' values are expanded (deferred) when used in the recipe below.

OUT_PREFIX = $(PFX_SPLIT_UNIQ_KMERS_FILE)$*_
GUIDE_FILE = $(if $(findstring 1,$*),$(EMPTY),$(PFX_KMERS_DATA_FILE)1.isect)

# SPLIT_KMER_EMPTY_TARGET is multiple targets, one per genome.
# Here, % is a genome number.

$(SPLIT_KMER_EMPTY_TARGET) : $(PFX_SPLIT_UNIQ_KMERS_FILE)%.isect.split : $(PFX_KMERS_DATA_FILE)%.isect | \
        $(DIR_SPLIT_KMERS) $(PATH_RSCRIPT) $(PATH_SPLIT_KMERS)
	@echo
	@echo "*** splitKmers PARAMS=$(PARAMS) $(CLEAN) GENOME=$* ***"
	@echo "Splitting $(K)-mer position file $< into files $(PFX_SPLIT_UNIQ_KMERS_FILE)$*_<refSeqID>.isect.split"
	$(TIME) $(PATH_RSCRIPT) $(PATH_SPLIT_KMERS) $(WD) $(OUT_PREFIX) $< $(GUIDE_FILE) $(REDIR)
	@touch $@
	@echo "Finished."

################################################################################
# findLCRs: Find locally conserved regions using common unique k-mers that have
# been split into separate files based on genome for which the k-mer position
# data applies and sequence ID within that genome.
################################################################################

# A list of all split k-mer empty-target files is already defined above: SPLIT_KMER_EMPTY_TARGET

# Phony target to make or clean PATH_LCR_FILE and PATH_BAD_KMERS_FILE files.
ifeq ($(CLEAN),)
findLCRs: $(PATH_LCR_FILE) $(PATH_BAD_KMERS_FILE)
	@echo
	@echo "findLCRs files are up to date."
else
findLCRs:
	@$(CMD_DELETE_WHEN_CLEANING) $(PATH_LCR_FILE) $(PATH_BAD_KMERS_FILE)
	@echo "findLCRs output files removed."
endif

# PATH_LCR_FILE and PATH_BAD_KMERS_FILE targets are built at the same time.
# Define target names with % in them, so we can create a pattern target.  The %
# tells make that all the target files are made at one time by the recipe.  Here
# the only thing in common between the names is "$(SFX_LCR_FILE).tsv", and
# $(SFX_LCR_FILE) might be empty, but % at least matches .tsv.
PTN_LCR_FILE := $(DIR_MAIN_OUTPUT)/$(PFX_LCR_FILE)%
PTN_BAD_KMERS_FILE := $(DIR_MAIN_OUTPUT)/$(PFX_BAD_KMERS_FILE)%

# Use the patterns in a pattern target.
$(PTN_LCR_FILE) $(PTN_BAD_KMERS_FILE) : $(SPLIT_KMER_EMPTY_TARGET) | \
        $(DIR_SCARF_OUT) $(PATH_RSCRIPT) $(PATH_FIND_LCRS)
	@echo
	@echo "*** findLCRs PARAMS=$(PARAMS) $(CLEAN) ***"
	@echo "Finding locally conserved regions using common unique $(K)-mers from .isect files into $@"
	$(TIME) $(PATH_RSCRIPT) $(PATH_FIND_LCRS) $(WD) $(PFX_SPLIT_UNIQ_KMERS_FILE) \
	    $(GENOME_LETTERS_SQUISHED) $(KMIN) $(LMIN) $(DMIN) $(DMAX) \
	    $(PATH_LCR_FILE) \
	    $(PATH_BAD_KMERS_FILE) \
	    $(INVESTIGATE_FINDLCRS) $(REDIR)
	@echo "Finished."

################################################################################
# findINDELS: Analyze locally conserved regions to find insertions/deletions.
################################################################################

# A list of all .idlens files is already defined above: IDLEN_FILES

# Phony target to make or clean PATH_OVERLAPPING_INDELS_FILE and PATH_NONOVERLAPPING_INDELS_FILE files.
ifeq ($(CLEAN),)
findINDELs: $(PATH_OVERLAPPING_INDELS_FILE) $(PATH_NONOVERLAPPING_INDELS_FILE)
	@echo
	@echo "findINDELs files are up to date."
else
findINDELs:
	@$(CMD_DELETE_WHEN_CLEANING) $(PATH_OVERLAPPING_INDELS_FILE) $(PATH_NONOVERLAPPING_INDELS_FILE)
	@echo "findINDELs output files removed."
endif

# PATH_OVERLAPPING_INDELS_FILE and PATH_NONOVERLAPPING_INDELS_FILE targets are built at the same time.
# Define target names with % in them, so we can create a pattern target.  The %
# tells make that all the target files are made at one time by the recipe.  Here
# the only thing in common between the names is "$(SFX_INDELS_FILE).tsv", and
# $(SFX_INDELS_FILE) might be empty, but % at least matches .tsv.
PTN_OVERLAPPING_INDELS_FILE := $(DIR_MAIN_OUTPUT)/$(PFX_OVERLAPPING_INDELS_FILE)%
PTN_NONOVERLAPPING_INDELS_FILE := $(DIR_MAIN_OUTPUT)/$(PFX_NONOVERLAPPING_INDELS_FILE)%

# Use the patterns in a pattern target.
$(PTN_OVERLAPPING_INDELS_FILE) $(PTN_NONOVERLAPPING_INDELS_FILE) : $(PATH_LCR_FILE) $(IDLEN_FILES) | \
        $(DIR_SCARF_OUT) $(PATH_RSCRIPT) $(PATH_FIND_INDELS)
	@echo
	@echo "*** findINDELs PARAMS=$(PARAMS) $(CLEAN) ***"
	@echo "Find indels using locally conserved regions in $< and write them to two output files"
	$(TIME) $(PATH_RSCRIPT) $(PATH_FIND_INDELS) $(WD) \
	    $(PATH_LCR_FILE) \
	    $(AMIN) $(AMAX) $(ADMIN) $(ADMAX) $(NDAMIN) $(MINFLANK) $(OVERLAP_REMOVAL) \
	    $(PATH_OVERLAPPING_INDELS_FILE) $(PATH_NONOVERLAPPING_INDELS_FILE) \
	    $(GENOME_LETTERS_SQUISHED) $(INVESTIGATE_ANALYZELCRS) \
	    $(IDLEN_FILES) $(REDIR)
	@echo "Finished."

################################################################################
# getDNAseqs: Extract DNA sequence around indels, for making primers.
################################################################################

# A list of all FASTA files already defined above: FASTA_FILES

# Likewise, a list of all .contigs files is already defined above: CONTIG_FILES

# A list of all .dnaseqs files to be produced.
DNA_SEQ_FILES := $(foreach G,$(GENOME_NUMBERS),$(PFX_GENOME_DATA_FILE)$(G).dnaseqs)

# Target .dnaseqs file(s) for genome GENOME.
ifeq ($(GENOME),ALL)
TARGET_DNASEQ := $(DNA_SEQ_FILES)
else
TARGET_DNASEQ := $(PFX_GENOME_DATA_FILE)$(GENOME).dnaseqs
endif

# Phony target to make or clean TARGET_DNASEQ file(s).
ifeq ($(CLEAN),)
getDNAseqs: $(TARGET_DNASEQ)
	@echo
	@echo "getDNAseqs files for genome(s) $(GENOME) are up to date."
else
getDNAseqs:
	@$(CMD_DELETE_WHEN_CLEANING) $(TARGET_DNASEQ) $(DIR_GENOME_OUT_DATA)/extract*.txt $(DIR_GENOME_OUT_DATA)/seqs*.txt
	@echo "getDNAseqs output files removed."
endif

# DNA_SEQ_FILES is multiple targets, one per genome.
# Here, % is a genome number.

$(DNA_SEQ_FILES) : $(PFX_GENOME_DATA_FILE)%.dnaseqs : $$(PATH_GENOME_FASTA_$$*) \
        $(PFX_GENOME_DATA_FILE)$$*.contigs $(PATH_OVERLAPPING_INDELS_FILE) | \
        $(DIR_SCARF_OUT) $(DIR_GENOME_OUT_DATA) \
        $(PATH_RSCRIPT) $(PATH_GET_DNA_SEQS) $(PATH_GET_SEQS_FASTA)
	@echo
	@echo "*** getDNAseqs PARAMS=$(PARAMS) $(CLEAN) GENOME=$* ***"
	@echo "Extract DNA sequence around indels and write to $@"
	$(TIME) $(PATH_RSCRIPT) $(PATH_GET_DNA_SEQS) $(WD) \
	    $(PATH_OVERLAPPING_INDELS_FILE) $* \
	    $@ \
		$(DIR_GENOME_OUT_DATA) $(EXTENSION_LEN) \
		$(PATH_PERL) $(PATH_GET_SEQS_FASTA) \
		$(PATH_GENOME_FASTA_$*) \
		$(PFX_GENOME_DATA_FILE)$*.contigs $(INVESTIGATE_GETDNASEQS) $(REDIR)
	@echo "Finished."

################################################################################
# findPrimers: Run primer3 to search for primers around indels.
################################################################################

# A list of all .dnaseqs files is already defined above: DNA_SEQ_FILES

# Phony target to make or clean PATH_MARKER_DATA_FILE file.
ifeq ($(CLEAN),)
findPrimers: $(PATH_MARKER_DATA_FILE)
	@echo
	@echo "findPrimers files are up to date."
else
findPrimers:
	@$(CMD_DELETE_WHEN_CLEANING) $(PATH_MARKER_DATA_FILE) $(PATH_PRIMER3_DATA) $(PATH_PRIMER3_OUT)
	@echo "findPrimers output files removed."
endif

# PATH_MARKER_DATA_FILE target.

$(PATH_MARKER_DATA_FILE) : $(PATH_OVERLAPPING_INDELS_FILE) $(DNA_SEQ_FILES) | $(PATH_RSCRIPT) \
        $(DIR_PRIMER_DATA) $(PATH_FIND_PRIMERS) $(PATH_PRIMER3CORE) $(PATH_PRIMER3_SETTINGS)
	@echo
	@echo "*** findPrimers PARAMS=$(PARAMS) $(CLEAN) ***"
	@echo "Find primers around indels in $< and write to $@"
	$(TIME) $(PATH_RSCRIPT) $(PATH_FIND_PRIMERS) $(WD) \
		$(PATH_OVERLAPPING_INDELS_FILE) \
		$(PFX_GENOME_DATA_FILE) \
		$(PATH_MARKER_DATA_FILE) \
		$(PATH_PRIMER3CORE) $(PATH_PRIMER3_SETTINGS) \
		$(PATH_PRIMER3_DATA) $(PATH_PRIMER3_OUT) \
		$(INVESTIGATE_FINDPRIMERS) $(REDIR)
	@echo "Finished."

################################################################################
# ePCRtesting: Run e-PCR to test all primer pairs for proper amplicon size.
################################################################################

# A list of all FASTA files already defined above: FASTA_FILES

# A list of all .bad.tsv files to be produced.
BAD_MARKER_FILES := $(foreach G,$(GENOME_LETTERS),$(PFX_BAD_MARKER_ERROR_PATH)_$(G).bad.tsv)

# Target .bad.tsv file(s) for genome GENOME.
ifeq ($(GENOME),ALL)
TARGET_BAD_MARKER := $(BAD_MARKER_FILES)
else
TARGET_BAD_MARKER := $(PFX_BAD_MARKER_ERROR_PATH)_$(word $(GENOME),$(GENOME_LETTERS)).bad.tsv
endif

# Phony target to make or clean TARGET_BAD_MARKER file(s).
ifeq ($(CLEAN),)
ePCRtesting: $(TARGET_BAD_MARKER)
	@echo
	@echo "ePCRtesting files for genome(s) $(GENOME) are up to date."
else
ePCRtesting:
	@$(CMD_DELETE_WHEN_CLEANING) $(TARGET_BAD_MARKER) $(DIR_GENOME_OUT_DATA)/*.sts $(DIR_GENOME_OUT_DATA)/*.epcr.out
	@echo "ePCRtesting output files removed."
endif

# BAD_MARKER_FILES is multiple targets, one per genome.
# Here, % is a genome letter, which is a pain in the ass to convert back to a genome number.

$(BAD_MARKER_FILES) : $(PFX_BAD_MARKER_ERROR_PATH)_%.bad.tsv : $$(PATH_GENOME_FASTA_$$(GENOME_$$*)) \
        $(PATH_MARKER_DATA_FILE) | $(DIR_SCARF_OUT) $(DIR_GENOME_OUT_DATA) $(DIR_PRIMER_DATA) \
        $(PATH_RSCRIPT) $(PATH_EPCR_TESTING) $(PATH_EPCR)
	@echo
	@echo "*** ePCRtesting PARAMS=$(PARAMS) $(CLEAN) GENOME=$* ***"
	@echo "Use e-PCR to test marker primer pairs and write errors to $@"
	$(TIME) $(PATH_RSCRIPT) $(PATH_EPCR_TESTING) $(WD) \
	    $(PATH_MARKER_DATA_FILE) $(GENOME_$*) \
	    $@ \
	    $(DIR_GENOME_OUT_DATA) $(PATH_EPCR) \
	    $(EPCR_MAX_DEV) $(EPCR_WORD_SIZE) $(EPCR_MAX_MISMATCH) $(EPCR_MAX_GAPS) \
		$(PATH_GENOME_FASTA_$(GENOME_$*)) $(INVESTIGATE_EPCRTESTING) $(REDIR)
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
	@echo "removeBadMarkers output files removed."
endif

# PATH_OVERLAPPING_MARKERS_FILE and PATH_NONOVERLAPPING_MARKERS_FILE targets are built at the same time.
# Define target names with % in them, so we can create a pattern target.  The %
# tells make that all the target files are made at one time by the recipe.  Here
# the only thing in common between the names is "$(SFX_GOOD_MARKER_FILE).tsv", and
# $(SFX_GOOD_MARKER_FILE) might be empty, but % at least matches .tsv.
PTN_OVERLAPPING_MARKERS_FILE := $(DIR_MAIN_OUTPUT)/$(PFX_OVERLAPPING_MARKERS_FILE)%
PTN_NONOVERLAPPING_MARKERS_FILE := $(DIR_MAIN_OUTPUT)/$(PFX_NONOVERLAPPING_MARKERS_FILE)%

# Use the patterns in a pattern target.
$(PTN_OVERLAPPING_MARKERS_FILE) $(PTN_NONOVERLAPPING_MARKERS_FILE) : $(PATH_MARKER_DATA_FILE) $(BAD_MARKER_FILES) | \
        $(PATH_RSCRIPT) $(PATH_RMV_BAD_MARKERS)
	@echo
	@echo "*** removeBadMarkers PARAMS=$(PARAMS) $(CLEAN) ***"
	@echo "Remove markers identified by e-PCR as bad from $< and write good ones to two output files"
	$(TIME) $(PATH_RSCRIPT) $(PATH_RMV_BAD_MARKERS) $(WD) $(OVERLAP_REMOVAL) \
	    $(PATH_MARKER_DATA_FILE) \
	    $(PFX_BAD_MARKER_ERROR_PATH) \
		$(PATH_OVERLAPPING_MARKERS_FILE) \
		$(PATH_NONOVERLAPPING_MARKERS_FILE) \
		$(INVESTIGATE_RMVBADMARKERS) $(REDIR)
	@echo "Finished."

################################################################################
# plotMarkers: Make density plots of final "good" candidate SCAR markers, both
# overlapping and non-overlapping.
################################################################################

# A list of all .idlens files is already defined above: IDLEN_FILES

# The full name of the counts plot output file.
MARKER_COUNTS_FILE := $(PFX_MARKER_COUNTS_PATH).pdf

# A list of all .png files to be produced.
MARKER_DENSITY_FILES := $(foreach G,$(GENOME_LETTERS),$(PFX_MARKER_DENSITY_PATH)_$(G).png)

# Phony target to make or clean MARKER_COUNTS_FILE and MARKER_DENSITY_FILES files.
ifeq ($(CLEAN),)
plotMarkers: $(MARKER_COUNTS_FILE) $(MARKER_DENSITY_FILES)
	@echo
	@echo "plotMarkers files are up to date."
else
plotMarkers:
	@$(CMD_DELETE_WHEN_CLEANING) $(MARKER_COUNTS_FILE) $(MARKER_DENSITY_FILES)
	@echo "plotMarkers output files removed."
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
	@echo "Make density plots of 'good' candidate SCAR markers to file $@"
	$(TIME) $(PATH_RSCRIPT) $(PATH_PLOT_MARKERS) $(WD) $(PLOT_NDAMIN) $(PLOT_ALPHA) \
	    $(PATH_OVERLAPPING_MARKERS_FILE) \
	    $(PATH_NONOVERLAPPING_MARKERS_FILE) \
		$(PFX_MARKER_COUNTS_PATH) \
		$(PFX_MARKER_DENSITY_PATH) \
		$(IDLEN_FILES) $(REDIR)
	@echo "Finished."

################################################################################
# End of file.
################################################################################
