# Copy final marker .tsv, .plot.pdf, and .plot.png files from directory given by
# the parameter file given by argument PARAMS= to the saveOutputFiles directory,
# same subdirectory name as given by DIR_SCARF_OUT in the PARAMS= file, but with
# simpler file names.

# Delete created files if error occurs creating them.
.DELETE_ON_ERROR:

# Run all commands of a target build using a single shell invocation.
.ONESHELL:

# We will use secondary expansion.
.SECONDEXPANSION:

# Use this shell.
SHELL=/bin/sh

# Some useful variables.
EMPTY :=
SPACE := $(EMPTY) $(EMPTY)
INDENT := $(EMPTY)    $(EMPTY)

# Cancel C-compiler implicit rules.
% : %.cpp
%.o : %.cpp
% : %.o

# Include parameter definition file.
include $(PARAMS)

# Get suffix of parameter definition file, which becomes part of the name of the
# copied files.
PARAM_SFX := $(patsubst allParameters.%,%,$(PARAMS))

# If variable PARAMS is not defined, show basic usage info, else do the copy.
ifeq ($(PARAMS),)
all:
	@echo "For basic usage, use this make command for example:"
	@echo
	@echo "$(INDENT)make -f CopyFiles.mak PARAMS=allParameters.test"
	@echo 
else

# The full name of the counts plot output file.
MARKER_COUNTS_FILE := $(PFX_MARKER_COUNTS_PATH).plot.pdf

# Names of all density plot files and all copied-to names.
DENSITY_PLOT_FILES := $(foreach G,$(GENOME_LETTERS),$(PFX_MARKER_DENSITY_PATH)_$(G).plot.png)
COPIED_DENSITY_PLOT_FILES := $(foreach G,$(GENOME_LETTERS),MarkerDensity_$(G).$(PARAM_SFX).plot.png)

# Define variables DENSITY_x := #, where x = a genome letter, equal to that density plot .plot.png file.
$(foreach G,$(GENOME_LETTERS),$(eval DENSITY_$(G) := $(PFX_MARKER_DENSITY_PATH)_$(G).plot.png))

# The output directory to copy files to.
COPY_TO := saveOutputFiles/$(DIR_SCARF_OUT)

all: $(COPIED_DENSITY_PLOT_FILES) $(PATH_OVERLAPPING_MARKERS_FILE) \
        $(PATH_NONOVERLAPPING_MARKERS_FILE) $(MARKER_COUNTS_FILE) | $(COPY_TO)
	@cp $(PATH_OVERLAPPING_MARKERS_FILE) $(COPY_TO)/MarkersOverlapping.$(PARAM_SFX).tsv
	@cp $(PATH_NONOVERLAPPING_MARKERS_FILE) $(COPY_TO)/MarkersNonoverlapping.$(PARAM_SFX).tsv
	@cp $(MARKER_COUNTS_FILE) $(COPY_TO)/MarkerCounts.$(PARAM_SFX).plot.pdf

$(COPY_TO):
	@mkdir -p $(COPY_TO)

$(COPIED_DENSITY_PLOT_FILES) : MarkerDensity_%.$(PARAM_SFX).plot.png : $(DENSITY_PLOT_FILES)
	cp $(DENSITY_$*) $(COPY_TO)/$@
endif

# Remove the destination files.

clean:
	-rm $(COPY_TO)/MarkersOverlapping.$(PARAM_SFX).tsv
	-rm $(COPY_TO)/MarkersNonoverlapping.$(PARAM_SFX).tsv
	-rm $(COPY_TO)/MarkerCounts.$(PARAM_SFX).plot.pdf
	-rm $(COPY_TO)/$(COPIED_DENSITY_PLOT_FILES)
