# The sub-directories are built before the current directory.
# In order to change this, include "." in the list of SUBDIRS.
SUBDIRS = ape

############ HEADERS #################

# Define the headers that will be installed in $(includedir):
include_HEADERS=pil.h pil_error.h

install-data-local:
	mkdir -p @includedir@/ape
	mv @includedir@/ape_*.h @includedir@/ape/
