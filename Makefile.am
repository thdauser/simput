dist_bin_SCRIPTS =simput-install.sh
dist_bin_SCRIPTS+=simput-install.csh

# The sub-directories are built before the current directory.
# In order to change this, include "." in the list of SUBDIRS.
SUBDIRS=extlib libsimput tools sample

dist_noinst_DATA=INSTALL.txt mac_solve_simput_dependencies.sh
EXTRA_DIST=tests

# Define the headers that will be installed in $(includedir):
include_HEADERS=simputconfig.h


# make sure all library files are removed
dir_cfitsio=$(distdir)/extlib/cfitsio
dir_wcslib=$(distdir)/extlib/wcslib
dist-hook:
	chmod u+w $(dir_cfitsio) $(dir_wcslib)
	cd $(dir_cfitsio) && make clean
	cd $(dir_wcslib) && make distclean

simputdistcheck: dist
	case '$(DIST_ARCHIVES)' in \
	*.tar.gz*) \
	  GZIP=$(GZIP_ENV) gzip -dc $(distdir).tar.gz | $(am__untar) ;;\
	*.tar.bz2*) \
	  bzip2 -dc $(distdir).tar.bz2 | $(am__untar) ;;\
	*.tar.lz*) \
	  lzip -dc $(distdir).tar.lz | $(am__untar) ;;\
	*.tar.xz*) \
	  xz -dc $(distdir).tar.xz | $(am__untar) ;;\
	*.tar.Z*) \
	  uncompress -c $(distdir).tar.Z | $(am__untar) ;;\
	*.shar.gz*) \
	  GZIP=$(GZIP_ENV) gzip -dc $(distdir).shar.gz | unshar ;;\
	*.zip*) \
	  unzip $(distdir).zip ;;\
	esac
	mkdir $(distdir)/_build $(distdir)/_inst
	test -d $(distdir) || exit 0; \
	dc_install_base=`$(am__cd) $(distdir)/_inst && pwd | sed -e 's,^[^:\\/]:[\\/],/,'` \
	&& echo ${dc_install_base} \
	  && dc_destdir="$${TMPDIR-/tmp}/am-dc-$$$$/" \
	  && am__cwd=`pwd` \
	  && $(am__cd) $(distdir) \
	  && ./configure \
	    $(AM_DISTCHECK_CONFIGURE_FLAGS) \
	    $(DISTCHECK_CONFIGURE_FLAGS) \
	    --srcdir=. --prefix="$$dc_install_base" \
	    --without-pgplot \
	  && $(MAKE) $(AM_MAKEFLAGS) \
	  && $(MAKE) $(AM_MAKEFLAGS) dvi \
	  && $(MAKE) $(AM_MAKEFLAGS) check \
	  && $(MAKE) $(AM_MAKEFLAGS) test \
	  && $(MAKE) $(AM_MAKEFLAGS) install \
	  && $(MAKE) $(AM_MAKEFLAGS) installcheck \
	  && $(MAKE) $(AM_MAKEFLAGS) uninstall \
	  && rm -rf "$$dc_destdir" \
	  && $(MAKE) $(AM_MAKEFLAGS) dist \
	  && rm -rf $(DIST_ARCHIVES) \
	  && cd "$$am__cwd" \
	  || exit 1
	$(am__post_remove_distdir)
	@(echo "$(distdir) archives ready for distribution: "; \
	  list='$(DIST_ARCHIVES)'; for i in $$list; do echo $$i; done) | \
	  sed -e 1h -e 1s/./=/g -e 1p -e 1x -e '$$p' -e '$$x'


.PHONY: test
test:
	cd tests/exec/ && make test
	cd tests/e2e/ && ./run_e2e_test.sh

crit_mac_version_major=10
crit_mac_version_minor=10
install-exec-hook:
	@if ([ 0$(OSX_VERSION_MAJOR) -eq $(crit_mac_version_major) ] && [ 0$(OSX_VERSION_MINOR) -gt $(crit_mac_version_minor) ]) || \
		[ 0$(OSX_VERSION_MAJOR) -gt $(crit_mac_version_major) ]; then\
		echo "\n  ***  Mac Version $(OSX_VERSION_MAJOR).$(OSX_VERSION_MINOR) detected *** ";\
		echo "\n   ... trying to solve known problems with this installation with this script:\n\n ";\
		./mac_solve_simput_dependencies.sh @prefix@;\
	fi
