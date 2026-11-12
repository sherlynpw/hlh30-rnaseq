# Top-level wrapper for the posthoc/rz/* figure subprojects.
# Delegates to each child Makefile (pca-stage, go-plot, sample-corr-heatmap,
# volcano) so you can build/upload everything at once.
#
#   make             # render every subproject
#   make upload      # sync every subproject's out/ back to S3
#   make clean       # rm out/ in every subproject
#   make distclean   # also clear download caches
#   make print       # list subprojects + their dirs

SUBDIRS ?= pca-stage pca-genotype go-plot sample-corr-heatmap volcano cluster-profiles

.PHONY: all upload clean distclean print $(SUBDIRS) \
        $(addprefix upload-,$(SUBDIRS)) \
        $(addprefix clean-,$(SUBDIRS)) \
        $(addprefix distclean-,$(SUBDIRS))
.DEFAULT_GOAL := all

all: $(SUBDIRS)

$(SUBDIRS):
	$(MAKE) -C $@

upload: $(addprefix upload-,$(SUBDIRS))
$(addprefix upload-,$(SUBDIRS)): upload-%:
	$(MAKE) -C $* upload

clean: $(addprefix clean-,$(SUBDIRS))
$(addprefix clean-,$(SUBDIRS)): clean-%:
	$(MAKE) -C $* clean

distclean: $(addprefix distclean-,$(SUBDIRS))
$(addprefix distclean-,$(SUBDIRS)): distclean-%:
	$(MAKE) -C $* distclean

print:
	@echo SUBDIRS: $(SUBDIRS)
	@for d in $(SUBDIRS); do echo "--- $$d ---"; $(MAKE) -C $$d print; done
