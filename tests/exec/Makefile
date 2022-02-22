# -*- mode: Make -*-

scripts := $(wildcard *.sh)
 .PHONY: $(scripts)	


.PHONY: test
test:
	@echo " "
	@echo " *************************** "
	@echo "   testing the sixte tools   "
	@echo " *************************** "
	@echo " "
	make $(scripts)	
	make clean

$(scripts): 
	./$(basename $@).sh > log/$(basename $@).log
	@echo " "


make clean:
	rm -f output/*.fits
