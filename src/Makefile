# This Makefile calls other Makefiles within the /src/ subdirectories.
# The other Makefiles should mostly provide installing, cleaning for themselves.
# The exception is for the Cmake generated files, they are cleaned inline of 'all'

.PHONY: all install clean distclean

all:
	cd build_pip; make all
	cd build_myddpost; make all

install: all
	cd build_pip; make install
	cd build_myddpost; make install
	install --mode 755 -D toolmenu ../bin/.

clean:
	cd build_myddpost; make clean
	cd build_pip; make clean

distclean: clean
	rm -f ../bin/ConvertFile
	rm -f ../bin/pip
	rm -f ../bin/toolmenu
	rm -f ../bin/myddpostprocess

