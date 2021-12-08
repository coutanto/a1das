
clean:
	cd ../A1-XCORPY; make clean; cd ../A1-SOSFILT; make clean; #cd A1-RAW2STRAINPY; make clean

python:
	cd ../A1-XCORPY; make python; cd ../A1-SOSFILT; make python; #cd A1-RAW2STRAINPY; make python; 

test:
	cd A1-TEST; make -f Makefile.test

doc:
	pdoc --html -force -o documentation a1das
