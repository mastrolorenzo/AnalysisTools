CC = g++
BDTLIB = bdt/src/bdtlib/
XMLLIB = bdt/src/tinyxml2/
CFLAGS = -I${BDTLIB} -I${XMLLIB} -c -O3 -g `root-config --cflags` -fPIC
LIBFLAGS = `root-config --libs` -g -O3 -lXMLIO -lMLP -lMinuit

libAutocategorizer.so: ${BDTLIB}CategoryReader.o ${BDTLIB}Node.o ${BDTLIB}Tree.o ${XMLLIB}tinyxml2.o
	$(CC) $(LIBFLAGS) -shared -o $@ ${BDTLIB}CategoryReader.o ${BDTLIB}Node.o ${BDTLIB}Tree.o ${XMLLIB}tinyxml2.o

${BDTLIB}CategoryReader.o: ${BDTLIB}CategoryReader.cxx
	$(CC) $(CFLAGS) ${BDTLIB}CategoryReader.cxx -o $@

${BDTLIB}Node.o: ${BDTLIB}Node.cxx
	$(CC) $(CFLAGS) ${BDTLIB}Node.cxx -o $@

${BDTLIB}Tree.o: ${BDTLIB}Tree.cxx
	$(CC) $(CFLAGS) ${BDTLIB}Tree.cxx -o $@

${XMLLIB}tinyxml2.o: ${XMLLIB}tinyxml2.cpp
	$(CC) $(CFLAGS) ${XMLLIB}tinyxml2.cpp -o $@

clean:
	rm ${BDTLIB}*.o
	rm ${XMLLIB}*.o
	rm libAutocategorizer.so

