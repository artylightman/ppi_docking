CXX = g++ -Wall -Wextra -O3 -std=c++1y 
CXXVORO = g++ -std=c++1y -g -O3
LUAFLAG = -DUSELUA

all:  LINK_luafree

GENERIC:  LINK 

obj/voro.o: lib/voro++/src/voro++.*
	$(CXXVORO) -c -o obj/voro.o lib/voro++/src/voro++.cc

obj/fileloader.o: src/fileloader.*
	$(CXX) -c -o obj/fileloader.o src/fileloader.cpp

obj/pointpattern.o: src/pointpattern.*
	$(CXX) -c -o obj/pointpattern.o src/pointpattern.cpp

obj/testduplicationremover.o: src/pointpattern.* src/unittests/testduplicationremover.cpp
	$(CXX) -c -o obj/testduplicationremover.o src/unittests/testduplicationremover.cpp

obj/testpointpattern.o: src/pointpattern.* src/unittests/testpointpattern.cpp
	$(CXX) -c -o obj/testpointpattern.o src/unittests/testpointpattern.cpp

obj/main_luafree.o: src/main.cpp  src/triangle.hpp src/duplicationremover.hpp src/writerpoly.hpp src/writeroff.hpp src/IWriter.hpp src/postprocessing.hpp src/parsexyz.hpp src/parseellipsoids.hpp src/parsesphcyl.hpp  src/GenericMatrix.h src/parsexyzr.hpp src/parsetetra.hpp src/output.hpp src/colorTable.hpp src/parsetetra_blunt.hpp src/tetrahedra.hpp src/cmdlparser.hpp
	mkdir -p obj
	mkdir -p bin
	$(CXX) -c -o obj/main_luafree.o src/main.cpp 

obj/main.o: src/main.cpp src/triangle.hpp src/duplicationremover.hpp src/writerpoly.hpp src/writeroff.hpp src/IWriter.hpp src/postprocessing.hpp src/parsexyz.hpp src/parsexyzr.hpp src/parseellipsoids.hpp src/parsesphcyl.hpp src/GenericMatrix.h src/parsetetra.hpp  src/output.hpp src/colorTable.hpp src/parsetetra_blunt.hpp src/tetrahedra.hpp src/cmdlparser.hpp
	mkdir -p obj
	mkdir -p bin
	$(CXX) -c -o obj/main.o src/main.cpp -I/ysm-gpfs/apps/software/Lua/5.3.4-GCCcore-6.4.0/include $(LUAFLAG)

LINK_luafree: obj/main_luafree.o obj/voro.o obj/fileloader.o obj/pointpattern.o
	$(CXX) obj/main_luafree.o obj/voro.o obj/fileloader.o obj/pointpattern.o -o bin/pomelo


LINK: obj/main.o obj/voro.o obj/fileloader.o obj/pointpattern.o
	$(CXX) obj/main.o obj/voro.o obj/fileloader.o obj/pointpattern.o -o bin/pomelo -I/ysm-gpfs/apps/software/Lua/5.3.4-GCCcore-6.4.0/include -L/ysm-gpfs/apps/software/Lua/5.3.4-GCCcore-6.4.0/lib/ -ldl -llua $(LUAFLAG) 


test: obj/pointpattern.o src/pointpattern.hpp src/pointpattern.cpp  obj/testpointpattern.o obj/testduplicationremover.o
	$(CXX) obj/pointpattern.o obj/testpointpattern.o -o bin/testpointpattern
	$(CXX) obj/pointpattern.o obj/testduplicationremover.o -o bin/testduplicationremover
	$(CXX) src/unittests/testcmdlparser.cpp -o bin/testcmdlparser
	$(CXX) src/unittests/testsplitstring.cpp -o bin/testsplitstring
	$(CXX) src/unittests/testtriangle.cpp obj/pointpattern.o -o bin/testtriangle
	bin/testpointpattern
	bin/testduplicationremover
	bin/testcmdlparser
	bin/testsplitstring
	bin/testtriangle

clean:
	rm obj/*
	rm bin/pomelo

