all: poa htslib/lib/libhts.a

htslib/lib/libhts.a:
	cd htslib && autoheader && autoconf && ./configure --disable-s3 --disable-lzma --disable-bz2 --prefix=$(PWD)/htslib/ && make -j 4 && make install

OPTS=-O0
poa: poa.cpp graph.cpp naiveSDP.cpp interval.cpp alignedRead.cpp htslib/lib/libhts.a
	g++ $(OPTS) $^ -o $@ -g -I htslib/ -L htslib/lib -lhts -Wl,-rpath,$(PWD)/htslib/lib -lz -lcurl -lpthread
#poa: poa.cpp graph.cpp naiveSDP.cpp interval.cpp alignedRead.cpp htslib/lib/libhts.a
#	g++ $(OPTS) $^ -o $@ -g -I htslib/ -L htslib/lib -lhts -Wl,-rpath,$(PWD)/htslib/lib -lz -lcurl -lpthread

clean: 
	rm poa
	echo Clean done
