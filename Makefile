# # Might need to play around with the Windows flags
# ifeq ($(OS),Windows_NT)
# 	CFLAGS := -std=c++17
# 	LDFLAGS := -Wl,-Bstatic
# 	STATIC_LIBFILE := bin/libjoemat.lib
# 	SHARED_LIBFILE := bin/libjoemat.dll
# else
# 	UNAME_S := $(shell uname -s)
# # Need to compile for x86_64 because MATLAB runs on it ...
# # Otherwise, MATLAB says that it's not compatible because
# # it's a 32-bit library (which is at least half right, haha).
# 	ifeq ($(UNAME_S),Darwin)
# 		CFLAGS := -std=c++17 -arch x86_64
# 		LDFLAGS := -Wl,-v
# 		STATIC_LIBFILE := bin/libjoemat.a
# 		SHARED_LIBFILE := bin/libjoemat.dylib
# 	else
# 		CFLAGS := -std=c++17
# 		LDFLAGS := -Wl,-Bstatic
# 		STATIC_LIBFILE := bin/libjoemat.a
# 		SHARED_LIBFILE := bin/libjoemat.so
# 	endif
# endif


# # Note that macOS's `ld` has no -Bstatic option,
# # so the only way to force a static link is by
# # removing all the shared versions of the libraries,
# # wherever they might appear in the PATH
# $(SHARED_LIBFILE): bin/shared/interface.o bin/shared/compat.o bin/shared/libjoemat.a
# 	$(CXX) $(CFLAGS) -shared bin/shared/interface.o bin/shared/compat.o -o $@ $(LDFLAGS) -Lbin/shared -ljoemat -lginac -lcln -lgmp

# bin/shared/%.o: src/%.cpp
# 	$(CXX) $(CFLAGS) -fPIC -c -o $@ $^ $(LDFLAGS)

# bin/shared/libjoemat.a:
# 	make -C lib/joemat
# 	cp lib/joemat/out/library.a bin/shared/libjoemat.a

# bin/testshared.o: test/test.cpp
# 	$(CXX) $(CFLAGS) -fPIC -c -o $@ $^ $(LDFLAGS)

# bin/testshared: bin/testshared.o
# 	$(CXX) $(CFLAGS) -o $@ -Lbin -ljoemat $^

# $(STATIC_LIBFILE): bin/static/interface.o bin/static/library.a
# 	ar rcs $@ $^

# bin/static/interface.o: src/interface.cpp
# 	$(CXX) $(CFLAGS) -c -o $@ $^ $(LDFLAGS)

# build_library:
# 	make all -C lib/joemat
# 	cp lib/joemat/out/*.o bin/static
# 	cp lib/joemat/out/library.a bin/static/

# bin/teststatic.o: test/test.cpp
# 	$(CXX) $(CFLAGS) -c -o $@ $^ $(LDFLAGS)

# bin/teststatic: bin/teststatic.o
# 	$(CXX) $(CFLAGS) -o $@ $^ $(LDFLAGS) -Lbin -ljoemat -lginac -lcln -lgmp 

# bin/joemat.o: src/main.cpp $(STATIC_LIBFILE)
# 	$(CXX) $(CFLAGS) -o $@ $^ $(LDFLAGS) -lginac -lcln -lgmp


# static: $(STATIC_LIBFILE)
# shared: $(SHARED_LIBFILE)
# teststatic: $(STATIC_LIBFILE) bin/teststatic
# testshared: $(SHARED_LIBFILE) bin/testshared
# all: $(STATIC_LIBFILE) $(SHARED_LIBFILE) teststatic testshared

# exec: build_library bin/joemat.o

# clean:
# 	rm -r bin/*
# 	mkdir -p bin/static
# 	mkdir -p bin/shared
# 	rm lib/joemat/out/*


CFLAGS :=  -std=c++17
LDFLAGS := -Wl,-Bstatic -lginac -lcln -lgmp -static-libstdc++

bin/joemat.o: src/main.cpp bin/library.a
	$(CXX) $(CFLAGS) -o $@ $^ $(LDFLAGS)

build_library:
	make all -C lib
	echo "library compiled"
	cp lib/out/library.a bin/library.a

mkdirBin:
	mkdir -p bin

exec: mkdirBin build_library bin/joemat.o

clean:
	rm -r bin/*
	rm lib/out/*