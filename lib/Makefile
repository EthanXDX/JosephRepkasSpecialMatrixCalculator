SOURCES=$(src/)
HEADERS=$(src/headers/)

UNAME_S := $(shell uname -s)

ifeq ($(UNAME_S),Darwin)
	CFLAGS :=  -std=c++20 -arch x86_64
	LDFLAGS := -Wl,-v -lginac -lcln -lgmp -static-libstdc++
else
	CFLAGS :=  -std=c++20
	LDFLAGS := -Wl,-Bstatic -lginac -lcln -lgmp -static-libstdc++ # this is more portable, but idk if it will work for all installs
endif

#Makefile Targets
out/%.o: src/%.cpp
	$(CXX) $(CFLAGS) -c -o $@ $^ $(LDFLAGS)

out/library.a: out/utils.o out/lin_alg.o out/lie_algebra.o out/compat.o out/interface.o
# $(CXX) $(CFLAGS) -shared -o $@ $^ # $(LDFLAGS)
	ar rcs $@ $^

out/test.o: src/test.cpp out/library.a
	$(CXX) $(CFLAGS) -o $@ $^ $(LDFLAGS)

out/example.o: src/example.cpp out/library.a #out/library.a
	$(CXX) $(CFLAGS) -o $@ $^ $(LDFLAGS)
#mv $@ out/$@

mkdirOut:
	mkdir -p out

.PHONY: all test clear
all: mkdirOut out/utils.o out/lin_alg.o out/lie_algebra.o out/library.a
test: out/test.o
example: out/example.o
clear:
	rm out/*
#LieAlgebra.a: LieAlgebra.o
#	$(AR) $(ARFLAGS) $@ $^
# If need static, change .so to .a and remove -shared flag
