CXX := g++
CXXFLAGS := -O3 -std=c++17 -pthread
LDFLAGS := -lm

nanore: src/nanore.cpp
	$(CXX) $(CXXFLAGS) -o $@ $< $(LDFLAGS)

clean:
	rm -f nanore

.PHONY: clean
