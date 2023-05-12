CXX=g++
CFLAGS=-Wall -Wextra -lm

all: encode_graph kmers_ascii kmers_encoded kmers_ascii_bit_operation

test: all
	@echo Comparing map hashing against bit operation hashing...
	@echo -n "Map "
	@./kmers_ascii | grep "Hashing took"
	@echo -n "Bit Operation "
	@./kmers_ascii_bit_operation | grep "Hashing took"
	@echo .
	@echo "Comparing result hashing against graph hashing + k-mer construction..."
	@echo -n "Result Collection and "
	@./kmers_ascii | grep "Hashing took"
	@echo -n "Graph "
	@./encode_graph | grep "Hashing took"
	@echo -n "Pre-Hashed K-mer "
	@./kmers_encoded | grep "Collection took"

test-bit-hash: all
	@echo Comparing map hashing against bit operation hashing...
	@echo -n "Map "
	@./kmers_ascii | grep "Hashing took"
	@echo -n "Bit Operation "
	@./kmers_ascii_bit_operation | grep "Hashing took"

test-graph-hash: all
	@echo "Comparing result hashing against graph hashing + k-mer construction..."
	@echo -n "Result Collection and "
	@./kmers_ascii | grep "Hashing took"
	@echo -n "Graph "
	@./encode_graph | grep "Hashing took"
	@echo -n "Pre-Hashed K-mer "
	@./kmers_encoded | grep "Collection took"

encode_graph: encode_graph.cpp hashing.cpp hashing.hpp
	$(CXX) $(CFLAGS) -o $@ $< hashing.cpp hashing.hpp

kmers_ascii: kmers_ascii.cpp hashing.cpp hashing.hpp
	$(CXX) $(CFLAGS) -o $@ $< hashing.cpp hashing.hpp

kmers_encoded: kmers_encoded.cpp hashing.cpp hashing.hpp
	$(CXX) $(CFLAGS) -o $@ $< hashing.cpp hashing.hpp

kmers_ascii_bit_operation: kmers_ascii_bit_operation.cpp hashing.cpp hashing.hpp
	$(CXX) $(CFLAGS) -o $@ $< hashing.cpp hashing.hpp
