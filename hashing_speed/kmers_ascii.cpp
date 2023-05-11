#include <stdio.h>
#include <iostream>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <cstdlib>
#include <chrono>
#include "hashing.hpp"

int main() {
	uint64_t length = 100000000;
	uint8_t k = 31;
	
	printf("Preparing sequence string of length %lu...\n", length);
	char *bases = (char *) malloc(sizeof(char) * length);
	memset(bases, 0, sizeof(char) * length);

	for (uint64_t tests = 0; tests < 1; tests++) {
		printf("Randomizing...\n");
		for (uint64_t i = 0; i < length; i++) {
			int rng = rand() % 4;
			bases[i] = ("ACTG")[rng];
		}

		printf("Preparing hashing map...\n");
		uint8_t map[256];
		fill_map_by_encoding(map, "ACTG");
	
		printf("Hashing...\n");
		auto start = std::chrono::high_resolution_clock::now();

		uint64_t max = length + 1 - k;
		uint64_t kmer;
		for (uint64_t i = 0; i < max; i++) {
			kmer = hash_min_kmer_by_map(bases + i, k, map);
		}

		auto stop = std::chrono::high_resolution_clock::now();
		auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(stop - start);
		std::cout << "Hashing took " << duration.count() << "ms" << std::endl;
	}

	printf("Cleaning up...\n");
	free(bases);
}
