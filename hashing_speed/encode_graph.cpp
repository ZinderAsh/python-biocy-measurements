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
	
	printf("Preparing sequence string of length %lu...\n", length);
	char *bases = (char *) malloc(sizeof(char) * length);
	memset(bases, 0, sizeof(char) * length);

	printf("Preparing encoded array...\n");
	uint64_t encoded_length = (length + 31) / 32;
	uint64_t *encoded = (uint64_t *) malloc(sizeof(uint64_t) * encoded_length);
	memset(encoded, 0, sizeof(uint64_t) * encoded_length);

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

		for (uint64_t i = 0; i < encoded_length - 1; i++) {
			encoded[i] = hash_max_kmer_by_map(bases + i * 32, 32, map);
		}
		encoded[encoded_length - 1] = hash_max_kmer_by_map(bases + (encoded_length - 1) * 32, length % 32, map);

		auto stop = std::chrono::high_resolution_clock::now();
		auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(stop - start);
		std::cout << "Hashing took " << duration.count() << "ms" << std::endl;
	}

	printf("Cleaning up...\n");
	free(bases);
	free(encoded);
}
