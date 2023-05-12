
test: test-vg test-odgi test-kivs

test-vg: vg data/yeast.vg
	@echo "Testing runtime for vg"
	@printf "Run\tCPU\tSYS\tMEM (kb)\n"
	@/usr/bin/time -f "1\t%U\t%S\t%M" ./vg kmers -k 31 -t 1 data/yeast.vg > /dev/null
	@/usr/bin/time -f "2\t%U\t%S\t%M" ./vg kmers -k 31 -t 1 data/yeast.vg > /dev/null
	@/usr/bin/time -f "3\t%U\t%S\t%M" ./vg kmers -k 31 -t 1 data/yeast.vg > /dev/null
	@/usr/bin/time -f "4\t%U\t%S\t%M" ./vg kmers -k 31 -t 1 data/yeast.vg > /dev/null
	@/usr/bin/time -f "5\t%U\t%S\t%M" ./vg kmers -k 31 -t 1 data/yeast.vg > /dev/null
	@/usr/bin/time -f "6\t%U\t%S\t%M" ./vg kmers -k 31 -t 1 data/yeast.vg > /dev/null
	@/usr/bin/time -f "7\t%U\t%S\t%M" ./vg kmers -k 31 -t 1 data/yeast.vg > /dev/null

test-odgi: odgi data/yeast.og
	@echo "Testing runtime for odgi"
	@printf "Run\tCPU\tSYS\tMEM (kb)\n"
	@/usr/bin/time -f "1\t%U\t%S\t%M" ./odgi kmers -i data/yeast.og -k 31 -t 1 > /dev/null
	@/usr/bin/time -f "2\t%U\t%S\t%M" ./odgi kmers -i data/yeast.og -k 31 -t 1 > /dev/null
	@/usr/bin/time -f "3\t%U\t%S\t%M" ./odgi kmers -i data/yeast.og -k 31 -t 1 > /dev/null
	@/usr/bin/time -f "4\t%U\t%S\t%M" ./odgi kmers -i data/yeast.og -k 31 -t 1 > /dev/null
	@/usr/bin/time -f "5\t%U\t%S\t%M" ./odgi kmers -i data/yeast.og -k 31 -t 1 > /dev/null
	@/usr/bin/time -f "6\t%U\t%S\t%M" ./odgi kmers -i data/yeast.og -k 31 -t 1 > /dev/null
	@/usr/bin/time -f "7\t%U\t%S\t%M" ./odgi kmers -i data/yeast.og -k 31 -t 1 > /dev/null

test-kivs: data/yeast.kivs
	@echo "Testing runtime for KIVS"
	@printf "Run\tCPU\tSYS\tMEM (kb)\n"
	@/usr/bin/time -f "1\t%U\t%S\t%M" python test_kmer_index_speed.py data/yeast.kivs > /dev/null
	@/usr/bin/time -f "2\t%U\t%S\t%M" python test_kmer_index_speed.py data/yeast.kivs > /dev/null
	@/usr/bin/time -f "3\t%U\t%S\t%M" python test_kmer_index_speed.py data/yeast.kivs > /dev/null
	@/usr/bin/time -f "4\t%U\t%S\t%M" python test_kmer_index_speed.py data/yeast.kivs > /dev/null
	@/usr/bin/time -f "5\t%U\t%S\t%M" python test_kmer_index_speed.py data/yeast.kivs > /dev/null
	@/usr/bin/time -f "6\t%U\t%S\t%M" python test_kmer_index_speed.py data/yeast.kivs > /dev/null
	@/usr/bin/time -f "7\t%U\t%S\t%M" python test_kmer_index_speed.py data/yeast.kivs > /dev/null

test-kivs-stdout: data/yeast.kivs
	@echo "Testing runtime for KIVS"
	@printf "Run\tCPU\tSYS\tMEM (kb)\n"
	@/usr/bin/time -f "1\t%U\t%S\t%M" python test_kmer_index_speed.py data/yeast.kivs stdout > /dev/null
	@/usr/bin/time -f "2\t%U\t%S\t%M" python test_kmer_index_speed.py data/yeast.kivs stdout > /dev/null
	@/usr/bin/time -f "3\t%U\t%S\t%M" python test_kmer_index_speed.py data/yeast.kivs stdout > /dev/null
	@/usr/bin/time -f "4\t%U\t%S\t%M" python test_kmer_index_speed.py data/yeast.kivs stdout > /dev/null
	@/usr/bin/time -f "5\t%U\t%S\t%M" python test_kmer_index_speed.py data/yeast.kivs stdout > /dev/null
	@/usr/bin/time -f "6\t%U\t%S\t%M" python test_kmer_index_speed.py data/yeast.kivs stdout > /dev/null
	@/usr/bin/time -f "7\t%U\t%S\t%M" python test_kmer_index_speed.py data/yeast.kivs stdout > /dev/null

test-kivs-full: data/yeast.kivs
	@echo "Testing runtime for KIVS"
	@printf "Run\tCPU\tSYS\tMEM (kb)\n"
	@/usr/bin/time -f "1\t%U\t%S\t%M" python test_kmer_index_speed.py data/yeast.kivs full > /dev/null
	@/usr/bin/time -f "2\t%U\t%S\t%M" python test_kmer_index_speed.py data/yeast.kivs full > /dev/null
	@/usr/bin/time -f "3\t%U\t%S\t%M" python test_kmer_index_speed.py data/yeast.kivs full > /dev/null
	@/usr/bin/time -f "4\t%U\t%S\t%M" python test_kmer_index_speed.py data/yeast.kivs full > /dev/null
	@/usr/bin/time -f "5\t%U\t%S\t%M" python test_kmer_index_speed.py data/yeast.kivs full > /dev/null
	@/usr/bin/time -f "6\t%U\t%S\t%M" python test_kmer_index_speed.py data/yeast.kivs full > /dev/null
	@/usr/bin/time -f "7\t%U\t%S\t%M" python test_kmer_index_speed.py data/yeast.kivs full > /dev/null

data/yeast.vg:
	./vg construct -r data/yeast.fa -v data/yeast.vcf > data/yeast.vg

data/yeast.gfa: data/yeast.vg
	./vg view data/yeast.vg > data/yeast.gfa

data/yeast.og: data/yeast.gfa
	./odgi build -g data/yeast.gfa -o data/yeast.og --optimize

data/yeast.kivs: data/yeast.gfa
	python gfa_to_bcg.py data/yeast.gfa data/yeast.kivs
