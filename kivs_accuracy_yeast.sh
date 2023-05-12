mkdir -p "kivs_accuracy/yeast/"

snakemake --use-conda --config --cores 1 --forceall test_yeast_full > "kivs_accuracy/yeast/kage_index_kage_signature.txt"
snakemake --use-conda --config use_kivs=True kivs_align_windows=True kivs_minimize_overlaps=True --cores 1 --forceall test_yeast_full > "kivs_accuracy/yeast/kage_index_kivs_align_minoverlap_signature.txt"
snakemake --use-conda --config use_kivs=True kivs_align_windows=True kivs_minimize_overlaps=False --cores 1 --forceall test_yeast_full > "kivs_accuracy/yeast/kage_index_kivs_align_signature.txt"
snakemake --use-conda --config use_kivs=True kivs_align_windows=False kivs_minimize_overlaps=True --cores 1 --forceall test_yeast_full > "kivs_accuracy/yeast/kage_index_kivs_minoverlap_signature.txt"
snakemake --use-conda --config use_kivs=True kivs_align_windows=False kivs_minimize_overlaps=False --cores 1 --forceall test_yeast_full > "kivs_accuracy/yeast/kage_index_kivs_signature.txt"
snakemake --use-conda --config use_kivs=True kivs_kmer_index=True kivs_align_windows=True kivs_minimize_overlaps=True --cores 1 --forceall test_yeast_full > "kivs_accuracy/yeast/kivs_index_kivs_align_minoverlap_signature.txt"
snakemake --use-conda --config use_kivs=True kivs_kmer_index=True kivs_align_windows=True kivs_minimize_overlaps=False --cores 1 --forceall test_yeast_full > "kivs_accuracy/yeast/kivs_index_kivs_align_signature.txt"
snakemake --use-conda --config use_kivs=True kivs_kmer_index=True kivs_align_windows=False kivs_minimize_overlaps=True --cores 1 --forceall test_yeast_full > "kivs_accuracy/yeast/kivs_index_kivs_minoverlap_signature.txt"
snakemake --use-conda --config use_kivs=True kivs_kmer_index=True kivs_align_windows=False kivs_minimize_overlaps=False --cores 1 --forceall test_yeast_full > "kivs_accuracy/yeast/kivs_index_kivs_signature.txt"

echo "--------------"
echo "KAGE Index + KAGE Signature + Align Windows + Minimize Overlaps"
grep "^[02]    " "kivs_accuracy/yeast/kage_index_kage_signature.txt" | sed -e "s/0   /Indel/" | sed -e "s/2 /SNP/"

echo "--------------"
echo "KAGE Index + KIVS Signature + Align Windows + Minimize Overlaps"
grep "^[02]    " "kivs_accuracy/yeast/kage_index_kivs_align_minoverlap_signature.txt" | sed -e "s/0   /Indel/" | sed -e "s/2 /SNP/"

echo "--------------"
echo "KAGE Index + KIVS Signature + Align Windows"
grep "^[02]    " "kivs_accuracy/yeast/kage_index_kivs_align_signature.txt" | sed -e "s/0   /Indel/" | sed -e "s/2 /SNP/"

echo "--------------"
echo "KAGE Index + KIVS Signature + Minimize Overlaps"
grep "^[02]    " "kivs_accuracy/yeast/kage_index_kivs_minoverlap_signature.txt" | sed -e "s/0   /Indel/" | sed -e "s/2 /SNP/"

echo "--------------"
echo "KAGE Index + KIVS Signature"
grep "^[02]    " "kivs_accuracy/yeast/kage_index_kivs_signature.txt" | sed -e "s/0   /Indel/" | sed -e "s/2 /SNP/"

echo "--------------"
echo "KIVS Index + KIVS Signature + Align Windows + Minimize Overlaps"
grep "^[02]    " "kivs_accuracy/yeast/kivs_index_kivs_align_minoverlap_signature.txt" | sed -e "s/0   /Indel/" | sed -e "s/2 /SNP/"

echo "--------------"
echo "KIVS Index + KIVS Signature + Align Windows"
grep "^[02]    " "kivs_accuracy/yeast/kivs_index_kivs_align_signature.txt" | sed -e "s/0   /Indel/" | sed -e "s/2 /SNP/"

echo "--------------"
echo "KIVS Index + KIVS Signature + Minimize Overlaps"
grep "^[02]    " "kivs_accuracy/yeast/kivs_index_kivs_minoverlap_signature.txt" | sed -e "s/0   /Indel/" | sed -e "s/2 /SNP/"

echo "--------------"
echo "KIVS Index + KIVS Signature"
grep "^[02]    " "kivs_accuracy/yeast/kivs_index_kivs_signature.txt" | sed -e "s/0   /Indel/" | sed -e "s/2 /SNP/"
