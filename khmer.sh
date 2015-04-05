python ~ubuntu/khmer/scripts/load-into-counting.py --ksize 17 --n_tables 4 --min-tablesize 1e9 --report-total-kmers SRR172903.fastq.ht.k17 SRR172903.fastq
python ~ubuntu/khmer/sandbox/hi-lo-abundance-by-position.py SRR172903.fastq.ht.k17 SRR172903.fastq
for x in {17..31..2}; do python unique-kmers.py -R unique_count2 -k $x SRR172903.fastq; done
