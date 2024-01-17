Build this program:
1. cmake .
2. make -j4

Feature:
1. Delta encodes similar chunks (4KB ~ 64KB) in fine-grained deduplication to generate delta chunks.
2. Faster than Xdelta, Zdelta, Ddelta, and Edelta.
3. Only removes inter-redundancy between similar chunks. Thus, it's necessary to use general compression (e.g., ZSTD) to further compress the delta chunks.
