Build

cmake .

make -j4

Feature:
- Delta encodes similar chunks (4KB ~ 64KB) in fine-grained deduplication to generate delta chunks.
- Faster than Xdelta, Zdelta, Ddelta, and Edelta.
- Only removes inter-redundancy between similar chunks. Thus, it's necessary to use general compression (e.g., ZSTD) to further compress the delta chunks.
