#include <cstdlib>
#include <cstring>
#include <cstdint>
#include <algorithm>
using namespace std;

#include "gdelta.h"
#include "gear_matrix.h"
//#include "jemalloc/jemalloc.h"

#pragma pack(push, 1)
/*
 * ABI:
 *
 * VarInt<N>: more 1 | pval N [more| VarInt<7>]
 * DeltaHead: flag 1 | VarInt<6>
 * DeltaUnit: DeltaHead [DeltaHead.flag| VarInt<7>]
 *
 * VarInt <- Val, Offset = Val | VarInt[i].pval << Offset, Offset + VarInt[i]::N
 */
template<uint8_t FLAGLEN>
struct _DeltaHead {
    uint8_t flag: FLAGLEN;
    uint8_t more: 1;
    uint8_t length: (7 - FLAGLEN);
    const static uint8_t lenbits = FLAGLEN;
//    const static uint8_t lenbits = (7 - FLAGLEN);
};

typedef _DeltaHead<1> DeltaHeadUnit;

typedef struct _VarIntPart {
    uint8_t more: 1;
    uint8_t subint: 7;
    const static uint8_t lenbits = 7;
} VarIntPart;

#pragma pack(pop)

typedef struct {
    uint8_t flag;
    uint64_t length;
    uint64_t offset;
} DeltaUnitMem;

// DeltaUnit/FlaggedVarInt: flag: 1, more: 1, len: 6
// VarInt: more: 1, len: 7
static_assert(sizeof(DeltaHeadUnit) == 1, "Expected DeltaHeads to be 1 byte");
static_assert(sizeof(VarIntPart) == 1, "Expected VarInt to be 1 byte");


typedef struct {
    uint8_t *buf;
    uint64_t cursor;
    uint64_t length;
} BufferStreamDescriptor;

void ensure_stream_length(BufferStreamDescriptor &stream, size_t length) {
    if (length > stream.length) {
        stream.buf = (uint8_t *) realloc(stream.buf, length);
        stream.length = length;
    }
}

template<typename T>
void write_field(BufferStreamDescriptor &buffer, const T &field) {
    ensure_stream_length(buffer, buffer.cursor + sizeof(T));
    memcpy(buffer.buf + buffer.cursor, &field, sizeof(T));
    buffer.cursor += sizeof(T);
    // TODO: check bounds (buffer->length)?
}


template<typename T>
void read_field(BufferStreamDescriptor &buffer, T &field) {
    memcpy(&field, buffer.buf + buffer.cursor, sizeof(T));
    buffer.cursor += sizeof(T);
    // TODO: check bounds (buffer->length)?
}

inline
void stream_into(BufferStreamDescriptor &dest, BufferStreamDescriptor &src, size_t length) {
    ensure_stream_length(dest, dest.cursor + length);
    memcpy(dest.buf + dest.cursor, src.buf + src.cursor, length);
    dest.cursor += length;
    src.cursor += length;
}

inline
void stream_from(BufferStreamDescriptor &dest, const BufferStreamDescriptor &src, size_t src_cursor, size_t length) {
    ensure_stream_length(dest, dest.cursor + length);
    memcpy(dest.buf + dest.cursor, src.buf + src_cursor, length);
    dest.cursor += length;
}

inline
void write_concat_buffer(BufferStreamDescriptor &dest, const BufferStreamDescriptor &src) {
    ensure_stream_length(dest, dest.cursor + src.cursor + 1);
    memcpy(dest.buf + dest.cursor, src.buf, src.cursor);
    dest.cursor += src.cursor;
}

inline
uint64_t read_varint(BufferStreamDescriptor &buffer) {
    VarIntPart vi;
    uint64_t val = 0;
    uint8_t offset = 0;
    do {
        read_field(buffer, vi);
        val |= vi.subint << offset;
        offset += VarIntPart::lenbits;
    } while (vi.more);
    return val;
}

inline
void read_unit(BufferStreamDescriptor &buffer, DeltaUnitMem &unit) {
    DeltaHeadUnit head;
    read_field(buffer, head);

    unit.flag = head.flag;
    unit.length = head.length;
    if (head.more) {
        unit.length = read_varint(buffer) << DeltaHeadUnit::lenbits | unit.length;
    }
    if (head.flag) {
        unit.offset = read_varint(buffer);
    }
#if DEBUG_UNITS
    fprintf(stderr, "Reading unit %d %zu %zu\n", unit.flag, unit.length, unit.offset);
#endif
}

const uint8_t varint_mask = ((1 << VarIntPart::lenbits) - 1);
const uint8_t head_varint_mask = ((1 << DeltaHeadUnit::lenbits) - 1);

void write_varint(BufferStreamDescriptor &buffer, uint64_t val) {
    VarIntPart vi;
    do {
        vi.subint = val & varint_mask;
        val >>= VarIntPart::lenbits;
        if (val == 0) {
            vi.more = 0;
            write_field(buffer, vi);
            break;
        }
        vi.more = 1;
        write_field(buffer, vi);
    } while (1);
}

void write_unit(BufferStreamDescriptor &buffer, const DeltaUnitMem &unit) {
    // TODO: Abort if length 0?
#if DEBUG_UNITS
    fprintf(stderr, "Writing unit %d %zu %zu\n", unit.flag, unit.length, unit.offset);
#endif

    DeltaHeadUnit head = {unit.flag, unit.length > head_varint_mask, (uint8_t) (unit.length & head_varint_mask)};
    write_field(buffer, head);
//  cout<<head_varint_mask<<endl;
    uint64_t remaining_length = unit.length >> DeltaHeadUnit::lenbits;
    if (remaining_length)
        write_varint(buffer, remaining_length);
//  write_varint(buffer, remaining_length);
    if (unit.flag) {
        write_varint(buffer, unit.offset);
    }
}


void GFixSizeChunking(unsigned char *data, int len, int begflag, int begsize,
                      uint32_t *hash_table, int mask, uint64_t hashMask) {
    if (len < WordSize)
        return;

    int i = 0;
    int movebitlength = sizeof(FPTYPE) * 8 / WordSize;
    if (sizeof(FPTYPE) * 8 % WordSize != 0)
        movebitlength++;
    FPTYPE fingerprint = 0;

    /** GEAR **/
    for (; i < WordSize; i++) {
        fingerprint = (fingerprint << (movebitlength)) + GEARmx[data[i]];
    }

    i -= WordSize;
    FPTYPE index = 0;
    int numChunks = len - WordSize + 1;

    int _begsize = begflag ? begsize : 0;
    int indexMoveLength = (sizeof(FPTYPE) * 8 - mask);

    while (i < numChunks) {
        index = (fingerprint) >> indexMoveLength;
        hash_table[index] = i + _begsize;
        fingerprint = (fingerprint << (movebitlength)) + GEARmx[data[i + WordSize]];
        i++;

    }

    return;
}


void GFixSizeChunking2(unsigned char *data, int len, int begflag, int begsize,
                       uint32_t *hash_table, int mask, uint64_t hashMask) {
    if (len < WordSize)
        return;

    int i = 0;
    int movebitlength = sizeof(FPTYPE) * 8 / WordSize;
    if (sizeof(FPTYPE) * 8 % WordSize != 0)
        movebitlength++;
    FPTYPE fingerprint = 0;

    /** GEAR **/
    for (; i < WordSize; i++) {
        fingerprint = (fingerprint << (movebitlength)) + GEARmx[data[i]];
    }

    i -= WordSize;
    FPTYPE index = 0;
    int numChunks = len - WordSize + 1;

    int _begsize = begflag ? begsize : 0;
    int indexMoveLength = (sizeof(FPTYPE) * 8 - mask);

    while (i < numChunks) {
        index = (fingerprint) >> indexMoveLength;
        hash_table[index] = i + _begsize;
        fingerprint = (fingerprint << 2) + Gearmx_l[data[i + WordSize]] + GEARmx[data[i + WordSize + 1]];
        i+=2;
    }
}

void GFixSizeChunking_3(unsigned char *data, int len, int begflag, int begsize,
                      uint32_t *hash_table, int mask, uint64_t hashMask) {
    if (len < WordSize)
        return;

    int i = 0;
    int movebitlength = sizeof(FPTYPE) * 8 / WordSize;
    if (sizeof(FPTYPE) * 8 % WordSize != 0)
        movebitlength++;
    FPTYPE fingerprint = 0;

    /** GEAR **/
    for (; i < WordSize; i++) {
        fingerprint = (fingerprint << (movebitlength)) + GEARmx[data[i]];
    }

    i -= WordSize;
    FPTYPE index = 0;
    int numChunks = len - WordSize + 1;


    int _begsize = begflag ? begsize : 0;
    int indexMoveLength = (sizeof(FPTYPE) * 8 - mask);

    while (i < numChunks) {
        index = (fingerprint) >> indexMoveLength;
        hash_table[index] = i + _begsize;
        fingerprint = (fingerprint << (movebitlength)) + GEARmx[data[i + WordSize]];
        fingerprint = (fingerprint << (movebitlength)) + GEARmx[data[i + WordSize + 1]];
        fingerprint = (fingerprint << (movebitlength)) + GEARmx[data[i + WordSize + 2]];
        i+=3;
    }
}

void GFixSizeChunking_4(unsigned char *data, int len, int begflag, int begsize,
                        uint32_t *hash_table, int mask, uint64_t hashMask) {
    if (len < WordSize)
        return;

    int i = 0;
    int movebitlength = sizeof(FPTYPE) * 8 / WordSize;
    if (sizeof(FPTYPE) * 8 % WordSize != 0)
        movebitlength++;
    FPTYPE fingerprint = 0;

    /** GEAR **/
    for (; i < WordSize; i++) {
        fingerprint = (fingerprint << (movebitlength)) + GEARmx[data[i]];
    }

    i -= WordSize;
    FPTYPE index = 0;
    int numChunks = len - WordSize + 1;

    int flag = 0;
    int _begsize = begflag ? begsize : 0;
    int indexMoveLength = (sizeof(FPTYPE) * 8 - mask);

    while (i < numChunks) {
        index = (fingerprint) >> indexMoveLength;
        hash_table[index] = i + _begsize;
        fingerprint = (fingerprint << (movebitlength)) + GEARmx[data[i + WordSize]];
        fingerprint = (fingerprint << (movebitlength)) + GEARmx[data[i + WordSize + 1]];
        fingerprint = (fingerprint << (movebitlength)) + GEARmx[data[i + WordSize + 2]];
        fingerprint = (fingerprint << (movebitlength)) + GEARmx[data[i + WordSize + 3]];
        i+=4;
    }

}


int gencode(uint8_t *newBuf, uint32_t newSize, uint8_t *baseBuf,
            uint32_t baseSize, uint8_t **deltaBuf, uint32_t *deltaSize) {
#if PRINT_PERF
    struct timespec tf0, tf1;
    clock_gettime(CLOCK_MONOTONIC, &tf0);
#endif

    /* detect the head and tail of one chunk */
    uint32_t beg = 0, end = 0, begSize = 0, endSize = 0;

    if (*deltaBuf == nullptr) {
        *deltaBuf = (uint8_t *) malloc(INIT_BUFFER_SIZE);
    }


    uint8_t *databuf = (uint8_t *) malloc(INIT_BUFFER_SIZE);
    uint8_t *instbuf = (uint8_t *) malloc(INIT_BUFFER_SIZE);


    // Find first difference
    // First in 8 byte blocks and then in 1 byte blocks for speed
    while (begSize + sizeof(uint64_t) <= baseSize &&
           begSize + sizeof(uint64_t) <= newSize &&
           *(uint64_t *) (baseBuf + begSize) == *(uint64_t *) (newBuf + begSize)) {
        begSize += sizeof(uint64_t);
    }

    while (begSize < baseSize &&
           begSize < newSize &&
           baseBuf[begSize] == newBuf[begSize]) {
        begSize++;
    }

    if (begSize > 16)
        beg = 1;
    else
        begSize = 0;

    // Find first difference (from the end)
    while (endSize + sizeof(uint64_t) <= baseSize &&
           endSize + sizeof(uint64_t) <= newSize &&
           *(uint64_t *) (baseBuf + baseSize - endSize - sizeof(uint64_t)) ==
           *(uint64_t *) (newBuf + newSize - endSize - sizeof(uint64_t))) {
        endSize += sizeof(uint64_t);
    }

    while (endSize < baseSize &&
           endSize < newSize &&
           baseBuf[baseSize - endSize - 1] == newBuf[newSize - endSize - 1]) {
        endSize++;
    }

    if (begSize + endSize > newSize)
        endSize = newSize - begSize;

    if (endSize > 16)
        end = 1;
    else
        endSize = 0;
    /* end of detect */

    BufferStreamDescriptor deltaStream = {*deltaBuf, 0, ChunkSize};
    BufferStreamDescriptor instStream = {instbuf, 0, INIT_BUFFER_SIZE}; // Instruction stream
    BufferStreamDescriptor dataStream = {databuf, 0, INIT_BUFFER_SIZE};
    BufferStreamDescriptor newStream = {newBuf, begSize, newSize};
    DeltaUnitMem unit = {}; // In-memory represtation of current working unit

    if (begSize + endSize >= baseSize) { // TODO: test this path
        if (beg) {
            // Data at start is from the original file, write instruction to copy from base
            unit.flag = true;
            unit.offset = 0;
            unit.length = begSize;
            write_unit(instStream, unit);
        }
        if (newSize - begSize - endSize > 0) {
            int32_t litlen = newSize - begSize - endSize;
            unit.flag = false;
            unit.length = litlen;
            write_unit(instStream, unit);
            stream_into(dataStream, newStream, litlen);
        }
        if (end) {
            int32_t matchlen = endSize;
            int32_t offset = baseSize - endSize;
            unit.flag = true;
            unit.offset = offset;
            unit.length = matchlen;
            write_unit(instStream, unit);
        }

        write_varint(deltaStream, instStream.cursor);
        write_concat_buffer(deltaStream, instStream);
        write_concat_buffer(deltaStream, dataStream);

        *deltaSize = deltaStream.cursor;
        *deltaBuf = deltaStream.buf;

#if PRINT_PERF
        clock_gettime(CLOCK_MONOTONIC, &tf1);
        fprintf(stderr, "gencode took: %zdns\n", (tf1.tv_sec - tf0.tv_sec) * 1000000000 + tf1.tv_nsec - tf0.tv_nsec);
#endif
        free(dataStream.buf);
        free(instStream.buf);
        return deltaStream.cursor;
    }

    /* chunk the baseFile */
    int32_t tmp = (baseSize - begSize - endSize) + 10;

    int32_t bit = 0;
    for (bit = 0; tmp; bit++)
        tmp >>= 1;

    uint64_t hashMask = 0XFFFFFFFFFFFFFFFF >> (64 - bit); // mask
    uint32_t handleBytes = begSize;
    FPTYPE fingerprint = 0;
    uint32_t inputPos = begSize;
    uint32_t cursor = 0;
    uint32_t matchNum = 0;
    int32_t moveBitLength = 0;
    uint32_t lastMatchPos = inputPos;
    uint32_t lastLiteral_begin = inputPos;
    uint32_t *hash_table = nullptr;
    uint32_t baseoffset = 0;
    uint32_t hash_size = hashMask + 1;
    bool isFindMatch = false;

    if (beg) {
        // Data at start is from the original file, write instruction to copy from base
        unit.flag = true;
        unit.offset = 0;
        unit.length = begSize;
        write_unit(instStream, unit);
        unit.length = 0; // Mark as written
    }

    uint32_t moveindex = (sizeof(FPTYPE) * 8 - bit);


    isFindMatch = true;
    hash_table = (uint32_t *) malloc(hash_size * sizeof(uint32_t));
    memset(hash_table, 0, sizeof(uint32_t) * hash_size);


#if PRINT_PERF
    struct timespec t0, t1;
    clock_gettime(CLOCK_MONOTONIC, &t0);
#endif


#ifdef BaseSampleRate
    if(BaseSampleRate == 2 && WordSize == 64)
    {
        GFixSizeChunking2(baseBuf + begSize, baseSize - begSize - endSize, beg, begSize, hash_table, bit, hashMask);
    }else if(BaseSampleRate == 3)
    {
        GFixSizeChunking_3(baseBuf + begSize, baseSize - begSize - endSize, beg, begSize, hash_table, bit, hashMask);
    }else if(BaseSampleRate == 4)
    {
        GFixSizeChunking_4(baseBuf + begSize, baseSize - begSize - endSize, beg, begSize, hash_table, bit, hashMask);
    } else
    {
        GFixSizeChunking(baseBuf + begSize, baseSize - begSize - endSize, beg, begSize, hash_table, bit, hashMask);
    }

#else
    GFixSizeChunking(baseBuf + begSize, baseSize - begSize - endSize, beg, begSize, hash_table, bit, hashMask);
#endif


#if PRINT_PERF
    clock_gettime(CLOCK_MONOTONIC, &t1);

    fprintf(stderr, "size:%d\n", baseSize - begSize - endSize);
    fprintf(stderr, "hash size:%d\n", hash_size);
    fprintf(stderr, "rolling hash:%.3fMB/s\n",
            (double)(baseSize - begSize - endSize) / 1024 / 1024 /
                ((t1.tv_sec - t0.tv_sec) * 1000000000 + t1.tv_nsec - t0.tv_nsec) *
                1000000000);
    fprintf(stderr, "rolling hash:%zd\n",
            (t1.tv_sec - t0.tv_sec) * 1000000000 + t1.tv_nsec - t0.tv_nsec);

    clock_gettime(CLOCK_MONOTONIC, &t0);
    fprintf(stderr, "hash table :%zd\n",
            (t0.tv_sec - t1.tv_sec) * 1000000000 + t0.tv_nsec - t1.tv_nsec);
#endif
    /* end of inserting */


    if (sizeof(FPTYPE) * 8 % WordSize == 0)
        moveBitLength = sizeof(FPTYPE) * 8 / WordSize;
    else
        moveBitLength = sizeof(FPTYPE) * 8 / WordSize + 1;


    for (uint32_t i = 0; i < WordSize && i < newSize - endSize - inputPos; i++) {
        fingerprint = (fingerprint << (moveBitLength)) + GEARmx[(newBuf + inputPos)[i]];
    }

    while (inputPos + WordSize <= newSize - endSize) {
        uint32_t length;
        bool matchflag = false;
        cursor = inputPos + WordSize;
        length = WordSize;

        int32_t index = (fingerprint) >> moveindex;

        uint32_t offset = 0;

        baseoffset = hash_table[index];

        if (baseoffset != 0 && memcmp(newBuf + inputPos, baseBuf + baseoffset, length) == 0) {
            matchflag = true;
            offset = baseoffset;
        }

#ifdef ReverseMatch
        if (baseoffset != 0 && !matchflag) {

            uint32_t i = length;
            uint8_t *basepos_end = baseBuf + baseoffset + length - 1;
            uint8_t *inputpos_end = newBuf + inputPos + length - 1;
            for (; i > 0 && (*basepos_end == *inputpos_end); i--, basepos_end--, inputpos_end--) {
            }

            uint32_t matchlen_end = length - i;

            if (matchlen_end > WordSize / 2) {

                int j = 0;
                {
                    uint8_t *p1 = baseBuf + baseoffset + length;
                    uint8_t *p2 = newBuf + inputPos + length;
                    uint32_t *basepos_end1 = (uint32_t *) p1;
                    uint32_t *inputpos_end1 = (uint32_t *) p2;

                    for (; p1 + 3 + j < baseBuf + baseSize &&
                           p2 + 3 + j < newBuf + newSize - endSize &&
                           (*basepos_end1 == *inputpos_end1); j += 4, basepos_end1++, inputpos_end1++) {
                    }

                }
                unit.flag = false;
                unit.length += i;
                stream_from(dataStream, newStream, lastLiteral_begin, unit.length);
                write_unit(instStream, unit);

                unit.flag = true;
                unit.offset = (basepos_end - baseBuf + 1);
                unit.length = matchlen_end + j;

                write_unit(instStream, unit);
                unit.length = 0; // Mark written
                handleBytes += (length + j);
                inputPos += (length + j);
                lastLiteral_begin = inputPos;
                for (int k = 0; k < WordSize && (inputPos + k < newSize - endSize); k++) {
                    fingerprint = (fingerprint << (moveBitLength)) + GEARmx[newBuf[inputPos + k]];
                }

                lastMatchPos = inputPos;
                continue;
            }
        }
#endif


        /* New data match found in hashtable/base data; attempt to create copy instruction*/
        if (matchflag) {
            matchNum++;
            // Check how much is possible to copy
            int32_t j = 0;
#if 1 /* 8-bytes optimization */
            while (offset + length + j + 7 < baseSize - endSize &&
                   cursor + j + 7 < newSize - endSize &&
                   *(uint64_t *) (baseBuf + offset + length + j) == *(uint64_t *) (newBuf + cursor + j)) {
                j += sizeof(uint64_t);
            }
            while (offset + length + j < baseSize - endSize &&
                   cursor + j < newSize - endSize &&
                   baseBuf[offset + length + j] == newBuf[cursor + j]) {
                j++;
            }
#endif
            cursor += j;

            int32_t matchlen = cursor - inputPos;
            handleBytes += cursor - inputPos;
            uint64_t _offset = offset;

            // Check if switching modes Literal -> Copy, and dump instruction if available
            if (!unit.flag && unit.length) {
                /* Detect if end of previous literal could have been a partial copy*/
                uint32_t k = 0;
                while (k + 1 <= offset && k + 1 <= unit.length) {
                    if (baseBuf[offset - (k + 1)] == newBuf[inputPos - (k + 1)])
                        k++;
                    else
                        break;
                }

                if (k > 0) {
                    // Reduce literal by the amount covered by the copy
                    unit.length -= k;
                    // Set up adjusted copy parameters
                    matchlen += k;
                    _offset -= k;
                }

                if (unit.length != 0)
                {
                    write_unit(instStream, unit);
                    stream_from(dataStream, newStream, lastLiteral_begin, unit.length);
                }

                unit.length = 0; // Mark written
            }

            unit.flag = true;
            unit.offset = _offset;
            unit.length = matchlen;
            write_unit(instStream, unit);
            unit.length = 0; // Mark written

            // Update cursor (inputPos) and fingerprint
            for (uint32_t k = cursor; k < cursor + WordSize && cursor + WordSize < newSize - endSize; k++) {
                fingerprint = (fingerprint << (moveBitLength)) + GEARmx[newBuf[k]];
            }

            inputPos = cursor;
            lastMatchPos = cursor;
            lastLiteral_begin = cursor;
        } else { // No match, need to write additional (literal) data
            /*
             * Accumulate length one byte at a time (as literal) in unit while no match is found
             * Pre-emptively write to datastream
             */
            unit.flag = false;
            unit.length += 1;

            handleBytes += 1;
            // Update cursor (inputPos) and fingerprint
            if (inputPos + WordSize < newSize - endSize) {
                fingerprint = (fingerprint << (moveBitLength)) + GEARmx[newBuf[inputPos + WordSize]];
            }
            inputPos++;
#ifdef SkipOn
            int step = ((inputPos - lastMatchPos) >> SkipStep) ;

            if(step <= WordSize)
            {
                for(int i = 0; i < step && (inputPos + WordSize < newSize - endSize); i++,inputPos++)
                {
                    fingerprint = (fingerprint << (moveBitLength)) + GEARmx[newBuf[inputPos + WordSize]];
                    handleBytes += 1;
                    unit.length += 1;
                }
            }
            else
            {
                fingerprint = 0;
                int cursor = inputPos + step;
                int len = 0;
                for(int i = 0; i < WordSize && (cursor + i < newSize - endSize); i++)
                {
                    fingerprint = (fingerprint << (moveBitLength)) + GEARmx[newBuf[cursor + i]];
                    len++;
                }
                int l = min(newSize - endSize, inputPos + step);
                int realStep = l - inputPos;
                handleBytes += realStep;
                unit.length += realStep;
                inputPos += realStep;
            }
#endif
        }
    }


#if PRINT_PERF
    clock_gettime(CLOCK_MONOTONIC, &t1);
    fprintf(stderr, "look up:%zd\n",
            (t1.tv_sec - t0.tv_sec) * 1000000000 + t1.tv_nsec - t0.tv_nsec);
    fprintf(stderr, "look up:%.3fMB/s\n",
            (double)(baseSize - begSize - endSize) / 1024 / 1024 /
                ((t1.tv_sec - t0.tv_sec) * 1000000000 + t1.tv_nsec - t0.tv_nsec) *
                1000000000);
#endif

    endloop:
    // If last unit was unwritten literal, update it to use the rest of the data
    if ((!unit.flag && unit.length) || !isFindMatch) {
        newStream.cursor = lastLiteral_begin;
        stream_into(dataStream, newStream, newSize - endSize - lastLiteral_begin);

        unit.length = (newSize - endSize - lastLiteral_begin);
        write_unit(instStream, unit);
        unit.length = 0;
    } else { // Last unit was Copy, need new instruction
        if (newSize - endSize - handleBytes) {
            newStream.cursor = lastLiteral_begin;
            stream_into(dataStream, newStream, newSize - endSize - lastLiteral_begin);

            unit.flag = false;
            unit.length = newSize - endSize - lastLiteral_begin;
            write_unit(instStream, unit);
            unit.length = 0;
        }
    }

    if (end) {
        int32_t matchLen = endSize;
        int32_t offset = baseSize - endSize;

        unit.flag = true;
        unit.offset = offset;
        unit.length = matchLen;
        write_unit(instStream, unit);
        unit.length = 0;
    }

    deltaStream.cursor = 0;
    write_varint(deltaStream, instStream.cursor);
    write_concat_buffer(deltaStream, instStream);
    write_concat_buffer(deltaStream, dataStream);
    *deltaSize = deltaStream.cursor;
    *deltaBuf = deltaStream.buf;

#if PRINT_PERF
    clock_gettime(CLOCK_MONOTONIC, &tf1);
    fprintf(stderr, "gencode took: %zdns\n", (tf1.tv_sec - tf0.tv_sec) * 1000000000 + tf1.tv_nsec - tf0.tv_nsec);
#endif

    free(dataStream.buf);
    free(instStream.buf);
    if (hash_table)
        free(hash_table);


    return deltaStream.cursor;
}



int gdecode(uint8_t *deltaBuf, uint32_t deltaSize, uint8_t *baseBuf, uint32_t baseSize,
            uint8_t **outBuf, uint32_t *outSize) {

    if (*outBuf == nullptr) {
        *outBuf = (uint8_t *) malloc(INIT_BUFFER_SIZE);
    }

#if PRINT_PERF
    struct timespec tf0, tf1;
    clock_gettime(CLOCK_MONOTONIC, &tf0);
#endif
    BufferStreamDescriptor deltaStream = {deltaBuf, 0, deltaSize}; // Instructions
    const uint64_t instructionLength = read_varint(deltaStream);
    const uint64_t instOffset = deltaStream.cursor;
    BufferStreamDescriptor addDeltaStream = {deltaBuf, deltaStream.cursor + instructionLength, deltaSize};
    BufferStreamDescriptor outStream = {*outBuf, 0, ChunkSize};   // Data out
    BufferStreamDescriptor baseStream = {baseBuf, 0, baseSize}; // Data in
    DeltaUnitMem unit = {};

    while (deltaStream.cursor < instructionLength + instOffset) {
        read_unit(deltaStream, unit);
        if (unit.flag) // Read from original file using offset
            stream_from(outStream, baseStream, unit.offset, unit.length);
        else          // Read from delta file at current cursor
            stream_into(outStream, addDeltaStream, unit.length);
    }

    *outSize = outStream.cursor;
    *outBuf = outStream.buf;
#if PRINT_PERF
    clock_gettime(CLOCK_MONOTONIC, &tf1);
    fprintf(stderr, "gdecode took: %zdns\n", (tf1.tv_sec - tf0.tv_sec) * 1000000000 + tf1.tv_nsec - tf0.tv_nsec);
#endif
    return outStream.cursor;
}
