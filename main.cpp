#include <iostream>
#include <cstdint>
#include "gdelta.h"
#include "cstring"

#define Chunk (10*1024*10)
#define BufferSize (Chunk + 20*1024)
#include "FileOperator.h"
int main() {

    uint8_t *newBuffer = (uint8_t*) malloc(BufferSize);
    uint8_t *baseBuffer = (uint8_t*) malloc(BufferSize);
    uint8_t *deltaBuffer = (uint8_t*) malloc(BufferSize);
    uint8_t *outBuffer = (uint8_t*) malloc(BufferSize);

    for(int i = 0; i < Chunk; i++)
    {
        if(i % 488 != 0)
        {
            newBuffer[i] = baseBuffer[i] = random();
        }
        else
        {
            newBuffer[i] = 2;
            baseBuffer[i] = 3;
        }
    }

    uint32_t dSize = 0;



    gencode(newBuffer, Chunk,
            baseBuffer, Chunk,
            &deltaBuffer, &dSize);

    uint32_t outSize = 0;

    gdecode(deltaBuffer, dSize,
            baseBuffer, Chunk,
            &outBuffer, &outSize);

    if(outSize != Chunk || memcmp(newBuffer, outBuffer, Chunk) != 0)
        std::cout<< "decode error" <<std::endl;

    std::cout<< ">>>>>>>>Success, Delta size: " << dSize << "<<<<<<<<" <<std::endl;
    std::cout<< ">>>>>>>>Success, new size: " << Chunk << "<<<<<<<<" <<std::endl;
    std::cout<< ">>>>>>>>Success, Delta ratio: " << (Chunk/(double)dSize) << "<<<<<<<<" <<std::endl;

    return 0;
}
