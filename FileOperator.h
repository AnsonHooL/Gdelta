//
// Created by BorelsetR on 2019/7/17.
//

#ifndef REDUNDANCY_DETECTION_FILEOPERATOR_H
#define REDUNDANCY_DETECTION_FILEOPERATOR_H

#include <sys/stat.h>
#include <unistd.h>
#include <string>
#include <cstring>
#include <cassert>

enum class FileOpenType {
    Read,
    Write,
    ReadWrite,
    Append,
};

struct OutputHead {
    uint64_t baseFid;
    uint64_t fileFid;
};

uint64_t fileCounter = 0;

class FileOperator {
public:
    FileOperator(char *path, FileOpenType fileOpenType) {
        switch (fileOpenType) {
            case FileOpenType::Write :
                file = fopen(path, "wb");
                break;
            case FileOpenType::ReadWrite :
                file = fopen(path, "wb+");
                break;
            case FileOpenType::Append :
                file = fopen(path, "a");
                break;
            case FileOpenType::Read :
            default:
                file = fopen(path, "rb");
                break;
        }
        if (!file) {
            printf("Open: %s\n", path);
            printf("error: %s\n", strerror(errno));
            printf("OpenFile:%lu\n", fileCounter);
            exit(0);
        }
        fileCounter++;
        //printf("open:%d\n", file->_fileno);
    }

    ~FileOperator() {
        if (file != NULL) {
            //printf("close:%d\n", file->_fileno);
            fclose(file);
            fileCounter--;
        }
    }

    uint64_t read(uint8_t *buffer, uint64_t length) {
        return fread(buffer, 1, length, file);
    }

    uint64_t write(uint8_t *buffer, uint64_t length) {
        return fwrite(buffer, 1, length, file);
    }

    int seek(uint64_t offset) {
        return fseek(file, offset, SEEK_SET);
    }

    static uint64_t size(const std::string &path) {
        struct stat statBuffer;
        stat(path.c_str(), &statBuffer);
        return statBuffer.st_size;
    }

    int fsync() {
        return fdatasync(file->_fileno);
    }

    int getFd() {
        return file->_fileno;
    }

private:
    FILE *file;
};


#endif //REDUNDANCY_DETECTION_FILEOPERATOR_H
