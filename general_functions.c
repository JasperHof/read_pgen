#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

///////////////////// reading varints //////////////////////

int read_varint(FILE *ptr){
    
    int continue_varint = 1, k = 0, varint = 0;
    uint8_t byte;

    while (continue_varint == 1){
        fread(&byte, sizeof(uint8_t), 1, ptr);

        uint8_t lower_7_bits = byte & 0x7F; // 0x7F is the bitmask with only the lower 7 bits set
        continue_varint = byte & 0x80;
        varint = varint * (1 << 7) + lower_7_bits;
    }

    return varint;
}