#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

int* find_location(char filename[], int chr, int pos_start, int pos_end);
uint8_t* find_types(FILE *ptr, int start_pos_types, int num_snps, int bits_per_record_type, int bytes_per_record_length);
uint8_t* find_lengths(FILE *ptr, int start_pos_types, int num_snps, int bits_per_record_type, int bytes_per_record_length);
int* find_snp_records(uint8_t* lengths, int start_pos_variant_records[], int num_snps);