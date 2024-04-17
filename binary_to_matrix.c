#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include "general_functions.h"
#include <math.h>

///////////////////////////////////////////////////////
// Before we start reading records, make code for interpreting difflists
///////////////////////////////////////////////////////

int read_difflist(FILE *ptr, int *snp_sample_values, int *snp_loc, uint8_t *types, uint8_t *lengths, int *snp_records_pos, int num_samples, int num_snps, int snp_index){
    int varint = read_varint(ptr);
    printf("Varint: %i\n", varint);

    if(varint == 0){
        return 0;
    }

    int n_groups = (varint - 1) / 64 + 1;
    int bytes_for_samples = (log2(num_samples) / 8) + 1;
    int sample_ids_start[n_groups];

    printf("Number of groups and bytes: %i, %i\n", n_groups, bytes_for_samples);

    for (int j = 0; j < n_groups; j++){
        fread(&sample_ids_start[j], bytes_for_samples, 1, ptr);
        printf("Sample ID start: %i\n", sample_ids_start[j]);
    }

    // Now read G-1 bytes, where byte j is (byte size of the j-th group) - 63.
    int group_byte_sizes[n_groups-1];

    for (int j = 1; j < n_groups; j++){
        fread(&group_byte_sizes[j-1], sizeof(int), 1, ptr);
        group_byte_sizes[j-1] = group_byte_sizes[j-1] + 63;
    }

    // If the 2bit genotype values are stored in the difflist, they come now in ceiling(L/4) bytes
    int num_snp_bytes = (varint + 3)/4;
    printf("Number of SNP bytes: %i\n", num_snp_bytes);

    if(snp_sample_values == NULL){
        printf("Memory allocation failed\n");
        return 0;
    }

    snp_sample_values[0] = varint;  // First value is the number of elements
    // First read the SNP values
    for(int i = 0; i < num_snp_bytes; i++){
        uint8_t diff_byte;
        fread(&diff_byte, sizeof(uint8_t), 1, ptr);

        uint8_t bottom_bits = 0x03; // Mask to extract 2 bits
        uint8_t snp1 = diff_byte & bottom_bits;
        uint8_t snp2 = (diff_byte >> 2) & bottom_bits;
        uint8_t snp3 = (diff_byte >> 4) & bottom_bits;
        uint8_t snp4 = (diff_byte >> 6) & bottom_bits;

        snp_sample_values[i*4+1] = snp1;
        snp_sample_values[i*4+2] = snp2;
        snp_sample_values[i*4+3] = snp3;
        snp_sample_values[i*4+4] = snp4;
    }

    // !!! Are the bit sizes correct?
    for(int i = 0; i < num_snp_bytes; i++){
        snp_sample_values[i*64 + varint + 1] = sample_ids_start[i];

        for (int j = 1; j < 64 && i*64 + j < varint; j++){
            int diff = read_varint(ptr);
            snp_sample_values[i*64+j + varint + 1] = snp_sample_values[i*64+j-1 + varint + 1] + diff;
        }
    }

    // printf("SNP and sample values: %i, %i, %i, %i, %i, %i\n", snp_sample_values[0], snp_sample_values[1], snp_sample_values[2], snp_sample_values[3], snp_sample_values[4], snp_sample_values[5]);

    // make a pointer of length 2*varint, with the pairs (sample, genotype)?
    return(1);
}


///////////////////////////////////////////////////////
// First option we consider: type = 0, no compression
///////////////////////////////////////////////////////

void read_type0_record(FILE *ptr, int *snp_data, int *snp_loc, uint8_t *types, uint8_t *lengths, int *snp_records_pos, int num_samples, int num_snps, int snp_index){

    fseek(ptr, snp_records_pos[snp_index], SEEK_SET);
    // two-bit PLINK encoding: 0 = HOM-REF; 1 = HET; 2 = double-ALT; 3 = missing

    for(int i = 0; i < lengths[snp_index]; i++){
            uint8_t byte;
            fread(&byte, sizeof(uint8_t), 1, ptr);

            uint8_t mask = 0x03; // Mask to extract 2 bits
            int piece1 = byte & mask;
            int piece2 = (byte >> 2) & mask;
            int piece3 = (byte >> 4) & mask;
            int piece4 = (byte >> 6) & mask;

            snp_data[i*4] = piece1;
            snp_data[i*4+1] = piece2;
            snp_data[i*4+2] = piece3;
            snp_data[i*4+3] = piece4;
        }
}

///////////////////////////////////////////////////////
// Read SNPs of type = 1; '1-bit representation'
///////////////////////////////////////////////////////

void read_type1_record(FILE *ptr, int *snp_data, int *snp_loc, uint8_t *types, uint8_t *lengths, int *snp_records_pos, int num_samples, int num_snps, int snp_index){

    fseek(ptr, snp_records_pos[snp_index], SEEK_SET);
    // two-bit PLINK encoding: 0 = HOM-REF; 1 = HET; 2 = double-ALT; 3 = missing

    int common, less_common;
    uint8_t byte;

    fread(&byte, sizeof(uint8_t), 1, ptr);

    // Define common / less common categories
    if(byte == 1 | byte == 2 | byte == 3){
        less_common = 0;
    } else if (byte == 5 | byte == 6) {
        less_common = 1;
    } else {less_common = 2;}

    if(byte == 1){
        common = 1;
    } else if (byte == 2 | byte == 5) {
        common = 2;
    } else {common = 3;}

    int bytes_for_snpinfo = (num_samples + 7)/8;

    for(int j = 0; j < bytes_for_snpinfo; j++){
        fread(&byte, sizeof(uint8_t), 1, ptr);

        for (int i = 0; i < 8; i++) {
            uint8_t mask = 1 << i;
            snp_data[j*8+i] = (byte & mask) ? common : less_common;
        }
    }

    // NEED TO WRITE THIS CODE: EXTRACT DIFFLIST AND REMAINING GENOTYPES
    //
    //
    // 

}

///////////////////////////////////////////////////////
// Read SNPs of type = 6; Difflist with samples not 2
///////////////////////////////////////////////////////

void read_type6_record(FILE *ptr, int *snp_data, int *snp_loc, uint8_t *types, uint8_t *lengths, int *snp_records_pos, int num_samples, int num_snps, int snp_index){
    
    printf("Reading num samps, pos; %i, %i; ", num_samples, snp_records_pos[snp_index]);
    
    fseek(ptr, snp_records_pos[snp_index], SEEK_SET);

    // compile with gdp (install)
    // valgrind --leak-check=full -g ./a.out

    printf("Allocating difflist ");
    int *snp_sample_values = malloc(1000 * sizeof(int)); // First SNP values, then sample values, 1000 should be enough to test
    printf("Difflist allocated \n");
    // read_difflist(ptr, snp_sample_values, snp_loc, types, lengths, snp_records_pos, num_samples, num_snps, snp_index);

    if(snp_sample_values == NULL){
        printf("Memory allocation failed \n");
    } 

    printf("Allocating SNP values; ");
    for (int i = 0; i < num_samples; i++){
        snp_data[i] = 2;
    }

    if (snp_sample_values == NULL){
        printf("Snp sample values NOT in difflist\n");
    } else {
        printf("Snp sample values in difflist\n");
    }

    printf("Free SNPsamp values ");
    free(snp_sample_values);
    snp_sample_values = NULL;
    printf("SNPsamp values freed ");
}


// write a function that converts the binary file to a matrix

int bin_to_mat(FILE *ptr, int *snp_loc, uint8_t *types, uint8_t *lengths, int *snp_records_pos, int num_samples, int num_snps){
    int n_blocks = (num_snps-1)/65536 + 1;      // Number of blocks
    int res = (num_snps-1) % 65536 + 1;         // Number of snps in the last block
    int index, count = 0;

    for (int i = 0; i < 10000; i++) {                 // This is a bit ugly, but sizeof(types) / sizeof(uint8_t) was not working
        if (snp_loc[i] != 0) {
            count++;
        }
    }

    int *snp_data;

    // for (int i = 0; i < count; i++) {
    for (int i = 0; i < 30; i++) {   
        index = snp_loc[i];        // Index of the SNP (i.e., row in the .bim file to be analysed)
        
        printf("SNP, Index: %i, %i\n", i, index);

        if(types[i] == 0){
            printf("allocate type 0 ");
            snp_data = malloc(num_samples * sizeof(int)); // Allocate memory
            printf("it is allocated");

            if(snp_data == NULL){
                printf("Memory allocation failed\n");
            }

            printf(" Snp data of SNP %i:\n", i);
            read_type0_record(ptr, snp_data, snp_loc, types, lengths, snp_records_pos, num_samples, num_snps, index);

            for (int j = 0; j < num_samples; j++){
                printf("%i ", snp_data[j]);
            }

            printf("free mem ");
            free(snp_data);
            snp_data = NULL;
            printf("mem free ");
            
            printf("\n");
        }
        if(types[i] == 1){
            printf("allocate type 1 ");
            snp_data = malloc(num_samples * sizeof(int)); // Allocate memory
            printf("it is allocated");

            if(snp_data == NULL){
                printf("Memory allocation failed\n");
            }

            printf("Snp data of SNP %i:\n", i);
            read_type1_record(ptr, snp_data, snp_loc, types, lengths, snp_records_pos, num_samples, num_snps, index);
            
            for (int j = 0; j < num_samples; j++){
                 printf("%i ", snp_data[j]);
            }
    
            printf("free mem ");
            free(snp_data);
            snp_data = NULL;
            printf("mem free ");

            printf("\n");
        }
        if(types[i] == 6){
            printf("allocate type 6 ");
            snp_data = malloc(num_samples * sizeof(int)); // Allocate memory
            printf("it is allocated");

            if(snp_data == NULL){
                printf("Memory allocation failed\n");
            }

            printf("Snp data of SNP %i:\n", i);
            read_type6_record(ptr, snp_data, snp_loc, types, lengths, snp_records_pos, num_samples, num_snps, index);
            
            for (int j = 0; j < num_samples; j++){
                 printf("%i ", snp_data[j]);
            }

            printf("free mem ");
            free(snp_data);
            snp_data = NULL;
            printf("mem free ");

            printf("\n");
        }

        
    }
    printf("Test\n");
    printf("Count: %i\n", count);
}