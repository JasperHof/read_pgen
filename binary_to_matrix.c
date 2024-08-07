#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include "general_functions.h"
#include <math.h>

///////////////////////////////////////////////////////
// Before we start reading records, make code for interpreting difflists
///////////////////////////////////////////////////////

// Define the custom 24-bit integer type
typedef struct {
    uint8_t bytes[3];
} uint24_t;

int read_difflist(FILE *ptr, double *snp_sample_values, int *snp_loc, int *types, int *lengths, int *snp_records_pos, int num_samples, int num_snps, int snp_index, int varint){
    
    /*
    int varint = read_varint(ptr);

    if(varint == 0){
        snp_sample_values[0] = 0;
        return 0;
    }
    */

    int n_groups = (varint - 1) / 64 + 1;
    int bytes_for_samples = (log2(num_samples) / 8) + 1;   

    if(bytes_for_samples == 1){
        // This variable needs to be dependent on the size of the number of samples!!!
        uint8_t sample_ids_start[n_groups];

        printf("Bytes for samples %i\n", bytes_for_samples);
        printf("Number groups %i\n", n_groups);
        
        for (int j = 0; j < n_groups; j++){
            fread(&sample_ids_start[j], bytes_for_samples, 1, ptr);

            printf("Group start %i\n", sample_ids_start[j]);
        }

        // Now read G-1 bytes, where byte j is (byte size of the j-th group) - 63.
        int group_byte_sizes[n_groups-1];

        for (int j = 1; j < n_groups; j++){
            fread(&group_byte_sizes[j-1], sizeof(uint8_t), 1, ptr);
            group_byte_sizes[j-1] = group_byte_sizes[j-1] + 63;
        }

        // If the 2bit genotype values are stored in the difflist, they come now in ceiling(L/4) bytes
        int num_snp_bytes = (varint + 3)/4;

        int snp_values[4];
        snp_sample_values[0] = (double) varint;  // First value is the number of elements

        // First read the SNP values
        printf("Read SNP values\n");
        for(int i = 0; i < num_snp_bytes; i++){
            uint8_t diff_byte;
            fread(&diff_byte, sizeof(uint8_t), 1, ptr);

            uint8_t bottom_bits = 0x03; // Mask to extract 2 bits
            snp_values[0] = diff_byte & bottom_bits;
            snp_values[1] = (diff_byte >> 2) & bottom_bits;
            snp_values[2] = (diff_byte >> 4) & bottom_bits;
            snp_values[3] = (diff_byte >> 6) & bottom_bits;

            for (int j = 0; j < 4 & i*4 + j < varint; j++){
                snp_sample_values[i*4 + j + 1] = (double)snp_values[j];
            }
        }

        // !!! Are the bit sizes correct now assuming some value?
        for(int i = 0; i < n_groups; i++){
            snp_sample_values[i*64 + varint + 1] = (double)sample_ids_start[i];

            for (int j = 1; j < 64 && i*64 + j < varint; j++){
                int diff = read_varint(ptr);

                snp_sample_values[i*64 + j + varint + 1] = (double)snp_sample_values[i*64 + j - 1 + varint + 1] + diff;
            }
        }

        return(1);
    } else if (bytes_for_samples == 2){
        // This variable needs to be dependent on the size of the number of samples!!!
        uint16_t sample_ids_start[n_groups];

        printf("Bytes for samples %i\n", bytes_for_samples);
        printf("Number groups %i\n", n_groups);
        
        for (int j = 0; j < n_groups; j++){
            fread(&sample_ids_start[j], bytes_for_samples, 1, ptr);

            printf("Group start %i\n", sample_ids_start[j]);
        }

        // Now read G-1 bytes, where byte j is (byte size of the j-th group) - 63.
        // How is this used??
        int group_byte_sizes[n_groups-1];

        for (int j = 1; j < n_groups; j++){
            fread(&group_byte_sizes[j-1], sizeof(uint8_t), 1, ptr);
            group_byte_sizes[j-1] = group_byte_sizes[j-1] + 63;
        }

        // If the 2bit genotype values are stored in the difflist, they come now in ceiling(L/4) bytes
        int num_snp_bytes = (varint + 3)/4;

        int snp_values[4];
        snp_sample_values[0] = (double) varint;  // First value is the number of elements

        // First read the SNP values
        printf("Read SNP values\n");
        for(int i = 0; i < num_snp_bytes; i++){
            uint8_t diff_byte;
            fread(&diff_byte, sizeof(uint8_t), 1, ptr);

            uint8_t bottom_bits = 0x03; // Mask to extract 2 bits
            snp_values[0] = diff_byte & bottom_bits;
            snp_values[1] = (diff_byte >> 2) & bottom_bits;
            snp_values[2] = (diff_byte >> 4) & bottom_bits;
            snp_values[3] = (diff_byte >> 6) & bottom_bits;

            for (int j = 0; j < 4 & i*4 + j < varint; j++){
                snp_sample_values[i*4 + j + 1] = (double)snp_values[j];
            }
        }

        // !!! Are the bit sizes correct now assuming some value?
        for(int i = 0; i < n_groups; i++){
            snp_sample_values[i*64 + varint + 1] = (double)sample_ids_start[i];

            for (int j = 1; j < 64 && i*64 + j < varint; j++){
                int diff = read_varint(ptr);

                snp_sample_values[i*64 + j + varint + 1] = (double)snp_sample_values[i*64 + j - 1 + varint + 1] + diff;
            }
        }

        return(1);
        
    } else if (bytes_for_samples == 3){

        // !!! This should be tested at some point
        // This variable needs to be dependent on the size of the number of samples!!!
        
        uint32_t sample_ids_start[n_groups];

        printf("Bytes for samples %i\n", bytes_for_samples);
        printf("Number groups %i\n", n_groups);
        
        for (int j = 0; j < n_groups; j++){
            fread(&sample_ids_start[j], bytes_for_samples, 1, ptr);

            printf("Group start %i\n", sample_ids_start[j]);
        }

        // Now read G-1 bytes, where byte j is (byte size of the j-th group) - 63.
        int group_byte_sizes[n_groups-1];

        for (int j = 1; j < n_groups; j++){
            fread(&group_byte_sizes[j-1], sizeof(uint8_t), 1, ptr);
            group_byte_sizes[j-1] = group_byte_sizes[j-1] + 63;
        }

        // If the 2bit genotype values are stored in the difflist, they come now in ceiling(L/4) bytes
        int num_snp_bytes = (varint + 3)/4;

        int snp_values[4];
        snp_sample_values[0] = (double) varint;  // First value is the number of elements

        // First read the SNP values
        printf("Read SNP values\n");
        for(int i = 0; i < num_snp_bytes; i++){
            uint8_t diff_byte;
            fread(&diff_byte, sizeof(uint8_t), 1, ptr);

            uint8_t bottom_bits = 0x03; // Mask to extract 2 bits
            snp_values[0] = diff_byte & bottom_bits;
            snp_values[1] = (diff_byte >> 2) & bottom_bits;
            snp_values[2] = (diff_byte >> 4) & bottom_bits;
            snp_values[3] = (diff_byte >> 6) & bottom_bits;

            for (int j = 0; j < 4 & i*4 + j < varint; j++){
                snp_sample_values[i*4 + j + 1] = (double)snp_values[j];
            }
        }

        // !!! Are the bit sizes correct now assuming some value?
        for(int i = 0; i < n_groups; i++){
            snp_sample_values[i*64 + varint + 1] = (double)sample_ids_start[i];

            for (int j = 1; j < 64 && i*64 + j < varint; j++){
                int diff = read_varint(ptr);
                snp_sample_values[i*64 + j + varint + 1] = (double)snp_sample_values[i*64 + j - 1 + varint + 1] + diff;
            }
        }

        return(1);
    }
}


///////////////////////////////////////////////////////
// First option we consider: type = 0, no compression
///////////////////////////////////////////////////////

void read_type0_record(FILE *ptr, double *snp_data, double *snp_data_all, int *snp_loc, int *types, int *lengths, int *snp_records_pos, int num_samples, int num_snps, int snp_index, int entry){

    fseek(ptr, snp_records_pos[snp_index], SEEK_SET);
    // two-bit PLINK encoding: 0 = HOM-REF; 1 = HET; 2 = double-ALT; 3 = missing
    double pieces[4];

    for(int i = 0; i < lengths[snp_index]; i++){
        uint8_t byte;
        fread(&byte, sizeof(uint8_t), 1, ptr);

        uint8_t mask = 0x03; // Mask to extract 2 bits
        pieces[0] = byte & mask;
        pieces[1] = (byte >> 2) & mask;
        pieces[2] = (byte >> 4) & mask;
        pieces[3] = (byte >> 6) & mask;

        for (int j = 0; j < 4 & i*4 + j < num_samples; j++){
            snp_data[i*4+j] = pieces[j];
            snp_data_all[entry * num_samples + i*4+j] = pieces[j];
        }
    }
}

///////////////////////////////////////////////////////
// Read SNPs of type = 1; '1-bit representation'
///////////////////////////////////////////////////////

void read_type1_record(FILE *ptr, double *snp_data, double *snp_data_all, int *snp_loc, int *types, int *lengths, int *snp_records_pos, int num_samples, int num_snps, int snp_index, int entry){

    fseek(ptr, snp_records_pos[snp_index], SEEK_SET);
    // two-bit PLINK encoding: 0 = HOM-REF; 1 = HET; 2 = double-ALT; 3 = missing

    double common, less_common;
    uint8_t byte;

    fread(&byte, sizeof(uint8_t), 1, ptr);

    // Define common / less common categories
    // However, the some people seem to be of third category, which is not well included!

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
    printf("Bytes for snpinfo: %i\n", bytes_for_snpinfo);

    for(int j = 0; j < bytes_for_snpinfo; j++){
        fread(&byte, sizeof(uint8_t), 1, ptr);

        for (int i = 0; i < 8 & j*8+i < num_samples; i++) {
            uint8_t mask = 1 << i;
            snp_data[j*8+i] = (byte & mask) ? common : less_common;
            snp_data_all[entry * num_samples + j*8+i] = (byte & mask) ? common : less_common;
        }
    }

    // ADD VALUES FROM THE DIFFLIST!!! //
    int varint = read_varint(ptr);

    printf("Varint %i \n", varint);

    if(!varint == 0){

        double *snp_sample_values = malloc((2*varint + 1) * sizeof(double)); // First Varint, then SNP values, then sample values
    
        if(snp_sample_values == NULL){
            printf("Memory allocation failed \n");
        }

        snp_sample_values[0] = (double) varint;

        read_difflist(ptr, snp_sample_values, snp_loc, types, lengths, snp_records_pos, num_samples, num_snps, snp_index, varint); 

        if (!varint == 0){    
            for (int i = 1; i <= varint; i++){
                printf("geno %f, sample %i ", snp_sample_values[i], (int)snp_sample_values[i + varint]);
                snp_data[(int)snp_sample_values[i + varint]] = (double)snp_sample_values[i];
                snp_data_all[entry * num_samples + (int)snp_sample_values[i + varint]] = (double)snp_sample_values[i];
            }        
        }

        free(snp_sample_values);
        snp_sample_values = NULL;
    }
}

///////////////////////////////////////////////////////
// Read SNPs of type = 2/3; 'LD-compressed'. Genotypes are the same as previous, 
// EXCEPT for the difflist
///////////////////////////////////////////////////////

void read_type23_record(FILE *ptr, double *snp_data, double *snp_data_all, int *snp_loc, int *types, int *lengths, int *snp_records_pos, int num_samples, int num_snps, int snp_index, int entry){
    
    // Must find the previous non-LD compressed snp.. so must go back in types to find SNP with not type 2 or 3
    int prev_snp = snp_index;
    int diff = 0;

    while(types[prev_snp] == 2 || types[prev_snp] == 3){
        prev_snp = prev_snp - 1;
        diff += 1;
    }
    printf("Previous SNP: %i\n", prev_snp);

    for (int i = 0; i < num_samples; i++){
        if (types[snp_index] == 2) {
            snp_data_all[entry * num_samples + i] = snp_data_all[(entry - diff) * num_samples + i];
        } else if (types[snp_index] == 3) {   // The NA should stay the same
            if(snp_data_all[(entry - diff) * num_samples + i] == 3){
                snp_data_all[entry * num_samples + i] = 3;
            } else {
                snp_data_all[entry * num_samples + i] = 2 - snp_data_all[(entry - diff) * num_samples + i];
            }
        }
        snp_data[i] = snp_data_all[entry * num_samples + i];
    }

    // ADD VALUES FROM THE DIFFLIST!!! //
    int varint = read_varint(ptr);

    if(!varint == 0){
        double *snp_sample_values = malloc((2*varint + 1) * sizeof(double)); // First SNP values, then sample values, 1000 should be enough to test
    
        if(snp_sample_values == NULL){
            printf("Memory allocation failed \n");
        }

        snp_sample_values[0] = varint;

        read_difflist(ptr, snp_sample_values, snp_loc, types, lengths, snp_records_pos, num_samples, num_snps, snp_index, varint); 

        int group_size = snp_sample_values[0];

        if (!group_size == 0){    
            for (int i = 1; i <= group_size; i++){
                //printf("geno %f, sample %i ", snp_sample_values[i], (int)snp_sample_values[i + group_size]);
                // Different for types 2 and 3
                if (types[snp_index] == 2) {
                    snp_data_all[entry * num_samples + (int)snp_sample_values[i + group_size]] = snp_sample_values[i];
                    snp_data[(int)snp_sample_values[i + group_size]] = (int)snp_sample_values[i];
                } else if (types[snp_index] == 3) {
                    if(snp_sample_values[i] == 3){
                        snp_data_all[entry * num_samples + (int)snp_sample_values[i + group_size]] = 3;
                        snp_data[(int)snp_sample_values[i + group_size]] = 3;
                    } else {
                        snp_data_all[entry * num_samples + (int)snp_sample_values[i + group_size]] = 2 - snp_sample_values[i];
                        snp_data[(int)snp_sample_values[i + group_size]] = 2 - (int)snp_sample_values[i];
                    }
                }

                //snp_data[(int)snp_sample_values[i + group_size]] = snp_sample_values[i];
                //snp_data_all[entry * num_samples + (int)snp_sample_values[i + group_size]] = snp_sample_values[i];
            }        
        }

        free(snp_sample_values);
        snp_sample_values = NULL;
    }
}

///////////////////////////////////////////////////////
// Read SNPs of type = 4/6/7; Difflist with samples not 0/2/3
///////////////////////////////////////////////////////

void read_type467_record(FILE *ptr, double *snp_data, double *snp_data_all, int *snp_loc, int *types, int *lengths, int *snp_records_pos, int num_samples, int num_snps, int snp_index, int default_value, int entry){
    fseek(ptr, snp_records_pos[snp_index], SEEK_SET);

    for (int i = 0; i < num_samples; i++){
        snp_data[i] = default_value - 4;
        snp_data_all[entry * num_samples + i] = default_value - 4;
    }

    // ADD VALUES FROM THE DIFFLIST!!! //
    int varint = read_varint(ptr);

    if(!varint == 0){
        double *snp_sample_values = malloc((2*varint + 1) * sizeof(double)); // First SNP values, then sample values, 1000 should be enough to test
    
        if(snp_sample_values == NULL){
            printf("Memory allocation failed \n");
        }

        snp_sample_values[0] = varint;

        read_difflist(ptr, snp_sample_values, snp_loc, types, lengths, snp_records_pos, num_samples, num_snps, snp_index, varint); 

        printf("Varint %i\n", varint);

        if (!varint == 0){    
            for (int i = 1; i <= varint; i++){
                //printf("geno %f, sample %i\n", snp_sample_values[i], (int)snp_sample_values[i + varint]);
                snp_data[(int)snp_sample_values[i + varint]] = snp_sample_values[i];
                snp_data_all[entry * num_samples + (int)snp_sample_values[i + varint]] = snp_sample_values[i];
            }        
        }

        free(snp_sample_values);
        snp_sample_values = NULL;
    }
}


// write a function that converts the binary file to a matrix

int bin_to_mat(FILE *ptr, double *snp_data_all, int *snp_loc, int *types, int *lengths, int *snp_records_pos, int num_samples, int num_snps, int incl){
    int n_blocks = (num_snps-1)/65536 + 1;      // Number of blocks
    int res = (num_snps-1) % 65536 + 1;         // Number of snps in the last block
    int index, default_value;
    double *snp_data;
    
    // for (int i = 0; i < count; i++) {
    for (int i = 0; i < incl - 1; i++) {   
        index = snp_loc[i] - 1;        // Index of the SNP (i.e., row in the .bim file to be analysed)

        printf("SNP %i, type %i, length %i, position %i\n", index, types[index], lengths[index], snp_records_pos[index]);

        if(types[index] == 0){
            snp_data = malloc(num_samples * sizeof(double)); // Allocate memory

            if(snp_data == NULL){
                printf("Memory allocation failed\n");
            }

            printf(" Snp data of SNP %i:\n", index);
            read_type0_record(ptr, snp_data, snp_data_all, snp_loc, types, lengths, snp_records_pos, num_samples, num_snps, index, i);
            
        }
        if(types[index] == 1){
            snp_data = malloc(num_samples * sizeof(double)); // Allocate memory

            if(snp_data == NULL){
                printf("Memory allocation failed\n");
            }

            printf("Snp data of SNP %i:\n", index);
            read_type1_record(ptr, snp_data, snp_data_all, snp_loc, types, lengths, snp_records_pos, num_samples, num_snps, index, i);
    
        }
        if(types[index] == 2 || types[index] == 3){
            snp_data = malloc(num_samples * sizeof(double)); // Allocate memory
            
            if(snp_data == NULL){
                printf("Memory allocation failed\n");
            }

            printf("Snp data of SNP %i:\n", index);
            read_type23_record(ptr, snp_data, snp_data_all, snp_loc, types, lengths, snp_records_pos, num_samples, num_snps, index, i);

        }   
        if(types[index] == 4 || types[index] == 6 || types[index] == 7){
            snp_data = malloc(num_samples * sizeof(double)); // Allocate memory
            double type = types[index];

            if(snp_data == NULL){
                printf("Memory allocation failed\n");
            }

            printf("Snp data of SNP %i:\n", index);
            read_type467_record(ptr, snp_data, snp_data_all, snp_loc, types, lengths, snp_records_pos, num_samples, num_snps, index, type, i);
            
        } 

        for (int j = 0; j < num_samples & j < 10; j++){
            printf("%f ", snp_data[j]);
        }

        free(snp_data);
        snp_data = NULL;
            
        printf("\n");
    }

    printf("printing all SNP values \n");
    for(int i = 0; i < (incl-1) * num_samples; i++){
        //printf("%f ", snp_data_all[i]);
    }

    printf("Test\n");
    return 0;
}