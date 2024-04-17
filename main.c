#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

#include "find_snp_loc.c"
#include "binary_to_matrix.c"
#include "general_functions.c"


// include the other functions

int read_pgen(FILE *pgen, char filename[], int chr, int pos_start, int pos_end);

// int* find_location(char filename[], int chr, int pos_start, int pos_end);
// int* find_snp_records(uint8_t* lengths, int start_pos_types[], int num_snps);
// int bin_to_mat(FILE *ptr, int *snp_loc, uint8_t *types, uint8_t *lengths, int *snp_records_pos, int num_samples, int num_snps);

// Functions for calling the different types of records
// int read_varint(FILE *ptr);
// int* read_difflist(FILE *ptr, int *snp_sample_values, int *snp_loc, uint8_t *types, uint8_t *lengths, int *snp_records_pos, int num_samples, int num_snps, int snp_index);
// void read_type0_record(FILE *ptr, int *snp_data, int *snp_loc, uint8_t *types, uint8_t *lengths, int *snp_records_pos, int num_snps, int snp_index);
// void read_type1_record(FILE *ptr, int *snp_data, int *snp_loc, uint8_t *types, uint8_t *lengths, int *snp_records_pos, int num_snps, int snp_index);
void read_type6_record(FILE *ptr, int *snp_data, int *snp_loc, uint8_t *types, uint8_t *lengths, int *snp_records_pos, int num_samples, int num_snps, int snp_index);

// uint8_t* find_types(FILE *ptr, int start_pos_types, int num_snps, int bits_per_record_type, int bytes_per_record_length);
// uint8_t* find_lengths(FILE *ptr, int start_pos_types, int num_snps, int bits_per_record_type, int bytes_per_record_length);

int main(int argc, char *argv[]){

    int chr, pos_start, pos_end;

    char filename[100];
    strcpy(filename, argv[1]);
    chr = atoi(argv[2]);
    pos_start = atoi(argv[3]);
    pos_end = atoi(argv[4]);

    int length = strlen(filename); // size of file name
    FILE *ptr = fopen(filename, "rb");

    printf("File name: %s\n", filename);

    // find_location(filename, chr, pos_start, pos_end);

    if (length >= 5 && strcmp(filename + length - 5, ".pgen") == 0){ // check the .pgen extension
        printf("Reading a .pgen file \n");
        
        FILE *ptr = fopen(filename, "rb");
        if (ptr == NULL) {
            printf("Error opening file\n");
            return 1;
        }

        read_pgen(ptr, filename, chr, pos_start, pos_end);
        fclose(ptr);

    } else {
        printf("File does not have the .pgen extension \n");
        fclose(ptr);    
    }
    
    return 0;
}

int read_pgen(FILE *ptr, char filename[], int chr, int pos_start, int pos_end){ // this function reads the .pgen file
    
    int num_snps, num_samples, bits_per_record_type, bytes_per_record_length;
    uint8_t info_snps_samples[4];
    size_t elements_read;
    unsigned char magic_numbers[2], filetype[1];

    // Read the 'magic' numbers
    elements_read = fread(magic_numbers, sizeof(unsigned char), 2, ptr); 
    if(!magic_numbers[0] == 128 || !magic_numbers[1] == 27){
        printf("This is not a .pgen file \n");
    } else {printf("This is a .pgen file \n");}

    // Check the file type, described by the third byte
    elements_read = fread(filetype, sizeof(unsigned char), 1, ptr);
    if(!filetype[0] == 16){
        printf("This is not a standard PLINK2 format \n");
        exit(1);
    } else {printf("This a standard PLINK2 format \n");}

    // Read the number of SNPs and samples
    fseek(ptr, 3, SEEK_SET);
    elements_read = fread(info_snps_samples, sizeof(uint8_t), 4, ptr);
    num_snps = (info_snps_samples[0]) | (info_snps_samples[1] << 8) | (info_snps_samples[2] << 16) | (info_snps_samples[3] << 24);

    fseek(ptr, 7, SEEK_SET);  // Line is redundant, but kept for clarity
    elements_read = fread(info_snps_samples, sizeof(uint8_t), 4, ptr);
    num_samples = (info_snps_samples[0]) | (info_snps_samples[1] << 8) | (info_snps_samples[2] << 16) | (info_snps_samples[3] << 24);
    printf("Number of SNPs and samples: %u and %u \n", num_samples, num_snps);

    // Read the number of bits per record type and bytes per record length
    fseek(ptr, 11, SEEK_SET);
    elements_read = fread(info_snps_samples, sizeof(uint8_t), 1, ptr);
    unsigned char SNP_info = info_snps_samples[0] & 0x0F;

    if(SNP_info >= 4){
        bits_per_record_type = 8; 
    } else {bits_per_record_type = 4;}
    bytes_per_record_length = (SNP_info % 4) + 1;

    printf("Bits per record: %u, bytes per record length: %u \n", bits_per_record_type, bytes_per_record_length);

    // Read bytes per allele count, and bytes per REF-allele bitarray (not too important for now)
    uint8_t bits_5_6 = (info_snps_samples[0] >> 4) & 0x03;
    uint8_t bits_7_8 = (info_snps_samples[0] >> 6) & 0x03;
    printf("Allele count: %u, REF bitarray: %u \n", bits_5_6, bits_7_8);

    // Find the starting positions of the variant records
    int n_blocks = (num_snps-1)/65536 + 1;
    int start_pos_variant_records[n_blocks], block_sizes[n_blocks];
    uint8_t vector[8];

    for (int i = 1; i <= n_blocks; i++) {
        fseek(ptr, 12 + (i-1)*8, SEEK_SET);
        elements_read = fread(vector, sizeof(uint8_t), 8, ptr);
        block_sizes[i] = (i == n_blocks) ? num_snps - 65536*(n_blocks-1) : 65536;
        start_pos_variant_records[i] = ((uint64_t)vector[0]) |
                ((uint64_t)vector[1] << 8) |
                ((uint64_t)vector[2] << 16) |
                ((uint64_t)vector[3] << 24) |
                ((uint64_t)vector[4] << 32) |
                ((uint64_t)vector[5] << 40) |
                ((uint64_t)vector[6] << 48) |
                ((uint64_t)vector[7] << 56);
        printf("SNPs in this block: %i \n", block_sizes[i]);
        printf("Starts at: %i \n", start_pos_variant_records[i]);
    }

    int start_pos_types = 12 + n_blocks*8;

    // Find the starting positions of the variant record types (either 1 or 2 stored, depending on bits_per_record_type)
    for (int i = 1; i <= n_blocks; i++) {
        // we expect ceiling(block_sizes[i] * bits_per_record_type/8) bytes
        int n_bytes_rec_types = (block_sizes[i] * bits_per_record_type + 7) / 8;
        int n_bytes_reclength_types = block_sizes[i] * bytes_per_record_length;
        // printf("Expecting %i bytes for variant record types \n", n_bytes_rec_types);
        // printf("Expecting %i bytes for variant lengths \n", n_bytes_reclength_types);
    }

    printf("Input %s %i %i %i \n", filename, chr, pos_start, pos_end);

    // Location, types, and lengths of the SNP records
    int *snp_loc = find_location(filename, chr, pos_start, pos_end);
    uint8_t *types = find_types(ptr, start_pos_types, num_snps, bits_per_record_type, bytes_per_record_length);
    uint8_t *lengths = find_lengths(ptr, start_pos_types, num_snps, bits_per_record_type, bytes_per_record_length);
    int *snp_records_pos = find_snp_records(lengths, start_pos_variant_records, num_snps);

    // information of the first ten SNPs:
    for (int j = 0; j < 15; j++){
        printf("Location; %i, type: %i, length: %i, position: %i \n", snp_loc[j], types[j], lengths[j], snp_records_pos[j]);

        fseek(ptr, snp_records_pos[j], SEEK_SET);
        for(int i = 0; i < lengths[j]; i++){
            uint8_t byte;
            fseek(ptr, snp_records_pos[j] + i, SEEK_SET);
            fread(&byte, sizeof(uint8_t), 1, ptr);

            printf("Byte %i, ", byte);
            for (int k = 7; k >= 0; k--) {
                printf("%d", (byte >> k) & 1); // Extract k-th bit and 
            }
            printf("  ");
        }
        printf("\n\n");
    }

    bin_to_mat(ptr, snp_loc, types, lengths, snp_records_pos, num_samples, num_snps);

    fclose(ptr);
    free(snp_loc);

    return(0);
}

