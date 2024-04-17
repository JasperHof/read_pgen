#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

int* find_location(char filename[], int chr, int pos_start, int pos_end){
    printf("test\n");
    
    char var_file[256], line[100]; // Assuming a maximum filename length of 255 characters
    
    FILE *var;

    // Copy filename to var_file
    strcpy(var_file, filename);
    strcpy(strstr(var_file, ".pgen"), ".pvar");

    printf("File name: %s\n", var_file);
    printf("Searching on chr %i between pos %i and %i \n", chr, pos_start, pos_end);
    
    var = fopen(var_file, "r");

    int i = 1;
    int j = 1;                                      // Index for the snps in the set
    int *snp_loc = malloc(1000000 * sizeof(int));   // Large vector, possible to look into 1M variants

    // Read each line from the file
    while (fgets(line, sizeof(line), var) != NULL) {

        // Parse the line to extract values in the first and third columns
        int col1, col2;

        // printf("Line: %s\n", line);

        if (sscanf(line, "%d %d %*s %*s %*s", &col1, &col2) == 2) {
            // Check if the conditions are met
            if (col1 == chr && col2 > pos_start && col2 < pos_end) {
                // Print the row number
                // printf("Row number: %ld\n", i);
                snp_loc[j] = i-1;                       // The header doesn't count, so substract 1
                j++;
            }
        }
        i++;
    }

    printf("Checked: %i\n", i);

    return snp_loc;
}

uint8_t* find_types(FILE *ptr, int start_pos_types, int num_snps, int bits_per_record_type, int bytes_per_record_length){

    uint8_t *types = malloc(num_snps * sizeof(uint8_t));
    int n_blocks = (num_snps-1)/65536 + 1;      // Number of blocks
    int res = (num_snps-1) % 65536 + 1;         // Number of snps in the last block
    int j = 0, block_bytes = (65536 * bits_per_record_type + 7) / 8 + 65536*bytes_per_record_length;

    printf("Block bytes: %i\n", block_bytes);

    for (int i = 1; i <= n_blocks; i++){

        fseek(ptr, start_pos_types, SEEK_SET);

        if(i == n_blocks){ // Reading the final block
            int n_bytes_rec_types = (res * bits_per_record_type + 7) / 8;
            
            for (int i = 0; i < n_bytes_rec_types; i++) {
                uint8_t byte;

                fseek(ptr, start_pos_types+i, SEEK_SET);
                fread(&byte, sizeof(uint8_t), 1, ptr);

                // Extract lower and higher 4 bits
                types[j] = byte & 0x0F;
                types[j+1] = (byte >> 4) & 0x0F;
                j += 2;

                if(i == n_bytes_rec_types-1){
                    printf("Last byte: %i\n", start_pos_types+i);
                }
            }

        } else {
            int n_bytes_rec_types = (65536 * bits_per_record_type + 7) / 8;
            
            for (int i = 0; i < n_bytes_rec_types; i++) {
                uint8_t byte;

                fseek(ptr, start_pos_types+i, SEEK_SET);
                fread(&byte, sizeof(uint8_t), 1, ptr);

                // Extract lower and higher 4 bits
                types[j] = byte & 0x0F;
                types[j+1] = (byte >> 4) & 0x0F;
                j += 2;
            }

            start_pos_types = start_pos_types + block_bytes;
        }
    }
    printf("Count: %i\n", j);

    return types;

}

uint8_t* find_lengths(FILE *ptr, int start_pos_types, int num_snps, int bits_per_record_type, int bytes_per_record_length){

    uint8_t *lengths = malloc(num_snps * sizeof(uint8_t));
    int n_blocks = (num_snps-1)/65536 + 1;      // Number of blocks
    int res = (num_snps-1) % 65536 + 1;         // Number of snps in the last block

    int skip = (65536 * bits_per_record_type + 7) / 8;
    int j = 0, block_bytes = (65536 * bits_per_record_type + 7) / 8 + 65536*bytes_per_record_length;

    printf("Block bytes: %i\n", block_bytes);

    for (int i = 1; i <= n_blocks; i++){

        if(i == n_blocks){ // Reading the final block
          
            for (int i = 0; i < res; i++) {
                uint8_t byte;
                fseek(ptr, start_pos_types+skip+i*bytes_per_record_length, SEEK_SET);
        
                fread(&byte, sizeof(uint8_t), 1, ptr);

                lengths[j] = byte;
                j += 1;
            }

        } else {

            for (int i = 0; i < 65536; i++) {
                uint8_t byte;
                fseek(ptr, start_pos_types+skip+i*bytes_per_record_length, SEEK_SET);
        
                fread(&byte, sizeof(uint8_t), 1, ptr);

                lengths[j] = byte;
                j += 1;
            }

            start_pos_types = start_pos_types + block_bytes;
        }
    }

    return lengths;
}

int* find_snp_records(uint8_t* lengths, int start_pos_variant_records[], int num_snps){
    
    int n_blocks = (num_snps-1)/65536 + 1;      // Number of blocks
    int res = (num_snps-1) % 65536 + 1;         // Number of snps in the last block
    int j = 0;
    // Compute bytes positions (might as well do already)

    int *byte_pos = malloc(num_snps * sizeof(int));

    for (int k = 1; k <= n_blocks; k++){
        byte_pos[j] = start_pos_variant_records[k];
        if(k == n_blocks){ // Reading the final block
            byte_pos[j] = start_pos_variant_records[k];
            j+=1;

            for (int i = 1; i < res; i++) {
                byte_pos[j] = byte_pos[j-1] + lengths[j-1];
                j += 1;
            }
        } else {
            byte_pos[j] = start_pos_variant_records[k];
            j+=1;

            for (int i = 1; i < 65536; i++) {
                byte_pos[j] = byte_pos[j-1] + lengths[j-1];
                j += 1;
            }
        }
    }

    return byte_pos;
}