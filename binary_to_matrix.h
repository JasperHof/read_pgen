#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

int read_difflist(FILE *ptr, int *snp_sample_values, int *snp_loc, uint8_t *types, uint8_t *lengths, int *snp_records_pos, int num_samples, int num_snps, int snp_index);
void read_type0_record(FILE *ptr, int *snp_data, int *snp_loc, uint8_t *types, uint8_t *lengths, int *snp_records_pos, int num_samples, int num_snps, int snp_index);
void read_type1_record(FILE *ptr, int *snp_data, int *snp_loc, uint8_t *types, uint8_t *lengths, int *snp_records_pos, int num_samples, int num_snps, int snp_index);
void read_type6_record(FILE *ptr, int *snp_data, int *snp_loc, uint8_t *types, uint8_t *lengths, int *snp_records_pos, int num_samples, int num_snps, int snp_index);
int bin_to_mat(FILE *ptr, int *snp_loc, uint8_t *types, uint8_t *lengths, int *snp_records_pos, int num_samples, int num_snps);
