#==============================================================================#
# Options to give in parameter file
#==============================================================================#

#------------------------------------------------------------------------------#
# Files
#------------------------------------------------------------------------------#

# set genotype, map, and pedigree file names
#genotype_file = "ALGP2_Cycle_1_5_Geno_BLUPF90.txt_clean"
genotype_file = "genotypes.txt"
map_file      = "ChrInfoBLUPF90.txt"
ped_file      = "dam_1_5.ped"

#------------------------------------------------------------------------------#
# General Options
#------------------------------------------------------------------------------#





#------------------------------------------------------------------------------#
# Map File Options
#------------------------------------------------------------------------------#

# Map file options
map_file_delimiter = ' '      # [char] needs to be a character with ' not "
map_file_header    = false    # [logical] is there a header on the map file
map_file_start_row = 1        # [Int] Line to start reading from



#------------------------------------------------------------------------------#
# Genotype File Options
#------------------------------------------------------------------------------#

# SNP format:
# "dense" =      ID1     0501102020102052
# "space" =      ID1     0 5 0 1 1 0 2 0 2 0 1 0 2 0 5 2
SNP_format = "dense"

# headers on map file or genotype file
genotype_file_delimiter   = ' '        # [char] set what separates in the genotype file
genotype_file_start_row   = 1          # [Int] row to start reading
genotype_missing_value    = 5          # [Int] set missing value for genotypes (will take anything not 0,1,2)

# NOT USED
n_skip_genotype_file      = 0          # [Int] number of lines to skip in the genotype file
genotype_header           = false      # [logical] ask if there is a header
genotype_start_col        = 2          # [Int] 

# number of SNPs to sample
#n_sample_M = 0              # Select all SNPs
n_sample_M = 50000                  # [Int] number to sample from genotype file




#------------------------------------------------------------------------------#
# QC Options
#------------------------------------------------------------------------------#

# Call Rates (minimum proportion of non-missing genotypes for each)
call_rate_animals  = 0.90
call_rate_SNPs = 0.90

# MAF threshold
MAF_threshold = 0.05

# Exclude IDs
exclude_animals = []
exclude_SNPs    = []
exclude_Chr     = []
exclude_SNP     = []

#----------------------------------------#
# PCA Options
#----------------------------------------#

# set number of Principal components to use and plot bivariate scatterplots
n_PCs = 3

#----------------------------------------#
# Save Options
#----------------------------------------#

# save clean SNPs
save_clean = "yes"

# Save Raw Relationship Matrices
saveG = "no"
saveA = "no"
saveH = "no"

# Save Inverses of relationship matrices
saveGinv = "no"
saveAinv = "no"
saveHinv = "no"









