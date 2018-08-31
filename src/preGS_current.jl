#==============================================================================#
# preGS.jl
#==============================================================================#

# Author:    Austin Putz
# Created:   Sept 27, 2017
# Modified:  August 31, 2018
# License:   GPLv2

# set version number
preGS_VERSION = "0.0.3"

# developed in Julia version __
JULIA_VERSION = "0.7.0"





#==============================================================================#
# Usage
#==============================================================================#

# Usage:
#			$ julia preGS.jl input_file.jl





#==============================================================================#
# Will determine by itself
#==============================================================================#

# 1. Type of storage you have based on the first few lines it will read in






#==============================================================================#
# Laundry list of shit to do!
#==============================================================================#

# Stars link to importance, even if not in order

# *** Compile into a standalone program (people won't have to download Julia)

# ** Take missing values








#==============================================================================#
# 
# PRINT HEADER
#
#==============================================================================#

# Images are https://www.text-image.com
# OR: http://patorjk.com/software/taag/#p=display&f=Graffiti&t=Type%20Something%20
# Simply select an image and upload it and it will convert it to ASCII Text

println("")
println("--------------------------------------------------------------------------------") 
println("--                                                                            --")
println("--                                           `              -:.//-::.         --")
println("--                                   `.`:+/:o+//o+:.::.--oy+ss-o+/yss         --")
println("--       /:`                      .-:yss+so-:s:.oo:oysos+/:: / /oy//+         --")
println("--       sho+oo+/`             .-:yys-:/  +  +  + .oo:ossoso//oso+sss         --")
println("--       dmshysshy://:.      ..yys-`/` /` /` + ./oso/+/:o+ysss+/s/.           --")
println("--       -d-:oyyyoysoos/`..`-yso`.: `/  / -- + +sso+///::/+oyys+:             --")
println("--        +y   +/`+sssosso++/+.-- `:``: -- / +.ss+.         /ssyo-            --")
println("--          y/  `s`  -- /yyssso///:.`:.`: /`/`+ss-           .soyy:           --")
println("--            .h.  //  .: .:o+-+yso+++:./-./`/--os/             `::-          --")
println("--        `  /s   s`  :yhs/   :++/sho+//--/:-os+.                             --")
println("--        o   s:  -/`+oyo-        -osooy++/:+++:                              --")
println("--        y:  `y`  o/dyo.             /o++o+//:`                              --")
println("--        -s`  :o yhyh.                  .-.-`                                --")
println("--         /+  /hyhs/.                   _________     _ __                   --")
println("--           y+/ddyy.   ____   _______  / ___/ __/    ( ) /                   --")
println("--       `yhyyh-.      /  _ \\/  __/ -_) (_ / \\ \\_    / / /                    --")
println("--       hhdyo-       /  .__/ _/  \\__/\\___/____(_)_ / /_/                     --")
println("--       h/          /__/                        |___/                        --")
println("--                                                                            --")
println("--                                                                            --")
println("--                preGS.jl - A genomic pre-processor based on                 --")
println("--                        Ignacio Aguilar's 'preGSf90'                        --")
println("--                                                                            --")
println("--                               Austin Putz                                  --")
println("--                          Iowa State University                             --")
println("--                             License: GPLv2                                 --")
println("--                   Will be available on GitHub soon...                      --")
println("--                                                                            --")
println("--------------------------------------------------------------------------------") 

#                      ________________ ______________
#  ______________________  ____/_  ___/ ______(_)__  /
#  ___  __ \_  ___/  _ \  / __ _____ \  _____  /__  / 
#  __  /_/ /  /   /  __/ /_/ / ____/ /______  / _  /  
#  _  .___//_/    \___/\____/  /____/_(_)__  /  /_/   
#  /_/                                  /___/         


#                 _________    _ __
#   ___  _______ / ___/ __/   (_) /
#  / _ \/ __/ -_) (_ /\ \_   / / / 
# / .__/_/  \__/\___/___(_)_/ /_/  
#/_/                     |___/     


#println("||--------------------------------------------------------------------------------||") 

# the following was taken from the website: www.network-science.de/ascii/
# using the 'speed' Font

# Print Header
#println("|                                       _________________                        |")
#println("|                 _____________________/   ____/    ____/                        |")
#println("|                 ___  __ \\_ / ___/  _ \\  /  ___ \\____ \\                         |") 
#println("|                 __  /_/ / / /  /  __/  /_/ /  ____/  /                         |")
#println("|                 _  .___/ /_/   \\___/\\_____/  /______/                          |")
#println("|                 /_/                                                            |")







#==============================================================================#
# Packages and Functions
#==============================================================================#

# load Packages
using DataFrames
using UnicodePlots
using ProgressMeter
using CSV
using Printf
using Distributed
using Statistics
using DelimitedFiles
using LinearAlgebra
#using PyPlots









#==============================================================================#
# Read in input files to change defaults
#==============================================================================#

# set directory
current_directory = pwd()

# source input file
parameter_file = ARGS[1]

# read in file with inputs
include(string(current_directory, "/", parameter_file))















#==============================================================================#
# 
# Begin Main Body
#
#==============================================================================#







#==============================================================================#
# Read in Map File
#==============================================================================#

println("\n--------------------------------------------------------------------------------") 
println("---------- Read in Map File ----------")
println("--------------------------------------------------------------------------------\n") 

if map_file == ""

	println("--- No map file specified ----")

else 

println("---- Read in Map File ----")

# read in Map File (SNP name, Chromosome, Position)
@time user_map = CSV.read(map_file, 
	delim      = map_file_delimiter, 
	header     = ["SNP", "Chr", "Pos"], 
	datarow    = map_file_start_row,
	types      = [String, Int64, Int64])

println("\n---- Head of Map File ----")

# head map file
println(head(user_map))

# describe the Map File
println("\n---- Summary of the Map File ---- ")
@printf "-- %-30s: %d \n" "Chromosome Minimum" minimum(user_map[:, :Chr])
@printf "-- %-30s: %d \n" "Chromosome Maximum" maximum(user_map[:, :Chr])
@printf "-- %-30s: %d \n" "Number of SNPs"     length(user_map[:,  :Chr])

# describe the map file
#println(describe(user_map[:, :Chr]))

end








#==============================================================================#
# Read in genotypes
#==============================================================================#

println("\n--------------------------------------------------------------------------------") 
println("---------- Read in Genotypes ----------")
println("--------------------------------------------------------------------------------\n") 






#------------------------------------------------------------------------------#
# Dense Storage
#------------------------------------------------------------------------------#

# Begin Dense Extraction
if SNP_format == "dense"

println("---- Reading in Genotypes From Compact Format ----")

# read in Map File (SNP name, Chromosome, Position)
#@time geno = CSV.read(genotype_file, delim = genotype_file_delimiter, header = ["ID", "SNPs"], datarow = genotype_file_start_row, types = [String, String])
@time geno = CSV.read(genotype_file, delim = ' ', header = ["ID", "SNPs"], datarow = 1, types = [String, String])

# set number of records
n_IDs_orig = size(geno, 1)

# set number of SNPs
n_SNPs_orig = length(geno[1, 2])

# print number of SNPs
@printf "-- %-50s: %d \n" "Number of SNPs in Raw Genotype File" n_SNPs_orig

# sample SNPs
if n_sample_M == 0
	# set the number to sample equal to all SNPs
	x = 1:n_SNPs_orig
	n_SNPs_sample = n_SNPs_orig
else 
	# sample a number of SNPs from the total
	x = unique(rand(1:n_SNPs_orig, n_sample_M))
	n_SNPs_sample = length(x)
end

# print number of SNPs and IDs
@printf "-- %-50s: %d \n" "Number of SNPs (or sample) from Genotype File" n_SNPs_sample
@printf "-- %-50s: %d \n" "Number of IDs" n_IDs_orig


# function to fill the M matrix
function fillM(n_IDs::Int, n_SNPs::Int, geno)
	
	# initiate M with the size
	M = zeros(n_IDs, n_SNPs);

	# show progress meter (took out because it doesn't work well with @parallel)
	#n = n_IDs*n_SNPs
	#p = Progress(n, 1)

	# split string of genotypes into the M matrix
	@distributed for i = 1:n_IDs
		@distributed for j = 1:n_SNPs
			# test if you are sampling genotypes
			if n_sample_M == 0
				M[i, j] = parse(Int64, geno[i, 2][j])
			else
				M[i, j] = parse(Int64, geno[i, 2][x[j]])
			end
			# use progress bar (doesn't work well with parallel computing so I removed
			#next!(p)
		end
	end

	# convert M to an Array
	M = Array(M)

	# return M
	return M

end 


println("-- Loop to fill M Matrix --")
@time M = fillM(n_IDs_orig, n_SNPs_sample, geno)

# print messages of the number of rows and columns
@printf "-- %-50s: %d \n" "Number of rows in M"    size(M, 1)
@printf "-- %-50s: %d \n" "Number of columns in M" size(M, 2)

# calculate the number of SNPs
n_SNPs_current = size(M, 2);


end  # end of "if SNP_format == "dense"





#------------------------------------------------------------------------------#
# Plot M
#------------------------------------------------------------------------------#

println("-- Reshape M --")

# reshape M into a 1-D Vector to Plot!
@time Mvec = vec(reshape(M, (length(M), 1)))

# print histogram of elements of M
println("-- Plot of elements in M --")
println(histogram(Mvec, bins = 5, title = "Plot of the Elements in M (raw)"))
println("-- Done Getting M --")




#==============================================================================#
# Allele Frequencies
#==============================================================================#

println("\n--------------------------------------------------------------------------------") 
println("---------- Allele Frequencies (prior to cleaning) ----------")
println("--------------------------------------------------------------------------------\n") 

# Function to calculate allele frequencies with missing values
function calcAlleleFreqs(M::Array)
	
	# initiate vector for allele frequencies
	allele_freqs = zeros(size(M, 2))
	
	# calculate allele frequencies
	@distributed for i = 1:size(M,2)
		
		# subset each vector
		current_SNP = vec(M[:, i])
		
		# subset to 0, 1, 2
		current_SNP = current_SNP[current_SNP .< 2.0001]
		current_SNP = current_SNP[current_SNP .> -0.0001]
		
		# take mean and divide by 2
		allele_freqs[i] = Statistics.mean(current_SNP) / 2
	
	end
	
	# return allele freqs
	return allele_freqs
	
end

# Calculate allele freqs
allele_freqs_before_clean = calcAlleleFreqs(M)

# histogram of allele frequencies
println("\n-- Plot of Allele Frequenceis --\n")
println(histogram(vec(allele_freqs_before_clean), bins=20, title = "Allele Frequencies (prior to cleaning)"))







#==============================================================================#
# Call Rates
#==============================================================================#

println("\n--------------------------------------------------------------------------------") 
println("---------- Call Rates ----------")
println("--------------------------------------------------------------------------------\n") 

# Initalize Data Frame for IDs
data_IDs  = DataFrame(ID_Number = 1:n_IDs_orig, ID_Orig = geno[:, 1]);

# Initalize Data Frame for SNPs
data_SNPs = DataFrame(SNP_Number = 1:n_SNPs_current);

# add allele freqs
data_SNPs[:AlleleFreqsBeforeClean] = allele_freqs_before_clean

# set new column if map file exits 
#(TO-DO)

# set up vectors to fill
ID_call_rate_orig  = zeros(n_IDs_orig);
SNP_call_rate_orig = zeros(n_SNPs_current);

# loop to check ID call rate
MissingMat = M .== genotype_missing_value

# sum to across row to get animal call rate
data_IDs[:CallRate]  = 1 .- (vec(sum(MissingMat, dims=2)) / size(MissingMat, 2))

# sum down row and divide to get SNP call rate
data_SNPs[:CallRate] = 1 .- (vec(sum(MissingMat, dims=1)) / size(MissingMat, 1))

# set true/false call rates
data_IDs[:CallRateTF]  = data_IDs[:CallRate]  .> call_rate_animals
data_SNPs[:CallRateTF] = data_SNPs[:CallRate] .> call_rate_SNPs

# Print
println(head(data_IDs))
println(head(data_SNPs))

# plot Call Rate IDs
println("\n---- Call Rate IDs ----")

# convert call rate to a vector to plot
IDs_call_rate_vec = vec(data_IDs[:, :CallRate])

# plot histogram of the Animal Call Rates
#println(histogram(IDs_call_rate_vec, bins=5, title="Call Rate - Animals"))
println(histogram(IDs_call_rate_vec, bins=10, title="Call Rate - Animals"))

# print stats
Printf.@printf "-- %-40s: %0.3f, %0.3f, %0.3f, %0.3f \n" "Minimum, Maximum, Average, Median" minimum(IDs_call_rate_vec) maximum(IDs_call_rate_vec) mean(IDs_call_rate_vec) median(IDs_call_rate_vec)

# print summary statistics for Animal Call Rate
describe(data_IDs[:, :CallRate])

# plot Call Rate IDs
println("\n---- Call Rate SNPs ----")
SNPs_call_rate_vec = vec(data_SNPs[:, :CallRate])

# plot histogram of the Animal Call Rates
println(histogram(SNPs_call_rate_vec, bins=10, title="Call Rate - SNPs"))

# print stats
Printf.@printf "-- %-40s: %0.3f, %0.3f, %0.3f, %0.3f \n" "Minimum, Maximum, Average, Median" minimum(SNPs_call_rate_vec) maximum(SNPs_call_rate_vec) mean(SNPs_call_rate_vec) median(SNPs_call_rate_vec)

# print summary statistics for Animal Call Rate
describe(data_SNPs[:, :CallRate])





#==============================================================================#
# Subset M matrix
#==============================================================================#

# Remove Rows for Call Rate on IDs
M = M[vec(data_IDs[:CallRateTF]), :];

# Remove Columns for Call Rate on SNPs
M = M[:, vec(data_SNPs[:CallRateTF])];

# Print size of M after cleaning for Call Rate
println("-- Size of M after Call Rates Removed --")
Printf.@printf "-- %-40s: %d, %d \n" "Size of new M matrix after subsetting for outliers" size(M,1) size(M,2)





#==============================================================================#
# Create Main Matrices
#==============================================================================#

println("\n\n\n--------------------------------------------------------------------------------") 
println("---------- Create Matrices for G ----------")
println("--------------------------------------------------------------------------------\n") 

# check for anything put 0/1/2 genotypes

# ToDo: need to write check for 0/1/2 genotypes




#------------------------------------------------------------------------------#
# Calculate column means of M
#------------------------------------------------------------------------------#

println("---- Column Means of M ----")

# get column means of M matrix
@time col_means = vec(Statistics.mean(M, dims=1))




#------------------------------------------------------------------------------#
# Calculate P
#------------------------------------------------------------------------------#

println("---- Calculate P ----")

# repeat column means (2*p) for P matrix
@time P = reshape(repeat(col_means, inner=size(M, 1)), (size(M,1), size(M,2)));




#------------------------------------------------------------------------------#
# Calculate Z
#------------------------------------------------------------------------------#

println("---- Calculate Z ----")

# subtract off column means
@time Z = M - P;

# calculate the size of G
Gdim = size(Z, 1);



#------------------------------------------------------------------------------#
# Calculate G
#------------------------------------------------------------------------------#

#----------------------------------------#
# Calculate G
#----------------------------------------#

println("---- Calculate G (unscaled) ----")

# calculate G
@time G = Z*Z'

println("---- Size of G ----")

# get size of G
@printf "-- %-40s: %d \n" "Number of rows/columns in G" size(G, 2)



#------------------------------------------------------------------------------#
# Allele Freqs & Scale Factor
#------------------------------------------------------------------------------#

println("---- Calculate Allele Frequencies (p & q) ----")

# get allele frequencies
@time p = vec(col_means ./ 2)  # column means already calculated
@time q = vec(1 .- p)

println("---- Calculate Scale Factor ----")

# set sum 2pq
@time sum2pq = 2*sum(p.*q)

println("Sum2pq is: $sum2pq")



#------------------------------------------------------------------------------#
# Scale G
#------------------------------------------------------------------------------#

println("---- Calculate G (scaled) ----")

# calculate G
@time G = Array(G / sum2pq)



# write G matrix
open("G.txt", "w") do io
	writedlm(io, G)
end



println("")





#==============================================================================#
# Plot components
#==============================================================================#

println("\n--------------------------------------------------------------------------------") 
println("---------- Plot Elements of G ----------")
println("--------------------------------------------------------------------------------\n") 




#----------------------------------------#
# Diagonals of G
#----------------------------------------#

println("---- Print Diagonal Elements of G ----\n")
# get diags
@time Gdiags = vec(LinearAlgebra.diag(G))

# histogram of my diagonals
println(histogram(Gdiags, bins = 20, title = "G Matrix Diagonals"))




#----------------------------------------#
# Off-Diagonals of G
#----------------------------------------#

println("---- Print Off-Diagonal Elements of G ----\n")

function extractOffdiags(G::Array)

	# calculate the size of G
	sizeG = size(G, 1)

	# initiate off-diagonal Vector
	# offdiagsG = [];
	#offdiagsG = zeros((((1+sizeG)*sizeG)/2), 1);  # equation for upper triangular part
	#offdiagsG = zeros(Float64, convert(Float64, (((1+sizeG)*sizeG)/2), 1));
	offdiag_length = convert(BigInt, ((1+sizeG)*sizeG) / 2)
	println("\nOffdiagonal length is : $offdiag_length\n")
	#offdiagsG = ones(offdiag_length, 1)
	offdiagsG = Array{Float64,1}(undef, offdiag_length)

	# set count
	count = 1;

	# start progress bar
	n = length(offdiagsG)
	p = Progress(n, 1)

	# loop to fill off-diag vector
	for i = 1:sizeG
		for j = i+1:sizeG
	
			# add off-diag to list
			if i < j 
				offdiagsG[count, 1] = G[i, j]
			end
	
			# advance count
			count = count + 1

			# advance progress bar
			next!(p)
	
		end
	end

	# convert to vector
	offdiagsG = vec(offdiagsG)

	# return off-diagonals
	return offdiagsG

end # END function

# extract the off-diagonals of G
@time offdiagsG = extractOffdiags(G)

# print offdiags of G
println(histogram(offdiagsG, bins = 20, title = "G Matrix Off-Diagonals"))


#----------------------------------------#
# Relationships of G
#----------------------------------------#

println("---- Print Relationships of G ----\n")

function extractRelationships(G::Array)

	# size of G
	sizeG = size(G, 1)

	# initalize the vector
	Grels = Float64[];

	# loop through to standardize the G matrix to relationships
	for i = 1:sizeG
		for j = i+1:sizeG
	
			# Skip if it's a diagonal or lower triangular
			if i == j || i > j
				continue
			end
			
			# fill G relationship matrix
			push!(Grels, (G[i, j] / sqrt( G[i,i] * G[j,j] )))
			
		end
	end

	# convert to vector
	Grels = vec(Grels)

	# return relationships
	return Grels

end # END Function

@time Grels = extractRelationships(G)

# print histogram of relationships
println(histogram(Grels, bins = 20, title = "G relationships"))




#==============================================================================#
# PCA Plot
#==============================================================================#

println("\n--------------------------------------------------------------------------------") 
println("---------- Plot PCA's ----------")
println("--------------------------------------------------------------------------------\n") 


println("---- Singular Value Decomposition for PCA's ----")
# Singular Value Decomposition
@time U, S, V = LinearAlgebra.svd(G)

# sum of s^2 values is the total variance
S2 = S.^2

# get each percent explained
S_pct_exp = vec((S2 ./ sum(S2)).*100)

# pull out the first few PCAs
S_1 = S_pct_exp[1]
S_2 = S_pct_exp[2]
S_3 = S_pct_exp[3]
S_4 = S_pct_exp[4]
S_5 = S_pct_exp[5]

# print how much they explain
@printf "\n-- %-40s: %0.2f, %0.2f, %0.2f, %0.2f, %0.2f \n " "Variances Explained of 1st 5 PCAs" S_1 S_2 S_3 S_4 S_5

# extract first 10 variance explained
S_pct_exp_first_10 = S_pct_exp[1:10]

# get vector list 1 to 10
first_10 = collect(1:10)

println("\n")

# plot
println(lineplot(first_10, S_pct_exp_first_10, title = "Scree Plot of First 10 PCAs (var exp)"))

println("\n\n")

# if you set only one PC, set to 2 to plot between 1 and 2
if n_PCs == 1
	n_PCs = 2
end

	# Loop through to make Plots
	for i = 1:n_PCs
		for j = 2:n_PCs

		if i == j || j < i || i == n_PCs
			continue
		end
	
		# plot the 1,2 principal components
		println(scatterplot(vec(U[:, i]), vec(U[:, j]), 
			color = :red, 
			title = "PC $i vs PC $j"))
	
		end
	end # end loop for PCA plots

# print space
println("\n\n")


















