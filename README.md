# preGS.jl

preGS.jl is a pre-genomic processor for analysis in animal breeding. 

preGS.jl was based on `preGSf90` by Ignacio Aguilar in the [BLUPF90](http://nce.ads.uga.edu/wiki/doku.php) 
family of programs. You can find the documentation for preGSf90 [here](http://nce.ads.uga.edu/wiki/doku.php?id=readme.pregsf90). 

Ignacio is a very talented programmer and scientist and this is my poor starting attempt 
at copying his program in Julia from scratch. I'd like to write it in a way that is more flexible
from the beginning, such as input formats and output formats. Ignacio has incorporated some of
this into his package. 

-------

I'm not making this a package through the package manager quite yet, 
but in time I hope I can make it available through the package manager. 
To make it easier. For now, you can simply copy the function to your
own computer. 

# How to run?

Currently I use this format

```julia
julia preGS.jl input.jl
```

Where the `input.jl` file will read in all of the inputs you need like a 
parameter file. I don't know how to write a simpler one at this time. 

# Formats

## Genotype File

Currently accepts only a simple format like `preGSf90`. Only the animal ID in column 1 and the string of genotypes in 

```
SIRE_100     201200201020100202
SIRE_200     201200201010221222
...
DAM_001      101200021002020102
```

## Map File

Simply the SNP name, the chromosome (1,2,...,n) in integer form, and position on chromosome (integer form). 

```
SNP-5051042 1 10232
SNP-0510203 1 10400
...
SNP-2004134 18 1401020
```

# Options

I'm trying to allow for many options to be fairly diverse in it's ability to do things. 
Right now I have you read in a 'parameter file' of sorts to set the options. 
Hopefully it's easier than setting 100 options through the function. 

# Limitations

I have only tested this on Ubuntu 16.04 (Linux). I will continue to test on other 
computers as I have time. 

# Comments, Questions, and Suggestions

Always welcome. Please drop me an email at putz.austin@gmail.com

Thanks for visiting my page. Please see others. 
