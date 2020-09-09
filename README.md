# mCSAC

Multi-Chromosome Constrained Self-Avoiding Chromatin

## Installation

1- Code requires installation of cmake and boost and eigen3 libraries.

2- After libraries are installed, make sure to change the paths in CMakeLists.txt file.

3- Create a folder called build in the directory where you have the code 
  ```
  cd build 
  cmake .. 
  make 
  ```
  
Above should compile the code. Hopefully no error.

## Running

4- to test, I have a test configuration file in the folder. You can run it as 
```
../build/bin/chromatin.sis.coarse -conf <<.ini file>> -prefix AA
```

The components
```
../build/bin/chromatin.sis.coarse  = executable created after building the code 
<<.ini file>> = configuration file needed to point to the input files 
AA = prefix for the output chain files. 
```
This helps to iterate over when creating multiple chains. 

## Brief explanation of the <<.ini>> file
```
out_path = output path
persistence_length = length of each unit in A (e.g. 30 nm chromatin fiber will have 30, each 10nm is a kb long) )
collision_length = 300 (in A, this is the diameter of each bead)
nucleus_sphere_diameter = 5200 (nucleus diameter in A)
number_nodes = number of beads in the chain (i.e., the total length of the chain is number of nodes x 3 kb) 
number_of_chains = number of chromosome
chain_length_i = number of nodes in the ith chromosome
number_sample_points = divides the sampling sphere around a monomer to 640 distrete points
```

# Citation

When referencing this code, please cite:

(1) Gamze Gursoy, Yun Xu, and Jie Liang. "Spatial organization of budding yeast genome from landmark constraints and predicting biological chromatin interactions from Chromosome Conformation Capture Data", PLOS Comp. Bio., 2017:13(7):e1005658

(2) Gamze Gursoy, Yun Xu and Jie Liang. "Computational predictions of structures of multichromosomes of budding yeast", Conf Proc IEEE Eng Med Biol Soc., 2014

# See also

CSAC

https://github.com/uic-lianglab/CSAC

nCSAC

https://github.com/uic-lianglab/nCSAC
