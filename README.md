# mCSAC

Multi-Chromosome Constrained Self-Avoiding Chromatin

# installation

1- Code requires installation of cmake and boost and eigen3 libraries.

2- After libraries are installed, make sure to change the paths in CMakeLists.txt file.

3- Create a folder called build in the directory where you have the code
cd build
cmake ..
make
Above should compile the code. Hopefully no error.

# running

4- to test, I have a test configuration file in the folder. You can run it as
../build/bin/ensemble -conf ./ensemble_try.ini -prefix AA -samplesize $SampleSize

ensemble_try.ini --> configuration file with the inputs needed for yeast genome
AA --> prefix of the output files
$SampleSize --> number of genomes to generate at one run

# Citation

When referencing this code, please cite:

(1) Gamze Gursoy, Yun Xu, and Jie Liang. "Spatial organization of budding yeast genome from landmark constraints and predicting biological chromatin interactions from Chromosome Conformation Capture Data", PLOS Comp. Bio., 2017:13(7):e1005658

(2) Gamze Gursoy, Yun Xu and Jie Liang. "Computational predictions of structures of multichromosomes of budding yeast", Conf Proc IEEE Eng Med Biol Soc., 2014

# See also

CSAC

https://github.com/uic-lianglab/CSAC

nCSAC

https://github.com/uic-lianglab/nCSAC
