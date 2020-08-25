# PolyFit: Clip-art vectorization
Code for the [SIGGRAPH 2020 paper](http://www.cs.ubc.ca/labs/imager/tr/2020/ClipArtVectorization/)

## Building
git clone --recursive https://github.com/dedoardo/polyfit 

Open the folder in Visual Studio 2017+ and build polyfit.exe in Release(WithDebInfo) mode, the output will be written to the cmake source directory.

## Running
polyfit.exe "data/binary_input/pear-32/input.png" "trained/random_forest_paper.txt" vectorization

The output will be written to vectorization/curves_closed

If you are interested in visualizing different stages of the algorithm, uncomment the respective lines in apps/polyfit.cpp and rebuild.

## Dependencies
The code depends on the following libraries, but they should be handled automatically through git submodule.
- [Eigen 3.3+](http://eigen.tuxfamily.org/index.php?title=Main_Page)
- [Andre's Random Forest](https://github.com/bjoern-andres/random-forest)
- [Cairo](https://cairographics.org/download/)
- [Sean Barret's stb](https://github.com/nothings/stb)

## Some of the known problems
- Some extremely complex pixel-art inputs will fail to be segmented and thus provoke a crash. (Our implementation of the depixelize segmentation needs more debugginb) 
- Visual Studio 2019 *might* not work, we encountered some problems with incorrect variable initialization.

## Immediate TODOs
- Automatically verify that this is the correct version of the code and classifier used to generate the results as in the paper
- Add the multicolor test data once that is done