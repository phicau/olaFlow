IHFOAM
======

# Description

IHFOAM is a set of solvers and boundary conditions to generate and absorb water waves actively at the boundaries and to simulate their interaction with porous coastal structures.

You can find all the information regarding the model at its web site and wiki:

http://ihfoam.ihcantabria.com/

http://openfoamwiki.net/index.php/Contrib/IHFOAM

# Download and compilation

## Basic download guide

To get a full copy of IHFOAM source and reference materials run in a terminal:

`git clone git://github.com/phicau/IHFOAM.git`

Code updates can be downloaded in the future from the IHFOAM folder as follows:

`git checkout`

`git pull`

## Basic compilation guide

First compile the boundary conditions:

`cd genAbs`

`./allMake`

Second compile the solvers (ihFoam and ihDyMFoam), selecting the correct version:

`cd solvers/ihFoamXXXXX`

`./allMake`

# Reference materials and tutorials

Reference materials and test cases are included in the IHFOAM download. Check the *reference.zip* file and the *tutorials* folder. For a full description on these materials see:

http://openfoamwiki.net/index.php/Contrib/IHFOAM#IHFOAM_Documentation_and_Tutorials

