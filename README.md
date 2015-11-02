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

# Reference materials

Reference materials and test cases are included in the IHFOAM download. Check the *reference.zip* file and the *tutorials* folder. For a full description on these materials see:

http://openfoamwiki.net/index.php/Contrib/IHFOAM#IHFOAM_Documentation_and_Tutorials

## Notes

Note that OpenFOAM 2.3.0 and 2.4.0 tutorials are now delivered. However, problems with pressure calculations (not linked in any way to IHFOAM implementation) are experienced in these versions and no corrections are expected from the OpenFOAM developers:

http://www.cfd-online.com/Forums/openfoam-bugs/139286-interfoam-pressure-miscalculation-2-3-wrt-previous-versions.html

# Who do I talk to?

https://openfoamwiki.net/index.php/Contrib/IHFOAM#Get_Connected
