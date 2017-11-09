OLAFOAM
======

# Description

OLAFOAM is a new numerical model initially branched from IHFOAM. This free and open source project is committed to bringing the latest advances in the simulation of wave dynamics to the OpenFOAM® and FOAM-extend communities.

OLAFOAM is a set of solvers and boundary conditions to generate and absorb water waves actively at the boundaries and to simulate their interaction with porous coastal structures.

You can find all the information regarding the model at its web site and wiki:

https://sites.google.com/view/olafoamcfd

https://openfoamwiki.net/index.php/Contrib/OLAFOAM

# Download and compilation

## Basic download guide

To get a full copy of OLAFOAM source and reference materials run in a terminal:

`git clone git://github.com/phicau/OLAFOAM.git`

Code updates can be downloaded in the future from the OLAFOAM folder as follows:

`git checkout`

`git pull`

## Basic compilation guide

Compilation has been simplified, simply run the following script from the OLAFOAM base folder:

`./allMake`

You can also compile the boundary conditions or the solvers (olaFoam and olaDyMFoam) independently:

`cd genAbs`

`./allMake`

`cd solvers/olaFoamXXXXX`

`./allMake`

In this case be sure to select the correct folder for your preferred OpenFOAM/FOAM-extend (OF/FE) version.

# Reference materials and tutorials

Reference materials and test cases are included in the OLAFOAM download. Check the *reference* and the *tutorials* folders. For a full description on these materials see:

https://sites.google.com/view/olafoamcfd/source-code/tutorials
https://openfoamwiki.net/index.php/Contrib/OLAFOAM#OLAFOAM_Documentation_and_Tutorials

----------------------------------------------------------
OPENFOAM®  is a registered trade mark of OpenCFD Limited, producer and distributor of the OpenFOAM software via wwww.openfoam.com.
