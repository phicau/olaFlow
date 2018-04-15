OLAFLOW
======

# Description

olaFlow (formerly known as olaFoam) is a numerical model initially branched from IHFOAM. This free and open source project is committed to bringing the latest advances in the simulation of wave dynamics to the OpenFOAM® and FOAM-extend communities.

olaFlow is a set of solvers and boundary conditions to generate and absorb water waves actively at the boundaries and to simulate their interaction with porous coastal structures.

You can find all the information regarding the model at its web site and wiki:

https://sites.google.com/view/olaflowcfd

https://openfoamwiki.net/index.php/Contrib/olaFlow

# Download and compilation

## Basic download guide

To get a full copy of olaFlow source and reference materials, run in a terminal:

`git clone git://github.com/phicau/olaFlow.git`

Code updates can be downloaded in the future from the olaFlow folder as follows:

`git checkout`

`git pull`

## Basic compilation guide

Compilation has been simplified, simply run the following script from the olaFlow base folder:

`./allMake`

You can also compile the boundary conditions or the solvers (olaFlow and olaDyMFlow) independently:

`cd genAbs`

`./allMake`

`cd solvers/olaFlowXXXXX`

`./allMake`

In this case be sure to select the correct folder for your preferred OpenFOAM/FOAM-extend (OF/FE) version.

# Reference materials and tutorials

Reference materials and test cases are included in the olaFlow download. Check the *reference* and the *tutorials* folders. For a full description on these materials see:

- https://sites.google.com/view/olaflowcfd/source-code/tutorials
- https://openfoamwiki.net/index.php/Contrib/olaFlow#olaFlow_Documentation_and_Tutorials

This repository only includes the tutorials for the most recent versions of OpenFOAM and OpenFOAM+.

The historic record of tutorials released along with olaFlow for the major past versions and FOAM-extend have been moved into their own repository, which you can find in: https://github.com/phicau/olaFlow_oldVersionTutorials

----------------------------------------------------------
OPENFOAM®  is a registered trade mark of OpenCFD Limited, producer and distributor of the OpenFOAM software via www.openfoam.com.
