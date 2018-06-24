OLAFLOW
======

# Description

[olaFlow](https://sites.google.com/view/olaflowcfd) (formerly known as olaFoam) is a numerical model conceived as a continuation of the work in [IHFOAM](https://github.com/phicau/IHFOAM). Its development has been continuous from ihFoam (Jul 8, 2014 - Feb 11, 2016) and olaFoam (Mar 2, 2016 - Nov 25, 2017).

This free and open source project is committed to bringing the latest advances in the simulation of wave dynamics to the OpenFOAM® and FOAM-extend communities.

olaFlow includes a set of solvers and boundary conditions to generate and absorb water waves actively at the boundaries and to simulate their interaction with porous coastal structures.

You can find all the information regarding the model at its web site and wiki:

- https://sites.google.com/view/olaflowcfd
- https://openfoamwiki.net/index.php/Contrib/olaFlow

# Citation

olaFlow now has a DOI that can be included in citations:

If you want to reference the model in your publications you can use the following phrase:

> olaFlow [DOI] is an open source project developed within the OpenFOAM® framework as a continuation of the work in Higuera et al. (xxxx). The numerical model enables simulating wave and porous structure interaction in the coastal and offshore fields.

Feel free to modify the phase and adapt it for your own needs.

You can also include any of the [following references](https://sites.google.com/view/olaflowcfd/numerical-model/references/references-internal) when citing the implementation, validation and applications.

# Download and compilation

## Basic download guide

To get a copy of olaFlow source and reference materials, run in a terminal:

`git clone git://github.com/phicau/olaFlow.git`

Code updates can be downloaded in the future from the olaFlow folder as follows:

`git checkout`  
`git pull`

If updating fails, try downloading the code again.

## Basic compilation guide

Compilation is straightforward, simply run the following script from the olaFlow base folder:

`./allMake`

# Reference materials and tutorials

Reference materials and test cases are included in the olaFlow download. Check the *reference* and the *tutorials* folders. For a full description on these materials see:

- https://sites.google.com/view/olaflowcfd/source-code/tutorials
- https://openfoamwiki.net/index.php/Contrib/olaFlow#olaFlow_Documentation_and_Tutorials

This repository only includes the tutorials for the most recent versions of OpenFOAM and OpenFOAM+. The historic record of tutorials released along with olaFlow for the major past versions and FOAM-extend have been moved into their own repository, which you can find in: https://github.com/phicau/olaFlow_oldVersionTutorials

Supplementary materials such as a solver coupled with [isoAdvector](https://github.com/isoAdvector/isoAdvector), turbulence models especially modified for wave simulations and additional tutorial cases can be found in: https://github.com/phicau/olaFlow_supplementary

----------------------------------------------------------
OPENFOAM®  is a registered trade mark of OpenCFD Limited, producer and distributor of the OpenFOAM software via www.openfoam.com.
