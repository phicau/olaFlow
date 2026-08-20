/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2014 OpenFOAM Foundation
    Copyright (C) 2016-2021 OpenCFD Ltd.
-------------------------------------------------------------------------------
License
    This file is a derivative work of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.
\*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*\
Application
    overOlaDyMFlow

Group
    grpMultiphaseSolvers grpMovingMeshSolvers

Description
    Solver for two incompressible, isothermal immiscible fluids using a VOF
    (volume of fluid) phase-fraction based interface capturing approach,
    with optional mesh motion and mesh topology changes including adaptive
    re-meshing.

    overOlaDyMFlow solves the VARANS equations.

\*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*\
| olaFlow Project                                       ll                    |
|                                                       l l                   |
|   Coder: Pablo Higuera Caubilla                 ooo   l l     aa            |
|   Bug reports: olaFlowCFD@gmail.com            o   o  l l    a  a           |
|                                                o   o  ll   l a  aa  aa      |
|                                                 ooo    llll   aa  aa        |
|                                                                             |
|                                                FFFFF L     OOOOO W   W      |
|                                                F     L     O   O W   W      |
|                                                FFFF  L     O   O W W W      |
|                                                F     L     O   O WW WW      |
|                                                F     LLLLL OOOOO W   W      |
|   -----------------------------------------------------------------------   |
| References:                                                                 |
|                                                                             |
| - Three-dimensional interaction of waves and porous coastal structures      |
|    using OpenFOAM. Part I: Formulation and validation.                      |
|    Higuera, P., Lara, J.L. and Losada, I.J. (2014)                          |
|    Coastal Engineering, Vol. 83, 243-258.                                   |
|    http://dx.doi.org/10.1016/j.coastaleng.2013.08.010                       |
|                                                                             |
| - Three-dimensional interaction of waves and porous coastal structures      |
|    using OpenFOAM. Part II: Application.                                    |
|    Higuera, P., Lara, J.L. and Losada, I.J. (2014)                          |
|    Coastal Engineering, Vol. 83, 259-270.                                   |
|    http://dx.doi.org/10.1016/j.coastaleng.2013.09.002                       |
|                                                                             |
| - Three-dimensional numerical wave generation with moving boundaries        |
|    Higuera, P., Lara, J.L. and Losada, I.J. (2015)                          |
|    Coastal Engineering, Vol. 101, 35-47.                                    |
|    https://doi.org/10.1016/j.coastaleng.2015.04.003                         |
|                                                                             |
\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "dynamicFvMesh.H"
#include "CMULES.H"
#include "EulerDdtScheme.H"
#include "localEulerDdtScheme.H"
#include "CrankNicolsonDdtScheme.H"
#include "subCycle.H"
#include "immiscibleIncompressibleTwoPhaseMixture.H"
#include "turbulentTransportModel.H"
#include "pimpleControl.H"
#include "fvOptions.H"
#include "fvcSmooth.H"
#include "cellCellStencilObject.H"
#include "localMin.H"
#include "oversetAdjustPhi.H"
#include "oversetPatchPhiErr.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::addNote
    (
        "Solver for two incompressible, isothermal immiscible fluids using"
        " VOF phase-fraction based interface capturing\n"
        "With optional mesh motion and mesh topology changes including"
        " adaptive re-meshing."
    );

    #include "postProcess.H"

    #include "setRootCaseLists.H"
    #include "createTime.H"
    #include "createDynamicFvMesh.H"
    #include "initContinuityErrs.H"

    #include "createDyMControls.H"
    #include "createFields.H"
    #include "createAlphaFluxes.H"
    #include "createFvOptions.H"

    volScalarField rAU
    (
        IOobject
        (
            "rAU",
            runTime.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("rAUf", dimTime/rho.dimensions(), 1.0)
    );

    #include "createUf.H"
    #include "createControls.H"

    #include "setCellMask.H"
    #include "setInterpolatedCells.H"

    turbulence->validate();

    if (!LTS)
    {
        #include "CourantNo.H"
        #include "setInitialDeltaT.H"
    }

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
    Info<< "\nStarting time loop\n" << endl;

    while (runTime.run())
    {
        #include "readDyMControls.H"
        #include "readOversetDyMControls.H"

        if (LTS)
        {
            #include "setRDeltaT.H"
        }
        else
        {
            #include "CourantNo.H"
            #include "alphaCourantNo.H"
            #include "setDeltaT.H"
        }

        ++runTime;

        Info<< "Time = " << runTime.timeName() << nl << endl;

        // --- Pressure-velocity PIMPLE corrector loop
        while (pimple.loop())
        {
            if (pimple.firstIter() || moveMeshOuterCorrectors)
            {
                scalar timeBeforeMeshUpdate = runTime.elapsedCpuTime();

                mesh.update();

                if (mesh.changing())
                {
                    Info<< "Execution time for mesh.update() = "
                        << runTime.elapsedCpuTime() - timeBeforeMeshUpdate
                        << " s" << endl;

                    // Do not apply previous time-step mesh compression flux
                    // if the mesh topology changed
                    if (mesh.topoChanging())
                    {
                        talphaPhi1Corr0.clear();
                    }

                    // Update cellMask field for blocking out hole cells
                    #include "setCellMask.H"
                    #include "setInterpolatedCells.H"
                    #include "correctPhiFaceMask.H"

                    gh = (g & mesh.C()) - ghRef;
                    ghf = (g & mesh.Cf()) - ghRef;


                    mixture.correct();

                    // Make the flux relative to the mesh motion
                    fvc::makeRelative(phi, U);

                }

                if (mesh.changing() && checkMeshCourantNo)
                {
                    #include "meshCourantNo.H"
                }
            }

            if (adjustFringe)
            {
                oversetAdjustPhi(phi, U, zoneIdMass);
            }

            #include "alphaControls.H"
            #include "alphaEqnSubCycle.H"

            rhoPhi *= faceMask;

            mixture.correct();

            #include "UEqn.H"

            // --- Pressure corrector loop
            while (pimple.correct())
            {
                #include "pEqn.H"
            }

            if (pimple.turbCorr())
            {
                turbulence->correct();
            }
        }

        runTime.write();
        // Write Porous Variables
        if( activePorosity && runTime.outputTime() )
        {
            porosity.write();
            porosityIndex.write();
        }

        runTime.printExecutionTime(Info);
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
