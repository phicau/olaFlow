/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2017 OpenFOAM Foundation
     \\/     M anipulation  |
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
    olaDyMFoam

Group
    grpMultiphaseSolvers grpMovingMeshSolvers

Description
    Solver for 2 incompressible, isothermal immiscible fluids using a VOF
    (volume of fluid) phase-fraction based interface capturing approach,
    with optional mesh motion and mesh topology changes including adaptive
    re-meshing.

    OLADYMFOAM solves the VARANS equations.

\*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*\
| olaFoam Project                                       ll                    |
|                                                       l l                   |
|   Coder: Pablo Higuera Caubilla                 ooo   l l     aa            |
|   Bug reports: phicau@gmail.com                o   o  l l    a  a           |
|                                                o   o  ll   l a  aa  aa      |
|                                                 ooo    llll   aa  aa        |
|                                                                             |
|                                                FFFFF OOOOO AAAAA M   M      |
|                                                F     O   O A   A MM MM      |
|                                                FFFF  O   O AAAAA M M M      |
|                                                F     O   O A   A M   M      |
|                                                F     OOOOO A   A M   M      |
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
|    Coastal Engineering, Vol. 83, 259–270.                                   |
|    http://dx.doi.org/10.1016/j.coastaleng.2013.09.002                       |
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
#include "CorrectPhi.H"
#include "fixedFluxPressureFvPatchScalarField.H"
#include "localEulerDdtScheme.H"
#include "fvcSmooth.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createDynamicFvMesh.H"
    #include "initContinuityErrs.H"

    pimpleControl pimple(mesh);

    #include "createControls.H"
    #include "createRDeltaT.H"
    #include "createFields.H"
    #include "createMRF.H"
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

    #include "correctPhi.H"
    #include "createUf.H"

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
        #include "readControls.H"

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

        runTime++;

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

                    gh = (g & mesh.C()) - ghRef;
                    ghf = (g & mesh.Cf()) - ghRef;
                }

                if (mesh.changing() && correctPhi)
                {
                    // Calculate absolute flux from the mapped surface velocity
                    phi = mesh.Sf() & Uf;

                    #include "correctPhi.H"

                    // Make the flux relative to the mesh motion
                    fvc::makeRelative(phi, U);

                    mixture.correct();
                }

                if (mesh.changing() && checkMeshCourantNo)
                {
                    #include "meshCourantNo.H"
                }
            }

            #include "alphaControls.H"
            #include "alphaEqnSubCycle.H"

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

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
