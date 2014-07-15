/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

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

Application
    ihDyMFoam

Description
    Solver for 2 incompressible, isothermal immiscible fluids using a VOF
    (volume of fluid) phase-fraction based interface capturing approach,
    with optional mesh motion and mesh topology changes including adaptive
    re-meshing.

    IHDYMFOAM solves the VARANS equations.

\*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*\
| Work developed at IH Cantabria       IIIII H   H FFFFF OOOOO AAAAA M   M    |
|                                        I   H   H F     O   O A   A MM MM    |
|   Coder: Pablo Higuera Caubilla        I   HHHHH FFFF  O   O AAAAA M M M    |
|   Bug reports: higuerap@unican.es      I   H   H F     O   O A   A M   M    |
|                                      IIIII H   H F     OOOOO A   A M   M    |
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
#include "MULES.H"
#include "subCycle.H"
#include "interfaceProperties.H"
#include "twoPhaseMixture.H"
#include "turbulenceModel.H"
#include "pimpleControl.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createDynamicFvMesh.H"
    #include "initContinuityErrs.H"
    #include "createFields.H"
    #include "readTimeControls.H"

    pimpleControl pimple(mesh);

    surfaceScalarField phiAbs("phiAbs", phi);
    fvc::makeAbsolute(phiAbs, U);

    #include "correctPhi.H"
    #include "CourantNo.H"
    #include "setInitialDeltaT.H"

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
    Info<< "\nStarting time loop\n" << endl;

    while (runTime.run())
    {
        #include "readControls.H"
        #include "alphaCourantNo.H"
        #include "CourantNo.H"

        #include "setDeltaT.H"

        runTime++;

        Info<< "Time = " << runTime.timeName() << nl << endl;

        scalar timeBeforeMeshUpdate = runTime.elapsedCpuTime();

        {
            // Calculate the relative velocity used to map the relative flux phi
            volVectorField Urel("Urel", U);

            if (mesh.moving())
            {
                Urel -= fvc::reconstruct(fvc::meshPhi(U));
            }

            // Do any mesh changes
            mesh.update();
        }

        if (mesh.changing())
        {
            Info<< "Execution time for mesh.update() = "
                << runTime.elapsedCpuTime() - timeBeforeMeshUpdate
                << " s" << endl;

            gh = g & mesh.C();
            ghf = g & mesh.Cf();
        }

        if (mesh.changing() && correctPhi)
        {
            #include "correctPhi.H"
        }

        if (mesh.changing() && checkMeshCourantNo)
        {
            #include "meshCourantNo.H"
        }

        twoPhaseProperties.correct();

        #include "alphaEqnSubCycle.H"

        // --- Pressure-velocity PIMPLE corrector loop
        while (pimple.loop())
        {
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
