/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2017 OpenFOAM Foundation
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
    setWaves

Description
    Applies wave models to the entire domain for case initialisation using
    level sets for second-order accuracy.

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
| Reference:                                                                  |
|                                                                             |
| - Enhancing active wave absorption in RANS models.                          |
|    Higuera, P. (2020)                                                       |
|    Applied Ocean Research, Vol. 94                                          |
|    https://doi.org/10.1016/j.apor.2019.102000                               |
|                                                                             |
\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "pointFields.H"
#include "timeSelector.H"
#include "wallDist.H"

#include "waveFun.H"
#include "memberFun.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    timeSelector::addOptions(false, false);

    Foam::argList::addOption
    (
        "U",
        "name",
        "name of the velocity field, default is \"U\""
    );

    Foam::argList::addOption
    (
        "alpha",
        "name",
        "name of the volume fraction field, default is \"" + alphaName() + "\""
    );

    #include "setRootCase.H"
    #include "createTime.H"

    instantList timeDirs = timeSelector::selectIfPresent(runTime, args);

    #include "createMesh.H"
    #include "readGravitationalAcceleration.H"
    #include "defaultVars.H"

    forAll(timeDirs, timeI)
    {
        runTime.setTime(timeDirs[timeI], timeI);

        Info<< "Time = " << runTime.timeName() << nl << endl;

        mesh.readUpdate();

        // Read the fields which are to be set
        volScalarField alpha
        (
            IOobject
            (
                args.optionFound("alpha") ? args["alpha"] : alphaName(),
                runTime.timeName(),
                mesh,
                IOobject::MUST_READ
            ),
            mesh
        );
        volVectorField U
        (
            IOobject
            (
                args.optionFound("U") ? args["U"] : "U",
                runTime.timeName(),
                mesh,
                IOobject::MUST_READ
            ),
            mesh
        );

        surfaceScalarField phi
        (
            IOobject
            (
                "phi",
                runTime.timeName(),
                mesh,
                IOobject::READ_IF_PRESENT,
                IOobject::NO_WRITE
            ),
            #if OFVERSION >= 400
            fvc::flux(U)
            #else
            fvc::interpolate(U) & mesh.Sf()
            #endif
        );

        Info << "Preparing wave field" << nl << endl;
        #include "prepareWaves.H"

        // Set initial water depth
        Info << "Setting initial water depth" << nl << endl;
        scalarList calculatedLevel (xComp.size(), waterDepth_);
        #include "alphaCycle.H"

        // Solving dummy equations to force wave BCs to work
        Info << "Loading boundary conditions" << nl << endl;
        fvScalarMatrix alphaEqn 
        ( 
            fvm::ddt(alpha) 
        );
        alphaEqn.solve();

        fvVectorMatrix UEqn 
        ( 
            fvm::ddt(U) 
        );
        UEqn.solve();

        // Correct BCs
        alpha.correctBoundaryConditions();
        U.correctBoundaryConditions();

        // Calculate wave field
        Info << "Calculating wave field" << nl << endl;
        if ( waveType_ == "regular" )
        {
            #include "calculatedLevelRegular.H"
        }
        else if ( waveType_ == "solitary" )
        {
            #include "calculatedLevelSolitary.H"
        }
        else if ( waveType_ == "irregular" )
        {
            #include "calculatedLevelIrregular.H"
        }
        else if ( waveType_ == "wavemaker" )
        {
            #include "calculatedLevelEta.H"
        }
        else if ( waveType_ == "current" )
        {
            #include "calculatedLevelCurrent.H"
        }

        // Alpha cycle
        #include "alphaCycle.H"

        // Velocity cycle
        #include "velocityCycle.H"

        // Correct BCs
        alpha.correctBoundaryConditions();
        U.correctBoundaryConditions();

        // Output
        Info<< "Writing " << alpha.name() << nl;
        alpha.write();
        Info<< "Writing " << U.name() << nl << endl;
        U.write();
    }

    Info<< "End\n" << endl;

    return 0;
}

// ************************************************************************* //
