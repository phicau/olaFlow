/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2016 OpenCFD Ltd.
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

Class
    Foam::shallowWaterAWA

Description

\*---------------------------------------------------------------------------*/

#ifndef activeWaveAbsorptionModels_shallowWaterAWA_H
#define activeWaveAbsorptionModels_shallowWaterAWA_H

#include "activeWaveAbsorptionModel.H"
#include "fvMesh.H"
#include "scalarField.H"
#include "vectorField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace activeWaveAbsorptionModels
{

/*---------------------------------------------------------------------------*\
                    Class shallowWaterAWA Declaration
\*---------------------------------------------------------------------------*/

class shallowWaterAWA
:
    public activeWaveAbsorptionModel
{
protected:

    // Protected Member Functions

        //- Set the alpha field
        virtual scalarField setAlpha();

        //- Set the velocity field
        virtual vectorField setVelocity
            (bool isPureAWA, scalarList refWaterDepths);

public:

    //- Runtime type information
    TypeName("shallowWaterAWA");

    //- Constructor
    shallowWaterAWA
    (
        const fvMesh& mesh,
        const polyPatch& patch,
        const dictionary& dict
    );

    //- Destructor
    virtual ~shallowWaterAWA();


    // Public Member Functions

        //- Read from dictionary
        virtual bool readDict(const dictionary& overrideDict);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace activeWaveAbsorptionModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
