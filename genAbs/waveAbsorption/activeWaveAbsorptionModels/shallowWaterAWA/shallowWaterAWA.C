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

\*---------------------------------------------------------------------------*/

#include "shallowWaterAWA.H"
#include "volFields.H"
#include "fvPatchFields.H"
#include "addToRunTimeSelectionTable.H"
#include "scalarList.H"

using namespace Foam::constant;
using namespace Foam;

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace activeWaveAbsorptionModels
{
    defineTypeNameAndDebug(shallowWaterAWA, 0);
    addToRunTimeSelectionTable(activeWaveAbsorptionModel, shallowWaterAWA, patch);
}
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

vectorField Foam::activeWaveAbsorptionModels::shallowWaterAWA::setVelocity
(
    bool isPureAWA,
    scalarList refWaterDepths
)
{
    // Get local U
    const volVectorField& U = mesh_.lookupObject<volVectorField>("U");
    vectorField UPatch = U.boundaryField()[patch_.index()].patchInternalField();
    const volScalarField& alpha = 
        mesh_.lookupObject<volScalarField>(alphaName());
    scalarField alphaPatch = 
        alpha.boundaryField()[patch_.index()].patchInternalField();

    // Auxiliar geometry
    const vector cMin = gMin(patch_.localPoints());
    const scalarField patchHeight = patch_.faceCentres().component(2) - cMin[2];

    // Correction velocity: U = -sqrt(g/h)*corrL
    scalarList corrLevels;
    if (isPureAWA)
    {
        corrLevels = waterDepths_ - initialWaterDepths_;
    }
    else
    {
        corrLevels = waterDepths_ - refWaterDepths;
    }
    
    inlinePrint( "Correction Levels", corrLevels );

    scalarList Uc = -sqrt(9.81/max(initialWaterDepths_,0.001))*corrLevels;

    // Limit in-flux
    for(label i = 0; i <= nPaddles_ - 1; i++)
    {
        if( (i <= nEdgeMin_ - 1) || (i > nPaddles_ - 1 - nEdgeMax_) )
        {
            if( pos(Uc[i]))
            {
                Uc[i] = 0.0;
            }
        }
    }

    forAll(UPatch, faceI)    
    {
        UPatch[faceI].component(0) =
        	Uc[faceGroup_[faceI] - 1]*cos(meanAngles_[faceGroup_[faceI] - 1]);
        UPatch[faceI].component(1) =
        	Uc[faceGroup_[faceI] - 1]*sin(meanAngles_[faceGroup_[faceI] - 1]);
    }

    // Clip based on VOF value
    UPatch *= pos(alphaPatch - 0.9);

    if (isPureAWA)
    {
        // Zero-gradient in Z and current velocity
        UPatch += uCurrent_*
            pos(alphaPatch - 0.5)*neg(patchHeight - min(waterDepths_));
    }
    else
    {
        // Just X and Y corrections
        #if OFVERSION >= 400
            UPatch.component(2).ref() = 0;
        #else
            UPatch.component(2) = 0*UPatch.component(2);
        #endif
    }

    return UPatch;
}


scalarField Foam::activeWaveAbsorptionModels::shallowWaterAWA::setAlpha()
{
    // Set alpha as zero-gradient
    const volScalarField& alpha =
        mesh_.lookupObject<volScalarField>(alphaName());
    scalarField alphaPatch = 
        alpha.boundaryField()[patch_.index()].patchInternalField();
    return alphaPatch;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::activeWaveAbsorptionModels::shallowWaterAWA::shallowWaterAWA
(
    const fvMesh& mesh,
    const polyPatch& patch,
    const dictionary& dict
)
:
    activeWaveAbsorptionModel(mesh, patch, dict)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::activeWaveAbsorptionModels::shallowWaterAWA::~shallowWaterAWA()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::activeWaveAbsorptionModels::shallowWaterAWA::readDict
(
    const dictionary& overrideDict
)
{
    return activeWaveAbsorptionModel::readDict(overrideDict);
}


// ************************************************************************* //
