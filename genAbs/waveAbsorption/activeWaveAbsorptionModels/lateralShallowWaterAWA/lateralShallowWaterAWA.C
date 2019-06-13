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

#include "lateralShallowWaterAWA.H"
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
    defineTypeNameAndDebug(lateralShallowWaterAWA, 0);
    addToRunTimeSelectionTable(activeWaveAbsorptionModel, lateralShallowWaterAWA, patch);
}
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

vectorField Foam::activeWaveAbsorptionModels::lateralShallowWaterAWA::setVelocity
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

    scalarList Ucalc = -sqrt(9.81/max(initialWaterDepths_, 0.001))*corrLevels;
    scalarList signoU(nPaddles_, 0.0); // Sign of Ucalc

    // Calculate mean velocities on each paddle
    scalarList groupUx (nPaddles_, 0.0);
    scalarList groupUy (nPaddles_, 0.0);
    calcUV
    (
        patch_,
        alphaPatch,
        faceGroup_,
        UPatch.component(0),
        &groupUx,
        UPatch.component(1),
        &groupUy
    );

    // Calculate the mean tangential velocity of water for each paddle
    const scalarList meanTgAngle = meanAngles_ + PI()/2.0;
    scalarList Utg (nPaddles_, 0.0); // Tg velocity
    scalarList Uc (nPaddles_, 0.0); // Correction velocity

    for (label i=0; i<=nPaddles_-1; i++)
    {
        Utg[i] = 
            groupUy[i]*cos(meanAngles_[i]) - groupUx[i]*sin(meanAngles_[i]);
        
        Uc[i] = sqr(Ucalc[i]) - sqr(Utg[i]);
        Uc[i] = sqrt( pos(Uc[i])*Uc[i] );

        signoU[i] = sign(Ucalc[i]);

        // Limit in-flux
        if( i <= nEdgeMin_-1 || i > nPaddles_-1-nEdgeMax_ )
        {
            Uc[i] = Ucalc[i];

            if( pos(Uc[i]) )
            {
                signoU[i] = 0.0;
            }
            else
            {
                signoU[i] = 1.0;
            }
        }
    }

    forAll(UPatch, faceI)    
    {
        UPatch[faceI].component(0) = signoU[faceGroup_[faceI] - 1]*
            Uc[faceGroup_[faceI] - 1]*alphaPatch[faceI]*
            cos(meanAngles_[faceGroup_[faceI] - 1]);
        UPatch[faceI].component(1) = signoU[faceGroup_[faceI] - 1]*
            Uc[faceGroup_[faceI] - 1]*alphaPatch[faceI]*
            sin(meanAngles_[faceGroup_[faceI] - 1]);
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


scalarField Foam::activeWaveAbsorptionModels::lateralShallowWaterAWA::setAlpha()
{
    // Set alpha as zero-gradient
    const volScalarField& alpha =
        mesh_.lookupObject<volScalarField>(alphaName());
    scalarField alphaPatch = 
        alpha.boundaryField()[patch_.index()].patchInternalField();
    return alphaPatch;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::activeWaveAbsorptionModels::lateralShallowWaterAWA::lateralShallowWaterAWA
(
    const fvMesh& mesh,
    const polyPatch& patch,
    const dictionary& dict
)
:
    activeWaveAbsorptionModel(mesh, patch, dict)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::activeWaveAbsorptionModels::lateralShallowWaterAWA::~lateralShallowWaterAWA()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::activeWaveAbsorptionModels::lateralShallowWaterAWA::readDict
(
    const dictionary& overrideDict
)
{
    return activeWaveAbsorptionModel::readDict(overrideDict);
}


// ************************************************************************* //
