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

#include "extendedRangeAWA.H"
#include "volFields.H"
#include "fvPatchFields.H"
#include "addToRunTimeSelectionTable.H"
#include "scalarList.H"

#include "waveFun.H"

using namespace Foam::constant;
using namespace Foam;

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace activeWaveAbsorptionModels
{
    defineTypeNameAndDebug(extendedRangeAWA, 0);
    addToRunTimeSelectionTable(activeWaveAbsorptionModel, extendedRangeAWA, patch);
}
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

vectorField Foam::activeWaveAbsorptionModels::extendedRangeAWA::setVelocity
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

    // Calculate wave parameters
    scalar waveK = 2.0*PI()/waveLength_;
    scalar waveOmega = 2.0*PI()/wavePeriod_;

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

    // Limit in-flux
    scalarList iF = -sqrt(9.81/max(initialWaterDepths_,0.001))*corrLevels;
    for(label i = 0; i <= nPaddles_ - 1; i++)
    {
        if( (i <= nEdgeMin_ - 1) || (i > nPaddles_ - 1 - nEdgeMax_) )
        {
            if( pos(iF[i]))
            {
                iF[i] = 0.0;
            }
        }
        else
        {
            iF[i] = 1.0;
        }
    }

    forAll(UPatch, faceI)    
    {
        scalar velAux = iF[faceGroup_[faceI] - 1]*
            StokesIFun::U
            (
                -2.*corrLevels[faceGroup_[faceI] - 1],
                initialWaterDepths_[faceGroup_[faceI] - 1], 
                waveK, 
                0., // x
                0., // ky
                0., // y
                waveOmega, 
                0., // t
                0., // phase
                patchHeight[faceI]
            );

        UPatch[faceI].component(0) =
        	velAux*cos(meanAngles_[faceGroup_[faceI] - 1]);
        UPatch[faceI].component(1) =
        	velAux*sin(meanAngles_[faceGroup_[faceI] - 1]);
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


scalarField Foam::activeWaveAbsorptionModels::extendedRangeAWA::setAlpha()
{
    // Set alpha as zero-gradient
    const volScalarField& alpha =
        mesh_.lookupObject<volScalarField>(alphaName());
    scalarField alphaPatch = 
        alpha.boundaryField()[patch_.index()].patchInternalField();
    return alphaPatch;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::activeWaveAbsorptionModels::extendedRangeAWA::extendedRangeAWA
(
    const fvMesh& mesh,
    const polyPatch& patch,
    const dictionary& dict
)
:
    activeWaveAbsorptionModel(mesh, patch, dict),
    wavePeriod_(-1.),
    waveLength_(-1.)
{
    wavePeriod_ = dict.lookupOrDefault<scalar>("wavePeriod", -1.);

    if (wavePeriod_ < 0.)
    {
        FatalError
            << "wavePeriod not defined. Required by extendedRangeAWA."
            << exit(FatalError);
    }
    
    // Calculate wave parameters
    waveLength_ = 
        StokesIFun::waveLength (max(initialWaterDepths_), wavePeriod_);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::activeWaveAbsorptionModels::extendedRangeAWA::~extendedRangeAWA()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::activeWaveAbsorptionModels::extendedRangeAWA::readDict
(
    const dictionary& overrideDict
)
{
    return activeWaveAbsorptionModel::readDict(overrideDict);
}

void Foam::activeWaveAbsorptionModels::extendedRangeAWA::info(Ostream& os) const
{
    os  << "AWA model: patch " << patch_.name() << nl
        << "    Theory: " << type() << nl
        << "    Number of paddles: " << nPaddles_ << nl
        << "    Reference water depth: " << initialWaterDepths_ << nl
        << "    Wave period: " << wavePeriod_ << nl;
}

// ************************************************************************* //
