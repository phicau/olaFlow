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

#include "activeWaveAbsorptionModel.H"
#include "fvMesh.H"
#include "polyPatch.H"
#include "uniformDimensionedFields.H"
#include "volFields.H"
#include "fvPatchFields.H"
#include "scalarList.H"

using namespace Foam::constant;
using namespace Foam;

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(activeWaveAbsorptionModel, 0);
    defineRunTimeSelectionTable(activeWaveAbsorptionModel, patch);
}

Foam::word Foam::activeWaveAbsorptionModel::modelID(const word& patchName)
{
    const word name = "AWA";
    return name + '.' + patchName;
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void Foam::activeWaveAbsorptionModel::initialise()
{
    // Read variables
    theory_ = dict_.lookupOrDefault<word>("theory", "");
    nPaddles_ = dict_.lookupOrDefault<label>("nPaddles", 1);
    initialWaterDepths_ = 
        dict_.lookupOrDefault<List<scalar> >
            ("initialWaterDepths", List<scalar> (1, -1.));
    absDir_ = dict_.lookupOrDefault<scalar>("absDir", 666.);
    uCurrent_ = dict_.lookupOrDefault<vector>("uCurrent", vector::zero);
    nEdgeMin_ = dict_.lookupOrDefault<label>("nEdgeMin", 0);
    nEdgeMax_ = dict_.lookupOrDefault<label>("nEdgeMax", 0);
    zSpanL_ = 
        dict_.lookupOrDefault<List<scalar> >("zSpanL", List<scalar> (1, -1.));
    meanAngles_ = 
        dict_.lookupOrDefault<List<scalar> >("meanAngles", List<scalar> (1, -1.));

    // Create faceGroup_
    const vector cMin = gMin(patch_.localPoints());
    const vector cMax = gMax(patch_.localPoints());
    const vector cSpan = cMax - cMin;

    scalar dMin = 0.0;
    scalar dSpan = 0.0;
    const scalarField patchD = patchDirection( patch_, cSpan, &dMin, &dSpan );

    const label nF = patch_.faceCells().size();
    faceGroup_ = Foam::labelList(nF, 0);
    scalarList dBreakPoints = Foam::scalarList(nPaddles_ + 1, dMin); 
    xPaddleGroup_ = Foam::scalarList(nPaddles_, 0.0);
    yPaddleGroup_ = Foam::scalarList(nPaddles_, 0.0);

    for (label i=0; i<nPaddles_; i++)
    {
        // Breakpoints, X & Y centre of the paddles
        dBreakPoints[i + 1] = dMin + dSpan/(nPaddles_)*(i + 1);
        xPaddleGroup_[i] = cMin[0] + cSpan[0]/(2.0*nPaddles_)
            + cSpan[0]/(nPaddles_)*i;
        yPaddleGroup_[i] = cMin[1] + cSpan[1]/(2.0*nPaddles_)
            + cSpan[1]/(nPaddles_)*i;
    }

    forAll(patchD, patchCells) 
    {
        for (label i=0; i<nPaddles_; i++)
        {
            if ( (patchD[patchCells]>=dBreakPoints[i]) && 
                (patchD[patchCells]<dBreakPoints[i + 1]) )
            {
                faceGroup_[patchCells] = i + 1; // Group of each face
                continue;
            }
        }      
    }
}


scalarList Foam::activeWaveAbsorptionModel::calcWaterLevels()
{
    scalarList groupTotalArea (nPaddles_, 0.0); // Paddle total area
    scalarList groupWaterArea (nPaddles_, 0.0); // Paddle wet area
    scalarList heights (nPaddles_, 0.0);
    
    // Patch variables
	const word& patchName = patch_.name();
	const label patchID = mesh_.boundaryMesh().findPatchID(patchName);

    // Fields
    const volScalarField& alpha =
        mesh_.lookupObject<volScalarField>(alphaName());
    const scalarField alphaCell =
    	alpha.boundaryField()[patchID].patchInternalField();

    const scalarField faceSurface = mag(patch_.faceAreas());

    forAll( faceSurface, faceI )
    {
        groupTotalArea[faceGroup_[faceI] - 1] += faceSurface[faceI];
        groupWaterArea[faceGroup_[faceI] - 1] += 
            faceSurface[faceI]*alphaCell[faceI];
    }

    for (label i = 0; i <= gMax(faceGroup_) - 1; i++)
    {
        reduce(groupTotalArea[i], sumOp<scalar>());
        reduce(groupWaterArea[i], sumOp<scalar>());
        // Free surface elevation at each paddle
        heights[i] = groupWaterArea[i]/groupTotalArea[i]*zSpanL_[i];
    }

    return heights;
}


scalarField Foam::activeWaveAbsorptionModel::setAlpha()
{
    // Patch variables
	const word& patchName = patch_.name();
	const label patchID = mesh_.boundaryMesh().findPatchID(patchName);

    // Fields
    const volScalarField& alpha =
        mesh_.lookupObject<volScalarField>(alphaName());
    const scalarField alphaCell =
    	alpha.boundaryField()[patchID].patchInternalField();
    return alphaCell;

    /*forAll(alpha_, facei)
    {
        const label paddlei = faceToPaddle_[facei];
        const scalar paddleCalc = level[paddlei];

        const scalar zMin0 = zMin_[facei] - zMin0_;
        const scalar zMax0 = zMax_[facei] - zMin0_;

        if (zMax0 < paddleCalc)
        {
            alpha_[facei] = 1.0;
        }
        else if (zMin0 > paddleCalc)
        {
            alpha_[facei] = 0.0;
        }
        else
        {
            scalar dz = paddleCalc - zMin0;
            alpha_[facei] = dz/(zMax0 - zMin0);
        }
    }*/
}


vectorField Foam::activeWaveAbsorptionModel::setVelocity
(
    bool isPureAWA,
    scalarList refWaterDepths
)
{
    // Patch variables
	const word& patchName = patch_.name();
	const label patchID = mesh_.boundaryMesh().findPatchID(patchName);

    // Fields
    const volVectorField& U = mesh_.lookupObject<volVectorField>("U");
    const vectorField UCell =
    	U.boundaryField()[patchID].patchInternalField();
    return UCell;
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::activeWaveAbsorptionModel::activeWaveAbsorptionModel
(
    const fvMesh& mesh,
    const polyPatch& patch,
    const dictionary& dict
)
:
    IOdictionary
    (
        IOobject
        (
            modelID(patch.name()),
            Time::timeName(mesh.time().startTime().value()),
            "uniform",
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        )
    ),
    mesh_(mesh),
    patch_(patch),
    dict_(dict),
    theory_(""),
    nPaddles_(1),
    nEdgeMin_(0),
    nEdgeMax_(0),
    xPaddleGroup_(List<scalar> (1, 0.0)),
    yPaddleGroup_(List<scalar> (1, 0.0)),
    absDir_(666),
    uCurrent_(vector::zero),
    faceGroup_(List<label> (1, -1)),
    meanAngles_(List<scalar> (1, 0.0)),
    zSpanL_(List<scalar> (1, 0.0)),
    initialWaterDepths_(List<scalar> (1, -1.0)),
    waterDepths_(List<scalar> (1, -1.0))
{
    initialise();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::activeWaveAbsorptionModel::~activeWaveAbsorptionModel()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::activeWaveAbsorptionModel::readDict(const dictionary& overrideDict)
{
    initialise();

    return true;
}


void Foam::activeWaveAbsorptionModel::correct()
{
    Info<< "Updating " << type() << " absorption model for patch "
        << patch_.name() << endl;

    // Update the calculated water levels
    waterDepths_ = calcWaterLevels();

    Info << "Water depths: " << waterDepths_ << endl;

}

const scalarList& Foam::activeWaveAbsorptionModel::xPaddleCentres() const
{
    return xPaddleGroup_;
}

const scalarList& Foam::activeWaveAbsorptionModel::yPaddleCentres() const
{
    return yPaddleGroup_;
}

void Foam::activeWaveAbsorptionModel::info(Ostream& os) const
{
    os  << "AWA model: patch " << patch_.name() << nl
        << "    Theory: " << type() << nl
        << "    Number of paddles: " << nPaddles_ << nl
        << "    Reference water depth: " << initialWaterDepths_ << nl;
}

// ************************************************************************* //
