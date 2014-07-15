/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2006-2010 OpenCFD Ltd.
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
/*---------------------------------------------------------------------------*\
| Work developed at IH Cantabria       IIIII H   H FFFFF OOOOO AAAAA M   M    |
|                                        I   H   H F     O   O A   A MM MM    |
|   Coder: Pablo Higuera Caubilla        I   HHHHH FFFF  O   O AAAAA M M M    |
|   Bug reports: higuerap@unican.es      I   H   H F     O   O A   A M   M    |
|                                      IIIII H   H F     OOOOO A   A M   M    |
|   -----------------------------------------------------------------------   |
| References:                                                                 |
|                                                                             |
| - Realistic wave generation and active wave absorption for Navier-Stokes    |
|    models: Application to OpenFOAM.                                         |
|    Higuera, P., Lara, J.L. and Losada, I.J. (2013)                          |
|    Coastal Engineering, Vol. 71, 102-118.                                   |
|    http://dx.doi.org/10.1016/j.coastaleng.2012.07.002                       |
|                                                                             |
| - Simulating coastal engineering processes with OpenFOAM                    |
|    Higuera, P., Lara, J.L. and Losada, I.J. (2013)                          |
|    Coastal Engineering, Vol. 71, 119-134.                                   |
|    http://dx.doi.org/10.1016/j.coastaleng.2012.06.002                       |
|                                                                             |
\*---------------------------------------------------------------------------*/

#include "IH_3D_2DAbsorption_InletVelocityFvPatchVectorField.H"
#include "volFields.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::
IH_3D_2DAbsorption_InletVelocityFvPatchVectorField::
IH_3D_2DAbsorption_InletVelocityFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchField<vector>(p, iF),
    nPaddles_(1),
    initialWaterDepths_(List<scalar> (1, -1.0)),
    absorptionDir_(400.0),
    meanAngles_(List<scalar> (1, -1.0)),
    zSpanL_(List<scalar> (1, -1.0)),
    nEdgeMin_(0),
    nEdgeMax_(0),
    allCheck_(false)
{}


Foam::
IH_3D_2DAbsorption_InletVelocityFvPatchVectorField::
IH_3D_2DAbsorption_InletVelocityFvPatchVectorField
(
    const IH_3D_2DAbsorption_InletVelocityFvPatchVectorField& ptf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchField<vector>(ptf, p, iF, mapper),
    nPaddles_(ptf.nPaddles_),
    initialWaterDepths_(ptf.initialWaterDepths_),
    absorptionDir_(ptf.absorptionDir_),
    meanAngles_(ptf.meanAngles_),
    zSpanL_(ptf.zSpanL_),
    nEdgeMin_(ptf.nEdgeMin_),
    nEdgeMax_(ptf.nEdgeMax_),
    allCheck_(ptf.allCheck_)
{}


Foam::
IH_3D_2DAbsorption_InletVelocityFvPatchVectorField::
IH_3D_2DAbsorption_InletVelocityFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchField<vector>(p, iF, dict),
    nPaddles_(dict.lookupOrDefault<label>("nPaddles", 1)),
    initialWaterDepths_( 
    	dict.lookupOrDefault("initialWaterDepths", List<scalar> (1, -1.0)) ),
    absorptionDir_(dict.lookupOrDefault<scalar>("absorptionDir", 400.0)),
    meanAngles_( dict.lookupOrDefault("meanAngles", List<scalar> (1, -1.0)) ),
    zSpanL_( dict.lookupOrDefault("zSpanL", List<scalar> (1, -1.0)) ),
    nEdgeMin_(dict.lookupOrDefault<label>("nEdgeMin", 0)),
    nEdgeMax_(dict.lookupOrDefault<label>("nEdgeMax", 0)),
    allCheck_(dict.lookupOrDefault<bool>("allCheck", false ))
{}


Foam::
IH_3D_2DAbsorption_InletVelocityFvPatchVectorField::
IH_3D_2DAbsorption_InletVelocityFvPatchVectorField
(
    const IH_3D_2DAbsorption_InletVelocityFvPatchVectorField& ptf
)
:
    fixedValueFvPatchField<vector>(ptf),
    nPaddles_(ptf.nPaddles_),
    initialWaterDepths_(ptf.initialWaterDepths_),
    absorptionDir_(ptf.absorptionDir_),
    meanAngles_(ptf.meanAngles_),
    zSpanL_(ptf.zSpanL_),
    nEdgeMin_(ptf.nEdgeMin_),
    nEdgeMax_(ptf.nEdgeMax_),
    allCheck_(ptf.allCheck_)
{}


Foam::
IH_3D_2DAbsorption_InletVelocityFvPatchVectorField::
IH_3D_2DAbsorption_InletVelocityFvPatchVectorField
(
    const IH_3D_2DAbsorption_InletVelocityFvPatchVectorField& ptf,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchField<vector>(ptf, iF),
    nPaddles_(ptf.nPaddles_),
    initialWaterDepths_(ptf.initialWaterDepths_),
    absorptionDir_(ptf.absorptionDir_),
    meanAngles_(ptf.meanAngles_),
    zSpanL_(ptf.zSpanL_),
    nEdgeMin_(ptf.nEdgeMin_),
    nEdgeMax_(ptf.nEdgeMax_),
    allCheck_(ptf.allCheck_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::IH_3D_2DAbsorption_InletVelocityFvPatchVectorField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    // PHC //
    Info << "3D_2D Absorption BC on patch " << this->patch().name() << endl;

    // 3D Variables
    const vector cMin = gMin(patch().patch().localPoints());
    const vector cMax = gMax(patch().patch().localPoints());
    const vector cSpan = cMax - cMin;

    scalar dMin = 0.0;
    scalar dSpan = 0.0;
    const scalarField patchD = patchDirection( cSpan, &dMin, &dSpan );

    // Variables & constants
    const fvMesh& mesh = dimensionedInternalField().mesh();
	const word& patchName = this->patch().name();
	const label patchID = mesh.boundaryMesh().findPatchID(patchName);
    const label nF = patch().faceCells().size();

    const volScalarField& alpha = 
    	db().lookupObject<volScalarField>(alphaName());
    const volVectorField& U = db().lookupObject<volVectorField>("U");

    const scalarField alphaCell = 
    	alpha.boundaryField()[patchID].patchInternalField();
    const vectorField UCell = U.boundaryField()[patchID].patchInternalField();
    scalarField UzCell = UCell.component(2);

    scalarField patchUc = Foam::scalarField(nF, 0.0); // Correction velocity
    scalarField patchVc = Foam::scalarField(nF, 0.0); // Correction velocity

    const scalar g = 9.81;

    // Check for errors - Just the first time
    if (!allCheck_)
    {
        // Check if the value of nPaddles is correct for the number of columns
        if (nPaddles_ < 1)
        {
            FatalError
                << "Check nPaddles value."
                << exit(FatalError);
        }

        if ( nPaddles_ > 1 )
        {
            nPaddles_ = decreaseNPaddles( nPaddles_, patchD, dMin, dSpan );
            reduce(nPaddles_, minOp<label>());
        }

        // Check nEdges
        if ( nEdgeMin_ < 0 || nEdgeMax_ < 0 )
        {
            FatalError
                << "Check nEdgeMin/Max value."
                << exit(FatalError); 
        }

        if ( nEdgeMin_ + nEdgeMax_ > nPaddles_ )
        {
            FatalError
                << "Check: nEdges > nPaddles."
                << exit(FatalError); 
        }
    }

    // Grouping part
    labelList faceGroup = Foam::labelList(nF, 0);
    scalarList dBreakPoints = Foam::scalarList(nPaddles_+1, dMin); 
    scalarList xGroup = Foam::scalarList(nPaddles_, 0.0);
    scalarList yGroup = Foam::scalarList(nPaddles_, 0.0);

    for (label i=0; i<nPaddles_; i++)
    {
    	// Breakpoints, X & Y centre of the paddles
        dBreakPoints[i+1] = dMin + dSpan/(nPaddles_)*(i+1);
        xGroup[i] = cMin[0] + cSpan[0]/(2.0*nPaddles_) + cSpan[0]/(nPaddles_)*i;
        yGroup[i] = cMin[1] + cSpan[1]/(2.0*nPaddles_) + cSpan[1]/(nPaddles_)*i;
    }

    forAll(patchD, patchCells) 
    {
        for (label i=0; i<nPaddles_; i++)
        {
            if ( (patchD[patchCells]>=dBreakPoints[i]) && 
            	(patchD[patchCells]<dBreakPoints[i+1]) )
            {
                faceGroup[patchCells] = i+1; // Group of each face
                continue;
            }
        }      
    }

    if (!allCheck_)
    {
	    // Calculate Z bounds of the faces
	    scalarField zSup, zInf;
		faceBoundsZ( &zSup, &zInf );
        // Z-span paddlewise
        zSpanL_ = zSpanList( faceGroup, zInf, zSup );
        //inlinePrint( "Z-span paddlewise", zSpanL_ );

        // Calculate initial water depth
        if ( gMin(initialWaterDepths_) < 0.0 )
        {
            initialWaterDepths_ = calcWL( alphaCell, faceGroup, zSpanL_ );
            inlinePrint( "Initial water depths for absorption", 
            	initialWaterDepths_ );
        }

        // Absorption direction part
        meanAngles_ = initialWaterDepths_;

        // Check if absorption is directional
        if ( absorptionDir_ > 360.0 ) // Automatic
        {
            meanAngles_ = meanPatchDirs( faceGroup );
        }
        else // Fixed
        {
            meanAngles_ = absorptionDir_*PI()/180.0;
        }
        //inlinePrint( "Paddle angles", meanAngles_ );

        allCheck_ = true;
    }

    // Calculate water measured levels
    scalarList measuredLevels = calcWL( alphaCell, faceGroup, zSpanL_ );

    // Correction velocity: U = -sqrt(g/h)*corrL
    scalarList corrLevels = measuredLevels - initialWaterDepths_;
    inlinePrint( "Correction Levels", corrLevels );

    scalarList Uc = -sqrt(g/max(initialWaterDepths_,0.1))*corrLevels;

    // Limit in-flux
    for(label i=0; i<=nPaddles_-1; i++)
    {
        if( (i <= nEdgeMin_-1) || (i > nPaddles_-1-nEdgeMax_) )
        {
            if( pos(Uc[i]) )
            {
                Uc[i] = 0.0;
            }
        }
    }

    forAll(patchUc, cellIndex)    
    {
        patchUc[cellIndex] = pos(alphaCell[cellIndex]-0.9)*
        	Uc[faceGroup[cellIndex]-1]*cos(meanAngles_[faceGroup[cellIndex]-1]);
        patchVc[cellIndex] = pos(alphaCell[cellIndex]-0.9)*
        	Uc[faceGroup[cellIndex]-1]*sin(meanAngles_[faceGroup[cellIndex]-1]);
        UzCell[cellIndex] = pos(alphaCell[cellIndex]-0.9)*UzCell[cellIndex];
    }

    const vectorField n1 = Foam::vectorField(nF, vector(1.0, 0.0, 0.0));
    const vectorField n2 = Foam::vectorField(nF, vector(0.0, 1.0, 0.0));
    const vectorField n3 = Foam::vectorField(nF, vector(0.0, 0.0, 1.0));

    operator==(n1*patchUc + n2*patchVc + n3*UzCell); 

    fixedValueFvPatchField<vector>::updateCoeffs();

}


void Foam::IH_3D_2DAbsorption_InletVelocityFvPatchVectorField::
write(Ostream& os) const
{
    fvPatchField<vector>::write(os);

    os.writeKeyword("absorptionDir") << absorptionDir_ << 
        token::END_STATEMENT << nl;
    os.writeKeyword("nPaddles") << nPaddles_ << token::END_STATEMENT << nl;

    initialWaterDepths_.writeEntry("initialWaterDepths", os);
    meanAngles_.writeEntry("meanAngles", os);
    zSpanL_.writeEntry("zSpanL", os);

    os.writeKeyword("nEdgeMin") << nEdgeMin_ << token::END_STATEMENT << nl;
    os.writeKeyword("nEdgeMax") << nEdgeMax_ << token::END_STATEMENT << nl;
    os.writeKeyword("allCheck") << allCheck_ << token::END_STATEMENT << nl;

    writeEntry("value", os);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
   makePatchTypeField
   (
       fvPatchVectorField,
       IH_3D_2DAbsorption_InletVelocityFvPatchVectorField
   );
}


// ************************************************************************* //
