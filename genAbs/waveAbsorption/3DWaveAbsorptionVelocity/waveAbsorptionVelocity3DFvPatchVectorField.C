/*---------------------------------------------------------------------------*\
License
    This file is part of olaFlow Project.

    olaFlow is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    olaFlow is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with olaFlow.  If not, see <http://www.gnu.org/licenses/>.

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

#include "waveAbsorptionVelocity3DFvPatchVectorField.H"
#include "volFields.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::
waveAbsorptionVelocity3DFvPatchVectorField::
waveAbsorptionVelocity3DFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchField<vector>(p, iF),
    nPaddles_(1),
    initialWaterDepths_(List<scalar> (1, -1.0)),
    meanAngles_(List<scalar> (1, -1.0)),
    zSpanL_(List<scalar> (1, -1.0)),
    nEdgeMin_(0),
    nEdgeMax_(0),
    allCheck_(false)
{}


Foam::
waveAbsorptionVelocity3DFvPatchVectorField::
waveAbsorptionVelocity3DFvPatchVectorField
(
    const waveAbsorptionVelocity3DFvPatchVectorField& ptf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchField<vector>(ptf, p, iF, mapper),
    nPaddles_(ptf.nPaddles_),
    initialWaterDepths_(ptf.initialWaterDepths_),
    meanAngles_(ptf.meanAngles_),
    zSpanL_(ptf.zSpanL_),
    nEdgeMin_(ptf.nEdgeMin_),
    nEdgeMax_(ptf.nEdgeMax_),
    allCheck_(ptf.allCheck_)
{}


Foam::
waveAbsorptionVelocity3DFvPatchVectorField::
waveAbsorptionVelocity3DFvPatchVectorField
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
    meanAngles_( dict.lookupOrDefault("meanAngles", List<scalar> (1, -1.0)) ),
    zSpanL_( dict.lookupOrDefault("zSpanL", List<scalar> (1, -1.0)) ),
    nEdgeMin_(dict.lookupOrDefault<label>("nEdgeMin", 0)),
    nEdgeMax_(dict.lookupOrDefault<label>("nEdgeMax", 0)),
    allCheck_(dict.lookupOrDefault<bool>("allCheck", false ))
{
    
    fixedValueFvPatchField<vector>::operator==(this->patchInternalField());
    fixedValueFvPatchField<vector>::updateCoeffs();
}


#if OFFLAVOUR == 3 && OFVERSION >= 900
#else
Foam::
waveAbsorptionVelocity3DFvPatchVectorField::
waveAbsorptionVelocity3DFvPatchVectorField
(
    const waveAbsorptionVelocity3DFvPatchVectorField& ptf
)
:
    fixedValueFvPatchField<vector>(ptf),
    nPaddles_(ptf.nPaddles_),
    initialWaterDepths_(ptf.initialWaterDepths_),
    meanAngles_(ptf.meanAngles_),
    zSpanL_(ptf.zSpanL_),
    nEdgeMin_(ptf.nEdgeMin_),
    nEdgeMax_(ptf.nEdgeMax_),
    allCheck_(ptf.allCheck_)
{}


Foam::
waveAbsorptionVelocity3DFvPatchVectorField::
waveAbsorptionVelocity3DFvPatchVectorField
(
    const waveAbsorptionVelocity3DFvPatchVectorField& ptf,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchField<vector>(ptf, iF),
    nPaddles_(ptf.nPaddles_),
    initialWaterDepths_(ptf.initialWaterDepths_),
    meanAngles_(ptf.meanAngles_),
    zSpanL_(ptf.zSpanL_),
    nEdgeMin_(ptf.nEdgeMin_),
    nEdgeMax_(ptf.nEdgeMax_),
    allCheck_(ptf.allCheck_)
{}
#endif


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::waveAbsorptionVelocity3DFvPatchVectorField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    // PHC //
    Info << "3D_3D Absorption BC on patch " << this->patch().name() << endl;

    // 3D Variables
    const vector cMin = gMin(patch().patch().localPoints());
    const vector cMax = gMax(patch().patch().localPoints());
    const vector cSpan = cMax - cMin;

    scalar dMin = 0.0;
    scalar dSpan = 0.0;
    const scalarField patchD = patchDirection( cSpan, &dMin, &dSpan );

    // Variables & constants
    const volScalarField& alpha = 
        db().lookupObject<volScalarField>(alphaName());
    const volVectorField& U = db().lookupObject<volVectorField>("U");

    const fvMesh& mesh = alpha.mesh();
	const word& patchName = this->patch().name();
	const label patchID = mesh.boundaryMesh().findPatchID(patchName);
    const label nF = patch().faceCells().size();

    const scalarField alphaCell = 
        alpha.boundaryField()[patchID].patchInternalField();
    const vectorField UCell = U.boundaryField()[patchID].patchInternalField();

    scalarField patchUc = Foam::scalarField(nF, 0.0); // Correction velocity
    scalarField patchVc = Foam::scalarField(nF, 0.0); // Correction velocity

    const scalarField UxCell = UCell.component(0);
    const scalarField UyCell = UCell.component(1);;
    scalarField UzCell = UCell.component(2);

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
        meanAngles_ = meanPatchDirs( faceGroup );
        // Info << "Paddle angles " << meanAngles_ << endl;

        allCheck_ = true;
    }

    // Calculate water measured levels
    scalarList measuredLevels = calcWL( alphaCell, faceGroup, zSpanL_ );

    // Correction velocity: Ucalc = -sqrt(g/h)*corrL
    scalarList corrLevels = measuredLevels - initialWaterDepths_;
    inlinePrint( "Correction Levels", corrLevels );

    scalarList Ucalc = -sqrt(g/max(initialWaterDepths_,0.1))*corrLevels;
    scalarList signoU(nPaddles_, 0.0); // Sign of Ucalc

    // Calculate mean velocities on each paddle
    scalarList groupUx (nPaddles_,0.0);
    scalarList groupUy (nPaddles_,0.0);
    calcUV( alphaCell, faceGroup, UxCell, &groupUx, UyCell, &groupUy );

    // Calculate the mean tangential velocity of water for each paddle
    const scalarList meanTgAngle = meanAngles_ + PI()/2.0;
    scalarList Utg (nPaddles_,0.0); // Tg velocity
    scalarList Uc (nPaddles_,0.0); // Correction velocity

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

    // Info << "Ucalc " << Ucalc << endl;
    // Info << "Utg " << Utg << endl;
    // Info << "Uc " << Uc << endl;
    // Info << "Signo " << signoU << endl;

    forAll(patchUc, cellIndex)    
    {
        patchUc[cellIndex] = signoU[faceGroup[cellIndex]-1]*
            Uc[faceGroup[cellIndex]-1]*alphaCell[cellIndex]*
            cos(meanAngles_[faceGroup[cellIndex]-1]);
        patchVc[cellIndex] = signoU[faceGroup[cellIndex]-1]*
            Uc[faceGroup[cellIndex]-1]*alphaCell[cellIndex]*
            sin(meanAngles_[faceGroup[cellIndex]-1]);
        UzCell[cellIndex] = (1.0-alphaCell[cellIndex])*UzCell[cellIndex];
    }

    const vectorField n1 = Foam::vectorField(nF, vector(1.0, 0.0, 0.0));
    const vectorField n2 = Foam::vectorField(nF, vector(0.0, 1.0, 0.0));
    const vectorField n3 = Foam::vectorField(nF, vector(0.0, 0.0, 1.0));

    operator==(n1*patchUc + n2*patchVc + n3*UzCell);

    fixedValueFvPatchField<vector>::updateCoeffs();

}


void Foam::waveAbsorptionVelocity3DFvPatchVectorField::
write(Ostream& os) const
{
    fvPatchField<vector>::write(os);
    #if OFFLAVOUR == 3 && OFVERSION >= 700
        #include "newWriting.H"
    #else
        #include "classicWriting.H"
    #endif
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
   makePatchTypeField
   (
       fvPatchVectorField,
       waveAbsorptionVelocity3DFvPatchVectorField
   );
}


// ************************************************************************* //
