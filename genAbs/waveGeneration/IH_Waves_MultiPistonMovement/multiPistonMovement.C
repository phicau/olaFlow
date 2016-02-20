/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2012 OpenFOAM Foundation
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
| olaFoam Project                                       ll                    |
|                                                       l l                   |
|   Coder: Pablo Higuera Caubilla                 ooo   l l     aa            |
|   Bug reports: phicau@gmail.com                o   o  l l    a  a           |
|                                                o   o  ll   l a  aa  aa      |
|                                                 ooo    llll   aa  aa        |
|                                                                             |
|                                                FFFFF OOOOO AAAAA M   M      |
|                                                F     O   O A   A MM MM      |
|  Formerly IHFOAM Project                       FFFF  O   O AAAAA M M M      |
|  Work initially developed at IH Cantabria      F     O   O A   A M   M      |
|                                                F     OOOOO A   A M   M      |
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
| - Three-dimensional numerical wave generation with moving boundaries        |
|    Higuera, P., Losada, I.J. and Lara, J.L. (2015)                          |
|    Coastal Engineering, Vol. 101, 35-47.                                    |
|    http://dx.doi.org/10.1016/j.coastaleng.2015.04.003                       |
|                                                                             |
\*---------------------------------------------------------------------------*/

#include "multiPistonMovement.H"
#include "addToRunTimeSelectionTable.H"
#include "pointPatchFieldMapper.H"
#include "pointFields.H"

#include "fvPatchFieldMapper.H"
#include "surfaceFields.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

makePointPatchTypeField
(
    pointPatchVectorField,
    IH_Waves_MultiPistonMovement
);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

IH_Waves_MultiPistonMovement::
IH_Waves_MultiPistonMovement
(
    const pointPatch& p,
    const DimensionedField<vector, pointMesh>& iF
)
:
    fixedValuePointPatchField<vector>(p, iF),
    multiPistonDictName_("IHMultiPistonMovementDict"),
    timeSeries_( List<scalar> (1, -1.0) ),
    nPaddles_(1),
    paddlePosition_( List<List<scalar> > (1,List<scalar> (1, -1.0)) ),
    paddleEta_( List<List<scalar> > (1,List<scalar> (1, -1.0)) ),
    initialWaterDepth_(-1),
    meanAngle_(vector (0,0,0)),
    genAbs_(false),
    DPS_(List<bool> (1, false)),
    maxStroke_(999.0),
    DPST_(25),
    DPSsign_(List<scalar> (nPaddles_, 1.0)),
    DPStIni_(List<scalar> (nPaddles_, -1.0)),
    instDPSCorrection_(List<scalar> (nPaddles_, 0.0)),
    cumDPSCorrection_(List<scalar> (nPaddles_, 0.0)),
    cumAbsCorrection_(List<scalar> (nPaddles_, 0.0)),
    tSmooth_(-1),
    tuningFactor_(1),
    allCheck_(false)
{}


IH_Waves_MultiPistonMovement::
IH_Waves_MultiPistonMovement
(
    const IH_Waves_MultiPistonMovement& ptf,
    const pointPatch& p,
    const DimensionedField<vector, pointMesh>& iF,
    const pointPatchFieldMapper& mapper
)
:
    fixedValuePointPatchField<vector>(ptf, p, iF, mapper),
    multiPistonDictName_(ptf.multiPistonDictName_),
    timeSeries_(ptf.timeSeries_),
    nPaddles_(ptf.nPaddles_),
    paddlePosition_(ptf.paddlePosition_),
    paddleEta_(ptf.paddleEta_),
    initialWaterDepth_(ptf.initialWaterDepth_),
    meanAngle_(ptf.meanAngle_),
    genAbs_(ptf.genAbs_),
    DPS_(ptf.DPS_),
    maxStroke_(ptf.maxStroke_),
    DPST_(ptf.DPST_),
    DPSsign_(ptf.DPSsign_),
    DPStIni_(ptf.DPStIni_),
    instDPSCorrection_(ptf.instDPSCorrection_),
    cumDPSCorrection_(ptf.cumDPSCorrection_),
    cumAbsCorrection_(ptf.cumAbsCorrection_),
    tSmooth_(ptf.tSmooth_),
    tuningFactor_(ptf.tuningFactor_),
    allCheck_(ptf.allCheck_)
{}


IH_Waves_MultiPistonMovement::
IH_Waves_MultiPistonMovement
(
    const pointPatch& p,
    const DimensionedField<vector, pointMesh>& iF,
    const dictionary& dict,
    const bool valueRequired
)
:
    fixedValuePointPatchField<vector>(p, iF, dict, valueRequired),
    multiPistonDictName_(dict.lookupOrDefault<word>("multiPistonDictName", "IHMultiPistonMovementDict")),
    timeSeries_( dict.lookupOrDefault("timeSeries", List<scalar> (1, -1.0)) ),
    nPaddles_(dict.lookupOrDefault<label>("nPaddles", 1)),
    paddlePosition_( dict.lookupOrDefault("paddlePosition", List<List<scalar> > (1,List<scalar> (1, -1.0)) )),
    paddleEta_( dict.lookupOrDefault("paddleEta", List<List<scalar> > (1,List<scalar> (1, -1.0)) )),
    initialWaterDepth_(dict.lookupOrDefault<scalar>("initialWaterDepth", -1 )),
    meanAngle_(dict.lookupOrDefault("meanAngle", vector (0,0,0) )),
    genAbs_(dict.lookupOrDefault<bool>("genAbs", false )),
    DPS_(dict.lookupOrDefault("DPS", List<bool> (1, false) )),
    maxStroke_(dict.lookupOrDefault<scalar>("maxStroke", 999.0)),
    DPST_(dict.lookupOrDefault<scalar>("DPST", 25)),
    DPSsign_(dict.lookupOrDefault("DPSsign", List<scalar> (nPaddles_, 0.0))),
    DPStIni_(dict.lookupOrDefault("DPStIni", List<scalar> (nPaddles_, -1.0))),
    instDPSCorrection_(dict.lookupOrDefault("instDPSCorrection", List<scalar> (nPaddles_, 0.0))),
    cumDPSCorrection_(dict.lookupOrDefault("cumDPSCorrection", List<scalar> (nPaddles_, 0.0))),
    cumAbsCorrection_(dict.lookupOrDefault("cumAbsCorrection", List<scalar> (nPaddles_, 0.0))),
    tSmooth_(dict.lookupOrDefault<scalar>("tSmooth", -1)),
    tuningFactor_(dict.lookupOrDefault<scalar>("tuningFactor", 1)),
    allCheck_(dict.lookupOrDefault<bool>("allCheck", false ))
{}


IH_Waves_MultiPistonMovement::
IH_Waves_MultiPistonMovement
(
    const IH_Waves_MultiPistonMovement& ptf
)
:
    fixedValuePointPatchField<vector>(ptf),
    multiPistonDictName_(ptf.multiPistonDictName_),
    timeSeries_(ptf.timeSeries_),
    nPaddles_(ptf.nPaddles_),
    paddlePosition_(ptf.paddlePosition_),
    paddleEta_(ptf.paddleEta_),
    initialWaterDepth_(ptf.initialWaterDepth_),
    meanAngle_(ptf.meanAngle_),
    genAbs_(ptf.genAbs_),
    DPS_(ptf.DPS_),
    maxStroke_(ptf.maxStroke_),
    DPST_(ptf.DPST_),
    DPSsign_(ptf.DPSsign_),
    DPStIni_(ptf.DPStIni_),
    instDPSCorrection_(ptf.instDPSCorrection_),
    cumDPSCorrection_(ptf.cumDPSCorrection_),
    cumAbsCorrection_(ptf.cumAbsCorrection_),
    tSmooth_(ptf.tSmooth_),
    tuningFactor_(ptf.tuningFactor_),
    allCheck_(ptf.allCheck_)
{}


IH_Waves_MultiPistonMovement::
IH_Waves_MultiPistonMovement
(
    const IH_Waves_MultiPistonMovement& ptf,
    const DimensionedField<vector, pointMesh>& iF
)
:
    fixedValuePointPatchField<vector>(ptf, iF),
    multiPistonDictName_(ptf.multiPistonDictName_),
    timeSeries_(ptf.timeSeries_),
    nPaddles_(ptf.nPaddles_),
    paddlePosition_(ptf.paddlePosition_),
    paddleEta_(ptf.paddleEta_),
    initialWaterDepth_(ptf.initialWaterDepth_),
    meanAngle_(ptf.meanAngle_),
    genAbs_(ptf.genAbs_),
    DPS_(ptf.DPS_),
    maxStroke_(ptf.maxStroke_),
    DPST_(ptf.DPST_),
    DPSsign_(ptf.DPSsign_),
    DPStIni_(ptf.DPStIni_),
    instDPSCorrection_(ptf.instDPSCorrection_),
    cumDPSCorrection_(ptf.cumDPSCorrection_),
    cumAbsCorrection_(ptf.cumAbsCorrection_),
    tSmooth_(ptf.tSmooth_),
    tuningFactor_(ptf.tuningFactor_),
    allCheck_(ptf.allCheck_)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

IH_Waves_MultiPistonMovement::
~IH_Waves_MultiPistonMovement()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void IH_Waves_MultiPistonMovement::updateCoeffs()
{
    if (this->updated())
    {
        return;
    }

    // PHC //
    Info << "Point displacement BC on patch " << this->patch().name() << endl;

    // Variables
    const scalar g = 9.81;

    const volScalarField& alpha = db().lookupObject<volScalarField>(alphaName());
    const fvMesh& mesh = alpha.mesh();

    label patchId = mesh.boundaryMesh().findPatchID(this->patch().name());

    const scalarField cellSurface = mesh.magSf().boundaryField()[patchId]; // Surface of the patch face
    const vectorField cellVector = mesh.Sf().boundaryField()[patchId]/cellSurface; // Unit vector of the patch face
    const scalarField alphaCell = alpha.boundaryField()[patchId]; // Alpha of the patch face


    scalarField patchD = mesh.Cf().boundaryField()[patchId].component(1);

    const scalar yMin = gMin(this->patch().localPoints().component(1)); // Min Y of the patch
    const scalar yMax = gMax(this->patch().localPoints().component(1)); // Max Y of the patch
    const scalar ySpan = yMax-yMin;

    const scalar zMin = gMin(this->patch().localPoints().component(2)); // Min Z of the patch
    const scalar zMax = gMax(this->patch().localPoints().component(2)); // Max Z of the patch
    const scalar zSpan = zMax-zMin;

    // Define dictionary
    IOdictionary IHMultiPistonMovementDict
    (
        IOobject
        (
            multiPistonDictName_,
            this->db().time().constant(),
            this->db(),
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );

    if (!allCheck_) // Just the first time
    {
        if ( initialWaterDepth_ <= 0.0 )
        {
            // Calculate water depth
            scalar groupTotalArea = gSum(cellSurface); // Total area
            scalar groupWaterArea = gSum(cellSurface*alphaCell); // Wet area

            initialWaterDepth_ = groupWaterArea/groupTotalArea*zSpan;

            Info << "Initial water depth Paddle_" 
                << this->patch().name() << " = " << initialWaterDepth_ << endl;
        }

        // Extracting values from dict
        timeSeries_ = (IHMultiPistonMovementDict.lookupOrDefault("timeSeries", List<scalar> (1, -1.0) ));

        paddlePosition_ = (IHMultiPistonMovementDict.lookupOrDefault("paddlePosition", List<List<scalar> > (1,List<scalar> (1, -1.0)) ));
        paddleEta_ = (IHMultiPistonMovementDict.lookupOrDefault("paddleEta", List<List<scalar> > (1,List<scalar> (1, -1.0)) ));

        genAbs_ = (IHMultiPistonMovementDict.lookupOrDefault<bool>("genAbs", false ));

        tSmooth_ = (IHMultiPistonMovementDict.lookupOrDefault<scalar>("tSmooth", -1.0 ));
        tuningFactor_ = (IHMultiPistonMovementDict.lookupOrDefault<scalar>("tuningFactor", 1.0 ));

        maxStroke_ = (IHMultiPistonMovementDict.lookupOrDefault<scalar>("maxStroke", 999.0 ));
        DPST_ = (IHMultiPistonMovementDict.lookupOrDefault<scalar>("DPST", 25 ));

        nPaddles_ = paddlePosition_.size();

        // Calculate mean horizontal angle for the boundary
        vector meanAngle = gSum(cellVector);
        meanAngle *= -1.0/returnReduce( cellSurface.size(), sumOp<label>() );
        meanAngle.component(2) = 0.0;
        meanAngle_ = meanAngle;

        // Checks
        // Check timeSeries
        label timePoints = timeSeries_.size();

        if( timePoints == 1 )
        {
            FatalError
                << "Check number of components of timeSeries (>1):\n"
                << "timeSeries_.size() = " << timePoints
                << exit(FatalError);
        }

        // Check paddleEta and nPaddles
        if ( genAbs_ && paddleEta_[0].size() == 1 )
        {
            WarningIn("IH_Waves_MultiPistonMovement::updateCoeffs()")
                << "No paddleEta provided. Assuming: paddleEta = paddlePosition"
                << endl;

            paddleEta_ = paddlePosition_;
        }

        if ( genAbs_ && nPaddles_ != paddleEta_.size() )
        {
            FatalError
                << "Check number of paddles "
                << "(elements in paddlePosition and paddleEta):\n"
                << "paddlePosition.size() = " << nPaddles_
                << "\npaddleEta.size() = " << paddleEta_.size()
                << exit(FatalError);
        }

        // Check number of elements in paddlePosition time series
        label minAuxPos =  999999, minAuxEta =  999999;
        label maxAuxPos = -999999, maxAuxEta = -999999;
        
        for(int i=0; i<nPaddles_; i++)
        {
            minAuxPos = min(minAuxPos, paddlePosition_[i].size());
            maxAuxPos = max(maxAuxPos, paddlePosition_[i].size());
            minAuxEta = min(minAuxEta, paddleEta_[i].size());
            maxAuxEta = max(maxAuxEta, paddleEta_[i].size());
        }

        if( timePoints != minAuxPos || minAuxPos != maxAuxPos )
        {
            FatalError
                << "Check number of components of each "
                << "of the series in paddlePosition:\n"
                << "Expected: " << timePoints
                << "; Min: " << minAuxPos << "; Max: " << maxAuxPos
                << exit(FatalError);
        }

        if( genAbs_ && (timePoints != minAuxEta || minAuxEta != maxAuxEta) )
        {
            FatalError
                << "Check number of components of each series of paddleEta:\n"
                << "Expected: " << timePoints
                << "; Min: " << minAuxEta << "; Max: " << maxAuxEta
                << exit(FatalError);
        }

        allCheck_ = true;

        // First time initializations
        DPS_ = List<bool>(nPaddles_, false);
        DPSsign_ = scalarList(nPaddles_, 1.0);
        DPStIni_ = scalarList(nPaddles_, -1.0);
        instDPSCorrection_ = scalarList(nPaddles_, 0.0);
        cumDPSCorrection_ = scalarList(nPaddles_, 0.0);
        cumAbsCorrection_ = scalarList(nPaddles_, 0.0);
    }

    // Time interpolation
    scalar currTime = this->db().time().value();
    scalar timeMult = tuningFactor_;

    if ( tSmooth_ > 0.0 && currTime < tSmooth_ )
    {
        timeMult = timeMult*currTime/tSmooth_;
    }

    if( currTime < min(timeSeries_) )
    {
        FatalError
            << "The first time in timeSeries should be <= startTime."
            << exit(FatalError);
    }
    else if( currTime > max(timeSeries_) )
    {
        FatalError
            << "The time series is not long enough."
            << exit(FatalError);
    }

    // Interpolate eta and displacement
    label indexF = 0;
    scalarList dispInterp = scalarList(nPaddles_, 0.0);
    scalarList etaInterp = scalarList(nPaddles_, 0.0);

    forAll( timeSeries_, contador )
    {
        if ( timeSeries_[contador] >= currTime )
        {
            indexF = contador;
            break;
        }
    }

    for(int i=0; i<nPaddles_; i++)
    {
        if ( indexF == 0)
        {
            dispInterp[i] = timeMult * paddlePosition_[i][0];
            if ( genAbs_ ){ etaInterp[i] = timeMult * paddleEta_[i][0]; }
        }
        else
        {
            dispInterp[i] = timeMult * linInterp(timeSeries_[indexF-1], timeSeries_[indexF], paddlePosition_[i][indexF-1], paddlePosition_[i][indexF], currTime);
            if ( genAbs_ )
            {
                etaInterp[i] = timeMult * linInterp(timeSeries_[indexF-1], timeSeries_[indexF], paddleEta_[i][indexF-1], paddleEta_[i][indexF], currTime);
            }
        }
    }

    // Active absorption correction
    if ( genAbs_ )
    {
        // Measure water depth
        scalarList measuredWaterLevel = calcWL( alphaCell, patchD, cellSurface, yMin, ySpan, zSpan );

        // Calculate correction
        scalar deltaT = db().time().deltaTValue();

        for(int i=0; i<nPaddles_; i++)
        {
            scalar expectedWaterLevel = initialWaterDepth_ + etaInterp[i];
            cumAbsCorrection_[i] -= timeMult * deltaT * (measuredWaterLevel[i] - expectedWaterLevel) * sqrt(g/initialWaterDepth_); // !Angle correction
        }
    }

    // Drift Prevention System (DPS)
    if( maxStroke_ != 999.0 ) 
    {
        for(int i=0; i<nPaddles_; i++)
        {
            if(!DPS_[i]) 
            {
                // Test if need to connect
                scalar dispAux = dispInterp[i] + cumAbsCorrection_[i] + cumDPSCorrection_[i];
                if( dispAux > 0.8*maxStroke_ )
                {
                    DPS_[i] = true;
                    DPSsign_[i] = -1.0;
                    DPStIni_[i] = currTime;
                }
                else if( dispAux < -0.8*maxStroke_ )
                {
                    DPS_[i] = true;
                    DPSsign_[i] = 1.0;
                    DPStIni_[i] = currTime;
                }
            }

            if(DPS_[i]) // If active
            {
                if( currTime-DPStIni_[i] < DPST_ )
                {
                    instDPSCorrection_[i] = DPSsign_[i] * DPSramp(maxStroke_, currTime-DPStIni_[i], DPST_);
                }
                else
                {
                    instDPSCorrection_[i] = 0.0;
                    cumDPSCorrection_[i] += DPSsign_[i] * maxStroke_;
                    DPS_[i] = false;
                }
            }
        }
    }

    // Displacement of the paddles
    scalarList displacements = scalarList(nPaddles_, 0.0);
    scalarList paddleCenters = scalarList(nPaddles_, 0.0);
    scalar paddleWidth = ySpan/nPaddles_;

    for(int i=0; i<nPaddles_; i++)
    {
        paddleCenters[i] = yMin + (0.5 + i) * paddleWidth;
        displacements[i] = dispInterp[i] + cumAbsCorrection_[i] + instDPSCorrection_[i] + cumDPSCorrection_[i];
    }    

    Info << "Displacement Paddle_" << this->patch().name() << " => " << displacements << endl;

    // Interpolation to the points
    vectorField auxPoints = this->patch().localPoints();

    forAll( auxPoints, point )
    {
        // Interpolation procedure
        scalar pDisp = paddleCosInterp(paddleCenters, displacements, auxPoints[point].component(1));

        auxPoints[point] = pDisp * meanAngle_;
    }

    (*this) == (auxPoints);

    this->fixedValuePointPatchField<vector>::updateCoeffs();

}

void IH_Waves_MultiPistonMovement::write(Ostream& os) const
{
    fixedValuePointPatchField<vector>::write(os);
    os.writeKeyword("allCheck") << allCheck_ << token::END_STATEMENT << nl;
    os.writeKeyword("multiPistonDictName") << multiPistonDictName_ << token::END_STATEMENT << nl;

    os.writeKeyword("nPaddles") << nPaddles_ << token::END_STATEMENT << nl;
    os.writeKeyword("initialWaterDepth") << initialWaterDepth_ << token::END_STATEMENT << nl;
    os.writeKeyword("meanAngle") << meanAngle_ << token::END_STATEMENT << nl;

    timeSeries_.writeEntry("timeSeries", os);
    // paddlePosition_.writeEntry("paddlePosition", os);
    os.writeKeyword("paddlePosition") << paddlePosition_ << token::END_STATEMENT << nl; // May not work for binary encoding

    if ( tSmooth_ != -1.0 )
    {
        os.writeKeyword("tSmooth") << tSmooth_ << token::END_STATEMENT << nl;
    }

    if ( tuningFactor_ != 1.0 )
    {
        os.writeKeyword("tuningFactor") << tuningFactor_ << token::END_STATEMENT << nl;
    }

    if ( genAbs_ )
    {
        os.writeKeyword("genAbs") << genAbs_ << token::END_STATEMENT << nl;
        cumAbsCorrection_.writeEntry("cumAbsCorrection", os);
        // paddleEta_.writeEntry("paddleEta", os);
        os.writeKeyword("paddleEta") << paddleEta_ << token::END_STATEMENT << nl; // May not work for binary encoding
    }

    if ( maxStroke_ != 999.0 ) // DPS active
    {
        DPS_.writeEntry("DPS", os);
        os.writeKeyword("maxStroke") << maxStroke_ << token::END_STATEMENT << nl;
        os.writeKeyword("DPST") << DPST_ << token::END_STATEMENT << nl;
        cumDPSCorrection_.writeEntry("cumDPSCorrection", os);
        DPSsign_.writeEntry("DPSsign", os);
        DPStIni_.writeEntry("DPStIni", os);
        instDPSCorrection_.writeEntry("instDPSCorrection", os);
    }

    //writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //

