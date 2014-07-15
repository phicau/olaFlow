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

#include "IH_Waves_InletVelocityFvPatchVectorField.H"
#include "volFields.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "surfaceFields.H"

#include "waveFun.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::
IH_Waves_InletVelocityFvPatchVectorField::
IH_Waves_InletVelocityFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchField<vector>(p, iF),
    wavePeriod_(-1),
    wavePeriods_( List<scalar> (1, -1.0) ),
    waveHeight_(-1),
    waveHeights_( List<scalar> (1, -1.0) ),
    waveLength_(-1),
    waveLengths_( List<scalar> (1, -1.0) ),
    waterDepth_(-1),
    wavePhase_(3.0*PI()/2.0),
    wavePhases_( List<scalar> (1, -1.0) ),
    timeLag_(0),
    timeLags_( List<scalar> (1, 0.0) ),
    lambdaStokesV_(-1),
    mCnoidal_(-1),
    uMean_(-1),
    Bjs_( List<scalar> (1, -1.0) ),
    Ejs_( List<scalar> (1, -1.0) ),
    uCurrent_( vector(0., 0., 0.) ),
    genAbs_(false),
    secondOrder_(false),
    nPaddles_(1),
    tSmooth_(-1),
    tuningFactor_(1),
    nComp_(1),
    waveDictName_("IHWavesDict"),
    waveType_("aaa"),
    waveTheory_("aaa"),
    waveTheoryOrig_("aaa"),
    allCheck_(false),
    absDir_(400.0),
    waveDir_(0),
    waveDirs_( List<scalar> (1, -1.0) ),
    timeSeries_( List<scalar> (1, -1.0) ),
    paddlePosition_( List<scalar> (1, -1.0) ),
    paddleVelocity_( List<scalar> (1, -1.0) ),
    paddleEta_( List<scalar> (1, -1.0) )
{}


Foam::
IH_Waves_InletVelocityFvPatchVectorField::
IH_Waves_InletVelocityFvPatchVectorField
(
    const IH_Waves_InletVelocityFvPatchVectorField& ptf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchField<vector>(ptf, p, iF, mapper),
    wavePeriod_(ptf.wavePeriod_),
    wavePeriods_(ptf.wavePeriods_),
    waveHeight_(ptf.waveHeight_),
    waveHeights_(ptf.waveHeights_),
    waveLength_(ptf.waveLength_),
    waveLengths_(ptf.waveLengths_),
    waterDepth_(ptf.waterDepth_),
    wavePhase_(ptf.wavePhase_),
    wavePhases_(ptf.wavePhases_),
    timeLag_(ptf.timeLag_),
    timeLags_(ptf.timeLags_),
    lambdaStokesV_(ptf.lambdaStokesV_),
    mCnoidal_(ptf.mCnoidal_),
    uMean_(ptf.uMean_),
    Bjs_(ptf.Bjs_),
    Ejs_(ptf.Ejs_),
    uCurrent_(ptf.uCurrent_),
    genAbs_(ptf.genAbs_),
    secondOrder_(ptf.secondOrder_),
    nPaddles_(ptf.nPaddles_),
    tSmooth_(ptf.tSmooth_),
    tuningFactor_(ptf.tuningFactor_),
    nComp_(ptf.nComp_),
    waveDictName_(ptf.waveDictName_),
    waveType_(ptf.waveType_),
    waveTheory_(ptf.waveTheory_),
    waveTheoryOrig_(ptf.waveTheoryOrig_),
    allCheck_(ptf.allCheck_),
    absDir_(ptf.absDir_),
    waveDir_(ptf.waveDir_),
    waveDirs_(ptf.waveDirs_),
    timeSeries_(ptf.timeSeries_),
    paddlePosition_(ptf.paddlePosition_),
    paddleVelocity_(ptf.paddleVelocity_),
    paddleEta_(ptf.paddleEta_)
{}


Foam::
IH_Waves_InletVelocityFvPatchVectorField::
IH_Waves_InletVelocityFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchField<vector>(p, iF, dict),
    wavePeriod_(dict.lookupOrDefault<scalar>("wavePeriod", -1)),
    wavePeriods_( dict.lookupOrDefault("wavePeriods", List<scalar> (1, -1.0)) ),
    waveHeight_(dict.lookupOrDefault<scalar>("waveHeight", -1)),
    waveHeights_( dict.lookupOrDefault("waveHeights", List<scalar> (1, -1.0)) ),
    waveLength_(dict.lookupOrDefault<scalar>("waveLength", -1)),
    waveLengths_( dict.lookupOrDefault("waveLengths", List<scalar> (1, -1.0)) ),
    waterDepth_(dict.lookupOrDefault<scalar>("waterDepth", -1 )),
    wavePhase_(dict.lookupOrDefault<scalar>("wavePhase", 3.0*PI()/2.0 )),
    wavePhases_( dict.lookupOrDefault("wavePhases", List<scalar> (1, -1.0)) ),
    timeLag_(dict.lookupOrDefault<scalar>("timeLag", 0)),
    timeLags_( dict.lookupOrDefault("timeLags", List<scalar> (1, 0.0)) ),
    lambdaStokesV_(dict.lookupOrDefault<scalar>("lambdaStokesV", -1 )),
    mCnoidal_(dict.lookupOrDefault<scalar>("mCnoidal", -1 )),
    uMean_(dict.lookupOrDefault<scalar>("uMean", -1)),
    Bjs_( dict.lookupOrDefault("Bjs", List<scalar> (1, -1.0)) ),
    Ejs_( dict.lookupOrDefault("Ejs", List<scalar> (1, -1.0)) ),
    uCurrent_(dict.lookupOrDefault("uCurrent", vector(0., 0., 0.))),
    genAbs_(dict.lookupOrDefault<bool>("genAbs", false )),
    secondOrder_(dict.lookupOrDefault<bool>("secondOrder", false )),
    nPaddles_(dict.lookupOrDefault<label>("nPaddles", 1)),
    tSmooth_(dict.lookupOrDefault<scalar>("tSmooth", -1)),
    tuningFactor_(dict.lookupOrDefault<scalar>("tuningFactor", 1)),
    nComp_(dict.lookupOrDefault<label>("nComp", 1)),
    waveDictName_(dict.lookupOrDefault<word>("waveDict", "IHWavesDict")),
    waveType_(dict.lookupOrDefault<word>("waveType", "aaa")),
    waveTheory_(dict.lookupOrDefault<word>("waveTheory", "aaa")),
    waveTheoryOrig_(dict.lookupOrDefault<word>("waveTheoryOrig", "aaa")),
    allCheck_(dict.lookupOrDefault<bool>("allCheck", false )),
    absDir_(dict.lookupOrDefault<scalar>("absDir", 400.0)),
    waveDir_(dict.lookupOrDefault<scalar>("waveDir", 0)),
    waveDirs_( dict.lookupOrDefault("waveDirs", List<scalar> (1, -1.0)) ),
    timeSeries_( dict.lookupOrDefault("timeSeries", List<scalar> (1, -1.0)) ),
    paddlePosition_( 
        dict.lookupOrDefault("paddlePosition", List<scalar> (1, -1.0)) ),
    paddleVelocity_( 
        dict.lookupOrDefault("paddleVelocity", List<scalar> (1, -1.0)) ),
    paddleEta_( dict.lookupOrDefault("paddleEta", List<scalar> (1, -1.0)) )
{
    word dictName = dict.lookupOrDefault<word>("waveDict", "empty");
    if(dictName!="empty")
    {
        Warning << 
        "Keyword waveDict defined in boundary condition.\n" << 
        "The new expected keyword is: waveDictName\n" << 
        "Using waveDictName = IHWavesDict by default.\n" << endl;
    }
}


Foam::
IH_Waves_InletVelocityFvPatchVectorField::
IH_Waves_InletVelocityFvPatchVectorField
(
    const IH_Waves_InletVelocityFvPatchVectorField& ptf
)
:
    fixedValueFvPatchField<vector>(ptf),
    wavePeriod_(ptf.wavePeriod_),
    wavePeriods_(ptf.wavePeriods_),
    waveHeight_(ptf.waveHeight_),
    waveHeights_(ptf.waveHeights_),
    waveLength_(ptf.waveLength_),
    waveLengths_(ptf.waveLengths_),
    waterDepth_(ptf.waterDepth_),
    wavePhase_(ptf.wavePhase_),
    wavePhases_(ptf.wavePhases_),
    timeLag_(ptf.timeLag_),
    timeLags_(ptf.timeLags_),
    lambdaStokesV_(ptf.lambdaStokesV_),
    mCnoidal_(ptf.mCnoidal_),
    uMean_(ptf.uMean_),
    Bjs_(ptf.Bjs_),
    Ejs_(ptf.Ejs_),
    uCurrent_(ptf.uCurrent_),
    genAbs_(ptf.genAbs_),
    secondOrder_(ptf.secondOrder_),
    nPaddles_(ptf.nPaddles_),
    tSmooth_(ptf.tSmooth_),
    tuningFactor_(ptf.tuningFactor_),
    nComp_(ptf.nComp_),
    waveDictName_(ptf.waveDictName_),
    waveType_(ptf.waveType_),
    waveTheory_(ptf.waveTheory_),
    waveTheoryOrig_(ptf.waveTheoryOrig_),
    allCheck_(ptf.allCheck_),
    absDir_(ptf.absDir_),
    waveDir_(ptf.waveDir_),
    waveDirs_(ptf.waveDirs_),
    timeSeries_(ptf.timeSeries_),
    paddlePosition_(ptf.paddlePosition_),
    paddleVelocity_(ptf.paddleVelocity_),
    paddleEta_(ptf.paddleEta_)
{}


Foam::
IH_Waves_InletVelocityFvPatchVectorField::
IH_Waves_InletVelocityFvPatchVectorField
(
    const IH_Waves_InletVelocityFvPatchVectorField& ptf,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchField<vector>(ptf, iF),
    wavePeriod_(ptf.wavePeriod_),
    wavePeriods_(ptf.wavePeriods_),
    waveHeight_(ptf.waveHeight_),
    waveHeights_(ptf.waveHeights_),
    waveLength_(ptf.waveLength_),
    waveLengths_(ptf.waveLengths_),
    waterDepth_(ptf.waterDepth_),
    wavePhase_(ptf.wavePhase_),
    wavePhases_(ptf.wavePhases_),
    timeLag_(ptf.timeLag_),
    timeLags_(ptf.timeLags_),
    lambdaStokesV_(ptf.lambdaStokesV_),
    mCnoidal_(ptf.mCnoidal_),
    uMean_(ptf.uMean_),
    Bjs_(ptf.Bjs_),
    Ejs_(ptf.Ejs_),
    uCurrent_(ptf.uCurrent_),
    genAbs_(ptf.genAbs_),
    secondOrder_(ptf.secondOrder_),
    nPaddles_(ptf.nPaddles_),
    tSmooth_(ptf.tSmooth_),
    tuningFactor_(ptf.tuningFactor_),
    nComp_(ptf.nComp_),
    waveDictName_(ptf.waveDictName_),
    waveType_(ptf.waveType_),
    waveTheory_(ptf.waveTheory_),
    waveTheoryOrig_(ptf.waveTheoryOrig_),
    allCheck_(ptf.allCheck_),
    absDir_(ptf.absDir_),
    waveDir_(ptf.waveDir_),
    waveDirs_(ptf.waveDirs_),
    timeSeries_(ptf.timeSeries_),
    paddlePosition_(ptf.paddlePosition_),
    paddleVelocity_(ptf.paddleVelocity_),
    paddleEta_(ptf.paddleEta_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::IH_Waves_InletVelocityFvPatchVectorField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    // PHC //
    Info << "Velocity BC on patch " << this->patch().name() << endl;

    // Auxiliar variables
    scalar auxiliar = 0; 
    scalar auxiliarTotal = 0;
    scalar auxiliarSolit = 0;
    scalar auxiliarSO = 0;

    // Variables solitary
    scalar Csolitary = 0;
    scalar etaSolit = 0;
    scalar ts = 0;
    scalar Xa = 0;
    scalar X0 = 0;
    scalarField patchXsolit;

    // Variables stream function
    scalar celerity = 0;
    scalar faseTot;

    // Variables tveta
    scalar etaInterp = 0;
    scalar UInterp = 0;
    label indexF = 0;

    // 3D Variables
    const vector cMin = gMin(patch().patch().localPoints());
    const vector cMax = gMax(patch().patch().localPoints());
    const vector cSpan = cMax - cMin;
    const scalar zSpan = cSpan[2];

    scalar dMin = 0.0;
    scalar dSpan = 0.0;
    const scalarField patchD = patchDirection( cSpan, &dMin, &dSpan );

    // Variables & constants
    const fvMesh& mesh = dimensionedInternalField().mesh();
	const word& patchName = this->patch().name();
	const label patchID = mesh.boundaryMesh().findPatchID(patchName);
    const label nF = patch().faceCells().size();
    labelList cellGroup = Foam::labelList(nF, 1);

    const volScalarField& alpha = 
        db().lookupObject<volScalarField>(alphaName());
    const volVectorField& U = db().lookupObject<volVectorField>("U");

    const scalarField alphaCell = 
        alpha.boundaryField()[patchID].patchInternalField();
    const vectorField UCell = U.boundaryField()[patchID].patchInternalField();

    const vectorField nVecCell = patch().nf();

    scalarField patchU = Foam::scalarField(nF, 0.0);
    scalarField patchUABS = Foam::scalarField(nF, 0.0);
    scalarField patchV = Foam::scalarField(nF, 0.0);
    scalarField patchVABS = Foam::scalarField(nF, 0.0);
    scalarField patchW = Foam::scalarField(nF, 0.0);

    const labelList celdas = patch().faceCells();

    const scalarField patchHeight = patch().Cf().component(2)-cMin[2];

    const scalar g = 9.81;

    // Calculate Z bounds of the faces
    scalarField zSup, zInf;
    faceBoundsZ( &zSup, &zInf );
    // Change in reference level (-= zMin)
    zSup -= cMin[2];
    zInf -= cMin[2];

    // Waves variables
    scalar waveOmega;
    scalarList waveOmegas;
    scalar waveK;
    scalarList waveKs;
    scalar waveAngle;
    scalarList waveAngles;
    scalar waveKx;
    scalarList waveKxs;
    scalar waveKy;
    scalarList waveKys;

    // Define dictionary
    IOdictionary IHWavesDict
    (
        IOobject
        (
            waveDictName_,
            this->db().time().constant(),
            this->db(),
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        )
    );

    // Check for errors - Just the first time
    if (!allCheck_)
    {
        waveType_ = (IHWavesDict.lookupOrDefault<word>("waveType", "aaa")); 

        tSmooth_ = (IHWavesDict.lookupOrDefault<scalar>("tSmooth", -1.0 ));
        tuningFactor_ = 
            (IHWavesDict.lookupOrDefault<scalar>("tuningFactor", 1.0 ));

        uCurrent_ = 
            (IHWavesDict.lookupOrDefault<vector>("uCurrent", vector(0,0,0) ));

        if ( waveType_ == "aaa" )    // No target specified
        {
            FatalError
                << "Wave type not specified. Use:\n"
                << "regular, solitary, irregular, wavemaker."
                << exit(FatalError);
        }

        if ( waveType_ == "regular" )
        {
            #include "checkInputErrorsRegular.H"
        }
        else if ( waveType_ == "solitary" )
        {
            #include "checkInputErrorsSolitary.H"
        }
        else if ( waveType_ == "irregular" )
        {
            #include "checkInputErrorsIrregular.H"
        }
        else if ( waveType_ == "wavemaker" )
        {
            #include "checkInputErrorsWavemaker.H"
        }
        else if ( waveType_ == "current" )
        {
            #include "checkInputErrorsCurrent.H"
        }
        else
        {
            FatalError
                << "Wave type not supported, use:\n"
                << "regular, solitary, irregular, wavemaker."
                << exit(FatalError);
        }

        allCheck_ = true;
    } // End of allCheck

    scalar currTime = this->db().time().value();
    scalar timeMult = tuningFactor_;

    // Setting the rest of the wave variables
    if ( tSmooth_ > 0.0 && currTime < tSmooth_ )
    {
        timeMult = timeMult*currTime/tSmooth_;
    }

    if ( waveType_ == "regular" )
    {
        waveOmega = (2.0*PI())/wavePeriod_;
        waveK = 2.0*PI()/waveLength_;

        celerity = waveLength_/wavePeriod_;

        waveAngle = waveDir_*PI()/180.0;
        waveKx = waveK*cos(waveAngle);
        waveKy = waveK*sin(waveAngle);

        // Add lag to correct the phase
        if 
        ( 
            waveTheory_ == "StokesII" || 
            waveTheory_ == "StokesV" || 
            waveTheory_ == "cnoidal"
        )
        { 
            currTime += timeLag_;
        }
    }
    else if ( waveType_ == "solitary" )
    {
        waveAngle = waveDir_*PI()/180.0;
        patchXsolit = 
            patch().Cf().component(0)*cos(waveAngle) 
            + patch().Cf().component(1)*sin(waveAngle);
        X0 = gMin(patchXsolit);
    }
    else if ( waveType_ == "irregular" )
    {
        waveOmegas = (2.0*PI())/wavePeriods_;
        waveKs = 2.0*PI()/waveLengths_;

        waveAngles = waveDirs_*PI()/180.0;
        waveKxs = waveKs*cos(waveAngles);
        waveKys = waveKs*sin(waveAngles);
    }

    // Grouping part
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
                cellGroup[patchCells] = i+1; // Group of each face
                continue;
            }
        }      
    }

    // Absorption direction part
    scalarList meanAngle (nPaddles_, 0.0);
    bool absDireccional = false;

    // Check if absorption is directional
    if ( absDir_ > 360.0 ) // Automatic
    {
        absDireccional = true;
        meanAngle = meanPatchDirs( cellGroup );
    }
    else // Fixed
    {
        meanAngle = absDir_*PI()/180.0;
    }
    // Info << "Paddle angle " << meanAngle << endl;

    scalarList calculatedLevel (nPaddles_,0.0);

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

    // Calculate water measured levels
    scalarList measuredLevels = calcWL( alphaCell, cellGroup, zSpan );

    // Auxiliar variables
    auxiliarTotal = 0.0;
    auxiliar = 0.0;

    // Define heights as minimum of calculatedLevel and measuredLevels

    scalarList heights (nPaddles_,0.0);
    forAll(heights, iterMin)
    {
        heights[iterMin] = 
            min(calculatedLevel[iterMin],measuredLevels[iterMin]);
    }

    bool noEta = false;
    if ( waveTheoryOrig_ == "tx" || waveTheoryOrig_ == "tv" )
    {
        noEta = true;
    }

    // Velocity cycle
    scalar corrLevel = 0.0;
    forAll(patchHeight, cellIndex)    
    {
        #include "velocityProfile.H"

        vector cellV = 
            vector(patchU[cellIndex], patchV[cellIndex], patchW[cellIndex]);

        // MODIFICADO --- CHECK DIRECTION
        patchU[cellIndex] = patchU[cellIndex]*alphaCell[cellIndex];
        patchV[cellIndex] = patchV[cellIndex]*alphaCell[cellIndex];
        patchW[cellIndex] = patchW[cellIndex]*alphaCell[cellIndex];

        //Info<< "U " << patchU[cellIndex] << endl;
        //Info<< "W " << patchW[cellIndex] << endl;

        if ( genAbs_ ) // Velocity correction when target free surface specified
        {
            corrLevel = // Theoretical - Measured
                calculatedLevel[cellGroup[cellIndex]-1]
                - measuredLevels[cellGroup[cellIndex]-1];

            if (waveType_ == "irregular")
            {
                corrLevel = pos(corrLevel)*corrLevel;
            }

            scalar corrABSfactor = 1.0;
            if(corrLevel > 0.0) // Theoretical > Measured == Inflow
            {
                corrABSfactor = 
                    calculatedLevel[cellGroup[cellIndex]-1]
                    /measuredLevels[cellGroup[cellIndex]-1];

                if (zSup[cellIndex] <= measuredLevels[cellGroup[cellIndex]-1])
                {
                    patchUABS[cellIndex] = 
                        corrABSfactor*pos(alphaCell[cellIndex]-0.9)
                        *corrLevel*sqrt(g/waterDepth_)
                        *cos(meanAngle[cellGroup[cellIndex]-1]);
                    patchVABS[cellIndex] = 
                        corrABSfactor*pos(alphaCell[cellIndex]-0.9)
                        *corrLevel*sqrt(g/waterDepth_)
                        *sin(meanAngle[cellGroup[cellIndex]-1]);
                }
            }
            else // Theoretical < Measured == Outflow
            {
                corrABSfactor = 
                    measuredLevels[cellGroup[cellIndex]-1]
                    /calculatedLevel[cellGroup[cellIndex]-1];

                patchUABS[cellIndex] = 
                    corrABSfactor*alphaCell[cellIndex]
                    *corrLevel*sqrt(g/waterDepth_)
                    *cos(meanAngle[cellGroup[cellIndex]-1]);
                patchVABS[cellIndex] = 
                    corrABSfactor*alphaCell[cellIndex]
                    *corrLevel*sqrt(g/waterDepth_)
                    *sin(meanAngle[cellGroup[cellIndex]-1]);

                if (zInf[cellIndex] >= calculatedLevel[cellGroup[cellIndex]-1])
                {
                    patchUABS[cellIndex] *= 
                        1.0 + 2.0*alphaCell[cellIndex]*(zInf[cellIndex]
                        - calculatedLevel[cellGroup[cellIndex]-1])
                        /(cMax[2]-calculatedLevel[cellGroup[cellIndex]-1]);
                    patchVABS[cellIndex] *= 
                        1.0 + 2.0*alphaCell[cellIndex]*(zInf[cellIndex]
                        - calculatedLevel[cellGroup[cellIndex]-1])
                        /(cMax[2]-calculatedLevel[cellGroup[cellIndex]-1]);
                }
            } 
        }
    }
//    inlinePrint( "Theoretical Level ", calculatedLevel );
//    inlinePrint( "Measured Level ", measuredLevels );

    const vectorField n1 = Foam::vectorField(nF, vector(1.0, 0.0, 0.0));
    const vectorField n2 = Foam::vectorField(nF, vector(0.0, 1.0, 0.0));
    const vectorField n3 = Foam::vectorField(nF, vector(0.0, 0.0, 1.0));

    // Set Velocity
    operator==
        (n1*patchU*timeMult*pos(alphaCell-0.9) + n1*patchUABS*timeMult
        + n2*patchV*timeMult*pos(alphaCell-0.9) + n2*patchVABS*timeMult
        + n3*patchW*timeMult 
        + uCurrent_*pos(alphaCell-0.9)*neg(patchHeight - min(measuredLevels))
        *timeMult);

    fixedValueFvPatchField<vector>::updateCoeffs();
}


void Foam::IH_Waves_InletVelocityFvPatchVectorField::write(Ostream& os) const
{
    fvPatchField<vector>::write(os);
    os.writeKeyword("waveType") << waveType_ << token::END_STATEMENT << nl;
    os.writeKeyword("waterDepth") << waterDepth_ << token::END_STATEMENT << nl;
    os.writeKeyword("genAbs") << genAbs_ << token::END_STATEMENT << nl;
    os.writeKeyword("nPaddles") << nPaddles_ << token::END_STATEMENT << nl; 
    os.writeKeyword("allCheck") << allCheck_ << token::END_STATEMENT << nl;
    os.writeKeyword("waveDictName") << waveDictName_ << token::END_STATEMENT << nl;
    os.writeKeyword("uCurrent") << uCurrent_ << token::END_STATEMENT << nl;

    if ( tSmooth_ != -1.0 )
    {
        os.writeKeyword("tSmooth") << tSmooth_ << token::END_STATEMENT << nl;
    }

    if ( tuningFactor_ != 1.0 )
    {
        os.writeKeyword("tuningFactor") << 
            tuningFactor_ << token::END_STATEMENT << nl;
    }

    if (genAbs_)
    {
        os.writeKeyword("absDir") << absDir_ << token::END_STATEMENT << nl;
    }

    if ( waveType_ == "irregular" )
    {
        waveHeights_.writeEntry("waveHeights", os);
        wavePeriods_.writeEntry("wavePeriods", os);
        waveLengths_.writeEntry("waveLengths", os);
        wavePhases_.writeEntry("wavePhases", os);
        waveDirs_.writeEntry("waveDirs", os);
        timeLags_.writeEntry("timeLags", os);

        os.writeKeyword("nComp") << nComp_ << token::END_STATEMENT << nl; 

        if ( secondOrder_ )
        {
            os.writeKeyword("secondOrder") << 
                secondOrder_ << token::END_STATEMENT << nl; 
        }

    }
    else if ( waveType_ == "regular" )
    {
        os.writeKeyword("waveTheory") << 
            waveTheory_ << token::END_STATEMENT << nl;
        os.writeKeyword("waveHeight") << 
            waveHeight_ << token::END_STATEMENT << nl;
        os.writeKeyword("waveDir") << waveDir_ << token::END_STATEMENT << nl;
        os.writeKeyword("timeLag") << timeLag_ << token::END_STATEMENT << nl;


        if ( waveTheory_ == "StokesI" || waveTheory_ == "StokesII" )
        {
            os.writeKeyword("waveLength") << 
                waveLength_ << token::END_STATEMENT << nl;
            os.writeKeyword("wavePeriod") << 
                wavePeriod_ << token::END_STATEMENT << nl;
            os.writeKeyword("wavePhase") << 
                wavePhase_ << token::END_STATEMENT << nl;
        }
        else if ( waveTheory_ == "StokesV" )
        {
            os.writeKeyword("waveLength") << 
                waveLength_ << token::END_STATEMENT << nl;
            os.writeKeyword("wavePeriod") << 
                wavePeriod_ << token::END_STATEMENT << nl;
            os.writeKeyword("wavePhase") << 
                wavePhase_ << token::END_STATEMENT << nl;
            os.writeKeyword("lambdaStokesV") << 
                lambdaStokesV_ << token::END_STATEMENT << nl;
        }
        else if ( waveTheory_ == "cnoidal" )
        {
            os.writeKeyword("waveLength") << 
                waveLength_ << token::END_STATEMENT << nl;
            os.writeKeyword("wavePeriod") << 
                wavePeriod_ << token::END_STATEMENT << nl;
            os.writeKeyword("wavePhase") << 
                wavePhase_ << token::END_STATEMENT << nl;
            os.writeKeyword("mCnoidal") << 
                mCnoidal_ << token::END_STATEMENT << nl;
        }
        else if ( waveTheory_ == "streamFunction" )
        {
            os.writeKeyword("uMean") << uMean_ << token::END_STATEMENT << nl;
            Bjs_.writeEntry("Bjs", os);
            Ejs_.writeEntry("Ejs", os);
        }

    }
    else if ( waveType_ == "wavemaker" )
    {
        os.writeKeyword("waveTheory") << 
            waveTheory_ << token::END_STATEMENT << nl;
        timeSeries_.writeEntry("timeSeries", os);
        paddleVelocity_.writeEntry("paddleVelocity", os);
        paddleEta_.writeEntry("paddleEta", os);

        if ( waveTheoryOrig_ != "aaa" )
        {
            os.writeKeyword("waveTheoryOrig") << 
                waveTheoryOrig_ << token::END_STATEMENT << nl;
        }
    }
    else if ( waveType_ == "solitary" )
    {
        os.writeKeyword("waveTheory") << 
            waveTheory_ << token::END_STATEMENT << nl;
        os.writeKeyword("waveHeight") << 
            waveHeight_ << token::END_STATEMENT << nl;
        os.writeKeyword("waveDir") << waveDir_ << token::END_STATEMENT << nl;
    }

    writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
   makePatchTypeField
   (
       fvPatchVectorField,
       IH_Waves_InletVelocityFvPatchVectorField
   );
}


// ************************************************************************* //
