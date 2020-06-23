// Calculate Z bounds of the faces
scalarField zSup, zInf;
cellBoundsZ( mesh, &zSup, &zInf );

// Set reference distances
scalar X0 = gMin(mesh.points().component(0)); // Variable for solitary
scalar Z0 = gMin(mesh.points().component(2));

// Auxiliar variables
scalar auxiliar = 0; 
scalar auxiliarTotal = 0;
scalar auxiliarSO = 0;

// Variables solitary
scalar Mmc = 1;
scalar Nmc = 1;

// Variables stream function
scalar celerity = 0;
scalar faseTot;

// Variables tveta
scalar etaInterp = 0, UInterp = 0, WInterp = 0;
label indexF = 0;

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
IOdictionary waveDict
(
    IOobject
    (
        waveDictName_,
        runTime.constant(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    )
);

// Read dictionary
waveType_ = (waveDict.lookupOrDefault<word>("waveType", "aaa")); 

tSmooth_ = (waveDict.lookupOrDefault<scalar>("tSmooth", -1.0 ));
tuningFactor_ = 
    (waveDict.lookupOrDefault<scalar>("tuningFactor", 1.0 ));

uCurrent_ = 
    (waveDict.lookupOrDefault<vector>("uCurrent", vector(0,0,0) ));
uWind_ = 
    (waveDict.lookupOrDefault<vector>("uWind", vector(0,0,0) ));

if ( waveType_ == "aaa" )    // No target specified
{
    FatalError
        << "Wave type (waveType) not specified. Use:\n"
        << "regular, solitary, irregular, wavemaker."
        << exit(FatalError);
}

waterLevel_ = (waveDict.lookupOrDefault<scalar>("setWaterLevel", -1666));

if ( waterLevel_ <= -1665.0 )
{
    FatalError
        << "Water level (setWaterLevel) not specified in "
        << waveDictName_
        << exit(FatalError);
}

waterDepth_ = (waveDict.lookupOrDefault<scalar>("setWaterDepth", -1));

if ( waterDepth_ <= 0.0 )
{
    waterDepth_ = waterLevel_ - Z0;

    if ( waterDepth_ <= 0.0 )
    {
        FatalError
            << "Water depth is negative (" << waterDepth_
            << ") Please fix a correct water depth/level (setWaterDepth or setWaterLevel) values."
            << exit(FatalError);
    }
}

// Change in reference level (+= h - l)
Z0 = waterDepth_ - waterLevel_;
zSup += Z0;
zInf += Z0;

// Cell centres
const scalarField xComp = mesh.C().component(0);
const scalarField yComp = mesh.C().component(1);
const scalarField zComp = mesh.C().component(2) + dimensionedScalar("Z0", dimLength, Z0);

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

extrapolation_ = waveDict.lookupOrDefault<bool>("extrapolation", true);
if (!extrapolation_)
{
    Info
        << "Velocity profile extrapolation over SWL not connected."
        << endl;
}

bool noEta = false;
if ( waveTheoryOrig_ == "tx" || waveTheoryOrig_ == "tv" )
{
    noEta = true;
}

scalar currTime = runTime.value();
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
    waveKx = waveK*Foam::cos(waveAngle);
    waveKy = waveK*Foam::sin(waveAngle);

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
    if(waveTheory_ == "McCowan")
    {
        Mmc = McCowanFun::Mcalc(waveHeight_, waterDepth_);
        Nmc = McCowanFun::Ncalc(waveHeight_, waterDepth_, Mmc);
    }
    if(nSolitaryWaves_ > 1)
    {
        double tSeparation = 
            waveDict.lookupOrDefault<scalar>("tSeparation", wavePeriod_ );
        if (currTime > wavePeriod_/2. + tSeparation/2.)
        {
            if (currTime < wavePeriod_/2. + (nSolitaryWaves_ - 1)*tSeparation)
            {
                while (currTime >= wavePeriod_/2. + tSeparation/2.)
                {
                    currTime -= tSeparation;
                }
            }
            else
            {
                currTime -= (nSolitaryWaves_ - 1)*tSeparation;
            }
        }
    }
}
else if ( waveType_ == "irregular" )
{
    waveOmegas = (2.0*PI())/wavePeriods_;
    waveKs = 2.0*PI()/waveLengths_;

    waveAngles = waveDirs_*PI()/180.0;
    waveKxs = waveKs*Foam::cos(waveAngles);
    waveKys = waveKs*Foam::sin(waveAngles);
}

// Auxiliar variables for Stokes III
if ( waveTheory_ == "StokesIII" )
{
    bool iniOK =
        StokesIIIFun::initialise
        (
            waveHeight_,
            waterDepth_,
            wavePeriod_,
            &aE_SIII_,
            &klE_SIII_
        );

    if( !iniOK )
    {
        FatalError
            << "StokesIII initialisation did not converge,"
            << "check your wave conditions."
            << exit(FatalError);
    }

    waveLength_ = 
        StokesIIIFun::waveLength (wavePeriod_, aE_SIII_, klE_SIII_);
}
