if ( initialWaterDepth_ <= 0.0 )
{
    // Calculate water depth
    scalar groupTotalArea = gSum(cellSurface); // Total area
    scalar groupWaterArea = gSum(cellSurface*alphaCell); // Wet area

    initialWaterDepth_ = groupWaterArea/groupTotalArea*zSpan;

    Info << "Initial water depth Paddle_" 
        << this->patch().name() << " = " << initialWaterDepth_ << endl;
}

wavemakerType_ = 
    (wavemakerMovementDict.lookupOrDefault<word>("wavemakerType", "aaa"));

// Extracting values from dict
timeSeries_ = 
    (wavemakerMovementDict.lookupOrDefault
        ("timeSeries", List<scalar> (1, -1.0))
    );

paddlePosition_ = 
    (wavemakerMovementDict.lookupOrDefault
        ("paddlePosition", List<List<scalar> > (1,List<scalar> (1, -1.0)) )
    );
paddleEta_ = 
    (wavemakerMovementDict.lookupOrDefault
        ("paddleEta", List<List<scalar> > (1,List<scalar> (1, -1.0)) )
    );
paddleTilt_ = 
    (wavemakerMovementDict.lookupOrDefault
        ("paddleTilt", List<List<scalar> > (1,List<scalar> (1, -1.0)) )
    );

tSmooth_ = (wavemakerMovementDict.lookupOrDefault<scalar>("tSmooth", -1.0 ));
tuningFactor_ = 
    (wavemakerMovementDict.lookupOrDefault<scalar>("tuningFactor", 1.0 ));
genAbs_ = (wavemakerMovementDict.lookupOrDefault<bool>("genAbs", false ));
hingeHeight_ = 
    (wavemakerMovementDict.lookupOrDefault<scalar>("hingeHeight", 999.0 ));
hingeLocation_ = 
    (wavemakerMovementDict.lookupOrDefault<scalar>("hingeLocation", 999.0 ));
maxStroke_ = 
    (wavemakerMovementDict.lookupOrDefault<scalar>("maxStroke", 999.0 ));
DPST_ = (wavemakerMovementDict.lookupOrDefault<scalar>("DPST", 25 ));

if ( wavemakerType_ == "Piston" )
{
    nPaddles_ = paddlePosition_.size();
}
else if ( wavemakerType_ == "Flap" )
{
    if ( hingeHeight_ > zMax || hingeHeight_ < zMin )
    {
        hingeHeight_ = zMin;
    }

    hingeLocation_ = gMin(this->patch().localPoints().component(0));

    nPaddles_ = paddleTilt_.size();
    genAbs_ = false;
}
else if ( wavemakerType_ == "Mixed" )
{
    if ( hingeHeight_ > zMax || hingeHeight_ < zMin )
    {
        hingeHeight_ = zMin;
    }

    nPaddles_ = paddlePosition_.size();
    genAbs_ = false;
}
else
{
    FatalError
        << "Wavemaker type not specified. Use:\n"
        << "Piston, Flap, Mixed."
        << exit(FatalError);
}

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
if ( wavemakerType_ == "Piston" )
{
    if ( genAbs_ && paddleEta_[0].size() == 1 )
    {
        WarningIn("wavemakerMovement::updateCoeffs()")
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
}
else if ( wavemakerType_ == "Mixed" )
{
    if ( nPaddles_ != paddleTilt_.size() )
    {
        FatalError
            << "Check number of paddles "
            << "(elements in paddlePosition and paddleTilt):\n"
            << "paddlePosition.size() = " << nPaddles_
            << "\npaddleTilt.size() = " << paddleTilt_.size()
            << exit(FatalError);
    }
}

// Check number of elements in paddlePosition time series
label minAuxPos =  999999, minAuxEta =  999999, minAuxTil =  999999;
label maxAuxPos = -999999, maxAuxEta = -999999, maxAuxTil = -999999;

for(int i=0; i<paddlePosition_.size(); i++)
{
    minAuxPos = min(minAuxPos, paddlePosition_[i].size());
    maxAuxPos = max(maxAuxPos, paddlePosition_[i].size());
}

for(int i=0; i<paddleEta_.size(); i++)
{
    minAuxEta = min(minAuxEta, paddleEta_[i].size());
    maxAuxEta = max(maxAuxEta, paddleEta_[i].size());
}

for(int i=0; i<paddleTilt_.size(); i++)
{
    minAuxTil = min(minAuxTil, paddleTilt_[i].size());
    maxAuxTil = max(maxAuxTil, paddleTilt_[i].size());
}

if ( wavemakerType_ == "Piston" )
{
    if( timePoints != minAuxPos || minAuxPos != maxAuxPos )
    {
        FatalError
            << "Check number of components of each "
            << "of the series in paddlePosition:\n"
            << "Expected: " << timePoints
            << "; Min: " << minAuxPos << "; Max: " << maxAuxPos
            << exit(FatalError);
    }
}
else if ( wavemakerType_ == "Flap" )
{
    if( timePoints != minAuxTil || minAuxTil != maxAuxTil )
    {
        FatalError
            << "Check number of components of each "
            << "of the series in paddleTilt:\n"
            << "Expected: " << timePoints
            << "; Min: " << minAuxTil << "; Max: " << maxAuxTil
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
}
else if ( wavemakerType_ == "Mixed" )
{
    if( timePoints != minAuxPos || minAuxPos != maxAuxPos )
    {
        FatalError
            << "Check number of components of each "
            << "of the series in paddlePosition:\n"
            << "Expected: " << timePoints
            << "; Min: " << minAuxPos << "; Max: " << maxAuxPos
            << exit(FatalError);
    }

    if( timePoints != minAuxTil || minAuxTil != maxAuxTil )
    {
        FatalError
            << "Check number of components of each "
            << "of the series in paddleTilt:\n"
            << "Expected: " << timePoints
            << "; Min: " << minAuxTil << "; Max: " << maxAuxTil
            << exit(FatalError);
    }
}

// First time initializations
DPS_ = List<bool>(nPaddles_, false);
DPSsign_ = scalarList(nPaddles_, 1.0);
DPStIni_ = scalarList(nPaddles_, -1.0);
instDPSCorrection_ = scalarList(nPaddles_, 0.0);
cumDPSCorrection_ = scalarList(nPaddles_, 0.0);
cumAbsCorrection_ = scalarList(nPaddles_, 0.0);

tiltOld_ = scalarList(nPaddles_, 0.0);
