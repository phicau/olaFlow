    // Interpolate eta
    if( currTime > max(timeSeries_) )
    {
        FatalError
            << "Time series is not long enough."
            << exit(FatalError);
    }
    else if( currTime < min(timeSeries_) )
    {
        FatalError
            << "First time should be <= startTime."
            << exit(FatalError);
    }

    forAll( timeSeries_, contador1 )
    {
        if ( timeSeries_[contador1] >= currTime )
        {
            indexF = contador1;
            break;
        }
    }

    if ( indexF == 0)
    {
        etaInterp = paddleEta_[0];
    }
    else
    {
        etaInterp = 
            interpolation
            (
                timeSeries_[indexF-1], 
                timeSeries_[indexF], 
                paddleEta_[indexF-1], 
                paddleEta_[indexF], 
                currTime
            );
    }

    forAll(calculatedLevel, it5)
    {
        calculatedLevel[it5] = waterDepth_ + timeMult * etaInterp;
    }

