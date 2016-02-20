    if ( indexF == 0)
    {
        UInterp = paddleVelocity_[0];
    }
    else
    {
        UInterp = 
            interpolation
            (
                timeSeries_[indexF-1], 
                timeSeries_[indexF], 
                paddleVelocity_[indexF-1], 
                paddleVelocity_[indexF], 
                currTime
            );
    }

    patchU[cellIndex] = UInterp*alphaCell[cellIndex];
    patchV[cellIndex] = 0;
    patchW[cellIndex] = 0; 
