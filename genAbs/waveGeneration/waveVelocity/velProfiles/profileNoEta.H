    if ( indexF == 0)
    {
        UInterp = paddleVelocityU_[0];

        if (waveTheory_ == "tvetaw")
        {
            WInterp = paddleVelocityW_[0];
        }
        else
        {
            WInterp = 0;
        }
    }
    else
    {
        UInterp = 
            interpolation
            (
                timeSeries_[indexF-1], 
                timeSeries_[indexF], 
                paddleVelocityU_[indexF-1], 
                paddleVelocityU_[indexF], 
                currTime
            );

        if (waveTheory_ == "tvetaw")
        {
            WInterp = 
                interpolation
                (
                    timeSeries_[indexF-1], 
                    timeSeries_[indexF], 
                    paddleVelocityW_[indexF-1], 
                    paddleVelocityW_[indexF], 
                    currTime
                );
        }
        else
        {
            WInterp = 0;
        }
    }

    patchU[cellIndex] = UInterp*alphaCell[cellIndex];
    patchV[cellIndex] = WInterp*alphaCell[cellIndex];
    patchW[cellIndex] = 0; 
