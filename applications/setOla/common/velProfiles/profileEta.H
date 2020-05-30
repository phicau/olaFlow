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

    if (zSup[cellIndex] <= calculatedLevel[cellIndex]) 
    {// Full cell
        patchU[cellIndex] = UInterp;
        patchV[cellIndex] = WInterp*zExtra/waterDepth_;
        patchW[cellIndex] = 0;        
    }
    else if (zInf[cellIndex] < calculatedLevel[cellIndex]) 
    {// Partial cell
        auxiliar = calculatedLevel[cellIndex] - zInf[cellIndex];
        auxiliarTotal = zSup[cellIndex]-zInf[cellIndex];
        auxiliarTotal = auxiliar/auxiliarTotal; // VOF of the cell

        // Height corresponding to the partial cell
        auxiliar = zInf[cellIndex] + auxiliar/scalar(2);

        patchU[cellIndex] = auxiliarTotal*UInterp;
        patchV[cellIndex] = auxiliarTotal*WInterp*zExtra/waterDepth_;
        patchW[cellIndex] = 0;  
    }
    else // Empty cell
    {
        patchU[cellIndex] = uWind_.component(0);
        patchV[cellIndex] = uWind_.component(1);
        patchW[cellIndex] = uWind_.component(2);
    }
