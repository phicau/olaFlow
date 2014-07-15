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

    if (zSup[cellIndex] <= heights[cellGroup[cellIndex]-1]) 
    {// Full cell
        patchU[cellIndex] = UInterp;
        patchV[cellIndex] = 0;
        patchW[cellIndex] = 0;        
    }
    else if (zInf[cellIndex] < heights[cellGroup[cellIndex]-1]) 
    {// Partial cell
        auxiliar = heights[cellGroup[cellIndex]-1] - zInf[cellIndex];
        auxiliarTotal = zSup[cellIndex]-zInf[cellIndex];
        auxiliarTotal = auxiliar/auxiliarTotal; // VOF of the cell

        // Height corresponding to the partial cell
        auxiliar = zInf[cellIndex] + auxiliar/scalar(2);

        patchU[cellIndex] = auxiliarTotal*UInterp;
        patchV[cellIndex] = 0;
        patchW[cellIndex] = 0;  
    }
    else // Empty cell
    {
        patchU[cellIndex] = 0;
        patchV[cellIndex] = 0;
        patchW[cellIndex] = 0;
    }
