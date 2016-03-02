        if (zSup[cellIndex] <= heights[cellGroup[cellIndex]-1]) 
        {// Full cell
            patchU[cellIndex] = 
                cnoidalFun::U
                (
                    waveHeight_, 
                    waterDepth_, 
                    mCnoidal_, 
                    waveKx, 
                    waveKy, 
                    wavePeriod_, 
                    xGroup[cellGroup[cellIndex]-1], 
                    yGroup[cellGroup[cellIndex]-1], 
                    currTime, 
                    patchHeight[cellIndex]
                );

            patchV[cellIndex] = patchU[cellIndex]*sin(waveAngle);
            patchU[cellIndex] = patchU[cellIndex]*cos(waveAngle);

            patchW[cellIndex] = 
                cnoidalFun::W
                (
                    waveHeight_, 
                    waterDepth_, 
                    mCnoidal_, 
                    waveKx, 
                    waveKy, 
                    wavePeriod_, 
                    xGroup[cellGroup[cellIndex]-1], 
                    yGroup[cellGroup[cellIndex]-1], 
                    currTime, 
                    patchHeight[cellIndex]
                );

        }
        else if (zInf[cellIndex] < heights[cellGroup[cellIndex]-1]) 
        {// Partial cell
            auxiliar = heights[cellGroup[cellIndex]-1] - zInf[cellIndex];
            auxiliarTotal = zSup[cellIndex]-zInf[cellIndex];
            auxiliarTotal = auxiliar/auxiliarTotal; // VOF of the cell

            // Height corresponding to the partial cell
            auxiliar = zInf[cellIndex] + auxiliar/scalar(2);

            patchU[cellIndex] = auxiliarTotal*
                cnoidalFun::U
                (
                    waveHeight_,
                    waterDepth_, 
                    mCnoidal_, 
                    waveKx, 
                    waveKy, 
                    wavePeriod_, 
                    xGroup[cellGroup[cellIndex]-1], 
                    yGroup[cellGroup[cellIndex]-1], 
                    currTime, 
                    auxiliar
                );

            patchV[cellIndex] = patchU[cellIndex]*sin(waveAngle);
            patchU[cellIndex] = patchU[cellIndex]*cos(waveAngle);

            patchW[cellIndex] = auxiliarTotal*
                cnoidalFun::W
                (
                    waveHeight_, 
                    waterDepth_, 
                    mCnoidal_, 
                    waveKx, 
                    waveKy, 
                    wavePeriod_, 
                    xGroup[cellGroup[cellIndex]-1], 
                    yGroup[cellGroup[cellIndex]-1], 
                    currTime, 
                    auxiliar
                );
        }
        else // Empty cell
        {
            patchU[cellIndex] = 0;
            patchV[cellIndex] = 0;
            patchW[cellIndex] = 0;
        }

