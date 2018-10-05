        if (zSup[cellIndex] <= heights[cellGroup[cellIndex]-1]) 
        {// Full cell
            patchU[cellIndex] = 
                StokesIIIFun::U
                (
                    waveHeight_, 
                    waterDepth_, 
                    waveKx, 
                    xGroup[cellGroup[cellIndex]-1], 
                    waveKy, 
                    yGroup[cellGroup[cellIndex]-1], 
                    waveOmega, 
                    currTime, 
                    wavePhase_, 
                    zExtra,
                    aE_SIII_,
                    klE_SIII_
                );

            patchV[cellIndex] = patchU[cellIndex]*sin(waveAngle);
            patchU[cellIndex] = patchU[cellIndex]*cos(waveAngle);

            patchW[cellIndex] = 
                StokesIIIFun::W
                (
                    waveHeight_, 
                    waterDepth_, 
                    waveKx, 
                    xGroup[cellGroup[cellIndex]-1], 
                    waveKy, 
                    yGroup[cellGroup[cellIndex]-1], 
                    waveOmega, 
                    currTime, 
                    wavePhase_, 
                    zExtra,
                    aE_SIII_,
                    klE_SIII_
                ); 
        }
        else if (zInf[cellIndex] < heights[cellGroup[cellIndex]-1]) 
        {// Partial cell
            auxiliar = heights[cellGroup[cellIndex]-1] - zInf[cellIndex];
            auxiliarTotal = zSup[cellIndex]-zInf[cellIndex];
            auxiliarTotal = auxiliar/auxiliarTotal; // VOF of the cell

            // Height corresponding to the partial cell
            auxiliar = zInf[cellIndex] + auxiliar/scalar(2);
            if (!extrapolation_)
            {
                auxiliar = zExtra;
            }

            patchU[cellIndex] = auxiliarTotal*
                StokesIIIFun::U
                (
                    waveHeight_, 
                    waterDepth_, 
                    waveKx, 
                    xGroup[cellGroup[cellIndex]-1], 
                    waveKy, 
                    yGroup[cellGroup[cellIndex]-1], 
                    waveOmega, 
                    currTime, 
                    wavePhase_, 
                    auxiliar,
                    aE_SIII_,
                    klE_SIII_
                );

            patchV[cellIndex] = patchU[cellIndex]*sin(waveAngle);
            patchU[cellIndex] = patchU[cellIndex]*cos(waveAngle);

            patchW[cellIndex] = auxiliarTotal*
                StokesIIIFun::W
                (
                    waveHeight_, 
                    waterDepth_, 
                    waveKx, 
                    xGroup[cellGroup[cellIndex]-1], 
                    waveKy, 
                    yGroup[cellGroup[cellIndex]-1], 
                    waveOmega, 
                    currTime, 
                    wavePhase_, 
                    auxiliar,
                    aE_SIII_,
                    klE_SIII_
                ); 
        }
        else // Empty cell
        {
            patchU[cellIndex] = 0;
            patchV[cellIndex] = 0;
            patchW[cellIndex] = 0;
        }