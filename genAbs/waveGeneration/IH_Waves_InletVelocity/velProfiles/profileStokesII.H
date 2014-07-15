        if (zSup[cellIndex] <= heights[cellGroup[cellIndex]-1]) 
        {// Full cell
            patchU[cellIndex] = 
                StokesIIFun::U
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
                    patchHeight[cellIndex]
                );

            patchV[cellIndex] = patchU[cellIndex]*sin(waveAngle);
            patchU[cellIndex] = patchU[cellIndex]*cos(waveAngle);

            patchW[cellIndex] = 
                StokesIIFun::W
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
                StokesIIFun::U
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
                    auxiliar
                );

            patchV[cellIndex] = patchU[cellIndex]*sin(waveAngle);
            patchU[cellIndex] = patchU[cellIndex]*cos(waveAngle);

            patchW[cellIndex] = auxiliarTotal*
                StokesIIFun::W
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
                    auxiliar
                ); 
        }
        else // Empty cell
        {
            patchU[cellIndex] = 0;
            patchV[cellIndex] = 0;
            patchW[cellIndex] = 0;
        }