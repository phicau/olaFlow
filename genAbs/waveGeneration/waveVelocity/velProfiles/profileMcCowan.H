        if (zSup[cellIndex] <= calculatedLevel[cellGroup[cellIndex]-1]) 
        {// Full cell
            patchU[cellIndex] = 
                McCowanFun::U
                (
                    waveHeight_, 
                    waterDepth_, 
                    xGroup[cellGroup[cellIndex]-1], 
                    yGroup[cellGroup[cellIndex]-1], 
                    waveAngle, 
                    currTime, 
                    X0, 
                    zExtra,
                    Mmc,
                    Nmc
                );

            patchV[cellIndex] = patchU[cellIndex]*sin(waveAngle);
            patchU[cellIndex] = patchU[cellIndex]*cos(waveAngle);

            patchW[cellIndex] = 
                McCowanFun::W
                (
                    waveHeight_, 
                    waterDepth_, 
                    xGroup[cellGroup[cellIndex]-1], 
                    yGroup[cellGroup[cellIndex]-1], 
                    waveAngle, 
                    currTime, 
                    X0, 
                    zExtra,
                    Mmc,
                    Nmc
                );
        }
        else if (zInf[cellIndex] < calculatedLevel[cellGroup[cellIndex]-1])
        {// Partial cell
            auxiliar = 
                calculatedLevel[cellGroup[cellIndex]-1] - zInf[cellIndex];
            auxiliarTotal = zSup[cellIndex]-zInf[cellIndex];
            auxiliarTotal = auxiliar/auxiliarTotal; // VOF of the cell

            // Height corresponding to the partial cell
            auxiliar = zInf[cellIndex] + auxiliar/scalar(2);
            if (!extrapolation_)
            {
                auxiliar = zExtra;
            }

            patchU[cellIndex] = auxiliarTotal*
                McCowanFun::U
                (
                    waveHeight_, 
                    waterDepth_, 
                    xGroup[cellGroup[cellIndex]-1], 
                    yGroup[cellGroup[cellIndex]-1], 
                    waveAngle, 
                    currTime, 
                    X0, 
                    auxiliar,
                    Mmc,
                    Nmc
                );

            patchV[cellIndex] = patchU[cellIndex]*sin(waveAngle);
            patchU[cellIndex] = patchU[cellIndex]*cos(waveAngle);

            patchW[cellIndex] = auxiliarTotal*
                McCowanFun::W
                (
                    waveHeight_, 
                    waterDepth_, 
                    xGroup[cellGroup[cellIndex]-1], 
                    yGroup[cellGroup[cellIndex]-1], 
                    waveAngle, 
                    currTime, 
                    X0, 
                    auxiliar,
                    Mmc,
                    Nmc
                );
        }
        else // Empty cell
        {
            patchU[cellIndex] = 0;
            patchV[cellIndex] = 0;
            patchW[cellIndex] = 0;
        }
