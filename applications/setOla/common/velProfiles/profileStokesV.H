        if (zSup[cellIndex] <= calculatedLevel[cellIndex]) 
        {// Full cell
            patchU[cellIndex] = 
                stokesVFun::U
                (
                    waterDepth_, 
                    waveKx, 
                    waveKy, 
                    lambdaStokesV_, 
                    wavePeriod_, 
                    xComp[cellIndex], 
                    yComp[cellIndex], 
                    currTime, 
                    wavePhase_, 
                    zExtra
                );

            patchV[cellIndex] = patchU[cellIndex]*Foam::sin(waveAngle);
            patchU[cellIndex] = patchU[cellIndex]*Foam::cos(waveAngle);

            patchW[cellIndex] =
                stokesVFun::V
                (
                    waterDepth_, 
                    waveKx, 
                    waveKy, 
                    lambdaStokesV_, 
                    wavePeriod_, 
                    xComp[cellIndex], 
                    yComp[cellIndex], 
                    currTime, 
                    wavePhase_, 
                    zExtra
                );
        }
        else if (zInf[cellIndex] < calculatedLevel[cellIndex]) 
        {// Partial cell
            auxiliar = calculatedLevel[cellIndex] - zInf[cellIndex];
            auxiliarTotal = zSup[cellIndex]-zInf[cellIndex];
            auxiliarTotal = auxiliar/auxiliarTotal; // VOF of the cell

            // Height corresponding to the partial cell
            auxiliar = zInf[cellIndex] + auxiliar/scalar(2);
            if (!extrapolation_)
            {
                auxiliar = zExtra;
            }

            patchU[cellIndex] = auxiliarTotal*
                stokesVFun::U
                (
                    waterDepth_, 
                    waveKx, 
                    waveKy, 
                    lambdaStokesV_, 
                    wavePeriod_, 
                    xComp[cellIndex], 
                    yComp[cellIndex], 
                    currTime, 
                    wavePhase_, 
                    auxiliar
                );

            patchV[cellIndex] = patchU[cellIndex]*Foam::sin(waveAngle);
            patchU[cellIndex] = patchU[cellIndex]*Foam::cos(waveAngle);

            patchW[cellIndex] = auxiliarTotal*
                stokesVFun::V
                (
                    waterDepth_, 
                    waveKx, 
                    waveKy, 
                    lambdaStokesV_, 
                    wavePeriod_, 
                    xComp[cellIndex], 
                    yComp[cellIndex], 
                    currTime, 
                    wavePhase_, 
                    auxiliar
                );
        }
        else // Empty cell
        {
            patchU[cellIndex] = uWind_.component(0);
            patchV[cellIndex] = uWind_.component(1);
            patchW[cellIndex] = uWind_.component(2);
        }