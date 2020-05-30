        if (zSup[cellIndex] <= calculatedLevel[cellIndex]) 
        {// Full cell
            patchU[cellIndex] = 
                StokesIIFun::U
                (
                    waveHeight_, 
                    waterDepth_, 
                    waveKx, 
                    xComp[cellIndex], 
                    waveKy, 
                    yComp[cellIndex], 
                    waveOmega, 
                    currTime, 
                    wavePhase_, 
                    zExtra
                );

            patchV[cellIndex] = patchU[cellIndex]*Foam::sin(waveAngle);
            patchU[cellIndex] = patchU[cellIndex]*Foam::cos(waveAngle);

            patchW[cellIndex] = 
                StokesIIFun::W
                (
                    waveHeight_, 
                    waterDepth_, 
                    waveKx, 
                    xComp[cellIndex], 
                    waveKy, 
                    yComp[cellIndex], 
                    waveOmega, 
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
                StokesIIFun::U
                (
                    waveHeight_, 
                    waterDepth_, 
                    waveKx, 
                    xComp[cellIndex], 
                    waveKy, 
                    yComp[cellIndex], 
                    waveOmega, 
                    currTime, 
                    wavePhase_, 
                    auxiliar
                );

            patchV[cellIndex] = patchU[cellIndex]*Foam::sin(waveAngle);
            patchU[cellIndex] = patchU[cellIndex]*Foam::cos(waveAngle);

            patchW[cellIndex] = auxiliarTotal*
                StokesIIFun::W
                (
                    waveHeight_, 
                    waterDepth_, 
                    waveKx, 
                    xComp[cellIndex], 
                    waveKy, 
                    yComp[cellIndex], 
                    waveOmega, 
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