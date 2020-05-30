        if (zSup[cellIndex] <= calculatedLevel[cellIndex]) 
        {// Full cell
            patchU[cellIndex] = 
                BoussinesqFun::U
                (
                    waveHeight_, 
                    waterDepth_, 
                    xComp[cellIndex], 
                    yComp[cellIndex], 
                    waveAngle, 
                    currTime, 
                    X0, 
                    zExtra
                );

            patchV[cellIndex] = patchU[cellIndex]*Foam::sin(waveAngle);
            patchU[cellIndex] = patchU[cellIndex]*Foam::cos(waveAngle);

            patchW[cellIndex] = 
                BoussinesqFun::W
                (
                    waveHeight_, 
                    waterDepth_, 
                    xComp[cellIndex], 
                    yComp[cellIndex], 
                    waveAngle, 
                    currTime, 
                    X0, 
                    zExtra
                );
        }
        else if (zInf[cellIndex] < calculatedLevel[cellIndex])
        {// Partial cell
            auxiliar = 
                calculatedLevel[cellIndex] - zInf[cellIndex];
            auxiliarTotal = zSup[cellIndex]-zInf[cellIndex];
            auxiliarTotal = auxiliar/auxiliarTotal; // VOF of the cell

            // Height corresponding to the partial cell
            auxiliar = zInf[cellIndex] + auxiliar/scalar(2);
            if (!extrapolation_)
            {
                auxiliar = zExtra;
            }

            patchU[cellIndex] = auxiliarTotal*
                BoussinesqFun::U
                (
                    waveHeight_, 
                    waterDepth_, 
                    xComp[cellIndex], 
                    yComp[cellIndex], 
                    waveAngle, 
                    currTime, 
                    X0, 
                    auxiliar
                );

            patchV[cellIndex] = patchU[cellIndex]*Foam::sin(waveAngle);
            patchU[cellIndex] = patchU[cellIndex]*Foam::cos(waveAngle);

            patchW[cellIndex] = auxiliarTotal*
                BoussinesqFun::W
                (
                    waveHeight_, 
                    waterDepth_, 
                    xComp[cellIndex], 
                    yComp[cellIndex], 
                    waveAngle, 
                    currTime, 
                    X0, 
                    auxiliar
                );
        }
        else // Empty cell
        {
            patchU[cellIndex] = uWind_.component(0);
            patchV[cellIndex] = uWind_.component(1);
            patchW[cellIndex] = uWind_.component(2);
        }
