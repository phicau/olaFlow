        if (zSup[cellIndex] <= calculatedLevel[cellIndex]) 
        {// Full cell
            patchU[cellIndex] = 0;
            patchW[cellIndex] = 0;

            forAll(Bjs_, iterSF2)
            {
                patchU[cellIndex] += (iterSF2 + 1)*Bjs_[iterSF2]
                    *Foam::cosh((iterSF2 + 1)*waveK*zExtra) 
                    /Foam::cosh((iterSF2 + 1)*waveK*waterDepth_)
                    *Foam::cos((iterSF2 + 1)*faseTot);
                patchW[cellIndex] += (iterSF2 + 1)*Bjs_[iterSF2]
                    *Foam::sinh((iterSF2 + 1)*waveK*zExtra) 
                    /Foam::cosh((iterSF2 + 1)*waveK*waterDepth_)
                    *Foam::sin((iterSF2 + 1)*faseTot);
            }

            patchU[cellIndex] = celerity - uMean_ 
                + Foam::sqrt(g.value()&vector(0,0,-1)*waterDepth_)*waveK*waterDepth_*patchU[cellIndex];

            patchV[cellIndex] = patchU[cellIndex]*Foam::sin(waveAngle);
            patchU[cellIndex] = patchU[cellIndex]*Foam::cos(waveAngle);

            patchW[cellIndex] = 
                Foam::sqrt(g.value()&vector(0,0,-1)*waterDepth_)*waveK*waterDepth_*patchW[cellIndex];        
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

            patchU[cellIndex] = 0;
            patchW[cellIndex] = 0;

            forAll(Bjs_, iterSF3)
            {
                patchU[cellIndex] += (iterSF3 + 1)*Bjs_[iterSF3]
                    *Foam::cosh((iterSF3 + 1)*waveK*auxiliar)
                    /Foam::cosh((iterSF3 + 1)*waveK*waterDepth_)
                    *Foam::cos((iterSF3 + 1)*faseTot);
                patchW[cellIndex] += (iterSF3 + 1)*Bjs_[iterSF3]
                    *Foam::sinh((iterSF3 + 1)*waveK*auxiliar)
                    /Foam::cosh((iterSF3 + 1)*waveK*waterDepth_)
                    *Foam::sin((iterSF3 + 1)*faseTot);
            }

            patchU[cellIndex] = celerity - uMean_ 
                + Foam::sqrt(g.value()&vector(0,0,-1)*waterDepth_)*waveK*waterDepth_*patchU[cellIndex];

            patchV[cellIndex] = patchU[cellIndex]*Foam::sin(waveAngle);
            patchU[cellIndex] = patchU[cellIndex]*Foam::cos(waveAngle);

            patchW[cellIndex] = 
                Foam::sqrt(g.value()&vector(0,0,-1)*waterDepth_)*waveK*waterDepth_*patchW[cellIndex];  
        }
        else // Empty cell
        {
            patchU[cellIndex] = uWind_.component(0);
            patchV[cellIndex] = uWind_.component(1);
            patchW[cellIndex] = uWind_.component(2);
        }