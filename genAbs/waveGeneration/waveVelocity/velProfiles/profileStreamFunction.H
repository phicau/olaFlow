        if (zSup[cellIndex] <= heights[cellGroup[cellIndex]-1]) 
        {// Full cell
            patchU[cellIndex] = 0;
            patchW[cellIndex] = 0;

            forAll(Bjs_, iterSF2)
            {
                patchU[cellIndex] += (iterSF2 + 1)*Bjs_[iterSF2]
                    *cosh((iterSF2 + 1)*waveK*patchHeight[cellIndex]) 
                    /cosh((iterSF2 + 1)*waveK*waterDepth_)
                    *cos((iterSF2 + 1)*faseTot);
                patchW[cellIndex] += (iterSF2 + 1)*Bjs_[iterSF2]
                    *sinh((iterSF2 + 1)*waveK*patchHeight[cellIndex]) 
                    /cosh((iterSF2 + 1)*waveK*waterDepth_)
                    *sin((iterSF2 + 1)*faseTot);
            }

            patchU[cellIndex] = celerity - uMean_ 
                + sqrt(g*waterDepth_)*waveK*waterDepth_*patchU[cellIndex];

            patchV[cellIndex] = patchU[cellIndex]*sin(waveAngle);
            patchU[cellIndex] = patchU[cellIndex]*cos(waveAngle);

            patchW[cellIndex] = 
                sqrt(g*waterDepth_)*waveK*waterDepth_*patchW[cellIndex];        
        }
        else if (zInf[cellIndex] < heights[cellGroup[cellIndex]-1]) 
        {// Partial cell
            auxiliar = heights[cellGroup[cellIndex]-1] - zInf[cellIndex];
            auxiliarTotal = zSup[cellIndex]-zInf[cellIndex];
            auxiliarTotal = auxiliar/auxiliarTotal; // VOF of the cell

            // Height corresponding to the partial cell
            auxiliar = zInf[cellIndex] + auxiliar/scalar(2);

            patchU[cellIndex] = 0;
            patchW[cellIndex] = 0;

            forAll(Bjs_, iterSF3)
            {
                patchU[cellIndex] += (iterSF3 + 1)*Bjs_[iterSF3]
                    *cosh((iterSF3 + 1)*waveK*auxiliar)
                    /cosh((iterSF3 + 1)*waveK*waterDepth_)
                    *cos((iterSF3 + 1)*faseTot);
                patchW[cellIndex] += (iterSF3 + 1)*Bjs_[iterSF3]
                    *sinh((iterSF3 + 1)*waveK*auxiliar)
                    /cosh((iterSF3 + 1)*waveK*waterDepth_)
                    *sin((iterSF3 + 1)*faseTot);
            }

            patchU[cellIndex] = celerity - uMean_ 
                + sqrt(g*waterDepth_)*waveK*waterDepth_*patchU[cellIndex];

            patchV[cellIndex] = patchU[cellIndex]*sin(waveAngle);
            patchU[cellIndex] = patchU[cellIndex]*cos(waveAngle);

            patchW[cellIndex] = 
                sqrt(g*waterDepth_)*waveK*waterDepth_*patchW[cellIndex];  
        }
        else // Empty cell
        {
            patchU[cellIndex] = 0;
            patchV[cellIndex] = 0;
            patchW[cellIndex] = 0;
        }