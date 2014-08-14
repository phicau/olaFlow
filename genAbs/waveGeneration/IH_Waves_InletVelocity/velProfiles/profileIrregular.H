        if (zSup[cellIndex] <= heights[cellGroup[cellIndex]-1]) 
        {// Full cell
            forAll(waveHeights_, ln2)
            {
                faseTot = 
                    waveKxs[ln2]*xGroup[cellGroup[cellIndex]-1] 
                    + waveKys[ln2]*yGroup[cellGroup[cellIndex]-1] 
                    - waveOmegas[ln2]*currTime + wavePhases_[ln2];
                
                // Sin extrapolacion
                //auxiliar = 0.0;
                //if ( patchHeight[cellIndex] > waterDepth_ )
                //{
                //    auxiliar = waterDepth_;
                //}
                //else
                //{
                    auxiliar = patchHeight[cellIndex];
                //}

                patchU[cellIndex] +=
                    sqr(cos( limAngle(waveAngles[ln2] 
                    - meanAngle[cellGroup[cellIndex]-1]) ))* // spreading function
                    waveHeights_[ln2]*0.5*waveOmegas[ln2]*
                    cos(faseTot)*
                    cosh(waveKs[ln2]*auxiliar)/
                    sinh(waveKs[ln2]*waterDepth_)*cos(waveAngles[ln2]);

                patchV[cellIndex] +=
                    sqr(cos( limAngle(waveAngles[ln2] 
                    - meanAngle[cellGroup[cellIndex]-1]) ))* // spreading function
                    waveHeights_[ln2]*0.5*waveOmegas[ln2]*
                    cos(faseTot)*
                    cosh(waveKs[ln2]*auxiliar)/
                    sinh(waveKs[ln2]*waterDepth_)*sin(waveAngles[ln2]);

                patchW[cellIndex] +=
                    sqr(cos( limAngle(waveAngles[ln2] 
                    - meanAngle[cellGroup[cellIndex]-1]) ))* // spreading function
                    waveHeights_[ln2]*0.5*waveOmegas[ln2]*
                    sin(faseTot)*
                    sinh(waveKs[ln2]*auxiliar)/
                    sinh(waveKs[ln2]*waterDepth_); 
            }

            if ( secondOrder_ )
            {
                auxiliarSO = 0.0;

                for ( label so11 = 0; so11 <= nComp_-2; so11++ )
                {
                    for ( label so21 = so11+1; so21 <= nComp_-1; so21++ )
                    {
                        auxiliarSO += 
                            secondOrderFun::uSO
                            (
                                waveHeights_[so11], 
                                waveHeights_[so21], 
                                waveOmegas[so11], 
                                waveOmegas[so21], 
                                wavePhases_[so11], 
                                wavePhases_[so21], 
                                waveKxs[so11], 
                                waveKys[so11], 
                                waveKxs[so21], 
                                waveKys[so21], 
                                xGroup[cellGroup[cellIndex]-1], 
                                yGroup[cellGroup[cellIndex]-1], 
                                currTime, 
                                waterDepth_, 
                                auxiliar
                            );

                        patchW[cellIndex] += 
                            secondOrderFun::wSO
                            (
                                waveHeights_[so11], 
                                waveHeights_[so21], 
                                waveOmegas[so11], 
                                waveOmegas[so21], 
                                wavePhases_[so11], 
                                wavePhases_[so21], 
                                waveKxs[so11], 
                                waveKys[so11], 
                                waveKxs[so21], 
                                waveKys[so21], 
                                xGroup[cellGroup[cellIndex]-1], 
                                yGroup[cellGroup[cellIndex]-1], 
                                currTime, 
                                waterDepth_, 
                                auxiliar
                            );
                    }
                }

                patchU[cellIndex] += auxiliarSO * cos(waveAngles[0]);
                patchV[cellIndex] += auxiliarSO * sin(waveAngles[0]);
            }
        }
        else if (zInf[cellIndex] < heights[cellGroup[cellIndex]-1]) 
        {// Partial cell
            auxiliar = heights[cellGroup[cellIndex]-1] - zInf[cellIndex];
            auxiliarTotal = zSup[cellIndex]-zInf[cellIndex];
            auxiliarTotal = auxiliar/auxiliarTotal; // VOF of the cell

            // Height corresponding to the partial cell
            auxiliar = zInf[cellIndex] + auxiliar/2.0;

            // Sin extrapolacion
            //if ( auxiliar > waterDepth_ )
            //{
            //    auxiliar = waterDepth_;
            //}

            forAll(waveHeights_, ln3)
            {
                faseTot = 
                    waveKxs[ln3]*xGroup[cellGroup[cellIndex]-1] 
                    + waveKys[ln3]*yGroup[cellGroup[cellIndex]-1] 
                    - waveOmegas[ln3]*currTime + wavePhases_[ln3];

                patchU[cellIndex] += 
                    sqr(cos( limAngle(waveAngles[ln3] 
                    - meanAngle[cellGroup[cellIndex]-1]) ))* // spreading function
                    auxiliarTotal*
                    waveHeights_[ln3]*0.5*waveOmegas[ln3]*
                    cos(faseTot)*
                    cosh(waveKs[ln3]*auxiliar)/
                    sinh(waveKs[ln3]*waterDepth_)*cos(waveAngles[ln3]);

                patchV[cellIndex] += 
                    sqr(cos( limAngle(waveAngles[ln3] 
                    - meanAngle[cellGroup[cellIndex]-1]) ))* // spreading function
                    auxiliarTotal*
                    waveHeights_[ln3]*0.5*waveOmegas[ln3]*
                    cos(faseTot)*
                    cosh(waveKs[ln3]*auxiliar)/
                    sinh(waveKs[ln3]*waterDepth_)*sin(waveAngles[ln3]);

                patchW[cellIndex] += 
                    sqr(cos( limAngle(waveAngles[ln3] 
                    - meanAngle[cellGroup[cellIndex]-1]) ))* // spreading function
                    auxiliarTotal*
                    waveHeights_[ln3]*0.5*waveOmegas[ln3]*
                    sin(faseTot)*
                    sinh(waveKs[ln3]*auxiliar)/
                    sinh(waveKs[ln3]*waterDepth_);  
            }

            if ( secondOrder_ )
            {
                auxiliarSO = 0.0;

                for ( label so12 = 0; so12 <= nComp_-2; so12++ )
                {
                    for ( label so22 = so12+1; so22 <= nComp_-1; so22++ )
                    {
                        auxiliarSO += 
                            secondOrderFun::uSO
                            (
                                waveHeights_[so12], 
                                waveHeights_[so22], 
                                waveOmegas[so12], 
                                waveOmegas[so22], 
                                wavePhases_[so12], 
                                wavePhases_[so22], 
                                waveKxs[so12], 
                                waveKys[so12], 
                                waveKxs[so22], 
                                waveKys[so22], 
                                xGroup[cellGroup[cellIndex]-1], 
                                yGroup[cellGroup[cellIndex]-1], 
                                currTime, 
                                waterDepth_, 
                                auxiliar
                            );

                        patchW[cellIndex] += auxiliarTotal*
                            secondOrderFun::wSO
                            (
                                waveHeights_[so12], 
                                waveHeights_[so22], 
                                waveOmegas[so12], 
                                waveOmegas[so22], 
                                wavePhases_[so12], 
                                wavePhases_[so22], 
                                waveKxs[so12], 
                                waveKys[so12], 
                                waveKxs[so22], 
                                waveKys[so22], 
                                xGroup[cellGroup[cellIndex]-1], 
                                yGroup[cellGroup[cellIndex]-1], 
                                currTime, 
                                waterDepth_, 
                                auxiliar
                            );
                    }
                }
                patchU[cellIndex] += 
                    auxiliarTotal*auxiliarSO*cos(waveAngles[0]);
                patchV[cellIndex] += 
                    auxiliarTotal*auxiliarSO*sin(waveAngles[0]);
            }
        }
        else // Empty cell
        {
            patchU[cellIndex] = 0;
            patchV[cellIndex] = 0;
            patchW[cellIndex] = 0;
        }