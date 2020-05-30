        if (zSup[cellIndex] <= calculatedLevel[cellIndex]) 
        {// Full cell
            forAll(waveHeights_, ln2)
            {
                faseTot = 
                    waveKxs[ln2]*xComp[cellIndex] 
                    + waveKys[ln2]*yComp[cellIndex] 
                    - waveOmegas[ln2]*currTime + wavePhases_[ln2];
                
                patchU[cellIndex] +=
                    //Foam::sqr(Foam::cos( limAngle(waveAngles[ln2] 
                    //- meanAngle[cellIndex]) ))* // spreading function
                    waveHeights_[ln2]*0.5*waveOmegas[ln2]*
                    Foam::cos(faseTot)*
                    Foam::cosh(waveKs[ln2]*zExtra)/
                    Foam::sinh(waveKs[ln2]*waterDepth_)*Foam::cos(waveAngles[ln2]);

                patchV[cellIndex] +=
                    //Foam::sqr(Foam::cos( limAngle(waveAngles[ln2] 
                    //- meanAngle[cellIndex]) ))* // spreading function
                    waveHeights_[ln2]*0.5*waveOmegas[ln2]*
                    Foam::cos(faseTot)*
                    Foam::cosh(waveKs[ln2]*zExtra)/
                    Foam::sinh(waveKs[ln2]*waterDepth_)*Foam::sin(waveAngles[ln2]);

                patchW[cellIndex] +=
                    //Foam::sqr(Foam::cos( limAngle(waveAngles[ln2] 
                    //- meanAngle[cellIndex]) ))* // spreading function
                    waveHeights_[ln2]*0.5*waveOmegas[ln2]*
                    Foam::sin(faseTot)*
                    Foam::sinh(waveKs[ln2]*zExtra)/
                    Foam::sinh(waveKs[ln2]*waterDepth_); 
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
                                xComp[cellIndex], 
                                yComp[cellIndex], 
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
                                xComp[cellIndex], 
                                yComp[cellIndex], 
                                currTime, 
                                waterDepth_, 
                                auxiliar
                            );
                    }
                }

                patchU[cellIndex] += auxiliarSO * Foam::cos(waveAngles[0]);
                patchV[cellIndex] += auxiliarSO * Foam::sin(waveAngles[0]);
            }
        }
        else if (zInf[cellIndex] < calculatedLevel[cellIndex]) 
        {// Partial cell
            auxiliar = calculatedLevel[cellIndex] - zInf[cellIndex];
            auxiliarTotal = zSup[cellIndex]-zInf[cellIndex];
            auxiliarTotal = auxiliar/auxiliarTotal; // VOF of the cell

            // Height corresponding to the partial cell
            auxiliar = zInf[cellIndex] + auxiliar/2.0;
            if (!extrapolation_)
            {
                auxiliar = zExtra;
            }

            forAll(waveHeights_, ln3)
            {
                faseTot = 
                    waveKxs[ln3]*xComp[cellIndex] 
                    + waveKys[ln3]*yComp[cellIndex] 
                    - waveOmegas[ln3]*currTime + wavePhases_[ln3];

                patchU[cellIndex] += 
                    //Foam::sqr(Foam::cos( limAngle(waveAngles[ln3] 
                    //- meanAngle[cellIndex]) ))* // spreading function
                    auxiliarTotal*
                    waveHeights_[ln3]*0.5*waveOmegas[ln3]*
                    Foam::cos(faseTot)*
                    Foam::cosh(waveKs[ln3]*auxiliar)/
                    Foam::sinh(waveKs[ln3]*waterDepth_)*Foam::cos(waveAngles[ln3]);

                patchV[cellIndex] += 
                    //Foam::sqr(Foam::cos( limAngle(waveAngles[ln3] 
                    //- meanAngle[cellIndex]) ))* // spreading function
                    auxiliarTotal*
                    waveHeights_[ln3]*0.5*waveOmegas[ln3]*
                    Foam::cos(faseTot)*
                    Foam::cosh(waveKs[ln3]*auxiliar)/
                    Foam::sinh(waveKs[ln3]*waterDepth_)*Foam::sin(waveAngles[ln3]);

                patchW[cellIndex] += 
                    //Foam::sqr(Foam::cos( limAngle(waveAngles[ln3] 
                    //- meanAngle[cellIndex]) ))* // spreading function
                    auxiliarTotal*
                    waveHeights_[ln3]*0.5*waveOmegas[ln3]*
                    Foam::sin(faseTot)*
                    Foam::sinh(waveKs[ln3]*auxiliar)/
                    Foam::sinh(waveKs[ln3]*waterDepth_);  
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
                                xComp[cellIndex], 
                                yComp[cellIndex], 
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
                                xComp[cellIndex], 
                                yComp[cellIndex], 
                                currTime, 
                                waterDepth_, 
                                auxiliar
                            );
                    }
                }
                patchU[cellIndex] += 
                    auxiliarTotal*auxiliarSO*Foam::cos(waveAngles[0]);
                patchV[cellIndex] += 
                    auxiliarTotal*auxiliarSO*Foam::sin(waveAngles[0]);
            }
        }
        else // Empty cell
        {
            patchU[cellIndex] = uWind_.component(0);
            patchV[cellIndex] = uWind_.component(1);
            patchW[cellIndex] = uWind_.component(2);
        }