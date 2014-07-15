    forAll(calculatedLevel, indWL)
    {
        calculatedLevel[indWL] = waterDepth_;

        forAll(waveHeights_, lnd)
        {
            calculatedLevel[indWL] += 
                timeMult*
                StokesIFun::eta
                (
                    waveHeights_[lnd], 
                    waveKxs[lnd], 
                    xGroup[indWL], 
                    waveKys[lnd], 
                    yGroup[indWL], 
                    waveOmegas[lnd], 
                    currTime, 
                    wavePhases_[lnd]
                );
        }

        if ( secondOrder_ )
        {
            for ( label so1 = 0; so1 <= nComp_-2; so1++ )
            {
                for ( label so2 = so1+1; so2 <= nComp_-1; so2++ )
                {
                    calculatedLevel[indWL] += 
                        timeMult*
                        secondOrderFun::etaSO
                        (
                            waveHeights_[so1], 
                            waveHeights_[so2], 
                            waveOmegas[so1], 
                            waveOmegas[so2], 
                            wavePhases_[so1], 
                            wavePhases_[so2], 
                            waveKxs[so1], 
                            waveKys[so1], 
                            waveKxs[so2], 
                            waveKys[so2], 
                            xGroup[indWL], 
                            yGroup[indWL], 
                            currTime, 
                            waterDepth_
                        );
                }
            }
        }
    }
