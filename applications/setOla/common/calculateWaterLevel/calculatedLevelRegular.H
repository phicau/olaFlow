
    if ( waveTheory_ == "StokesI" )
    {
        forAll(calculatedLevel, itS1)
        {
            calculatedLevel[itS1] += timeMult*
                StokesIFun::eta
                (
                    waveHeight_, 
                    waveKx, 
                    xComp[itS1], 
                    waveKy, 
                    yComp[itS1], 
                    waveOmega, 
                    currTime, 
                    wavePhase_
                );
        }
    }
    else if ( waveTheory_ == "StokesII" )
    {
        forAll(calculatedLevel, itS2)
        {
            calculatedLevel[itS2] += timeMult*
                StokesIIFun::eta
                (
                    waveHeight_, 
                    waterDepth_, 
                    waveKx, 
                    xComp[itS2], 
                    waveKy, 
                    yComp[itS2], 
                    waveOmega, 
                    currTime, 
                    wavePhase_
                );
        }
    }
    else if ( waveTheory_ == "StokesIII" )
    {
        forAll(calculatedLevel, itS2)
        {
            calculatedLevel[itS2] += timeMult*
                StokesIIIFun::eta
                (
                    waveHeight_, 
                    waterDepth_, 
                    waveKx, 
                    xComp[itS2], 
                    waveKy, 
                    yComp[itS2], 
                    waveOmega, 
                    currTime, 
                    wavePhase_,
                    aE_SIII_,
                    klE_SIII_
                );
        }
    }
    else  if ( waveTheory_ == "StokesV" )
    {
        forAll(calculatedLevel, it1)
        {
            calculatedLevel[it1] += timeMult*
                stokesVFun::eta
                (
                    waterDepth_, 
                    waveKx, 
                    waveKy, 
                    lambdaStokesV_, 
                    wavePeriod_, 
                    xComp[it1], 
                    yComp[it1], 
                    currTime, 
                    wavePhase_
                );
        }
    }
    else if ( waveTheory_ == "cnoidal" )
    {
        forAll(calculatedLevel, it3)
        {
            calculatedLevel[it3] += timeMult*
                cnoidalFun::eta
                (
                    waveHeight_, 
                    mCnoidal_, 
                    waveKx, 
                    waveKy, 
                    wavePeriod_, 
                    xComp[it3], 
                    yComp[it3], 
                    currTime
                );
        }
    }
    else if ( waveTheory_ == "streamFunction" )
    {
        scalar strfnAux = 0.0;
        forAll(calculatedLevel, it4)
        {
            strfnAux = 0.0;

            forAll(Ejs_, iterSF)
            {
                strfnAux +=  Ejs_[iterSF]*Foam::cos((iterSF + 1)*(waveKx*xComp[it4] 
                    + waveKy*yComp[it4] - waveOmega*currTime + wavePhase_));
            }

            calculatedLevel[it4] += timeMult*waterDepth_*strfnAux;
        }
    }
    else 
    {
        FatalError
            << "Wave theory not supported, use:\n"
            << "StokesI, StokesII, StokesV, cnoidal, "
            << "streamFunction"
            << exit(FatalError);
    }
