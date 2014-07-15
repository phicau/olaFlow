
    if ( waveTheory_ == "StokesI" )
    {
        forAll(calculatedLevel, itS1)
        {
            calculatedLevel = waterDepth_ + timeMult*
                StokesIFun::eta
                (
                    waveHeight_, 
                    waveKx, 
                    xGroup[itS1], 
                    waveKy, 
                    yGroup[itS1], 
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
            calculatedLevel = waterDepth_ + timeMult*
                StokesIIFun::eta
                (
                    waveHeight_, 
                    waterDepth_, 
                    waveKx, 
                    xGroup[itS2], 
                    waveKy, 
                    yGroup[itS2], 
                    waveOmega, 
                    currTime, 
                    wavePhase_
                );
        }
    }
    else  if ( waveTheory_ == "StokesV" )
    {
        forAll(calculatedLevel, it1)
        {
            calculatedLevel[it1] = waterDepth_ + timeMult*
                stokesVFun::eta
                (
                    waterDepth_, 
                    waveKx, 
                    waveKy, 
                    lambdaStokesV_, 
                    wavePeriod_, 
                    xGroup[it1], 
                    yGroup[it1], 
                    currTime, 
                    wavePhase_
                );
        }
    }
    else if ( waveTheory_ == "cnoidal" )
    {
        forAll(calculatedLevel, it3)
        {
            calculatedLevel[it3] = waterDepth_ + timeMult*
                cnoidalFun::eta
                (
                    waveHeight_, 
                    mCnoidal_, 
                    waveKx, 
                    waveKy, 
                    wavePeriod_, 
                    xGroup[it3], 
                    yGroup[it3], 
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
                strfnAux +=  Ejs_[iterSF]*cos((iterSF + 1)*(waveKx*xGroup[it4] 
                    + waveKy*yGroup[it4] - waveOmega*currTime + wavePhase_));
            }

            calculatedLevel[it4] = waterDepth_ + timeMult*waterDepth_*strfnAux;
        }
    }
    else 
    {
        FatalError
            << "Wave theory not supported, use:\n"
            << "StokesI, StokesII, StokesV, cnoidal, "
            << "solitaryBoussinesq, streamFunction"
            << exit(FatalError);
    }
