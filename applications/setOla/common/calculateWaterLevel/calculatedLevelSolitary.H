    if ( waveTheory_ == "Boussinesq" )
    {
        forAll(calculatedLevel, it2)
        {
            calculatedLevel[it2] += 
                BoussinesqFun::eta
                (
                    waveHeight_, 
                    waterDepth_, 
                    xComp[it2], 
                    yComp[it2], 
                    waveAngle, 
                    currTime, 
                    X0
                );
        }
    }
    else if (waveTheory_ == "McCowan")
    {
        forAll(calculatedLevel, it2)
        {
            calculatedLevel[it2] += 
                McCowanFun::eta
                (
                    waveHeight_, 
                    waterDepth_, 
                    xComp[it2], 
                    yComp[it2], 
                    waveAngle, 
                    currTime, 
                    X0,
                    Mmc,
                    Nmc
                );
        }
    }
    else if (waveTheory_ == "Grimshaw")
    {
        forAll(calculatedLevel, it2)
        {
            calculatedLevel[it2] += 
                GrimshawFun::eta
                (
                    waveHeight_, 
                    waterDepth_, 
                    xComp[it2], 
                    yComp[it2], 
                    waveAngle, 
                    currTime, 
                    X0
                );
        }
    }
    else 
    {
        FatalError
            << "Wave theory not supported, use:\n"
            << "Boussinesq, McCowan, Grimshaw"
            << exit(FatalError);
    }
