    if ( waveTheory_ == "Boussinesq" )
    {
        forAll(calculatedLevel, it2)
        {
            calculatedLevel[it2] = waterDepth_ + 
                BoussinesqFun::eta
                (
                    waveHeight_, 
                    waterDepth_, 
                    xGroup[it2], 
                    yGroup[it2], 
                    waveAngle, 
                    currTime, 
                    X0
                );
        }
    }
    else if (waveTheory_ == "Grimshaw")
    {
        forAll(calculatedLevel, it2)
        {
            calculatedLevel[it2] = waterDepth_ + 
                GrimshawFun::eta
                (
                    waveHeight_, 
                    waterDepth_, 
                    xGroup[it2], 
                    yGroup[it2], 
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
            << "Boussinesq, Grimshaw"
            << exit(FatalError);
    }
