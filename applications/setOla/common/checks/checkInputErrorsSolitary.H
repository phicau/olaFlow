    // Extracting values from dict
    waveTheory_ = (waveDict.lookupOrDefault<word>("waveTheory", "aaa"));
    waveHeight_ = (waveDict.lookupOrDefault<scalar>("waveHeight", -1));
    waveDir_ = (waveDict.lookupOrDefault<scalar>("waveDir", 0));
    genAbs_ = (waveDict.lookupOrDefault<bool>("genAbs", false ));
    nPaddles_ = (waveDict.lookupOrDefault<label>("nPaddles", 1));
    nSolitaryWaves_ = (waveDict.lookupOrDefault<label>("nSolitaryWaves", 1));

    // Check if the value of nPaddles is correct for the number of columns
    if (nPaddles_ < 1)
    {
        FatalError
            << "Check nPaddles value."
            << exit(FatalError);
    }

    if ( waveHeight_ < 0.0 )
    {
        FatalError
            << "Check wave height value."
            << exit(FatalError);
    }

    if ( nSolitaryWaves_ < 1 )
    {
        FatalError
            << "Check nSolitaryWaves value."
            << exit(FatalError);
    }

    if ( waveTheory_ == "aaa" )    // No wave theory specified
    {
        FatalError
            << "Wave theory not specified, use:\n"
            << "Boussinesq, McCowan, Grimshaw"
            << exit(FatalError);
    }
    else if
    ( 
        waveTheory_ == "Boussinesq" ||
        waveTheory_ == "McCowan" ||
        waveTheory_ == "Grimshaw"
    )
    {
        if ( waveTheory_ == "Boussinesq" )
        {
            waveLength_ = BoussinesqFun::waveLength(waveHeight_, waterDepth_);
            wavePeriod_ = BoussinesqFun::wavePeriod(waveHeight_, waterDepth_);
        }
        else if ( waveTheory_ == "McCowan" )
        {
            scalar Mtemp = McCowanFun::Mcalc(waveHeight_, waterDepth_);
            waveLength_ = McCowanFun::waveLength(waveHeight_, waterDepth_);
            wavePeriod_ = McCowanFun::wavePeriod(waveHeight_, waterDepth_, Mtemp);
        }
        else if ( waveTheory_ == "Grimshaw" )
        {
            waveLength_ = GrimshawFun::waveLength(waveHeight_, waterDepth_);
            wavePeriod_ = GrimshawFun::wavePeriod(waveHeight_, waterDepth_);
        }

        Info << "\nWave Generation BC" << endl;
        Info << "Wave theory: " << waveTheory_ << endl;
        Info << "H: " << waveHeight_ << endl;
        Info << "h: " << waterDepth_ << endl;
        Info << "c: " << waveLength_/wavePeriod_ << endl;
        Info << "Equivalent T: " << wavePeriod_ << endl;
        Info << "Equivalent L: " << waveLength_ << endl;
        Info << "Direction: " << waveDir_ << "º" <<  "\n\n" << endl;
    }
    else
    {
        FatalError
            << "Wave theory not supported, use:\n"
            << "Boussinesq, McCowan, Grimshaw"
            << exit(FatalError);
    }
