    // Extracting values from dict
    waveTheory_ = (waveDict.lookupOrDefault<word>("waveTheory", "aaa"));
    waveHeight_ = (waveDict.lookupOrDefault<scalar>("waveHeight", -1));
    waveDir_ = (waveDict.lookupOrDefault<scalar>("waveDir", 0));
    genAbs_ = (waveDict.lookupOrDefault<bool>("genAbs", false ));
    nPaddles_ = (waveDict.lookupOrDefault<label>("nPaddles", 1));
    absDir_ = (waveDict.lookupOrDefault<scalar>("absDir", 0));

    // Check if the value of nPaddles is correct for the number of columns
    if (nPaddles_ < 1)
    {
        FatalError
            << "Check nPaddles value."
            << exit(FatalError);
    }

    if ( nPaddles_ > 1 )
    {
        nPaddles_ = decreaseNPaddles( nPaddles_, patchD, dMin, dSpan );
        reduce(nPaddles_, minOp<label>());
    }

    if ( waterDepth_ <= 0.0 )
    {
        waterDepth_ = calcWL(alphaCell, cellGroup, zSpan)[0];
    }

    if ( waveHeight_ < 0.0 )
    {
        FatalError
            << "Check wave height value."
            << exit(FatalError);
    }

    if ( waveTheory_ == "aaa" )    // No wave theory specified
    {
        FatalError
            << "Wave theory not specified."
            << exit(FatalError);
    }
    else if ( waveTheory_ == "Boussinesq" || waveTheory_ == "Grimshaw")
    {
        Info << "\nWave Generation BC" << endl;
        Info << "Wave theory: " << waveTheory_ << endl;
        Info << "H: " << waveHeight_ << endl;
        Info << "h: " << waterDepth_ << endl;
        Info << "Direction: " << waveDir_ << "º" <<  "\n\n" << endl;
    }
    else
    {
        FatalError
            << "Wave theory not supported, use:\n"
            << "Boussinesq, Grimshaw"
            << exit(FatalError);
    }
