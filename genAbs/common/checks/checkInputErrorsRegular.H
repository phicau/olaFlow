    // Extracting values from dict
    waveTheory_ = (waveDict.lookupOrDefault<word>("waveTheory", "aaa"));
    waveHeight_ = (waveDict.lookupOrDefault<scalar>("waveHeight", -1));
    wavePeriod_ = (waveDict.lookupOrDefault<scalar>("wavePeriod", -1));
    waveDir_ = (waveDict.lookupOrDefault<scalar>("waveDir", 0));
    genAbs_ = (waveDict.lookupOrDefault<bool>("genAbs", false ));
    wavePhase_ = 
        (waveDict.lookupOrDefault<scalar>("wavePhase", 3.0*PI()/2.0 ));
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
            << "Wave theory not specified, use:\n"
            << "StokesI, StokesII, StokesIII, StokesV, cnoidal, streamFunction"
            << exit(FatalError);
    }
    else if ( waveTheory_ == "StokesI" || waveTheory_ == "StokesII" )
    {
        if ( wavePeriod_ <= 0.0 )
        {
            FatalError
                << "Check wave period value."
                << exit(FatalError);
        }

        word generation = "a"; // Auxiliar

        waveLength_ = StokesIFun::waveLength (waterDepth_, wavePeriod_);

        waveK = 2.0*PI()/waveLength_;

        if ( waveK*waterDepth_ > PI() )
        {
            generation = "Deep";
        }
        else if ( waveK*waterDepth_ < PI()/10.0 )
        {
            generation = "Shallow";
        }
        else
        {
            generation = "Intermediate";
        }

        Info << "\nWave Generation BC" << endl;
        Info << "Wave theory: " << waveTheory_ << endl;
        Info << "H: " << waveHeight_ << endl;
        Info << "T: " << wavePeriod_ << endl;
        Info << "h: " << waterDepth_ << endl;
        Info << "L: " << waveLength_ << endl;
        Info << "Direction: " << waveDir_ << "º" << endl;
        Info << "Generation in: " << generation << " waters." << endl;
        Info << "Relative depth (kh): " << 
            waveK*waterDepth_ << "\n\n" << endl;

        if (waveTheory_ == "StokesII")
        {
            timeLag_ = 
                StokesIIFun::timeLag
                (
                    waveHeight_, 
                    waterDepth_, 
                    waveK, 
                    0, 0, 0, 
                    wavePeriod_, 
                    wavePhase_
                );
        }
    }
    else if ( waveTheory_ == "StokesIII" )
    {
        if ( wavePeriod_ <= 0.0 )
        {
            FatalError
                << "Check wave period value."
                << exit(FatalError);
        }

        word generation = "a"; // Auxiliar

        bool iniOK =
            StokesIIIFun::initialise
            (
                waveHeight_,
                waterDepth_,
                wavePeriod_,
                &aE_SIII_,
                &klE_SIII_
            );

        if( !iniOK )
        {
            FatalError
                << "StokesIII initialisation did not converge,"
                << "check your wave conditions."
                << exit(FatalError);
        }

        waveLength_ = 
            StokesIIIFun::waveLength (wavePeriod_, aE_SIII_, klE_SIII_);


        waveK = 2.0*PI()/waveLength_;

        if ( waveK*waterDepth_ > PI() )
        {
            generation = "Deep";
        }
        else if ( waveK*waterDepth_ < PI()/10.0 )
        {
            generation = "Shallow";
        }
        else
        {
            generation = "Intermediate";
        }

        Info << "\nWave Generation BC" << endl;
        Info << "Wave theory: " << waveTheory_ << endl;
        Info << "H: " << waveHeight_ << endl;
        Info << "T: " << wavePeriod_ << endl;
        Info << "h: " << waterDepth_ << endl;
        Info << "L: " << waveLength_ << endl;
        Info << "Direction: " << waveDir_ << "º" << endl;
        Info << "Generation in: " << generation << " waters." << endl;
        Info << "Relative depth (kh): " << 
            waveK*waterDepth_ << "\n\n" << endl;
    }
    else if ( waveTheory_ == "StokesV" )
    {
        if ( wavePeriod_ <= 0.0 )
        {
            FatalError
                << "Check wave period value."
                << exit(FatalError);
        }

        word generation = "a"; // Auxiliar

        scalar f1;
        scalar f2;

        stokesVFun::StokesVNR
        (
            waveHeight_, 
            waterDepth_, 
            wavePeriod_, 
            &waveK, 
            &lambdaStokesV_, 
            &f1, 
            &f2 
        );

        Info << "f1 residual " << f1 << endl;
        Info << "f2 residual " << f2 << endl;

        if ( f1 > 0.001 || f2 > 0.001 )
        {
	        FatalError
	            << "No convergence for Stokes V wave theory.\n"
                << exit(FatalError);
        }

        waveLength_ = 2.0*PI()/waveK;

        timeLag_ = 
            stokesVFun::timeLag
            (
                waterDepth_, 
                waveK, 
                0, 
                lambdaStokesV_, 
                wavePeriod_, 
                0, 0, 
                wavePhase_
            );

        if ( waveK*waterDepth_ > PI() )
        {
            generation = "Deep";
        }
        else if ( waveK*waterDepth_ < PI()/10.0 )
        {
            generation = "Shallow";
        }
        else
        {
            generation = "Intermediate";
        }

        Info << "\nWave Generation BC" << endl;
        Info << "Wave theory: " << waveTheory_ << endl;
        Info << "H: " << waveHeight_ << endl;
        Info << "T: " << wavePeriod_ << endl;
        Info << "h: " << waterDepth_ << endl;
        Info << "L: " << waveLength_ << endl;
        Info << "Lambda: " << lambdaStokesV_ << endl;
        Info << "Direction: " << waveDir_ << "º"  << endl;
        Info << "Generation in: " << generation << " waters." << endl;
        Info << "Relative depth (kh): " << 
            waveK*waterDepth_ << "\n\n" << endl;
    }
    else if ( waveTheory_ == "cnoidal" )
    {
        if ( wavePeriod_ <= 0.0 )
        {
            FatalError
                << "Check wave period value."
                << exit(FatalError);
        }

        cnoidalFun::calculations
        (
            waveHeight_, 
            waterDepth_, 
            wavePeriod_, 
            &mCnoidal_, 
            &waveLength_
        );

        timeLag_ = 
            cnoidalFun::timeLag
            (
                waveHeight_, 
                mCnoidal_, 
                waveLength_, 
                0, 
                wavePeriod_, 
                0, 0, 
                wavePhase_
            );

        Info << "\nWave Generation BC" << endl;
        Info << "Wave theory: " << waveTheory_ << endl;
        Info << "H: " << waveHeight_ << endl;
        Info << "T: " << wavePeriod_ << endl;
        Info << "h: " << waterDepth_ << endl;
        Info << "L: " << waveLength_ << endl;
        Info << "m parameter: " << mCnoidal_ << endl;
        Info << "Direction: " << waveDir_ << "º" <<  "\n\n" << endl;
    }
    else if ( waveTheory_ == "streamFunction" )
    {
        waveLength_ = (waveDict.lookupOrDefault<scalar>("waveLength", -1));

        if ( wavePeriod_ <= 0.0 )
        {
            FatalError
                << "Check wave period value."
                << exit(FatalError);
        }

        if ( waveLength_ <= 0.0 )
        {
            FatalError
                << "Check waveLength value."
                << exit(FatalError);
        }

        uMean_ = (waveDict.lookupOrDefault<scalar>("uMean", -1));

        if ( uMean_ <= 0.0 )
        {
            FatalError
                << "Check uMean value."
                << exit(FatalError);
        }

        Bjs_ = (waveDict.lookupOrDefault("Bj", List<scalar> (1, -1.0) ));
        Ejs_ = (waveDict.lookupOrDefault("Ej", List<scalar> (1, -1.0) ));

        if ( Bjs_.size() != Ejs_.size() || Bjs_.size() == 1 )
        {
            FatalError
                << "Check the number of components for Ej and Bj."
                << exit(FatalError);
        }

        Info << "\nWave Generation BC" << endl;
        Info << "Wave theory: " << waveTheory_ << endl;
        Info << "H: " << waveHeight_ << endl;
        Info << "T: " << wavePeriod_ << endl;
        Info << "h: " << waterDepth_ << endl;
        Info << "L: " << waveLength_ << endl;
        Info << "Bj: " << Bjs_ << endl;
        Info << "Ej: " << Ejs_ << endl;
        Info << "Direction: " << waveDir_ << "º" << endl;

    }
    else
    {
        FatalError
            << "Wave theory not supported, use:\n"
            << "StokesI, StokesII, StokesIII, StokesV, cnoidal, streamFunction"
            << exit(FatalError);
    }
