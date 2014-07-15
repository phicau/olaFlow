    // Extracting values from dict
    waveHeights_ = 
        (IHWavesDict.lookupOrDefault("waveHeights", List<scalar> (1, -1.0) ));
    wavePeriods_ = 
        (IHWavesDict.lookupOrDefault("wavePeriods", List<scalar> (1, -1.0) ));
    waveDirs_ = 
        (IHWavesDict.lookupOrDefault("waveDirs", List<scalar> (1, -1.0) ));
    wavePhases_ = 
        (IHWavesDict.lookupOrDefault("wavePhases", List<scalar> (1, -1.0) ));

    genAbs_ = (IHWavesDict.lookupOrDefault<bool>("genAbs", false ));
    nPaddles_ = (IHWavesDict.lookupOrDefault<label>("nPaddles", 1));
    absDir_ = (IHWavesDict.lookupOrDefault<scalar>("absDir", 0));

    secondOrder_ = (IHWavesDict.lookupOrDefault<bool>("secondOrder", false ));

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

    if ( min(waveHeights_) < 0.0 )
    {
        FatalError
            << "Check wave height values."
            << exit(FatalError);
    }

    if ( min(wavePeriods_) <= 0.0 )
    {
        FatalError
            << "Check wave period values."
            << exit(FatalError);
    }

    nComp_ = waveHeights_.size();

    if ( nComp_ != wavePeriods_.size() || nComp_ != waveDirs_.size() || 
            nComp_ != wavePhases_.size() || nComp_ == 1 )
    {
        FatalError
            << "Check the number of components for "
            << "waveHeights, wavePeriods, wavePhases and waveDirs."
            << exit(FatalError);
    }

    if ( secondOrder_ )
    {
        if ( max(waveDirs_) != min(waveDirs_) || max(waveDirs_) != 0.0 )
        {
            FatalError
                << "Check waveDirs must be the same and equal to 0º"
                << " to apply second order formulation."
                << exit(FatalError); 
        }
    }

    waveLengths_ = 0.0*waveHeights_;

    forAll(waveLengths_, compInd)
    {
        waveLengths_[compInd] = 
            StokesIFun::waveLength (waterDepth_, wavePeriods_[compInd]);
    }
