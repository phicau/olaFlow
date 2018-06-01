    // Extracting values from dict
    waveTheory_ = (waveDict.lookupOrDefault<word>("waveTheory", "aaa"));
    timeSeries_ = 
        (waveDict.lookupOrDefault("timeSeries", List<scalar> (1, -1.0) ));
    paddlePosition_ = 
        (waveDict.lookupOrDefault("paddlePosition", List<scalar> (1, -1.0)));
    paddleVelocityU_ = 
        (waveDict.lookupOrDefault("paddleVelocity", List<scalar> (1, -1.0)));
    paddleVelocityW_ = 
        (waveDict.lookupOrDefault("paddleVelocityW", List<scalar> (1, -1.0)));
    paddleEta_ = 
        (waveDict.lookupOrDefault("paddleEta", List<scalar> (1, -1.0)));

    genAbs_ = (waveDict.lookupOrDefault<bool>("genAbs", false ));
    nPaddles_ = (waveDict.lookupOrDefault<label>("nPaddles", 1));
    absDir_ = (waveDict.lookupOrDefault<scalar>("absDir", 0.0));

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

    if ( waveTheory_ == "aaa" )    // No wave theory specified
    {
        FatalError
            << "Wave theory not specified for wavemaker, use:\n"
            << "tx, tv, tvw, txeta, tveta, tvetaw."
            << exit(FatalError);
    }

    if( waveTheory_ == "tx" )
    {
        // Check number of elements
        if( timeSeries_.size() != paddlePosition_.size() || 
            timeSeries_.size() == 1 )
        {
            FatalError
                << "Check number of components of all time series"
                << " or waveTheory."
                << exit(FatalError);
        }

        // Calculate v
        paddleVelocityU_ = 0.0*timeSeries_;

        for( auxiliar=0; auxiliar<timeSeries_.size()-1; auxiliar++ )
        {
            paddleVelocityU_[auxiliar] = 
                (paddlePosition_[auxiliar+1] - paddlePosition_[auxiliar]) 
                /(timeSeries_[auxiliar+1] - timeSeries_[auxiliar]);
        }

        // Calculate eta
        paddleEta_ = paddlePosition_ - paddlePosition_[0];
        
        // Write & update results
        waveTheoryOrig_ = waveTheory_;
        waveDict.set("waveTheoryOrig", waveTheory_);
        waveTheory_ = "tveta";
        waveDict.set("waveTheory", waveTheory_);
        waveDict.set("paddleVelocity", paddleVelocityU_);
        waveDict.set("paddleEta", paddleEta_);
        waveDict.Foam::regIOobject::write();

    }
    else if( waveTheory_ == "tv" )
    {
        // Check number of elements
        if( timeSeries_.size() != paddleVelocityU_.size() || 
            timeSeries_.size() == 1 )
        {
            FatalError
                << "Check number of components of all time series"
                << " or waveTheory."
                << exit(FatalError);
        }

        // Calculate eta
        paddleEta_ = 0.0*timeSeries_;

        for( auxiliar=1; auxiliar<timeSeries_.size(); auxiliar++ )
        {
            paddleEta_[auxiliar] = 
                paddleEta_[auxiliar-1] 
                + paddleVelocityU_[auxiliar-1]
                *(timeSeries_[auxiliar] - timeSeries_[auxiliar-1]);
        }

        // Write & update results
        waveTheoryOrig_ = waveTheory_;
        waveDict.set("waveTheoryOrig", waveTheory_);
        waveTheory_ = "tveta";
        waveDict.set("waveTheory", waveTheory_);
        waveDict.set("paddleEta", paddleEta_);
        waveDict.Foam::regIOobject::write();

    }
    else if( waveTheory_ == "tvw" )
    {
        // Check number of elements
        if( timeSeries_.size() != paddleVelocityU_.size() ||
            timeSeries_.size() != paddleVelocityW_.size() ||
            timeSeries_.size() == 1 )
        {
            FatalError
                << "Check number of components of all time series"
                << " or waveTheory."
                << exit(FatalError);
        }

        // Calculate eta
        paddleEta_ = 0.0*timeSeries_;

        for( auxiliar=1; auxiliar<timeSeries_.size(); auxiliar++ )
        {
            paddleEta_[auxiliar] = 
                paddleEta_[auxiliar-1] 
                + paddleVelocityU_[auxiliar-1]
                *(timeSeries_[auxiliar] - timeSeries_[auxiliar-1]);
        }

        // Write & update results
        waveTheoryOrig_ = waveTheory_;
        waveDict.set("waveTheoryOrig", waveTheory_);
        waveTheory_ = "tvetaw";
        waveDict.set("waveTheory", waveTheory_);
        waveDict.set("paddleEta", paddleEta_);
        waveDict.Foam::regIOobject::write();

    }
    else if( waveTheory_ == "txeta" )
    {
        // Check number of elements
        if( timeSeries_.size() != paddlePosition_.size() || 
            timeSeries_.size() != paddleEta_.size() || 
            timeSeries_.size() == 1 )
        {
            FatalError
                << "Check number of components of all time series"
                << " or waveTheory."
                << exit(FatalError);
        }

        // Calculate v
        paddleVelocityU_ = 0.0*timeSeries_;

        for( auxiliar=0; auxiliar<timeSeries_.size()-1; auxiliar++ )
        {
            paddleVelocityU_[auxiliar] = 
                (paddlePosition_[auxiliar+1] 
                - paddlePosition_[auxiliar])
                /(timeSeries_[auxiliar+1] - timeSeries_[auxiliar]);
        }

        // Write & update results
        waveTheoryOrig_ = waveTheory_;
        waveDict.set("waveTheoryOrig", waveTheory_);
        waveTheory_ = "tveta";
        waveDict.set("waveTheory", waveTheory_);
        waveDict.set("paddleVelocity", paddleVelocityU_);
        waveDict.Foam::regIOobject::write();

    }
    else if( waveTheory_ == "tveta" )
    {
        // Check number of elements
        if( timeSeries_.size() != paddleVelocityU_.size() || 
            timeSeries_.size() != paddleEta_.size() || 
            timeSeries_.size() == 1 )
        {
            FatalError
                << "Check number of components of all time series"
                << " or waveTheory."
                << exit(FatalError);
        }
        // All set, just interpolate
        waveTheoryOrig_ = 
            (waveDict.lookupOrDefault<word>("waveTheoryOrig", "aaa"));
    }
    else if( waveTheory_ == "tvetaw" )
    {
        // Check number of elements
        if( timeSeries_.size() != paddleVelocityU_.size() ||
            timeSeries_.size() != paddleVelocityW_.size() ||
            timeSeries_.size() != paddleEta_.size() ||
            timeSeries_.size() == 1 )
        {
            FatalError
                << "Check number of components of all time series"
                << " or waveTheory."
                << exit(FatalError);
        }
        // All set, just interpolate
        waveTheoryOrig_ = 
            (waveDict.lookupOrDefault<word>("waveTheoryOrig", "aaa"));
    }
    else
    {
        FatalError
            << "Wave theory not supported for wavemaker, use:\n"
            << "tx, tv, tvw, txeta, tveta, tvetaw."
            << exit(FatalError);
    }
