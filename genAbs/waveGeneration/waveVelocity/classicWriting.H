    os.writeKeyword("waveType") << waveType_ << token::END_STATEMENT << nl;
    os.writeKeyword("waterDepth") << waterDepth_ << token::END_STATEMENT << nl;
    os.writeKeyword("genAbs") << genAbs_ << token::END_STATEMENT << nl;
    os.writeKeyword("nPaddles") << nPaddles_ << token::END_STATEMENT << nl; 
    os.writeKeyword("allCheck") << allCheck_ << token::END_STATEMENT << nl;
    os.writeKeyword("waveDictName") << waveDictName_ << token::END_STATEMENT << nl;

    #if OFVERSION >= 1712
        os.writeEntryIfDifferent<vector>("uCurrent", vector::zero, uCurrent_);
        os.writeEntryIfDifferent<scalar>("tSmooth", -1.0, tSmooth_);
        os.writeEntryIfDifferent<scalar>("tuningFactor", 1.0, tuningFactor_);
        os.writeEntryIfDifferent<bool>("extrapolation", true, extrapolation_);
    #else
        writeEntryIfDifferent<vector>(os, "uCurrent", vector::zero, uCurrent_);
        writeEntryIfDifferent<scalar>(os, "tSmooth", -1.0, tSmooth_);
        writeEntryIfDifferent<scalar>(os, "tuningFactor", 1.0, tuningFactor_);
        writeEntryIfDifferent<scalar>(os, "extrapolation", true, extrapolation_);
    #endif

    if (genAbs_)
    {
        os.writeKeyword("absDir") << absDir_ << token::END_STATEMENT << nl;
    }

    if ( waveType_ == "irregular" )
    {
        waveHeights_.writeEntry("waveHeights", os);
        wavePeriods_.writeEntry("wavePeriods", os);
        waveLengths_.writeEntry("waveLengths", os);
        wavePhases_.writeEntry("wavePhases", os);
        waveDirs_.writeEntry("waveDirs", os);
        timeLags_.writeEntry("timeLags", os);

        os.writeKeyword("nComp") << nComp_ << token::END_STATEMENT << nl;

        #if OFVERSION >= 1712
            os.writeEntryIfDifferent<bool>("secondOrder", false, secondOrder_);
        #else
            writeEntryIfDifferent<bool>(os, "secondOrder", false, secondOrder_);
        #endif
    }
    else if ( waveType_ == "regular" )
    {
        os.writeKeyword("waveTheory") << 
            waveTheory_ << token::END_STATEMENT << nl;
        os.writeKeyword("waveHeight") << 
            waveHeight_ << token::END_STATEMENT << nl;
        os.writeKeyword("waveDir") << waveDir_ << token::END_STATEMENT << nl;
        os.writeKeyword("timeLag") << timeLag_ << token::END_STATEMENT << nl;

        if ( waveTheory_ == "streamFunction" )
        {
            os.writeKeyword("uMean") << uMean_ << token::END_STATEMENT << nl;
            Bjs_.writeEntry("Bjs", os);
            Ejs_.writeEntry("Ejs", os);
        }
        else
        {
            os.writeKeyword("waveLength") << 
                waveLength_ << token::END_STATEMENT << nl;
            os.writeKeyword("wavePeriod") << 
                wavePeriod_ << token::END_STATEMENT << nl;
            os.writeKeyword("wavePhase") << 
                wavePhase_ << token::END_STATEMENT << nl;

            if ( waveTheory_ == "StokesIII" )
            {
                #if OFVERSION >= 1712
                    os.writeEntryIfDifferent<scalar>("aE_SIII", -1, aE_SIII_);
                    os.writeEntryIfDifferent<scalar>("klE_SIII", -1, klE_SIII_);
                #else
                    writeEntryIfDifferent<scalar>(os, "aE_SIII", -1, aE_SIII_);
                    writeEntryIfDifferent<scalar>(os, "klE_SIII", -1, klE_SIII_);
                #endif
            }
            else if ( waveTheory_ == "StokesV" )
            {
                os.writeKeyword("lambdaStokesV") << 
                    lambdaStokesV_ << token::END_STATEMENT << nl;
            }
            else if ( waveTheory_ == "cnoidal" )
            {
                os.writeKeyword("mCnoidal") << 
                    mCnoidal_ << token::END_STATEMENT << nl;
            }
        }
    }
    else if ( waveType_ == "wavemaker" )
    {
        os.writeKeyword("waveTheory") << 
            waveTheory_ << token::END_STATEMENT << nl;
        timeSeries_.writeEntry("timeSeries", os);
        paddleVelocityU_.writeEntry("paddleVelocity", os);
        paddleVelocityW_.writeEntry("paddleVelocityW", os);
        paddleEta_.writeEntry("paddleEta", os);

        #if OFVERSION >= 1712
            os.writeEntryIfDifferent<word>("waveTheoryOrig", "aaa", waveTheoryOrig_);
        #else
            writeEntryIfDifferent<word>(os, "waveTheoryOrig", "aaa", waveTheoryOrig_);
        #endif
    }
    else if ( waveType_ == "solitary" )
    {
        os.writeKeyword("waveTheory") << 
            waveTheory_ << token::END_STATEMENT << nl;
        os.writeKeyword("waveHeight") << 
            waveHeight_ << token::END_STATEMENT << nl;
        os.writeKeyword("waveLength") << 
            waveLength_ << token::END_STATEMENT << nl;
        os.writeKeyword("wavePeriod") << 
            wavePeriod_ << token::END_STATEMENT << nl;
        os.writeKeyword("waveDir") << waveDir_ << token::END_STATEMENT << nl;

        #if OFVERSION >= 1712
            os.writeEntryIfDifferent<label>("nSolitaryWaves", 1, nSolitaryWaves_);
        #else
            writeEntryIfDifferent<label>(os, "nSolitaryWaves", 1, nSolitaryWaves_);
        #endif 
    }

    writeEntry("value", os);