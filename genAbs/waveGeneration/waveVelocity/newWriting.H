    os.writeKeyword("waveType") << waveType_ << token::END_STATEMENT << nl;
    os.writeKeyword("waterDepth") << waterDepth_ << token::END_STATEMENT << nl;
    os.writeKeyword("genAbs") << genAbs_ << token::END_STATEMENT << nl;
    os.writeKeyword("nPaddles") << nPaddles_ << token::END_STATEMENT << nl; 
    os.writeKeyword("allCheck") << allCheck_ << token::END_STATEMENT << nl;
    os.writeKeyword("waveDictName") << waveDictName_ << token::END_STATEMENT << nl;

    writeEntryIfDifferent<vector>(os, "uCurrent", vector::zero, uCurrent_);
    writeEntryIfDifferent<scalar>(os, "tSmooth", -1.0, tSmooth_);
    writeEntryIfDifferent<scalar>(os, "tuningFactor", 1.0, tuningFactor_);
    writeEntryIfDifferent<scalar>(os, "extrapolation", true, extrapolation_);

    if (genAbs_)
    {
        os.writeKeyword("absDir") << absDir_ << token::END_STATEMENT << nl;
    }

    if ( waveType_ == "irregular" )
    {
        writeEntry(os, "waveHeights", waveHeights_);
        writeEntry(os, "wavePeriods", wavePeriods_);
        writeEntry(os, "waveLengths", waveLengths_);
        writeEntry(os, "wavePhases", wavePhases_);
        writeEntry(os, "waveDirs", waveDirs_);
        writeEntry(os, "timeLags", timeLags_);

        os.writeKeyword("nComp") << nComp_ << token::END_STATEMENT << nl;

        writeEntryIfDifferent<bool>(os, "secondOrder", false, secondOrder_);
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
            writeEntry(os, "Bjs", Bjs_);
            writeEntry(os, "Ejs", Ejs_);
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
                writeEntryIfDifferent<scalar>(os, "aE_SIII", -1, aE_SIII_);
                writeEntryIfDifferent<scalar>(os, "klE_SIII", -1, klE_SIII_);
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
        writeEntry(os, "timeSeries", timeSeries_);
        writeEntry(os, "paddleVelocity", paddleVelocityU_);
        writeEntry(os, "paddleVelocityW", paddleVelocityW_);
        writeEntry(os, "paddleEta", paddleEta_);

        writeEntryIfDifferent<word>(os, "waveTheoryOrig", "aaa", waveTheoryOrig_);
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

        writeEntryIfDifferent<label>(os, "nSolitaryWaves", 1, nSolitaryWaves_);
    }

    writeEntry(os, "value", *this);