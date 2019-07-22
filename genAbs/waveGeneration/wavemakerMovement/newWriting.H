    os.writeKeyword("allCheck") << allCheck_ << token::END_STATEMENT << nl;
    os.writeKeyword("wavemakerDictName") << wavemakerDictName_
        << token::END_STATEMENT << nl;
    os.writeKeyword("wavemakerType") << wavemakerType_
        << token::END_STATEMENT << nl;

    os.writeKeyword("nPaddles") << nPaddles_ << token::END_STATEMENT << nl;
    os.writeKeyword("initialWaterDepth") << initialWaterDepth_
        << token::END_STATEMENT << nl;
    os.writeKeyword("meanAngle") << meanAngle_ << token::END_STATEMENT << nl;

    writeEntry(os, "timeSeries", timeSeries_);

    if ( tSmooth_ != -1.0 )
    {
        os.writeKeyword("tSmooth") << tSmooth_ << token::END_STATEMENT << nl;
    }

    if ( tuningFactor_ != 1.0 )
    {
        os.writeKeyword("tuningFactor") << tuningFactor_
            << token::END_STATEMENT << nl;
    }

    if ( wavemakerType_ == "Piston" )
    {
        // paddlePosition_.writeEntry("paddlePosition", os);
        os.writeKeyword("paddlePosition") << paddlePosition_
            << token::END_STATEMENT << nl; // May not work for binary encoding

        if ( genAbs_ )
        {
            os.writeKeyword("genAbs") << genAbs_ << token::END_STATEMENT << nl;
            writeEntry(os, "cumAbsCorrection", cumAbsCorrection_);
            os.writeKeyword("paddleEta") << paddleEta_
                << token::END_STATEMENT << nl; // May not work for binary encoding
        }
    }
    else if ( wavemakerType_ == "Flap" )
    {
        os.writeKeyword("paddleTilt") << paddleTilt_
            << token::END_STATEMENT << nl; // May not work for binary encoding
        os.writeKeyword("hingeHeight") << hingeHeight_
            << token::END_STATEMENT << nl;
        os.writeKeyword("hingeLocation") << hingeLocation_
            << token::END_STATEMENT << nl;

        writeEntry(os, "tiltOld", tiltOld_);
    }
    else if ( wavemakerType_ == "Mixed" )
    {
        os.writeKeyword("paddlePosition") << paddlePosition_
            << token::END_STATEMENT << nl; // May not work for binary encoding
        os.writeKeyword("paddleTilt") << paddleTilt_
            << token::END_STATEMENT << nl; // May not work for binary encoding

        os.writeKeyword("hingeHeight") << hingeHeight_
            << token::END_STATEMENT << nl;
        writeEntry(os, "tiltOld", tiltOld_);
    }

    if ( wavemakerType_ == "Piston" || wavemakerType_ == "Mixed" )
    {
        if ( maxStroke_ != 999.0 ) // DPS active
        {
            writeEntry(os, "DPS", DPS_);
            os.writeKeyword("maxStroke") << maxStroke_
                << token::END_STATEMENT << nl;
            os.writeKeyword("DPST") << DPST_ << token::END_STATEMENT << nl;
            writeEntry(os, "cumDPSCorrection", cumDPSCorrection_);
            writeEntry(os, "DPSsign", DPSsign_);
            writeEntry(os, "DPStIni", DPStIni_);
            writeEntry(os, "instDPSCorrection", instDPSCorrection_);
        }
    }