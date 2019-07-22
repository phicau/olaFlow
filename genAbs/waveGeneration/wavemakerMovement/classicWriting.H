    os.writeKeyword("allCheck") << allCheck_ << token::END_STATEMENT << nl;
    os.writeKeyword("wavemakerDictName") << wavemakerDictName_
        << token::END_STATEMENT << nl;
    os.writeKeyword("wavemakerType") << wavemakerType_
        << token::END_STATEMENT << nl;

    os.writeKeyword("nPaddles") << nPaddles_ << token::END_STATEMENT << nl;
    os.writeKeyword("initialWaterDepth") << initialWaterDepth_
        << token::END_STATEMENT << nl;
    os.writeKeyword("meanAngle") << meanAngle_ << token::END_STATEMENT << nl;

    timeSeries_.writeEntry("timeSeries", os);

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
            cumAbsCorrection_.writeEntry("cumAbsCorrection", os);
            // paddleEta_.writeEntry("paddleEta", os);
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

        tiltOld_.writeEntry("tiltOld", os);
    }
    else if ( wavemakerType_ == "Mixed" )
    {
        os.writeKeyword("paddlePosition") << paddlePosition_
            << token::END_STATEMENT << nl; // May not work for binary encoding
        os.writeKeyword("paddleTilt") << paddleTilt_
            << token::END_STATEMENT << nl; // May not work for binary encoding

        os.writeKeyword("hingeHeight") << hingeHeight_
            << token::END_STATEMENT << nl;
        tiltOld_.writeEntry("tiltOld", os);
    }

    if ( wavemakerType_ == "Piston" || wavemakerType_ == "Mixed" )
    {
        if ( maxStroke_ != 999.0 ) // DPS active
        {
            DPS_.writeEntry("DPS", os);
            os.writeKeyword("maxStroke") << maxStroke_
                << token::END_STATEMENT << nl;
            os.writeKeyword("DPST") << DPST_ << token::END_STATEMENT << nl;
            cumDPSCorrection_.writeEntry("cumDPSCorrection", os);
            DPSsign_.writeEntry("DPSsign", os);
            DPStIni_.writeEntry("DPStIni", os);
            instDPSCorrection_.writeEntry("instDPSCorrection", os);
        }
    }