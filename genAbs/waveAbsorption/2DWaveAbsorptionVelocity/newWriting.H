    os.writeKeyword("absorptionDir") << absorptionDir_ << 
        token::END_STATEMENT << nl;
    os.writeKeyword("nPaddles") << nPaddles_ << token::END_STATEMENT << nl;

    writeEntry(os, "initialWaterDepths", initialWaterDepths_);
    writeEntry(os, "meanAngles", meanAngles_);
    writeEntry(os, "zSpanL", zSpanL_);

    writeEntryIfDifferent<vector>(os, "uCurrent", vector::zero, uCurrent_);
    writeEntryIfDifferent<label>(os, "nEdgeMin", 0, nEdgeMin_);
    writeEntryIfDifferent<label>(os, "nEdgeMax", 0, nEdgeMax_);

    os.writeKeyword("allCheck") << allCheck_ << token::END_STATEMENT << nl;

    writeEntry(os, "value", *this);
