    if ( waveType_ == "regular" )
    {
        if ( waveTheory_ == "StokesI" )
        {
            faseTot = waveKx*xGroup[cellGroup[cellIndex]-1] 
                + waveKy*yGroup[cellGroup[cellIndex]-1] 
                - waveOmega*currTime + wavePhase_;
            #include "profileStokesI.H"
        }
        else if ( waveTheory_ == "StokesII" )
        {
            faseTot = waveKx*xGroup[cellGroup[cellIndex]-1] 
                + waveKy*yGroup[cellGroup[cellIndex]-1] 
                - waveOmega*currTime + wavePhase_;
            #include "profileStokesII.H"
        }
        else  if ( waveTheory_ == "StokesV" )
        {
            #include "profileStokesV.H"
        }
        else if ( waveTheory_ == "cnoidal" )
        {
            #include "profileCnoidal.H"
        }
        else if ( waveTheory_ == "streamFunction" )
        {
            faseTot = waveKx*xGroup[cellGroup[cellIndex]-1] 
                + waveKy*yGroup[cellGroup[cellIndex]-1] 
                - waveOmega*currTime + wavePhase_;
            #include "profileStreamFunction.H"
        }
    }
    else if ( waveType_ == "solitary" )
    {
        #include "profileBoussinesq.H"
    }
    else if ( waveType_ == "irregular" )
    {
        #include "profileIrregular.H"
    }
    else if ( waveType_ == "wavemaker" )
    {
        if ( noEta && genAbs_ == false )
        {
            #include "profileNoEta.H"
        }
        else
        {
            #include "profileEta.H"
        }
    }
    else if ( waveType_ == "current" )
    {
        #include "profileCurrent.H"
    }
