Deta1 = 
    8.0*auxiliarSolit*waterDepth_
    *exp(2.0*auxiliarSolit*Xa)
    *(1.0-exp(2.0*auxiliarSolit*Xa))
    /pow(1.0+exp(2.0*auxiliarSolit*Xa),3);

Deta2 = 
    16.0*sqr(auxiliarSolit)*waterDepth_
    *exp(2.0*auxiliarSolit*Xa)
    *(exp(4.0*auxiliarSolit*Xa)-
    4.0*exp(2.0*auxiliarSolit*Xa)+1.0)
    /pow(1.0+exp(2.0*auxiliarSolit*Xa),4);

Deta3 = 
    -32.0*pow(auxiliarSolit,3)*waterDepth_
    *exp(2.0*auxiliarSolit*Xa)
    *(exp(6.0*auxiliarSolit*Xa)-
    11.0*exp(4.0*auxiliarSolit*Xa)+
    11.0*exp(2.0*auxiliarSolit*Xa)-1.0)
    /pow(1.0+exp(2.0*auxiliarSolit*Xa),5);
