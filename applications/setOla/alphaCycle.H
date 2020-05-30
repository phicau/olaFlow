forAll(calculatedLevel, cellI)
{
    if(zSup[cellI] <= calculatedLevel[cellI])
    {
        alpha[cellI] = 1;
    }
    else if(zInf[cellI] <= calculatedLevel[cellI])
    {
        alpha[cellI] = 
            (calculatedLevel[cellI] - zInf[cellI])/
            (zSup[cellI] - zInf[cellI]);
    }
    else
    {
        alpha[cellI] = 0;
    }
}