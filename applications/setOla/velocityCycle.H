scalar zExtra = 0.0;

scalarList patchU (xComp.size(), 0.0);
scalarList patchV (xComp.size(), 0.0);
scalarList patchW (xComp.size(), 0.0);

forAll(calculatedLevel, cellIndex)
{
    if (zComp[cellIndex] > waterDepth_ && !extrapolation_)
    {
        zExtra = waterDepth_;
    }
    else
    {
        zExtra = zComp[cellIndex];
    }

    if (zComp[cellIndex] < 0.)
    {
        zExtra = 0.;
    }

    #include "velocityProfile.H"

    vector cellV = 
        vector(patchU[cellIndex], patchV[cellIndex], patchW[cellIndex])
        + uCurrent_*alpha[cellIndex];

    U[cellIndex] = cellV;
}
