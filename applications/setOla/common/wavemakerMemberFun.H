//- Return Pi
virtual scalar PI() const
{
    #if OFVERSION >= 200
        const scalar PI = constant::mathematical::pi;
    #else
        const scalar PI = mathematicalConstant::pi;
    #endif
    return PI;
}

//- Return alphaName
virtual word alphaName() const
{
    #if OFVERSION >= 230
        const word an = "alpha.water";
    #else
        const word an = "alpha1";
    #endif
    return an;
}

//- Degree to Rad conversion
virtual scalar deg2rad(scalar num) const
{
    return num*PI()/180.0;
}

virtual scalarList deg2rad(scalarList nums) const
{
    return nums*PI()/180.0;
}

//- Return prevalent direction - X or Y to sort groups of columns
virtual scalarField patchDirection (
    vector span, 
    scalar* dMin, 
    scalar* dSpan, 
    vectorField patchCf ) const
{
    // Check prevalent direction
    if (span[0] >= span[1]) // More perpendicular to Y plane
    {
        *dMin = gMin(this->patch().localPoints().component(0));
        *dSpan = gMax(this->patch().localPoints().component(0)) - *dMin;
        return patchCf.component(0); // patchX
    }
    else // More perpendicular to X plane
    {
        *dMin = gMin(this->patch().localPoints().component(1));
        *dSpan = gMax(this->patch().localPoints().component(1)) - *dMin;
        return patchCf.component(1);  // patchY
    }
}

//- Paddle index of a specific point
virtual label paddleIndex(scalar pointD, scalar yMin, scalar ySpan) const
{
    scalarList dBreakPoints = scalarList(nPaddles_+1, yMin); 

    for (label i=0; i<nPaddles_; i++)
    {
        dBreakPoints[i+1] = yMin + ySpan/nPaddles_*(i+1); // Breakpoints
    }

    // Number of the group which the cell belongs to
    for (label i=0; i<nPaddles_; i++)
    {
        if ( (pointD >= dBreakPoints[i]) 
            && (pointD <= dBreakPoints[i+1]) )
        {
            return i;
        }
    }      
    return nPaddles_-1;
}

//- Simple linear interpolation
virtual scalar linInterp (
    scalar x1, 
    scalar x2, 
    scalar y1, 
    scalar y2, 
    scalar xInt) const
{
    scalar yInt = y1 + (y2-y1)/(x2-x1)*(xInt-x1);
    return yInt;
}

//- Simple ramp cosine interpolation: f'= 0 at ends
virtual scalar DPSramp (scalar S, scalar t, scalar T) const
{
    scalar s = S * (1.0-cos( PI()*t/T ))/2.0;
    return s;
}

//- Simple linear interpolation for paddles
virtual scalar paddleLinInterp (
    scalarList space,
    scalarList disp,
    scalar spInt) const
{
    label length = space.size();

    if ( spInt <= space[0] )
    {
        return disp[0];
    }
    else if ( spInt >= space[length-1] )
    {
        return disp[length-1];
    }

    label indexF = 0;

    forAll( space, contador )
    {
        if ( space[contador] >= spInt )
        {
            indexF = contador;
            break;
        }
    }
    
    scalar deltaSpace = space[indexF] - space[indexF-1];
    scalar deltaDisp = disp[indexF] - disp[indexF-1];

    scalar rampRelWidth = min(abs(deltaDisp/deltaSpace), 1.0);

    scalar midPoint = (space[indexF] + space[indexF-1])/2.0;
    scalar rampWidth = (space[indexF] - space[indexF-1])*rampRelWidth;

    if ( spInt <= midPoint-rampWidth/2.0 )
    {
        return disp[indexF-1];
    }
    else if ( spInt >= midPoint+rampWidth/2.0 )
    {
        return disp[indexF];
    }
    else
    {
        return linInterp
        (
            midPoint-rampWidth/2.0,
            midPoint+rampWidth/2.0,
            disp[indexF-1],
            disp[indexF],
            spInt
        );
    }
}

//- Simple ramp cosine interpolation for paddles
virtual scalar paddleCosInterp (
    scalarList space, 
    scalarList disp, 
    scalar spInt) const
{
    label length = space.size();

    if ( spInt <= space[0] )
    {
        return disp[0];
    }
    else if ( spInt >= space[length-1] )
    {
        return disp[length-1];
    }

    label indexF = 0;

    forAll( space, contador )
    {
        if ( space[contador] >= spInt )
        {
            indexF = contador;
            break;
        }
    }

    scalar deltaSpace = space[indexF] - space[indexF-1];
    scalar deltaDisp = disp[indexF] - disp[indexF-1];

    scalar rampRelWidth = min(fabs(deltaDisp/deltaSpace), 1.0);

    scalar midPoint = (space[indexF] + space[indexF-1])/2.0;
    scalar rampWidth = (space[indexF] - space[indexF-1])*rampRelWidth;

    if ( spInt <= midPoint-rampWidth/2.0 )
    {
        return disp[indexF-1];
    }
    else if ( spInt >= midPoint+rampWidth/2.0 )
    {
        return disp[indexF];
    }
    else
    {
        scalar Saux = disp[indexF] - disp[indexF-1];
        scalar taux = spInt - (midPoint-rampWidth/2.0);
        scalar Taux = rampWidth;

        return disp[indexF-1] + DPSramp(Saux, taux, Taux);
    }
}

//- Simple ramp cosine interpolation for paddles tilt
virtual scalar paddleTiltCosInterp (
    scalarList space,
    scalarList tilt,
    scalarList disp,
    scalar spInt) const
{
    label length = space.size();

    if ( spInt <= space[0] )
    {
        return tilt[0];
    }
    else if ( spInt >= space[length-1] )
    {
        return tilt[length-1];
    }

    label indexF = 0;

    forAll( space, contador )
    {
        if ( space[contador] >= spInt )
        {
            indexF = contador;
            break;
        }
    }

    scalar deltaSpace = space[indexF] - space[indexF-1];
    scalar deltaDisp = disp[indexF] - disp[indexF-1];

    scalar rampRelWidth = min(fabs(deltaDisp/deltaSpace), 1.0);

    scalar midPoint = (space[indexF] + space[indexF-1])/2.0;
    scalar rampWidth = (space[indexF] - space[indexF-1])*rampRelWidth;

    if ( spInt <= midPoint-rampWidth/2.0 )
    {
        return tilt[indexF-1];
    }
    else if ( spInt >= midPoint+rampWidth/2.0 )
    {
        return tilt[indexF];
    }
    else
    {
        scalar Saux = tilt[indexF] - tilt[indexF-1];
        scalar taux = spInt - (midPoint-rampWidth/2.0);
        scalar Taux = rampWidth;

        return tilt[indexF-1] + DPSramp(Saux, taux, Taux);
    }
}

//- Calculate water levels on each paddle (same zSpan)
virtual scalarList calcWL ( 
    scalarField alphaCell, 
    scalarField patchD, 
    scalarField cellSurface, 
    scalar yMin, 
    scalar ySpan, 
    scalar zSpan ) const
{
    // Grouping part
    scalarList dBreakPoints = scalarList(nPaddles_+1, yMin); 

    for (label i=0; i<nPaddles_; i++)
    {
        dBreakPoints[i+1] = yMin + ySpan/nPaddles_*(i+1); // Breakpoints
    }

    // Number of the group which the cell belongs to
    labelList cellGroup = Foam::labelList(alphaCell.size(), 0); 

    forAll(cellGroup, patchCells) 
    {
        for (label i=0; i<nPaddles_; i++)
        {
            if ( (patchD[patchCells]>=dBreakPoints[i]) 
                && (patchD[patchCells]<dBreakPoints[i+1]) )
            {
                cellGroup[patchCells] = i+1; // Group of each face
                continue;
            }
        }      
    }

    // Free surface elevation measurement
    scalarList groupTotalArea (nPaddles_,0.0);
    scalarList groupWaterArea (nPaddles_,0.0);
    scalarList heights (nPaddles_,0.0);

    forAll( cellGroup, index ) 
    {
        groupTotalArea[cellGroup[index]-1] += cellSurface[index];
        groupWaterArea[cellGroup[index]-1] += 
            cellSurface[index] * alphaCell[index];
    } 

    for (label i=0; i<=gMax(cellGroup)-1; i++)
    {
        reduce(groupTotalArea[i], sumOp<scalar>());
        reduce(groupWaterArea[i], sumOp<scalar>());

        heights[i] = groupWaterArea[i]/groupTotalArea[i]*zSpan;
    }

    return heights;
}