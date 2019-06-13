// TODO Integrate with memberFun.H

        //- Return Pi
        scalar PI()
        {
            #if OFVERSION >= 200
                const scalar PI = constant::mathematical::pi;
            #else
                const scalar PI = mathematicalConstant::pi;
            #endif
            return PI;
        }

        //- Return alphaName
        word alphaName()
        {
            #if OFVERSION >= 230
                const word an = "alpha.water";
            #else
                const word an = "alpha1";
            #endif
            return an;
        }

        //- Return prevalent direction - X or Y to sort groups of columns
        scalarField patchDirection ( vector span, scalar* dMin, scalar* dSpan )
        {
            const pointField& points = patch().patch().localPoints();
            // Check prevalent direction
            if (span[0] >= span[1]) // More perpendicular to Y plane
            {
                *dMin = gMin(points.component(0));
                *dSpan = gMax(points.component(0)) - *dMin;
                return patch().Cf().component(0); // patchX
            }
            else // More perpendicular to X plane
            {
                *dMin = gMin(points.component(1));
                *dSpan = gMax(points.component(1)) - *dMin;
                return patch().Cf().component(1);  // patchY
            }
        }

        //- Return Zmax and Zmin of the patch faces
        void faceBoundsZ(scalarField* zSup, scalarField* zInf)
        {
            const label nF = patch().faceCells().size();
            scalarField zMax = Foam::scalarField(nF, -9999.0);
            scalarField zMin = Foam::scalarField(nF, 9999.0);
            
            const faceList& faces = this->patch().patch().localFaces();
            const List<vector>& points = this->patch().patch().localPoints();

            forAll( faces, faceI )
            {
                const face& f = faces[faceI];
                forAll(f, point)
                {
                    scalar auxiliar = points[f[point]].component(2);

                    zMax[faceI] = max(zMax[faceI], auxiliar);
                    zMin[faceI] = min(zMin[faceI], auxiliar);
                }
            }

            *zSup = zMax;
            *zInf = zMin;
        }

        //- Calculate water levels on each paddle (same zSpan)
        scalarList calcWL
        ( 
            scalarField alphaCell, 
            labelList cellGroup, 
            scalar zSpan, 
            label nPaddles
        )
        {
            scalarList groupTotalArea (nPaddles,0.0); // Paddle total area
            scalarList groupWaterArea (nPaddles,0.0); // Paddle wet area
            scalarList heights (nPaddles,0.0);

            const scalarField faceSurface = patch().magSf();

            forAll( patch().faceCells(), index ) 
            {
                groupTotalArea[cellGroup[index]-1] += faceSurface[index];
                groupWaterArea[cellGroup[index]-1] += 
                    faceSurface[index]*alphaCell[index];
            } 

            for (label i=0; i<=gMax(cellGroup)-1; i++)
            {
                reduce(groupTotalArea[i], sumOp<scalar>());
                reduce(groupWaterArea[i], sumOp<scalar>());
                // Free surface elevation at each paddle
                heights[i] = groupWaterArea[i]/groupTotalArea[i]*zSpan;
            }

            return heights;
        }

        //- Calculate water levels on each paddle (different zSpan)
        scalarList calcWL
        ( 
            scalarField alphaCell, 
            labelList cellGroup, 
            scalarList zSpanL, 
            label nPaddles
        )
        {
            scalarList groupTotalArea (nPaddles,0.0); // Paddle total area
            scalarList groupWaterArea (nPaddles,0.0); // Paddle wet area
            scalarList heights (nPaddles,0.0);

            const scalarField faceSurface = patch().magSf();

            forAll( patch().faceCells(), index ) 
            {
                groupTotalArea[cellGroup[index]-1] += faceSurface[index];
                groupWaterArea[cellGroup[index]-1] += 
                    faceSurface[index]*alphaCell[index];
            } 

            for (label i=0; i<=gMax(cellGroup)-1; i++)
            {
                reduce(groupTotalArea[i], sumOp<scalar>());
                reduce(groupWaterArea[i], sumOp<scalar>());
                // Free surface elevation at each paddle
                heights[i] = groupWaterArea[i]/groupTotalArea[i]*zSpanL[i];
            }

            return heights;
        }

        //- Reduce number of paddles if too high (paddles without cells)
        label decreaseNPaddles 
        (
            label np, 
            scalarField patchD, 
            scalar dMin, 
            scalar dSpan
        )
        {
            scalarList dBreakPoints(np+1, 0.0); // Paddles' start-end points

            while(1)
            {
                scalarList usedGroup (np,0.0); // If group is used = 1, else 0

                for (label i=0; i<=np; i++)
                {
                    dBreakPoints[i] = dMin + dSpan/(np)*i; // Breakpoint coords
                }

                forAll(patch().faceCells(), patchCells) 
                {
                    for (label j=0; j<=np; j++) // Grouping loop
                    {
                        if ( (patchD[patchCells]>=dBreakPoints[j]) && 
                            (patchD[patchCells]<dBreakPoints[j+1]) )
                        {
                            usedGroup[j] = 1.0;
                            continue;
                        }
                    }      
                }

                reduce(usedGroup, maxOp<scalarList>());

                if ( gMin(usedGroup) < 0.5 ) // Some group not used (=0)
                {
                    Info << 
                        "Reduced nPaddles from " << np << 
                        " to " << np-1 << endl;
                    np -= 1;
                }
                else // All groups are used
                {
                    break;
                }
            }
            return np;
        }
        
        //- Own convenient definition of atan
        scalar arcTan ( scalar x, scalar y )
        {
            scalar A = 0.0;

            if( x > 0.0 )
            {
                A = atan(y/x);
            }
            else if ( y >= 0.0 && x < 0.0 )
            {
                A = PI() + atan(y/x);
            }
            else if ( y < 0.0 && x < 0.0 )
            {
                A = -PI() + atan(y/x);
            }
            else if ( y > 0.0 && x == 0.0 )
            {
                A = PI()/2.0;
            }        
            else if ( y < 0.0 && x == 0.0 )
            {
                A = -PI()/2.0;
            }
            else if ( y == 0.0 && x == 0.0 )
            {
                A = 0;
            }

            return A;
        }

        //- Calculate mean horizontal angle for each paddle
        scalarList meanPatchDirs ( labelList cellGroup, label nPaddles )
        {
            const vectorField nVecCell = patch().nf();

            scalarList numAngles (nPaddles, 0.0);
            scalarList meanAngle (nPaddles, 0.0);

            forAll(patch().faceCells(), patchCells) 
            {
                meanAngle[cellGroup[patchCells]-1] += 
                    arcTan( nVecCell[patchCells].component(0), 
                        nVecCell[patchCells].component(1));
                numAngles[cellGroup[patchCells]-1] += 1.0;
            }

            for (label i=0; i<=nPaddles-1; i++)
            {
                reduce(meanAngle[i], sumOp<scalar>());
                reduce(numAngles[i], sumOp<scalar>());
                meanAngle[i] = meanAngle[i]/numAngles[i] + PI();
            }
            return meanAngle;
        }

        //- Calculate z-bounds for each paddle
        scalarList zSpanList 
        ( 
            labelList cellGroup, 
            scalarField zInf, 
            scalarField zSup, 
            label nPaddles
        )
        {
            scalarList zMaxL (nPaddles, -9999.0);
            scalarList zMinL (nPaddles,  9999.0);

            forAll(cellGroup, patchCells) 
            {
                zMaxL[cellGroup[patchCells]-1] = 
                    max(zMaxL[cellGroup[patchCells]-1], zSup[patchCells] );
                zMinL[cellGroup[patchCells]-1] = 
                    min(zMinL[cellGroup[patchCells]-1], zInf[patchCells] );
            }

            for (label i=0; i<=nPaddles-1; i++)
            {
                reduce(zMaxL[i], maxOp<scalar>());
                reduce(zMinL[i], minOp<scalar>());
            }
            return zMaxL-zMinL;
        }

        //- In-line print for scalarLists
        label inlinePrint ( std::string name, scalarList SL )
        {
            Info << name << " " << SL.size() << "( ";
            forAll(SL, i)
            {
                Info << SL[i] << " ";
            }
            Info << ")\n" << endl;
            return 0;
        }