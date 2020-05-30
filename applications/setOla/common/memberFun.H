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

        //- Return Zmax and Zmin of the cells
        void cellBoundsZ(fvMesh& mesh, scalarField* zSup, scalarField* zInf)
        {
            const label nF = mesh.C().size();
            scalarField zMax = Foam::scalarField(nF, -9999.0);
            scalarField zMin = Foam::scalarField(nF, 9999.0);

            const faceList& ff = mesh.faces();
            const pointField& pp = mesh.points();

            forAll(mesh.C(), cellI)
            {
                const cell& cc = mesh.cells()[cellI];
                labelList pLabels(cc.labels(ff));
                pointField pLocal(pLabels.size(), vector::zero);

                forAll (pLabels, pointi)
                {
                    scalar auxiliar = pp[pLabels[pointi]].component(2);

                    zMax[cellI] = max(zMax[cellI], auxiliar);
                    zMin[cellI] = min(zMin[cellI], auxiliar);
                }
            }
            *zSup = zMax;
            *zInf = zMin;
        }

        //- Simple linear interpolation
        scalar interpolation 
        (
            scalar x1, 
            scalar x2, 
            scalar y1, 
            scalar y2, 
            scalar xInt
        )
        {
            scalar yInt = y1 + (y2-y1)/(x2-x1)*(xInt-x1);
            return yInt;
        }

        //- Limit angle between -pi/2 and pi/2
        scalar limAngle (scalar ang)
        {
            ang = abs(ang);

            while (ang >= 2.0*PI())
            {
                ang -= 2.0*PI();
            }

            if ( ang >= PI()/2.0 && ang <= 3.0*PI()/2.0 )
            {
                return PI()/2.0;
            }
            else
            {
                return ang;
            }
        }
