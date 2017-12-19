    Info<< "Reading field porosityIndex\n" << endl;

    volScalarField porosityIndex
    (
        IOobject
        (
            "porosityIndex",
            runTime.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        mesh,
        0.0
    );

    scalar nPor = gMax(porosityIndex); // Number of porous media
    bool activePorosity = (nPor > 0.0);

    volScalarField porosity
    (
        IOobject
        (
            "porosity",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        1.0
    );

    surfaceScalarField porosityF = fvc::interpolate(porosity);

    // Porosity variables
    scalarList aPor, bPor, cPor, D50Por, phiPor; 
    bool debugPor = false;
    bool useTransient = false;
    scalar KC = 1.0;

    // Define fields related with porosity
    volScalarField aPorField = volScalarField("aPorField", 0.0 * porosityIndex);
    volScalarField bPorField = volScalarField("bPorField", aPorField);
    volScalarField cPorField = volScalarField("cPorField", aPorField);
    volScalarField useTransMask = volScalarField("useTrans", aPorField);
    volScalarField KCPorField = volScalarField("KCPorField", aPorField + 1.0);

    volScalarField D50Field
    (
        IOobject
        (
            "D50Field",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar( "D50", dimLength, 1.0 )
    );

    if( activePorosity ) 
    {
        Info << "Porosity activated\n"
            << "Reading porosity variables\n" << endl;

        // Define dictionary
        IOdictionary porosityDict
        (
            IOobject
            (
                "porosityDict",
                runTime.constant(),
                mesh,
                IOobject::MUST_READ,
                IOobject::NO_WRITE
            )
        );

        // Read porosity variables
        aPor = porosityDict.lookupOrDefault("a", List<scalar>(1, 0.0));
        bPor = porosityDict.lookupOrDefault("b", List<scalar>(1, 0.0));
        cPor = porosityDict.lookupOrDefault("c", List<scalar>(1, 0.0));
        D50Por  = porosityDict.lookupOrDefault("D50", List<scalar>(1, 1.0));
        phiPor  = porosityDict.lookupOrDefault("porosity", List<scalar>(1, 1.0));
        debugPor = porosityDict.lookupOrDefault<bool>("debugPor", false );
        useTransient = porosityDict.lookupOrDefault<bool>("useTransient", false);
        KC = porosityDict.lookupOrDefault<scalar>("KC", 1.0);

        Info << "Number of materials: " << nPor+1.0 << endl;
        Info << "a        = " << aPor << endl;
        Info << "b        = " << bPor << endl;
        Info << "c        = " << cPor << endl;
        Info << "D50      = " << D50Por << endl;
        Info << "porosity = " << phiPor << endl;
        Info << "Transient formulation is " << useTransient << endl;
        if ( useTransient )
        {
            Info << "KC is " << KC << endl;
        }
        
        if (aPor.size() != bPor.size() || aPor.size() != cPor.size() || 
            aPor.size() != phiPor.size() || aPor.size() != D50Por.size() || 
            aPor.size() != nPor+1)
        {
            FatalError
                << "Check the number of components for aPor, bPor, dPor,"
                << " phiPor and D50Por within porosityDict or the maximum"
                << " index within porosityIndex field"
            << exit(FatalError);
        }

        Info << "\nCreating fields:"
            << " porosity, aPor, bPor, cPor and D50Por\n" << endl;

        forAll(porosityIndex, item)
        {
            if( porosityIndex[item] > 0.0 )
            {
                aPorField[item] = aPor[porosityIndex[item]];
                bPorField[item] = bPor[porosityIndex[item]];
                cPorField[item] = cPor[porosityIndex[item]];
                KCPorField[item] = KC;
                D50Field[item]  = D50Por[porosityIndex[item]];
                porosity[item]  = phiPor[porosityIndex[item]];
                if ( useTransient )
                {
                    useTransMask[item] = 1.0;
                } 
            }
        }
        
        porosityF = fvc::interpolate(porosity);

        // Write out porosity
        porosity.write();

        // Write other fields (debug)
        if (debugPor)
        {
            aPorField.write();
            bPorField.write();
            cPorField.write();
            D50Field.write();
        }
    }
    else
    {
        Info << "Porosity NOT activated\n" << endl;
    }
