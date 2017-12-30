{
    if (mesh.changing())
    {
        forAll(U.boundaryField(), patchi)
        {
            if (U.boundaryField()[patchi].fixesValue())
            {
                U.boundaryField()[patchi].initEvaluate();
            }
        }

        forAll(U.boundaryField(), patchi)
        {
            if (U.boundaryField()[patchi].fixesValue())
            {
                U.boundaryField()[patchi].evaluate();

                phi.boundaryField()[patchi] =
                U.boundaryField()[patchi] & mesh.Sf().boundaryField()[patchi];
            }
        }
    }

#   include "continuityErrs.H"

    volScalarField pcorr
    (
        IOobject
        (
            "pcorr",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("pcorr", pd.dimensions(), 0.0),
        pcorrTypes
    );

    dimensionedScalar rAUf
    (
        "(1|A(U))",
        dimTime/rho.dimensions(),
        runTime.deltaT().value()
    );

    phi = (fvc::interpolate(U) & mesh.Sf());

    adjustPhi(phi, U, pcorr);

    for(int nonOrth=0; nonOrth<=nNonOrthCorr; nonOrth++)
    {
        fvScalarMatrix pcorrEqn
        (
            fvm::laplacian(rAUf, pcorr) == fvc::div(phi)
        );

        pcorrEqn.setReference(pdRefCell, pdRefValue);
        pcorrEqn.solve();

        if (nonOrth == nNonOrthCorr)
        {
            phi -= pcorrEqn.flux();
        }
    }

#   include "continuityErrs.H"
#   include "CourantNo.H"

    // Recalculate rhoPhi from rho
    rhoPhi = fvc::interpolate(rho)*phi;
}
