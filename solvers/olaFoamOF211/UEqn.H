    surfaceScalarField muEff
    (
        "muEff",
        twoPhaseProperties.muf()
      + fvc::interpolate(rho*turbulence->nut())
    );

    fvVectorMatrix UEqn
    (
        (1.0 + cPorField) / porosity * fvm::ddt(rho, U)
        + 1.0/porosity * fvm::div(rhoPhi/porosityF, U)
        - fvm::laplacian(muEff/porosityF , U)
        - 1.0/porosity * ( fvc::grad(U) & fvc::grad(muEff) )
        // Closure Terms
        + aPorField * pow(1.0 - porosity, 3) / pow(porosity,3)
            * twoPhaseProperties.mu() / pow(D50Field,2) * U
        + bPorField * rho * (1.0 - porosity) / pow(porosity,3) / D50Field
            * mag(U) * U * 
        // Transient formulation
        (1.0 + useTransMask * 7.5 / KCPorField)
    );

    UEqn.relax();

    if (pimple.momentumPredictor())
    {
        solve
        (
            UEqn
         ==
            fvc::reconstruct
            (
                (
                    fvc::interpolate(interface.sigmaK())*fvc::snGrad(alpha1)
                  - ghf*fvc::snGrad(rho)
                  - fvc::snGrad(p_rgh)
                ) * mesh.magSf()
            )
        );
    }
