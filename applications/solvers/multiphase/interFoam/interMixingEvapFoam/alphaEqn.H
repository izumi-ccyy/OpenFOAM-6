{
    word alphaScheme("div(phi,alpha)");
    word alpharScheme("div(phirb,alpha)");

    surfaceScalarField phir // phir = c_alpha * phi / Sf * nHatf
    (
        IOobject
        (
            "phir",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        thermo().cAlpha()*mag(phi/mesh.magSf())*thermo().nHatf()
    );

    for (int gCorr=0; gCorr<nAlphaCorr; gCorr++)
    {
        // Create the limiter to be used for all phase-fractions
        scalarField allLambda(mesh.nFaces(), 1.0);

        // Split the limiter into a surfaceScalarField
        slicedSurfaceScalarField lambda
        (
            IOobject
            (
                "lambda",
                mesh.time().timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            mesh,
            dimless,
            allLambda,
            false   // Use slices for the couples
        );

        // add mass flow rate
        Pair<tmp<volScalarField>> vDotAlphalA = mixture->mDot();

        const volScalarField& vDotcAlphalA = vDotAlphalA[0]();
        const volScalarField& vDotvAlphalA = vDotAlphalA[1]();
        const volScalarField vDotvmcAlphalA(vDotvAlphalA - vDotcAlphalA);

        Info<< "alpha Equation 1" << nl << endl;

        // Create the complete convection flux for alpha1
        surfaceScalarField alphaPhi1 // add compression
        (
            fvc::flux
            (
                phi,
                alpha1,
                alphaScheme
            )
          + fvc::flux // compression for alpha2
            (
                -fvc::flux(-phir, alpha2, alpharScheme),
                alpha1,
                alpharScheme
            )
          + fvc::flux // compression for alpha3
            (
                -fvc::flux(-phir, alpha3, alpharScheme),
                alpha1,
                alpharScheme
            )
          //+ fvc::interpolate(vDotcAlphalA/rho1)  // add mass loss
        );

        // Create the bounded (upwind) flux for alpha1
        surfaceScalarField alphaPhi1BD
        (
            upwind<scalar>(mesh, phi).flux(alpha1)
        );

        // Calculate the flux correction for alpha1
        alphaPhi1 -= alphaPhi1BD;

        // Calculate the limiter for alpha1
        if (LTS)
        {
            const volScalarField& rDeltaT =
                fv::localEulerDdt::localRDeltaT(mesh);

            MULES::limiter
            (
                allLambda,
                rDeltaT,
                geometricOneField(),
                alpha1,
                alphaPhi1BD,
                alphaPhi1,
                zeroField(),
                zeroField(),
                oneField(),
                zeroField()
            );
        }
        else
        {
            MULES::limiter
            (
                allLambda,
                1.0/runTime.deltaT().value(),
                geometricOneField(),
                alpha1,
                alphaPhi1BD,
                alphaPhi1,
                zeroField(),
                zeroField(),
                oneField(),
                zeroField()
            );
        }

        Info<< "alpha Equation 2" << nl << endl;

        // Create the complete flux for alpha2
        surfaceScalarField alphaPhi2 // interface only exsists between 1 and 2, 1 and 3
        (
            fvc::flux
            (
                phi,
                alpha2,
                alphaScheme
            )
          + fvc::flux // add compression between 1 and 2
            (
                -fvc::flux(phir, alpha1, alpharScheme),
                alpha2,
                alpharScheme
            )
          //+ fvc::interpolate(vDotcAlphalA/rho2) // add mass loss
        );

        // Create the bounded (upwind) flux for alpha2
        surfaceScalarField alphaPhi2BD
        (
            upwind<scalar>(mesh, phi).flux(alpha2)
        );

        // Calculate the flux correction for alpha2
        alphaPhi2 -= alphaPhi2BD;

        // Further limit the limiter for alpha2
        if (LTS)
        {
            const volScalarField& rDeltaT =
                fv::localEulerDdt::localRDeltaT(mesh);

            MULES::limiter
            (
                allLambda,
                rDeltaT,
                geometricOneField(),
                alpha2,
                alphaPhi2BD,
                alphaPhi2,
                zeroField(),
                zeroField(),
                oneField(),
                zeroField()
            );
        }
        else
        {
            MULES::limiter
            (
                allLambda,
                1.0/runTime.deltaT().value(),
                geometricOneField(),
                alpha2,
                alphaPhi2BD,
                alphaPhi2,
                zeroField(),
                zeroField(),
                oneField(),
                zeroField()
            );
        }

        // Construct the limited fluxes
        alphaPhi1 = alphaPhi1BD + lambda*alphaPhi1;
        alphaPhi2 = alphaPhi2BD + lambda*alphaPhi2;

        // Solve for alpha1
        solve(fvm::ddt(alpha1) + fvc::div(alphaPhi1));

        // Create the diffusion coefficients for alpha2<->alpha3
        volScalarField Dc23(D23*max(alpha3, scalar(0))*pos0(alpha2));
        volScalarField Dc32(D23*max(alpha2, scalar(0))*pos0(alpha3));

        // Add the diffusive flux for alpha3->alpha2
        alphaPhi2 -= fvc::interpolate(Dc32)*mesh.magSf()*fvc::snGrad(alpha1);

        // Solve for alpha2
        fvScalarMatrix alpha2Eqn
        (
            fvm::ddt(alpha2)
          + fvc::div(alphaPhi2)
          - fvm::laplacian(Dc23 + Dc32, alpha2)
        );
        alpha2Eqn.solve();

        // Construct the complete mass flux
        rhoPhi = // alphaPhi1 rho1 + alphaPhi2 rho2 + alphaPhi3 rho3
              alphaPhi1*(rho1 - rho3)
            + (alphaPhi2 + alpha2Eqn.flux())*(rho2 - rho3)
            + phi*rho3;

        alpha3 = 1.0 - alpha1 - alpha2;
    }

    Info<< "Air phase volume fraction = "
        << alpha1.weightedAverage(mesh.V()).value()
        << "  Min(" << alpha1.name() << ") = " << min(alpha1).value()
        << "  Max(" << alpha1.name() << ") = " << max(alpha1).value()
        << endl;

    Info<< "Liquid phase volume fraction = "
        << alpha2.weightedAverage(mesh.V()).value()
        << "  Min(" << alpha2.name() << ") = " << min(alpha2).value()
        << "  Max(" << alpha2.name() << ") = " << max(alpha2).value()
        << endl;
}
