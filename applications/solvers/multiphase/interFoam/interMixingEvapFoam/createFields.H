#include "createRDeltaT.H"

Info<< "Reading field p_rgh\n" << endl;
volScalarField p_rgh
(
    IOobject
    (
        "p_rgh",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info<< "Reading field U\n" << endl;
volVectorField U
(
    IOobject
    (
        "U",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

#include "createPhi.H"

// Creating e based thermo
autoPtr<threePhaseMixtureEThermo> thermo // thermo, defined by twoPhaseMixtureEThermo
(
    new threePhaseMixtureEThermo(U, phi)
);

// Create mixture and
Info<< "Creating temperaturePhaseChangeThreePhaseMixture\n" << endl;
autoPtr<temperaturePhaseChangeThreePhaseMixture> mixture =
    temperaturePhaseChangeThreePhaseMixture::New(thermo(), mesh); // mixture, defined by left


// three alpha
volScalarField& alpha1(thermo->alpha1());
volScalarField& alpha2(thermo->alpha2());
volScalarField& alpha3(thermo->alpha3());

// three rho
const dimensionedScalar& rho1 = thermo->rho1();
const dimensionedScalar& rho2 = thermo->rho2();
const dimensionedScalar& rho3 = thermo->rho3();

dimensionedScalar D23("D23", dimViscosity, mixture);

// Need to store rho for ddt(rho, U)
volScalarField rho
(
    IOobject
    (
        "rho",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT
    ),
    alpha1*rho1 + alpha2*rho2 + alpha3*rho3
);
rho.oldTime();


// Mass flux
// Initialisation does not matter because rhoPhi is reset after the
// alpha solution before it is used in the U equation.
surfaceScalarField rhoPhi
(
    IOobject
    (
        "rhoPhi",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    rho1*phi
);


// Construct incompressible turbulence model
autoPtr<incompressible::turbulenceModel> turbulence
(
    incompressible::turbulenceModel::New(U, phi, mixture)
);


#include "readGravitationalAcceleration.H"
#include "readhRef.H"
#include "gh.H"


volScalarField& p = thermo->p(); // define pressure from thermo

label pRefCell = 0;
scalar pRefValue = 0.0;
setRefCell
(
    p,
    p_rgh,
    pimple.dict(),
    pRefCell,
    pRefValue
);

if (p_rgh.needReference())
{
    p += dimensionedScalar
    (
        "p",
        p.dimensions(),
        pRefValue - getRefCellValue(p, pRefCell)
    );
    p_rgh = p - rho*gh;
}

mesh.setFluxRequired(p_rgh.name());
mesh.setFluxRequired(alpha2.name());

// Turbulent Prandtl number
dimensionedScalar Prt("Prt", dimless, thermo->transportPropertiesDict()); // define Prt from thermo

volScalarField kappaEff // define kappaEff by thermo
(
    IOobject
    (
        "kappaEff",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    thermo->kappa()
);

// Need to store rho for ddt(rhoCp, U)
volScalarField rhoCp // define rho * cp by thermo
(
    IOobject
    (
        "rhoCp",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    rho*thermo->Cp()
);

rhoCp.oldTime();

#include "createMRF.H"
#include "createFvOptions.H"
