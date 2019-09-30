/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2016-2019 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "threePhaseMixtureEThermo.H"

#include "zeroGradientFvPatchFields.H"
#include "fixedEnergyFvPatchScalarField.H"
#include "gradientEnergyFvPatchScalarField.H"
#include "mixedEnergyFvPatchScalarField.H"
#include "threePhaseMixtureEThermo.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(threePhaseMixtureEThermo, 0);
}

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void Foam::threePhaseMixtureEThermo::eBoundaryCorrection(volScalarField& h)
{
    volScalarField::Boundary& hbf = h.boundaryFieldRef();

    forAll(hbf, patchi)
    {
        if (isA<gradientEnergyFvPatchScalarField>(hbf[patchi]))
        {
            refCast<gradientEnergyFvPatchScalarField>(hbf[patchi]).gradient()
                = hbf[patchi].fvPatchField::snGrad();
        }
        else if (isA<mixedEnergyFvPatchScalarField>(hbf[patchi]))
        {
            refCast<mixedEnergyFvPatchScalarField>(hbf[patchi]).refGrad()
                = hbf[patchi].fvPatchField::snGrad();
        }
    }
}


void Foam::threePhaseMixtureEThermo::init()
{
    const volScalarField alpha1Rho1(alpha1()*rho1());
    const volScalarField alpha2Rho2(alpha2()*rho2());

    e_ =
        (
            (T_ - TSat_)*(alpha1Rho1*Cv1() + alpha2Rho2*Cv2())
          + (alpha1Rho1*Hf1() + alpha2Rho2*Hf2())
        )
       /(alpha1Rho1 + alpha2Rho2);

    e_.correctBoundaryConditions();
}


// * * * * * * * * * * * * * * * * Constructor * * * * * * * * * * * * * * * //

Foam::threePhaseMixtureEThermo::threePhaseMixtureEThermo
(
    const volVectorField& U,
    const surfaceScalarField& phi
)
:
    basicThermo(U.mesh(), word::null),
    thermoIncompressibleThreePhaseMixture(U, phi),

    e_
    (
        volScalarField
        (
            IOobject
            (
                "e",
                U.mesh().time().timeName(),
                U.mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            U.mesh(),
            dimensionedScalar(dimEnergy/dimMass, Zero),
            heBoundaryTypes()
        )
    ),

    TSat_("TSat", dimTemperature, static_cast<const basicThermo&>(*this)),

    pDivU_(basicThermo::lookupOrDefault<Switch>("pDivU", true))

{
    // Initialise e
    init();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::threePhaseMixtureEThermo::correct() // correct T by alpha1 and alpha2
{
    incompressibleTwoPhaseMixture::correct();

    const volScalarField alpha1Rho1(alpha1()*rho1());
    const volScalarField alpha2Rho2(alpha2()*rho2());

    T_ =
        (
            (e_*(alpha1Rho1 + alpha2Rho2))
         -  (alpha1Rho1*Hf1() + alpha2Rho2*Hf2())
        )
       /(alpha1Rho1*Cv1() + alpha2Rho2*Cv2())
       + TSat_;

    T().correctBoundaryConditions();
}


Foam::word Foam::threePhaseMixtureEThermo::thermoName() const
{
    NotImplemented;
    return word::null;
}


Foam::tmp<Foam::volScalarField> Foam::threePhaseMixtureEThermo::he
(
    const volScalarField& p,
    const volScalarField& T
) const
{
    const volScalarField alpha1Rho1(alpha1()*rho1());
    const volScalarField alpha2Rho2(alpha2()*rho2());

    return
    (
        (T - TSat_)*(alpha1Rho1*Cv1() + alpha2Rho2*Cv2())
        + (alpha1Rho1*Hf1() + alpha2Rho2*Hf2())
    )
    / (alpha1Rho1 + alpha2Rho2);
}


Foam::tmp<Foam::scalarField> Foam::threePhaseMixtureEThermo::he
(
    const scalarField& p,
    const scalarField& T,
    const labelList& cells
) const
{
    tmp<scalarField> the(new scalarField(T.size()));
    scalarField& he = the.ref();

    const volScalarField alpha1Rho1(alpha1()*rho1());
    const volScalarField alpha2Rho2(alpha2()*rho2());

    forAll(T, i)
    {
        const label celli = cells[i];
        he[i] =
            (
                (T[i] - TSat_.value())
               *(
                   alpha1Rho1[celli]*Cv1().value()
                 + alpha2Rho2[celli]*Cv2().value()
                )
              + (
                    alpha1Rho1[celli]*Hf1().value()
                  + alpha2Rho2[celli]*Hf2().value()
                )
            )
            / (alpha1Rho1[celli] + alpha2Rho2[celli]);
    }

    return the;
}


Foam::tmp<Foam::scalarField> Foam::threePhaseMixtureEThermo::he
(
    const scalarField& p,
    const scalarField& T,
    const label patchi
) const
{
    const scalarField& alpha1p = alpha1().boundaryField()[patchi];
    const scalarField& alpha2p = alpha2().boundaryField()[patchi];

    const scalarField& Tp = T_.boundaryField()[patchi];

    return
    (
        (
            (Tp - TSat_.value())
           *(
               alpha1p*rho1().value()*Cv1().value()
             + alpha2p*rho2().value()*Cv2().value()
            )
          + (
               alpha1p*rho1().value()*Hf1().value()
             + alpha2p*rho2().value()*Hf2().value()
            )
        )
        / (alpha1p*rho1().value() + alpha2p*rho2().value())
    );
}


Foam::tmp<Foam::volScalarField> Foam::threePhaseMixtureEThermo::hc() const // define hc
{
    const fvMesh& mesh = this->T_.mesh();

    return tmp<volScalarField>::New
    (
        IOobject
        (
            "hc",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("hc", Hf2() - Hf1())
    );
}


Foam::tmp<Foam::scalarField> Foam::threePhaseMixtureEThermo::THE
(
    const scalarField& h,
    const scalarField& p,
    const scalarField& T0,      // starting temperature
    const labelList& cells
) const
{
    NotImplemented;
    return nullptr;
}


Foam::tmp<Foam::scalarField> Foam::threePhaseMixtureEThermo::THE
(
    const scalarField& h,
    const scalarField& p,
    const scalarField& T0,      // starting temperature
    const label patchi
) const
{
    NotImplemented;
    return nullptr;
}


Foam::tmp<Foam::volScalarField> Foam::threePhaseMixtureEThermo::Cp() const
{
    const volScalarField limitedAlpha1
    (
        min(max(alpha1_, scalar(0)), scalar(1))
    );

    return tmp<volScalarField>
    (
        new volScalarField
        (
            "cp",
            limitedAlpha1*Cp1() + (scalar(1) - limitedAlpha1)*Cp2()
        )
    );
}


Foam::tmp<Foam::scalarField> Foam::threePhaseMixtureEThermo::Cp
(
    const scalarField& p,
    const scalarField& T,
    const label patchi
) const
{
    const volScalarField limitedAlpha1
    (
        min(max(alpha1_, scalar(0)), scalar(1))
    );

    const scalarField& alpha1p = limitedAlpha1.boundaryField()[patchi];

    return
    (
        alpha1p*Cp1().value() + (scalar(1) - alpha1p)*Cp2().value()
    );
}


Foam::tmp<Foam::volScalarField> Foam::threePhaseMixtureEThermo::rho() const
{
    const volScalarField limitedAlpha1
    (
        min(max(alpha1_, scalar(0)), scalar(1))
    );

    return tmp<volScalarField>
    (
        new volScalarField
        (
            "rho",
            limitedAlpha1*rho1().value()
          + (scalar(1) - limitedAlpha1)*rho2().value()
        )
    );
}


Foam::tmp<Foam::scalarField> Foam::threePhaseMixtureEThermo::rho
(
    const label patchi
) const
{
    const volScalarField limitedAlpha1
    (
        min(max(alpha1_, scalar(0)), scalar(1))
    );

    const scalarField& alpha1p = limitedAlpha1.boundaryField()[patchi];

    return
    (
        alpha1p*rho1().value() + (scalar(1) - alpha1p)*rho2().value()
    );
}


Foam::tmp<Foam::volScalarField> Foam::threePhaseMixtureEThermo::Cv() const
{
    const volScalarField limitedAlpha1
    (
        min(max(alpha1_, scalar(0)), scalar(1))
    );

    return tmp<volScalarField>
    (
        new volScalarField
        (
            "cv",
            limitedAlpha1*Cv1() + (scalar(1) - limitedAlpha1)*Cv2()
        )
    );
}


Foam::tmp<Foam::scalarField> Foam::threePhaseMixtureEThermo::Cv
(
    const scalarField& p,
    const scalarField& T,
    const label patchi
) const
{
    const volScalarField limitedAlpha1
    (
        min(max(alpha1_, scalar(0)), scalar(1))
    );

    const scalarField& alpha1p = limitedAlpha1.boundaryField()[patchi];

    return
    (
        alpha1p*Cv1().value() + (scalar(1) - alpha1p)*Cv2().value()
    );
}


Foam::tmp<Foam::volScalarField> Foam::threePhaseMixtureEThermo::gamma() const
{
    return tmp<volScalarField>
    (
        (alpha1_*Cp1() + alpha2_*Cp2())/(alpha1_*Cv1() + alpha2_*Cv2())
    );
}


Foam::tmp<Foam::scalarField> Foam::threePhaseMixtureEThermo::gamma
(
    const scalarField& p,
    const scalarField& T,
    const label patchi
) const
{
    return
    (
        gamma()().boundaryField()[patchi]
    );
}


Foam::tmp<Foam::volScalarField> Foam::threePhaseMixtureEThermo::Cpv() const
{
     // This is an e thermo (Cpv = Cv)
     return Cv();
}


Foam::tmp<Foam::scalarField> Foam::threePhaseMixtureEThermo::Cpv
(
    const scalarField& p,
    const scalarField& T,
    const label patchi
) const
{
    // This is an e thermo (Cpv = Cv)
    return Cv(p, T, patchi);
}


Foam::tmp<Foam::volScalarField> Foam::threePhaseMixtureEThermo::CpByCpv() const
{
    NotImplemented;
    return nullptr;
}


Foam::tmp<Foam::volScalarField> Foam::threePhaseMixtureEThermo::W() const
{
    NotImplemented;
    return nullptr;
}


Foam::tmp<Foam::scalarField> Foam::threePhaseMixtureEThermo::CpByCpv
(
    const scalarField& p,
    const scalarField& T,
    const label patchi
) const
{
    NotImplemented;
    return nullptr;
}


Foam::tmp<Foam::volScalarField> Foam::threePhaseMixtureEThermo::kappa() const
{
    const volScalarField limitedAlpha1
    (
        min(max(alpha1_, scalar(0)), scalar(1))
    );

    return tmp<volScalarField>
    (
        new volScalarField
        (
            "kappa",
            limitedAlpha1*kappa1() + (scalar(1) - limitedAlpha1)*kappa2()
        )
    );
}


Foam::tmp<Foam::scalarField> Foam::threePhaseMixtureEThermo::kappa
(
    const label patchi
) const
{
    const volScalarField limitedAlpha1
    (
        min(max(alpha1_, scalar(0)), scalar(1))
    );

    const scalarField& alpha1p = limitedAlpha1.boundaryField()[patchi];

    return (alpha1p*kappa1().value() + (1 - alpha1p)*kappa2().value());
}


Foam::tmp<Foam::volScalarField> Foam::threePhaseMixtureEThermo::alphahe() const
{
    NotImplemented;
    return nullptr;
}


Foam::tmp<Foam::scalarField> Foam::threePhaseMixtureEThermo::alphahe
(
    const label patchi
) const
{
    NotImplemented;
    return nullptr;
}

Foam::tmp<Foam::volScalarField> Foam::threePhaseMixtureEThermo::kappaEff
(
    const volScalarField& kappat
) const
{
    tmp<Foam::volScalarField> kappaEff(kappa() + kappat);
    kappaEff.ref().rename("kappaEff");
    return kappaEff;
}


Foam::tmp<Foam::scalarField> Foam::threePhaseMixtureEThermo::kappaEff
(
    const scalarField& kappat,
    const label patchi
) const
{
    const volScalarField limitedAlpha1
    (
        min(max(alpha1_, scalar(0)), scalar(1))
    );

    const scalarField& alpha1p = limitedAlpha1.boundaryField()[patchi];

    return
        (alpha1p*kappa1().value() + (1 - alpha1p)*kappa2().value()) + kappat;

}


Foam::tmp<Foam::volScalarField> Foam::threePhaseMixtureEThermo::alphaEff
(
    const volScalarField& alphat
) const
{
    const volScalarField rho
    (
        alpha1_*rho1() + (1.0 - alpha1_)*rho2()
    );

    return (kappa()/Cp()/rho + alphat);
}


Foam::tmp<Foam::scalarField> Foam::threePhaseMixtureEThermo::alphaEff
(
    const scalarField& alphat,
    const label patchi
) const
{
    const volScalarField limitedAlpha1
    (
        min(max(alpha1_, scalar(0)), scalar(1))
    );

    const scalarField& alpha1p = limitedAlpha1.boundaryField()[patchi];

    const scalarField rho
    (
        alpha1p*rho1().value() + (1.0 - alpha1p)*rho2().value()
    );

    const scalarField kappa
    (
        alpha1p*kappa1().value() + (1.0 - alpha1p)*kappa2().value()
    );

    const scalarField Cp
    (
        alpha1p*Cp1().value() + (1.0 - alpha1p)*Cp2().value()
    );

    return kappa/Cp/rho + alphat;
}


bool Foam::threePhaseMixtureEThermo::read()
{
    if (basicThermo::read() && thermoIncompressibleThreePhaseMixture::read())
    {
        basicThermo::readIfPresent("pDivU", pDivU_);
        basicThermo::readEntry("TSat", TSat_);
        return true;
    }

    return false;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
