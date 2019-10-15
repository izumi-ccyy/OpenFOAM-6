/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2016 OpenCFD Ltd.
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

#include "constant.H"
#include "addToRunTimeSelectionTable.H"
#include "fvcGrad.H"
#include "../threePhaseMixtureEThermo/threePhaseMixtureEThermo.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace temperaturePhaseChangeThreePhaseMixtures
{
    defineTypeNameAndDebug(constant, 0);
    addToRunTimeSelectionTable
    (
        temperaturePhaseChangeThreePhaseMixture,
        constant,
        components
    );
}
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::temperaturePhaseChangeThreePhaseMixtures::constant::constant
(
    const thermoIncompressibleThreePhaseMixture& mixture,
    const fvMesh& mesh
)
:
    temperaturePhaseChangeThreePhaseMixture(mixture, mesh),
    coeffC_
    (
        "coeffC", dimless/dimTime/dimTemperature, subDict(type() + "Coeffs")
    ),
    coeffE_
    (
        "coeffE", dimless/dimTime/dimTemperature, subDict(type() + "Coeffs")
    )
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::Pair<Foam::tmp<Foam::volScalarField>>
Foam::temperaturePhaseChangeThreePhaseMixtures::constant::mDotAlphal() const
{
    const volScalarField& T = mesh_.lookupObject<volScalarField>("T");

    const threePhaseMixtureEThermo& thermo =
        refCast<const threePhaseMixtureEThermo>
        (
            mesh_.lookupObject<basicThermo>(basicThermo::dictName)
        );

    const dimensionedScalar& TSat = thermo.TSat();

    const dimensionedScalar T0("T0", dimTemperature, Zero);

    return Pair<tmp<volScalarField>>
    (
        coeffC_*mixture_.rho2()*max(TSat - T.oldTime(), T0),
       -coeffE_*mixture_.rho1()*max(T.oldTime() - TSat, T0)
    );
}


Foam::Pair<Foam::tmp<Foam::volScalarField>>
Foam::temperaturePhaseChangeThreePhaseMixtures::constant::mDot() const
{

    volScalarField limitedAlpha1
    (
        min(max(mixture_.alpha1(), scalar(0)), scalar(1))
    );

    volScalarField limitedAlpha2
    (
        min(max(mixture_.alpha2(), scalar(0)), scalar(1))
    );

    const volScalarField& T = mesh_.lookupObject<volScalarField>("T");

    const threePhaseMixtureEThermo& thermo =
        refCast<const threePhaseMixtureEThermo>
        (
            mesh_.lookupObject<basicThermo>(basicThermo::dictName)
        );

    const dimensionedScalar& TSat = thermo.TSat();

    const dimensionedScalar T0("T0", dimTemperature, Zero);

    return Pair<tmp<volScalarField>>
    (
        coeffC_*mixture_.rho2()*limitedAlpha2*max(TSat - T.oldTime(), T0),
        coeffE_*mixture_.rho1()*limitedAlpha1*max(T.oldTime() - TSat, T0)
    );
}


Foam::Pair<Foam::tmp<Foam::volScalarField>>
Foam::temperaturePhaseChangeThreePhaseMixtures::constant::mDotDeltaT() const
{
    volScalarField limitedAlpha1
    (
        min(max(mixture_.alpha1(), scalar(0)), scalar(1))
    );

    volScalarField limitedAlpha2
    (
        min(max(mixture_.alpha2(), scalar(0)), scalar(1))
    );

    const volScalarField& T = mesh_.lookupObject<volScalarField>("T");

    const threePhaseMixtureEThermo& thermo =
        refCast<const threePhaseMixtureEThermo>
        (
            mesh_.lookupObject<basicThermo>(basicThermo::dictName)
        );

    const dimensionedScalar& TSat = thermo.TSat();

    return Pair<tmp<volScalarField>>
    (
        coeffC_*mixture_.rho2()*limitedAlpha2*pos(TSat - T.oldTime()),
        coeffE_*mixture_.rho1()*limitedAlpha1*pos(T.oldTime() - TSat)
    );
}


void Foam::temperaturePhaseChangeThreePhaseMixtures::constant::correct()
{
}


bool Foam::temperaturePhaseChangeThreePhaseMixtures::constant::read()
{
    if (temperaturePhaseChangeThreePhaseMixture::read())
    {
        /*
        subDict(type() + "Coeffs").readEntry("coeffC", coeffC_);
        subDict(type() + "Coeffs").readEntry("coeffE", coeffE_);
        */

        subDict(type() + "Coeffs").lookup("coeffC") >> coeffC_;
        subDict(type() + "Coeffs").lookup("coeffE") >> coeffE_;

        return true;
    }

    return false;
}


// ************************************************************************* //
