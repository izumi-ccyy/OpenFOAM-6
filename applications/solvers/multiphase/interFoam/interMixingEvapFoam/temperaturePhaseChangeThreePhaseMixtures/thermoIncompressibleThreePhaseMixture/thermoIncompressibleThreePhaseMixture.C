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

#include "./thermoIncompressibleThreePhaseMixture.H"


// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(thermoIncompressibleThreePhaseMixture, 0);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::thermoIncompressibleThreePhaseMixture::thermoIncompressibleThreePhaseMixture
(
    const volVectorField& U,
    const surfaceScalarField& phi
)
:
    immiscibleIncompressibleThreePhaseMixture(U, phi),

    // add variables for phase 3
    kappa1_
    (
        "kappa1",
        dimEnergy/dimTime/dimLength/dimTemperature,
        subDict(phase1Name_).lookup("kappa")
    ),

    kappa2_
    (
        "kappa2",
        kappa1_.dimensions(),
        subDict(phase2Name_).lookup("kappa")
    ),

    kappa3_
    (
        "kappa3",
        kappa1_.dimensions(),
        subDict(phase3Name_).lookup("kappa")
    ),

    Cp1_
    (
        "Cp1",
        dimEnergy/dimTemperature/dimMass,
        subDict(phase1Name_).lookup("Cp")
    ),

    Cp2_
    (
        "Cp2",
        dimEnergy/dimTemperature/dimMass,
        subDict(phase2Name_).lookup("Cp")
    ),

    Cp3_
    (
        "Cp3",
        dimEnergy/dimTemperature/dimMass,
        subDict(phase3Name_).lookup("Cp")
    ),

    Cv1_
    (
        "Cv1",
        dimEnergy/dimTemperature/dimMass,
        subDict(phase1Name_).lookup("Cv")
    ),

    Cv2_
    (
        "Cv2",
        dimEnergy/dimTemperature/dimMass,
        subDict(phase2Name_).lookup("Cv")
    ),

    Cv3_
    (
        "Cv3",
        dimEnergy/dimTemperature/dimMass,
        subDict(phase3Name_).lookup("Cv")
    ),

    Hf1_
    (
        "Hf1",
        dimEnergy/dimMass,
        subDict(phase1Name_).lookup("hf")
    ),

    Hf2_
    (
        "Hf2",
        dimEnergy/dimMass,
        subDict(phase2Name_).lookup("hf")
    ),

    Hf3_
    (
        "Hf3",
        dimEnergy/dimMass,
        subDict(phase3Name_).lookup("hf")
    )
{

}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

bool Foam::thermoIncompressibleThreePhaseMixture::read()
{
    if (immiscibleIncompressibleThreePhaseMixture::read())
    {
        // add functions for phase 3
        /*
        subDict(phase1Name_).readEntry("kappa", kappa1_);
        subDict(phase2Name_).readEntry("kappa", kappa2_);
        subDict(phase3Name_).readEntry("kappa", kappa3_);

        subDict(phase1Name_).readEntry("Cp", Cp1_);
        subDict(phase2Name_).readEntry("Cp", Cp2_);
        subDict(phase3Name_).readEntry("Cp", Cp3_);

        subDict(phase1Name_).readEntry("Cv", Cv1_);
        subDict(phase2Name_).readEntry("Cv", Cv2_);
        subDict(phase3Name_).readEntry("Cv", Cv3_);

        subDict(phase1Name_).readEntry("hf", Hf1_);
        subDict(phase2Name_).readEntry("hf", Hf2_);
        subDict(phase3Name_).readEntry("hf", Hf3_);
        */
        subDict(phase1Name_).lookup("kappa") >> kappa1_;
        subDict(phase2Name_).lookup("kappa") >> kappa2_;
        subDict(phase3Name_).lookup("kappa") >> kappa3_;

        subDict(phase1Name_).lookup("Cp") >> Cp1_;
        subDict(phase2Name_).lookup("Cp") >> Cp2_;
        subDict(phase3Name_).lookup("Cp") >> Cp3_;

        subDict(phase1Name_).lookup("Cv") >> Cv1_;
        subDict(phase2Name_).lookup("Cv") >> Cv2_;
        subDict(phase3Name_).lookup("Cv") >> Cv3_;

        subDict(phase1Name_).lookup("hf") >> Hf1_;
        subDict(phase2Name_).lookup("hf") >> Hf2_;
        subDict(phase3Name_).lookup("hf") >> Hf3_;

        return true;
    }

    return false;
}


// ************************************************************************* //
