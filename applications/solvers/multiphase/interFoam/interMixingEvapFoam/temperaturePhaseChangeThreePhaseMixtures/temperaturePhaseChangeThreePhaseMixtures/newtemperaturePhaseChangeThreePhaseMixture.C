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

#include "temperaturePhaseChangeThreePhaseMixture.H"
#include "basicThermo.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::temperaturePhaseChangeThreePhaseMixture>
Foam::temperaturePhaseChangeThreePhaseMixture::New
(
    const thermoIncompressibleThreePhaseMixture& thermo,
    const fvMesh& mesh
)
{
    IOdictionary phaseChangePropertiesDict
    (
        IOobject
        (
            "phaseChangeProperties",
            mesh.time().constant(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE,
            false
        )
    );

    const word modelType
    (
        phaseChangePropertiesDict.lookup("phaseChangeThreePhaseModel")
    );

    Info<< "Selecting phaseChange model " << modelType << endl;

    auto cstrIter = componentsConstructorTablePtr_->find(modelType);

    if (cstrIter == componentsConstructorTablePtr_->end())
    {
        FatalErrorInFunction
            << "Unknown temperaturePhaseChangeThreePhaseMixture type "
            << modelType << nl << nl
            << "Valid temperaturePhaseChangeThreePhaseMixture types :" << endl
            << componentsConstructorTablePtr_->sortedToc()
            << exit(FatalError);
    }

    return autoPtr<temperaturePhaseChangeThreePhaseMixture>
        (cstrIter()(thermo, mesh));
}


// ************************************************************************* //
