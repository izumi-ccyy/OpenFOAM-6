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

Class
    Foam::thermoIncompressibleThreePhaseMixture

Description
    A two-phase incompressible transportModel

SourceFiles
    thermoIncompressibleThreePhaseMixture.C

\*---------------------------------------------------------------------------*/

#ifndef thermoIncompressibleThreePhaseMixture_H
#define thermoIncompressibleThreePhaseMixture_H

// use immiscibleIncompressibleThreePhaseMixture.H
//#include "incompressibleTwoPhaseMixture.H"
#include "immiscibleIncompressibleThreePhaseMixture.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
              Class thermoIncompressibleThreePhaseMixture Declaration
\*---------------------------------------------------------------------------*/

class thermoIncompressibleThreePhaseMixture
:
    public immiscibleIncompressibleThreePhaseMixture
{
protected:

    // Protected data


        //- Thermal variables
        //- add Thermal variables for phase 3
        dimensionedScalar kappa1_;
        dimensionedScalar kappa2_;
        dimensionedScalar kappa3_;

        dimensionedScalar Cp1_;
        dimensionedScalar Cp2_;
        dimensionedScalar Cp3_;

        dimensionedScalar Cv1_;
        dimensionedScalar Cv2_;
        dimensionedScalar Cv3_;

        dimensionedScalar Hf1_;
        dimensionedScalar Hf2_;
        dimensionedScalar Hf3_;


public:

    TypeName("thermoIncompressibleThreePhaseMixture");


    // Constructors

        //- Construct from U and phi
        thermoIncompressibleThreePhaseMixture
        (
            const volVectorField& U,
            const surfaceScalarField& phi
        );


    //- Destructor
    virtual ~thermoIncompressibleThreePhaseMixture()
    {}


    // Access function
    // add function for phase 3

        //- Return const-access to phase1 kappa
        const dimensionedScalar& kappa1() const
        {
            return kappa1_;
        };

        //- Return const-access to phase2 kappa
        const dimensionedScalar& kappa2() const
        {
            return kappa2_;
        };

        //- Return const-access to phase3 kappa
        const dimensionedScalar& kappa3() const
        {
            return kappa3_;
        };

        //- Return const-access to phase1 Cp
        const dimensionedScalar& Cp1() const
        {
            return Cp1_;
        };

        //- Return const-access to phase2 Cp
        const dimensionedScalar& Cp2() const
        {
            return Cp2_;
        };

        //- Return const-access to phase3 Cp
        const dimensionedScalar& Cp3() const
        {
            return Cp3_;
        };

        //- Return const-access to phase1 Cv
        const dimensionedScalar& Cv1() const
        {
            return Cv1_;
        };

        //- Return const-access to phase2 Cv
        const dimensionedScalar& Cv2() const
        {
            return Cv2_;
        };

        //- Return const-access to phase3 Cv
        const dimensionedScalar& Cv3() const
        {
            return Cv3_;
        };

        //- Return latent heat for phase 1
        const dimensionedScalar& Hf1() const
        {
            return Hf1_;
        };

        //- Return latent heat for phase 2
        const dimensionedScalar& Hf2() const
        {
            return Hf2_;
        };

        //- Return latent heat for phase 3
        const dimensionedScalar& Hf3() const
        {
            return Hf3_;
        };

        //- Read base transportProperties dictionary
        virtual bool read();

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
