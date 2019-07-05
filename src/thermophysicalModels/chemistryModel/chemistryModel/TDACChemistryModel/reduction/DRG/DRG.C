/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2016-2018 OpenFOAM Foundation
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

#include "DRG.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CompType, class ThermoType>
Foam::chemistryReductionMethods::DRG<CompType, ThermoType>::DRG
(
    const IOdictionary& dict,
    TDACChemistryModel<CompType, ThermoType>& chemistry
)
:
    chemistryReductionMethod<CompType, ThermoType>(dict, chemistry),
    searchInitSet_(this->coeffsDict_.subDict("initialSet").size())
{
    label j=0;
    //define initial set with dictionary
    dictionary initSet = this->coeffsDict_.subDict("initialSet");
    //search species i of the intial set in the mechanism
    //if species i is found, j = j + 1
    for (label i=0; i<chemistry.nSpecie(); i++)
    {
        if (initSet.found(chemistry.Y()[i].member()))
        {
            searchInitSet_[j++] = i;
        }
    }
    //if j < the element number of intial set
    //species in the initial set is not in the mechanism
    if (j<searchInitSet_.size())
    {
        FatalErrorInFunction
            << searchInitSet_.size()-j
            << " species in the initial set is not in the mechanism "
            << initSet
            << exit(FatalError);
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class CompType, class ThermoType>
Foam::chemistryReductionMethods::DRG<CompType, ThermoType>::~DRG()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CompType, class ThermoType>
void Foam::chemistryReductionMethods::DRG<CompType, ThermoType>::reduceMechanism
(
    const scalarField &c,
    const scalar T,
    const scalar p
)
{
    //define c1 = number of species + 2
    //+ 2 means temperature T and pressure p
    scalarField c1(this->nSpecie_+2, 0.0);
    //c for every species, T and p are assigned to c1
    for(label i=0; i<this->nSpecie_; i++)
    {
        c1[i] = c[i];
    }

    c1[this->nSpecie_] = T;
    c1[this->nSpecie_+1] = p;

    // Compute the rAB matrix
    //define the matrix of rABNum as a RectangularMatrix
    //with nSpecie row and nSpecie column, and the values
    //are assigned as 0.0
    //rABNum means numerator of rAB
    RectangularMatrix<scalar> rABNum(this->nSpecie_,this->nSpecie_,0.0);
    //define a scalarField rABDen with length of nSpecie, and 
    //the values are assigned as 0.0
    //rABDen means denominator of rAB
    scalarField rABDen(this->nSpecie_,0.0);

    // Number of initialized rAB for each lines
    Field<label> NbrABInit(this->nSpecie_,0);

    // Position of the initialized rAB, -1 when not initialized
    RectangularMatrix<label> rABPos(this->nSpecie_, this->nSpecie_, -1);

    // Index of the other species involved in the rABNum
    RectangularMatrix<label> rABOtherSpec(this->nSpecie_, this->nSpecie_, -1);

    scalar pf, cf, pr, cr;
    label lRef, rRef;
    //for all reaction equations
    forAll(this->chemistry_.reactions(), i)
    {
        const Reaction<ThermoType>& R = this->chemistry_.reactions()[i];
        // For each reaction compute omegai
        scalar omegai = this->chemistry_.omega
        (
         R, c1, T, p, pf, cf, lRef, pr, cr, rRef
         );


        // Then for each pair of species composing this reaction,
        // compute the rAB matrix (separate the numerator and
        // denominator)
        //define a list named wA and a list named wAID
        //with the length of the number of species in the
        //left hand of the reaction equation plus the number 
        //of species in the right hand of the reaction equation
        //the size of a DynamicList is defined by the number of
        //species in the left hand and in the right hand, so for 
        //every reaction equations, the size of wA and wAID may
        //be different, and for every reaction equation, wA and
        //wAID are empty list
        //
        DynamicList<scalar> wA(R.lhs().size()+R.rhs().size());
        DynamicList<label> wAID(R.lhs().size()+R.rhs().size());

        //for every species in the left hand of the reaction equation
        forAll(R.lhs(), s)
        {
            //define the index of species s as ss
            label ss = R.lhs()[s].index;
            //define the stoichiometric coefficient of s as sl
            //namely, sl = nu'
            scalar sl = -R.lhs()[s].stoichCoeff;
            //define found as false
            bool found(false);
            //for every species in the list wA
            forAll(wAID, id)
            {
                //if species s is in the wA
                if (ss==wAID[id])
                {
                    //wA for species s = wA + sl * omegai
                    //namely, wA[s] = wA[s] + nu's * omegai
                    wA[id] += sl*omegai;
                    //set found as true
                    found = true;
                    break;
                }
            }
            //if species s is in wA, then it is considered in last if
            //so next if is avoided. if the species s is not in wA, so
            //it will be appended to wA in next if
            if (!found)
            {
                wA.append(sl*omegai);
                wAID.append(ss);
            }
        }
        //the same process for the species in the right hand 
        //of the reaction equation
        forAll(R.rhs(), s)
        {
            label ss = R.rhs()[s].index;
            scalar sl = R.rhs()[s].stoichCoeff;
            bool found(false);
            forAll(wAID, id)
            {
                if (ss==wAID[id])
                {
                    wA[id] += sl*omegai;
                    found = true;
                    break;
                }
            }
            if (!found)
            {
                wA.append(sl*omegai);
                wAID.append(ss);
            }
        }

        // Now that all nuAi*wi are computed, without counting twice species
        // present in both rhs and lhs, we can update the denominator and
        // numerator for the rAB
        //DynamicList wAID shrinks to become a compact list
        wAID.shrink();
        //for every species in wAID
        forAll(wAID, id)
        {
            //curID means index of A
            label curID = wAID[id];

            // Absolute value of aggregated value
            //make all wA positive
            //a guess: the coefficients of products(speices in the right 
            //hand of reaction equation) are negative
            scalar curwA = ((wA[id]>=0) ? wA[id] : -wA[id]);

            //define a list deltaBi with length of nSpecie
            List<bool> deltaBi(this->nSpecie_, false);
            //define a first in first out stack usedIndex
            FIFOStack<label> usedIndex;
            //for all species in left hand of reaction equation
            forAll(R.lhs(), j)
            {
                //index of species j is sj
                label sj = R.lhs()[j].index;
                //push sj into the stack usedIndex
                usedIndex.push(sj);
                //assign deltaBi[sj] as true
                deltaBi[sj] = true;
            }
            //for all species in right hand of reaction equation
            //do the same operation as last for loop
            forAll(R.rhs(), j)
            {
                label sj = R.rhs()[j].index;
                usedIndex.push(sj);
                deltaBi[sj] = true;
            }

            // Disable for self reference (by definition rAA=0)
            //so rAA do not need to be calculated
            deltaBi[curID] = false;
            //do while the stack is not empty
            while(!usedIndex.empty())
            {
                //pop a index from the stack
                //curIndex means index of B
                label curIndex = usedIndex.pop();

                //deltaBi is initialized as true
                //if deltabi is not assigned a value
                if (deltaBi[curIndex])
                {
                    // Disable to avoid counting it more than once
                    //set deltaBi as false in case count it twice or more
                    deltaBi[curIndex] = false;

                    // Test if this rAB is not initialized
                    if (rABPos(curID, curIndex)==-1)
                    {
                        //the initial value of NbrANInit is 0.0
                        //so rABPos(curID, curIndex) (rABPos[A,B]) changes 
                        //from -1 to 0.0
                        //the number of edge of A
                        //B appears for the first time
                        rABPos(curID, curIndex) = NbrABInit[curID];
                        //edge number of A add 1
                        NbrABInit[curID]++;
                        //the edge of numerator of rAB is curwA
                        rABNum(curID, rABPos(curID, curIndex)) = curwA;
                        //store the index of B in rABOtherSpec
                        //rABPos(curID, curIndex) means the order of B
                        //in the set of edge from A
                        rABOtherSpec(curID, rABPos(curID, curIndex)) = curIndex;
                    }
                    //B does not appear for the first time
                    else
                    {
                        rABNum(curID, rABPos(curID, curIndex)) += curwA;
                    }
                }
            }
            //assign the value of denominator of rAB
            if (rABDen[curID] == 0.0)
            {
                //curwA is positive
                rABDen[curID] = curwA;
            }
            else
            {
                rABDen[curID] +=curwA;
            }
        }
    }
    // rii = 0.0 by definition

    //define the number of active species
    label speciesNumber = 0;

    // Set all species to inactive and activate them according
    // to rAB and initial set
    for (label i=0; i<this->nSpecie_; i++)
    {
        this->activeSpecies_[i] = false;
    }

    //define a FIFO stack Q
    FIFOStack<label> Q;

    // Initialize the list of active species with the search initiating set
    // (SIS)
    for (label i=0; i<searchInitSet_.size(); i++)
    {
        //make species in SIS active
        label q = searchInitSet_[i];
        this->activeSpecies_[q] = true;
        //the number of active species plus 1
        speciesNumber++;
        //push that species into the stack Q
        Q.push(q);
    }

    // Breadth first search with rAB
    while (!Q.empty())
    {
        //pop a index from Q and assign it to u
        label u = Q.pop();
        //define Den as the denominator of rAB of index u (or species A)
        scalar Den = rABDen[u];

        //if Den is not too small
        if (Den > vSmall)
        {
            //NBrABInit is the number of edge from A
            for (label v=0; v<NbrABInit[u]; v++)
            {
                //get label of B
                label otherSpec = rABOtherSpec(u, v);
                //clculate rAB
                scalar rAB = rABNum(u, v)/Den;

                if (rAB > 1)
                {
                    rAB = 1;
                }

                // Include B only if rAB is above the tolerance and if the
                // species was not searched before
                if
                (
                    rAB >= this->tolerance()
                 && !this->activeSpecies_[otherSpec]
                )
                {
                    //add B to the set
                    //namely push B into the stack Q
                    //because it is a FIFO stack, so it will start from
                    //B for the next search
                    Q.push(otherSpec);
                    //set B as an active species
                    this->activeSpecies_[otherSpec] = true;
                    //the number of active species plus 1
                    speciesNumber++;
                }
            }
        }
    }
    //after the loop of the stack, all the species set is found out

    // Put a flag on the reactions containing at least one removed species
    forAll(this->chemistry_.reactions(), i)
    {
        //for every reaction in the full mechanism
        const Reaction<ThermoType>& R = this->chemistry_.reactions()[i];
        //first set the reaction equation as active
        this->chemistry_.reactionsDisabled()[i] = false;

        //for every species in the left hand of the reaction equation
        forAll(R.lhs(), s)
        {
            //define the index of species s as ss
            label ss = R.lhs()[s].index;

            // The species is inactive then the reaction is removed
            //if species s is active, then it will not excute
            //if species s is inactive, then the reaction equation
            //is set to inactive and break the for loop
            if (!this->activeSpecies_[ss])
            {
                // Flag the reaction to disable it
                this->chemistry_.reactionsDisabled()[i] = true;
                break;
            }
        }

        // If the reaction has not been disabled yet
        //if there is no inactive species in the left hand
        //of the reaction equation, then the species in the
        //right hand is checked
        if (!this->chemistry_.reactionsDisabled()[i])
        {
            forAll(R.rhs(), s)
            {
                label ss = R.rhs()[s].index;
                if (!this->activeSpecies_[ss])
                {
                    this->chemistry_.reactionsDisabled()[i] = true;
                    break;
                }
            }
        }
    }
    //after this for loop, all reaction equations containing inactive
    //species are disabled (set as inactive) 

    //define NsSimp as the number of species for simplified mechanism
    //and equals to the active species
    this->NsSimp_ = speciesNumber;
    //define simplified c with length of NsSimp + 2, for temperature 
    //and pressure
    this->chemistry_.simplifiedC().setSize(this->NsSimp_+2);
    this->chemistry_.simplifiedToCompleteIndex().setSize(this->NsSimp_);

    label j = 0;
    //for every species in the full mechanism 
    for (label i=0; i<this->nSpecie_; i++)
    {
        //if species i is active
        if (this->activeSpecies_[i])
        {
            //then assign the index of species i to the index of
            //the simplifed set
            this->chemistry_.simplifiedToCompleteIndex()[j] = i;
            //c is also assigned
            this->chemistry_.simplifiedC()[j] = c[i];
            //move j forward
            this->chemistry_.completeToSimplifiedIndex()[i] = j++;
            //if species i is not active, then set it to active
            if (!this->chemistry_.active(i))
            {
                this->chemistry_.setActive(i);
            }
        }
        //if species i is inaactive, the index set to -1
        else
        {
            this->chemistry_.completeToSimplifiedIndex()[i] = -1;
        }
    }

    //set the last two value of simplifiedc as temperature and pressure
    this->chemistry_.simplifiedC()[this->NsSimp_] = T;
    this->chemistry_.simplifiedC()[this->NsSimp_+1] = p;
    this->chemistry_.setNsDAC(this->NsSimp_);

    // Change temporary Ns in chemistryModel
    // to make the function nEqns working
    this->chemistry_.setNSpecie(this->NsSimp_);
}


// ************************************************************************* //
