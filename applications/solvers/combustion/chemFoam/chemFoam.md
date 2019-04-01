# chemFoam

## thermo

in createFields.H 

```cpp
autoPtr<BasicChemistryModel<rhoReactionThermo>> pChemistry
(
    BasicChemistryModel<rhoReactionThermo>::New(thermo)
);
```

rhoReactionThermo is used, details are in src\thermophysicalModels\chemistryModel\chemistryModel\StandardChemistryModel\StandardChemistryModel.C

### StandardChemistryModel.H

variables 

```cpp
//- Reference to the field of specie mass fractions
PtrList<volScalarField>& Y_;

//- Reactions
const PtrList<Reaction<ThermoType>>& reactions_;

//- Thermodynamic data of the species
const PtrList<ThermoType>& specieThermo_;

//- Number of species
label nSpecie_;

//- Number of reactions
label nReaction_;

//- Temperature below which the reaction rates are assumed 0
scalar Treact_;

//- List of reaction rate per specie [kg/m3/s]
PtrList<volScalarField::Internal> RR_;

//- Temporary concentration field
mutable scalarField c_;

//- Temporary rate-of-change of concentration field
mutable scalarField dcdt_;


// Protected Member Functions

//- Write access to chemical source terms
//  (e.g. for multi-chemistry model)
inline PtrList<volScalarField::Internal>& RR();
```

### StandardChemistryModel.C

#### calculateRR

##### calculate concentration $C_k$

the concentration of a species k:

$$
C_k = \frac{Y_k \rho}{(MW)_k}
$$

```cpp
const scalar rhoi = rho[celli];
const scalar Ti = T[celli];
const scalar pi = p[celli];

for (label i=0; i<nSpecie_; i++)
{
    const scalar Yi = Y_[i][celli];
    c_[i] = rhoi*Yi/specieThermo_[i].W();
}
```

compare the equation and the code, we can find that 

+ rhoi is total density ($\rho$)
+ pi is total pressure ($p$) 
+ i in code means species i (or k in equation)
+ Yi is the mass fraction of species k ($Y_k$)
+ specieThermo_[i].W() is the molecular weight of species i (or $(MW)_k$)

#####

```cpp
//get omega with 10 variables
const Reaction<ThermoType>& R = reactions_[ri];
const scalar omegai = R.omega(pi, Ti, c_, pf, cf, lRef, pr, cr, rRef);

//get left RR
forAll(R.lhs(), s)
{
    if (si == R.lhs()[s].index)
    {
        RR[celli] -= R.lhs()[s].stoichCoeff*omegai;
    }
}

//get right RR
forAll(R.rhs(), s)
{
    if (si == R.rhs()[s].index)
    {
        RR[celli] += R.rhs()[s].stoichCoeff*omegai;
    }
}

RR[celli] *= specieThermo_[si].W();
}

//the source term used in species transport equation
return tRR;
```

tRR is the source term used in species transport equation

$$
tRR = \sum^{m}_{j=1}\left((v''_{k} - v'_{k})(k_f \prod_{k}^{Nl} (C_{M_{kl}})^{v'_{kj}}  - k_r \prod_{k}^{Nr} (C_{M_{kr}})^{v''_{kj}})\right) (MW)_k
$$

### TDACChemistryModel.C

located in src\thermophysicalModels\chemistryModel\chemistryModel\TDACChemistryModel\TDACChemistryModel.C

#### omega with 4 variables

```cpp
template<class ReactionThermo, class ThermoType>
void Foam::TDACChemistryModel<ReactionThermo, ThermoType>::omega
(
    const scalarField& c, // Contains all species even when mechRed is active
    const scalar T,
    const scalar p,
    scalarField& dcdt
) const
```

```cpp
scalar omegai = R.omega
(
    p, T, c, pf, cf, lRef, pr, cr, rRef
);

forAll(R.lhs(), s)
{
    label si = R.lhs()[s].index;
    if (reduced)
    {
        si = completeToSimplifiedIndex_[si];
    }

    const scalar sl = R.lhs()[s].stoichCoeff;
    dcdt[si] -= sl*omegai;
}
forAll(R.rhs(), s)
{
    label si = R.rhs()[s].index;
    if (reduced)
    {
        si = completeToSimplifiedIndex_[si];
    }

    const scalar sr = R.rhs()[s].stoichCoeff;
    dcdt[si] += sr*omegai;
}
```



#### omega with 10 variables

```cpp
template<class ReactionThermo, class ThermoType>
Foam::scalar Foam::TDACChemistryModel<ReactionThermo, ThermoType>::omega
(
    const Reaction<ThermoType>& R,
    const scalarField& c, // Contains all species even when mechRed is active
    const scalar T,
    const scalar p,
    scalar& pf,
    scalar& cf,
    label& lRef,
    scalar& pr,
    scalar& cr,
    label& rRef
) const
```

##### pf

```cpp
pf = kf;
for (label s=1; s<Nl; s++)
{
    const label si = R.lhs()[s].index;

    if (c[si] < c[lRef])
    {
        const scalar exp = R.lhs()[slRef].exponent;
        pf *= pow(max(c[lRef], 0), exp);
        lRef = si;
        slRef = s;
    }
    else
    {
        const scalar exp = R.lhs()[s].exponent;
        pf *= pow(max(c[si], 0), exp);
    }
}
cf = max(c[lRef], 0);

{
    const scalar exp = R.lhs()[slRef].exponent;
    if (exp < 1)
    {
        if (cf > small)
        {
            pf *= pow(cf, exp - 1);
        }
        else
        {
            pf = 0;
        }
    }
    else
    {
        pf *= pow(cf, exp - 1);
    }
}
```

calculate from large concentration to small concentration with the equation:

$$
\mathbf{RR} = \frac{d(C_{product})}{dt} = k_f \prod_{k}^{N} (C_{M_k})^{v'_{kj}}
$$

in the code

$$
pf = k_f \prod_{k \not = least}^{Nl - 1} (C_{M_{kl}})^{v'_{kj}}
$$

for the species k with the smallest concentration $C_k$ (cf)

if exponent < 1 and cf is very small, then $pf = 0$, because pf multiply a very small value equals zero.

else:

$$
pf = k_f \prod_{k \not = least}^{Nl -1} (C_{M_{kl}})^{v'_{kj}} (C_{M_{least}l})^{v'_{leastj}-1} 
$$

##### pr

calculate pr

```cpp
label srRef = 0;
rRef = R.rhs()[srRef].index;

// Find the matrix element and element position for the rhs
pr = kr;
for (label s=1; s<Nr; s++)
{
    const label si = R.rhs()[s].index;
    if (c[si] < c[rRef])
    {
        const scalar exp = R.rhs()[srRef].exponent;
        pr *= pow(max(c[rRef], 0), exp);
        rRef = si;
        srRef = s;
    }
    else
    {
        const scalar exp = R.rhs()[s].exponent;
        pr *= pow(max(c[si], 0), exp);
    }
}
cr = max(c[rRef], 0);

{
    const scalar exp = R.rhs()[srRef].exponent;
    if (exp < 1)
    {
        if (cr > small)
        {
            pr *= pow(cr, exp - 1);
        }
        else
        {
            pr = 0;
        }
    }
    else
    {
        pr *= pow(cr, exp - 1);
    }
}
```

similar to pf

$$
pr = k_r \prod_{k \not = least}^{Nr -1} (C_{M_{kr}})^{v''_{kj}} (C_{M_{least}r})^{v''_{leastj}-1} 
$$

##### return omega

```cpp
return pf*cf - pr*cr;
```

because for the smallest concentration, the exponent is `exp-1`, so at the `return` the left 1 is multiplied, and we have

$$
omega = k_f \prod_{k \not = least}^{Nl -1} (C_{M_{kl}})^{v'_{kj}} (C_{M_{least}l})^{v'_{leastj}} - k_r \prod_{k \not = least}^{Nr -1} (C_{M_{kr}})^{v''_{kj}} (C_{M_{least}r})^{v''_{leastj}} 
$$

namely

$$
\omega = k_f \prod_{k}^{Nl} (C_{M_{kl}})^{v'_{kj}}  - k_r \prod_{k}^{Nr} (C_{M_{kr}})^{v''_{kj}}
$$

**NOTE: This omega does not contain $(v''_{k}-v'_{k})$**