# Preparation

+ copy and rename `interMixingFoam` as `interMixingEvapFoam`
+ rename `interMixingFoam.C` as `interMixingEvapFoam.C`
+ modify Make folder, modify `files`
+ copy `UEqn.H` and `pEqn.H` from interFoam to here, since `pEqn.H` should be modified later
+ copy `temperaturePhaseChangeTwoPhaseMixtures` folder and `Allwmake` and `Allwclean` from `interCondensationEvaporationFoam` of `OpenFOAM-PLUS`

## temperaturePhaseChangeTwoPhaseMixtures

+ rename `temperaturePhaseChangeTwoPhaseMixtures` as `temperaturePhaseChangeThreePhaseMixtures` and the following files
+ modify `options` and `files` in `Make`
+ modify files in `temperaturePhaseChangeThreePhaseMixture`
+ moditf files in `thermoIncompressibleThreePhaseMixture`
+ modify files in `threePhaseMixtureEThermo`

## thermoIncompressibleThreePhaseMixture

+ change from 'incompressibleTwoPhaseMixture.H' to `immiscibleIncompressibleThreePhaseMixture.H`
+ replace 'incompressibleTwoPhaseMixture' with `immiscibleIncompressibleThreePhaseMixture`

### thermoIncompressibleThreePhaseMixture.H

+ add variables and member functions for phase 3

```cpp
dimensionedScalar kappa3_;
dimensionedScalar Cp3_;
dimensionedScalar Cv3_;
dimensionedScalar Hf3_;
```

```cpp
//- Return const-access to phase3 kappa
const dimensionedScalar& kappa3() const
{
    return kappa3_;
};

//- Return const-access to phase3 Cp
const dimensionedScalar& Cp3() const
{
    return Cp3_;
};

//- Return const-access to phase3 Cv
const dimensionedScalar& Cv3() const
{
    return Cv3_;
};

//- Return latent heat for phase 3
const dimensionedScalar& Hf3() const
{
    return Hf3_;
};
```

### thermoIncompressibleThreePhaseMixture.C

In the defination of mixture

```cpp
kappa3_
(
    "kappa3",
    kappa1_.dimensions(),
    subDict(phase3Name_),
    "kappa"
),

Cp3_
(
    "Cp3",
    dimEnergy/dimTemperature/dimMass,
    subDict(phase3Name_),
    "Cp"
),

Cv3_
(
    "Cv3",
    dimEnergy/dimTemperature/dimMass,
    subDict(phase3Name_),
    "Cv"
),

Hf3_
(
    "Hf3",
    dimEnergy/dimMass,
    subDict(phase3Name_),
    "hf"
)
```

In the defination of member function

```cpp
subDict(phase2Name_).readEntry("kappa", kappa3_);
subDict(phase2Name_).readEntry("Cp", Cp3_);
subDict(phase2Name_).readEntry("Cv", Cv3_);
subDict(phase2Name_).readEntry("hf", Hf3_);
```

## threePhaseMixtureEThermo

### threePhaseMixtureEThermo.H

There is only declaration, so do not have to change.

### threePhaseMixtureEThermo.C

#### init()

$$
e_1 = C_{v1} ( T - T_{sat}) + H_{f1}
$$

$$
e_2 = C_{v2} ( T - T_{sat}) + H_{f2}
$$

$$
e_3 = C_{v3} ( T - T_{sat}) + H_{f3}
$$

$$
e = \frac{\alpha_1 \rho_1 e_1 + \alpha_2 \rho_2 e_2 + \alpha_3 \rho_3 e_3}{\alpha_1 \rho_1 + \alpha_2 \rho_2 + \alpha_3 \rho_3}
$$

$$
e = \frac{(\alpha_1 \rho_1 C_{v1} + \alpha_2 \rho_2 C_{v2} + \alpha_3 \rho_3 C_{v3})(T-T_{sat}) + (\alpha_1 \rho_1 H_{f1} + \alpha_2 \rho_2 H_{f2} + \alpha_3 \rho_3 H_{f3})}{\alpha_1 \rho_1 + \alpha_2 \rho_2 + \alpha_3 \rho_3}
$$

```cpp
void Foam::threePhaseMixtureEThermo::init()
{
    // add phase 3
    const volScalarField alpha1Rho1(alpha1()*rho1());
    const volScalarField alpha2Rho2(alpha2()*rho2());
    const volScalarField alpha3Rho3(alpha3()*rho3());

    e_ =
        (
            (T_ - TSat_)*(alpha1Rho1*Cv1() + alpha2Rho2*Cv2() + alpha3Rho3*Cv3())
          + (alpha1Rho1*Hf1() + alpha2Rho2*Hf2() + alpha3Rho3*Hf3())
        )
       /(alpha1Rho1 + alpha2Rho2 + alpha3Rho3);

    e_.correctBoundaryConditions();
}
```

#### correct()

$$
T = \frac{e (\alpha_1 \rho_1 + \alpha_2 \rho_2 + \alpha_3 \rho_3) - (\alpha_1 \rho_1 H_{f1} + \alpha_2 \rho_2 H_{f2} + \alpha_3 \rho_3 H_{f3})}{\alpha_1 \rho_1 C_{v1} + \alpha_2 \rho_2 C_{v2} + \alpha_3 \rho_3 C_{v3}} + T_{sat}
$$

```cpp
void Foam::threePhaseMixtureEThermo::correct() // correct T by alpha1 and alpha2
{
    incompressibleTwoPhaseMixture::correct();

    // add phase 3
    const volScalarField alpha1Rho1(alpha1()*rho1());
    const volScalarField alpha2Rho2(alpha2()*rho2());
    const volScalarField alpha3Rho3(alpha3()*rho3());

    T_ =
        (
            (e_*(alpha1Rho1 + alpha2Rho2 + alpha3Rho3))
         -  (alpha1Rho1*Hf1() + alpha2Rho2*Hf2() + alpha3Rho3*Hf3())
        )
       /(alpha1Rho1*Cv1() + alpha2Rho2*Cv2() + alpha3Rho3*Cv3())
       + TSat_;

    T().correctBoundaryConditions();
}
```

#### he(p, T)

```cpp
Foam::tmp<Foam::volScalarField> Foam::threePhaseMixtureEThermo::he
(
    const volScalarField& p,
    const volScalarField& T
) const
{
    // add phase 3
    const volScalarField alpha1Rho1(alpha1()*rho1());
    const volScalarField alpha2Rho2(alpha2()*rho2());
    const volScalarField alpha3Rho3(alpha3()*rho3());

    return
    (
        (T - TSat_)*(alpha1Rho1*Cv1() + alpha2Rho2*Cv2() + alpha3Rho3*Cv3())
        + (alpha1Rho1*Hf1() + alpha2Rho2*Hf2( + alpha3Rho3*Hf3()))
    )
    / (alpha1Rho1 + alpha2Rho2 + alpha3Rho3);
}
```

#### h(p, T, celli)

```cpp
Foam::tmp<Foam::scalarField> Foam::threePhaseMixtureEThermo::he
(
    const scalarField& p,
    const scalarField& T,
    const labelList& cells
) const
{
    tmp<scalarField> the(new scalarField(T.size()));
    scalarField& he = the.ref();

    // add phase 3
    const volScalarField alpha1Rho1(alpha1()*rho1());
    const volScalarField alpha2Rho2(alpha2()*rho2());
    const volScalarField alpha3Rho3(alpha3()*rho3());

    forAll(T, i)
    {
        const label celli = cells[i];
        he[i] =
            (
                (T[i] - TSat_.value())
               *(
                   alpha1Rho1[celli]*Cv1().value()
                 + alpha2Rho2[celli]*Cv2().value()
                 + alpha3Rho3[celli]*Cv3().value()
                )
              + (
                    alpha1Rho1[celli]*Hf1().value()
                  + alpha2Rho2[celli]*Hf2().value()
                  + alpha3Rho3[celli]*Hf3().value()
                )
            )
            / (alpha1Rho1[celli] + alpha2Rho2[celli] + alpha3Rho3[celli]);
    }

    return the;
}
```

#### h(p, T, patchi)

```cpp
Foam::tmp<Foam::scalarField> Foam::threePhaseMixtureEThermo::he
(
    const scalarField& p,
    const scalarField& T,
    const label patchi
) const
{
    // add phase 3
    const scalarField& alpha1p = alpha1().boundaryField()[patchi];
    const scalarField& alpha2p = alpha2().boundaryField()[patchi];
    const scalarField& alpha3p = alpha3().boundaryField()[patchi];

    const scalarField& Tp = T_.boundaryField()[patchi];

    return
    (
        (
            (Tp - TSat_.value())
           *(
               alpha1p*rho1().value()*Cv1().value()
             + alpha2p*rho2().value()*Cv2().value()
             + alpha3p*rho3().value()*Cv3().value()
            )
          + (
               alpha1p*rho1().value()*Hf1().value()
             + alpha2p*rho2().value()*Hf2().value()
             + alpha3p*rho3().value()*Hf3().value()
            )
        )
        / (alpha1p*rho1().value() + alpha2p*rho2().value() + alpha3p*rho3().value())
    );
}
```

#### Cp()

mass average:

$$
C_p = \frac{\alpha_1 \rho_1 C_{p1} + \alpha_2 \rho_2 C_{p2} + \alpha_3 \rho_3 C_{p3}}{\alpha_1 \rho_1 + \alpha_2 \rho_2 + \alpha_3 \rho_3}
$$

```cpp
Foam::tmp<Foam::volScalarField> Foam::threePhaseMixtureEThermo::Cp() const
{
    // use mass average
    // add phase 3
    const volScalarField alpha1Rho1(alpha1()*rho1());
    const volScalarField alpha2Rho2(alpha2()*rho2());
    const volScalarField alpha3Rho3(alpha3()*rho3());

    return tmp<volScalarField>
    (
        new volScalarField
        (
            "cp",
            (
                (alpha1Rho1*Cp1() + alpha2Rho2*Cp2() + alpha3Rho3*Cp3()) 
                / (alpha1Rho1 + alpha2Rho2 + alpha3Rho3)
            )
        )
    );
}
```

#### Cp(p, T, patchi)

```cpp
Foam::tmp<Foam::scalarField> Foam::threePhaseMixtureEThermo::Cp
(
    const scalarField& p,
    const scalarField& T,
    const label patchi
) const
{
    // use mass average
    // add phase 3    
    const scalarField& alpha1p = alpha1().boundaryField()[patchi];
    const scalarField& alpha2p = alpha2().boundaryField()[patchi];
    const scalarField& alpha3p = alpha3().boundaryField()[patchi];

    return
    (
        (
            alpha1p*rho1().value()*Cp1().value()
          + alpha2p*rho2().value()*Cp2().value()
          + alpha3p*rho3().value()*Cp3().value()
        )
        /  (alpha1p*rho1().value() + alpha2p*rho2().value() + alpha3p*rho3().value())
    );
}
```

#### rho()

add limitation to $alpha_2$ as:

$$
0 \leqslant \alpha_2 \leqslant 1 - \alpha_1
$$

```cpp
Foam::tmp<Foam::volScalarField> Foam::threePhaseMixtureEThermo::rho() const
{
    const volScalarField limitedAlpha1
    (
        min(max(alpha1_, scalar(0)), scalar(1))
    );
    // add limitation to alpha2
    const volScalarField limitedAlpha2
    (
        min(max(alpha2_, scalar(0)), scalar(1) - limitedAlpha1)
    );

    return tmp<volScalarField>
    (
        new volScalarField
        (
            "rho",
            (
                limitedAlpha1*rho1().value()
                limitedAlpha2*rho2().value()
              + (scalar(1) - limitedAlpha1 - limitedAlpha3)*rho3().value()
            )
        )
    );
}
```

#### rho(patchi)

```cpp
Foam::tmp<Foam::scalarField> Foam::threePhaseMixtureEThermo::rho
(
    const label patchi
) const
{
    const volScalarField limitedAlpha1
    (
        min(max(alpha1_, scalar(0)), scalar(1))
    );
    // add limitation to alpha2
    const volScalarField limitedAlpha2
    (
        min(max(alpha2_, scalar(0)), scalar(1) - limitedAlpha1)
    );

    const scalarField& alpha1p = limitedAlpha1.boundaryField()[patchi];
    const scalarField& alpha2p = limitedAlpha2.boundaryField()[patchi];

    return
    (
        alpha1p*rho1().value() + alpha2p*rho2().value()
      + (scalar(1) - alpha1p - alpha2p)*rho3().value()
    );
}
```

#### Cv()

mass average like Cp()

```cpp
Foam::tmp<Foam::volScalarField> Foam::threePhaseMixtureEThermo::Cv() const
{
    // use mass average
    // add phase 3
    const volScalarField alpha1Rho1(alpha1()*rho1());
    const volScalarField alpha2Rho2(alpha2()*rho2());
    const volScalarField alpha3Rho3(alpha3()*rho3());

    return tmp<volScalarField>
    (
        new volScalarField
        (
            "cv",
            (
                (alpha1Rho1*Cv1() + alpha2Rho2*Cv2() + alpha3Rho3*Cv3()) 
                / (alpha1Rho1 + alpha2Rho2 + alpha3Rho3)
            )
        )
    );
}
```

#### Cv(p, T, patchi)

```cpp
Foam::tmp<Foam::scalarField> Foam::threePhaseMixtureEThermo::Cv
(
    const scalarField& p,
    const scalarField& T,
    const label patchi
) const
{
    // use mass average
    // add phase 3    
    const scalarField& alpha1p = alpha1().boundaryField()[patchi];
    const scalarField& alpha2p = alpha2().boundaryField()[patchi];
    const scalarField& alpha3p = alpha3().boundaryField()[patchi];

    return
    (
        (
            alpha1p*rho1().value()*Cv1().value()
          + alpha2p*rho2().value()*Cv2().value()
          + alpha3p*rho3().value()*Cv3().value()
        )
        /  (alpha1p*rho1().value() + alpha2p*rho2().value() + alpha3p*rho3().value())
    );
}
```

#### gamma()

```cpp
Foam::tmp<Foam::volScalarField> Foam::threePhaseMixtureEThermo::gamma() const
{
    // add phase 3
    return tmp<volScalarField>
    (
        (alpha1_*Cp1() + alpha2_*Cp2() + alpha3_*Cp3())
      / (alpha1_*Cv1() + alpha2_*Cv2() + alpha3_*Cv3())
    );
}
```

#### kappa()

```cpp
Foam::tmp<Foam::volScalarField> Foam::threePhaseMixtureEThermo::kappa() const
{
    const volScalarField limitedAlpha1
    (
        min(max(alpha1_, scalar(0)), scalar(1))
    );
    // add limitation to alpha2
    const volScalarField limitedAlpha2
    (
        min(max(alpha2_, scalar(0)), scalar(1) - limitedAlpha1)
    );

    return tmp<volScalarField>
    (
        new volScalarField
        (
            "kappa",
            (
                limitedAlpha1*kappa1() + limitedAlpha2*kappa2() 
              + (scalar(1) - limitedAlpha1  - limitedAlpha2)*kappa3()
            )
        )
    );
}
```

#### kappa(patchi)

```cpp
Foam::tmp<Foam::scalarField> Foam::threePhaseMixtureEThermo::kappa
(
    const label patchi
) const
{
    const volScalarField limitedAlpha1
    (
        min(max(alpha1_, scalar(0)), scalar(1))
    );
    // add limitation to alpha2
    const volScalarField limitedAlpha2
    (
        min(max(alpha2_, scalar(0)), scalar(1) - limitedAlpha1)
    );

    const scalarField& alpha1p = limitedAlpha1.boundaryField()[patchi];
    const scalarField& alpha2p = limitedAlpha2.boundaryField()[patchi];

    return 
    (
        alpha1p*kappa1().value() + alpha2p*kappa2().value() 
      + (1 - alpha1p - alpha2p)*kappa3().value()
    );
}
```

### kappaEff(kappat, patchi)

```cpp
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
    // add limitation to alpha2
    const volScalarField limitedAlpha2
    (
        min(max(alpha2_, scalar(0)), scalar(1) - limitedAlpha1)
    );

    const scalarField& alpha1p = limitedAlpha1.boundaryField()[patchi];
    const scalarField& alpha2p = limitedAlpha2.boundaryField()[patchi];

    return
    (    
        (alpha1p*kappa1().value() + alpha2p*kappa2().value() 
      + (1 - alpha1p - alpha2p)*kappa3().value()) + kappat
    );

}
```

#### alphaEff(alphat)

```cpp
Foam::tmp<Foam::volScalarField> Foam::threePhaseMixtureEThermo::alphaEff
(
    const volScalarField& alphat
) const
{
    // add phase 3
    const volScalarField rho
    (
        alpha1_*rho1() + alpha2_*rho2() + (1.0 - alpha1_ - alpha2_)*rho3()
    );

    return (kappa()/Cp()/rho + alphat);
}
```

#### alphaEff(alphat, patchi)

```cpp
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
    // add limitation to alpha2
    const volScalarField limitedAlpha2
    (
        min(max(alpha2_, scalar(0)), scalar(1) - limitedAlpha1)
    );

    const scalarField& alpha1p = limitedAlpha1.boundaryField()[patchi];
    const scalarField& alpha2p = limitedAlpha2.boundaryField()[patchi];

    const scalarField rho
    (
        alpha1p*rho1().value() + alpha2p*rho2().value() + (1.0 - alpha1p - alpha2p)*rho3().value()
    );

    const scalarField kappa
    (
        alpha1p*kappa1().value() + alpha2p*kappa2().value() + (1.0 - alpha1p - alpha2p)*kappa3().value()
    );

    const scalarField Cp
    (
        alpha1p*Cp1().value() + alpha2p*Cp2().value() + (1.0 - alpha1p - alpha2p)*Cp3().value()
    );

    return kappa/Cp/rho + alphat;
}
```


## interMixingEvapFoam.C

+ replace all `interMixingFoam` with `interMixingEvapFoam`

# Theory

## Definition

The interMixingFoam solver traces three fluids 1, 2 and 3, whereof the second and third are miscible. For droplet evaporation, liquid and gas phases need to be traced, the gas phase consisting of both bulk gas (air) and vapor. **Thus, phase 1 is considered to be liquid, phase 2 vapor and phase 3 bulk gas.**

## Tracing variable

$$
\alpha_i = \frac{V_i}{V}
$$

where $V_i$ is traced volume, $V$ is cell volume

## Evaporation Rate

In the following, physical properties of the different fluids are denoted with subscripts $l$ for liquid, $g$ for gas and $v$ for vapor, respectively.

The mass flow rate $\dot m$ due to evaporation is given as

$$
\dot{m} = A_\Gamma \frac{D_{vg}  \rho_{gp}}{1 - X_v}\langle\nabla X_v, \hat {\mathbf{n}}_\Gamma \rangle
$$

where $A_\Gamma$ is interface area, $D_{vg}$ is diffusion coefficient of vapor in bulk gas, $\rho_{gp}$ is gas phase density.

$X_v$ is the vaper mass fraction and can be calculated as:

$$
X_v = \frac{\alpha_2}{1-\alpha_1 + \epsilon} \frac{\rho_v}{\rho_{gp}}
$$

where, $\alpha_1$ is $\alpha_l$ (for liquid), $\alpha_2$ is $\alpha_v$ (for vapor), $1 - \alpha_1 = \alpha_2 + \alpha_3 = \alpha_v + \alpha_g$, $\epsilon$ is a very small value (VSMALL of OpenFOAM, avoiding denominator is zero)

$$
\rho_{gp} = \frac{\rho_g(1 - \alpha_1 - \alpha_2) + \rho_v \alpha_2 + 0.5 \epsilon (\rho_v + \rho_g)}{1 - \alpha_1 + \epsilon}
$$

$$
\hat \mathbf{n} _\Gamma = \frac{\nabla \alpha_l}{|\alpha_l| + \epsilon_\Gamma}
$$

Express the evaporation rate as a volumetric mass source:

$$
\dot m ''' = \frac{\dot m}{V}
$$

### Interface Area Estimation

Variables denoted with subscript $f$ are calculated on cell faces, others in the cell center.

The interface area on cell faces $A_{\Gamma, f}$ is computed as

$$
A_{\Gamma, f} = \Delta \alpha_f \langle \mathbf{S}_f, \hat \mathbf n _\Gamma \rangle
$$

where $\Delta \alpha_f = \alpha_{1, nei} - \alpha_{1, own}$ is the $\alpha$ difference between neighbor and owner cell of the face, $\mathbf S_f$ is normal face area, $\hat \mathbf n _ \Gamma$ is interface normal.

This interface area then needs to be mapped back to the cell centers of both cells adjacent to the face, using the following interpolation:

$$
A_{\Gamma, nei} = (1-\omega) |A_{\Gamma, f}|
$$

$$
A_{\Gamma, own} = \omega |A_{\Gamma, f}|
$$

where $\omega$ is weight, calculated as:

$$
\omega = 0.5 (1 - erf(\log(r)))
$$

$$
r = \frac{|\alpha_{1, nei} - 0.5|}{|\alpha_{1, own} - 0.5| + \epsilon}
$$

### Vapor Concentration Gradient Calculation

$$
\nabla X_v = \nabla(\alpha_1 X_{v, sat} + (1-\alpha_1)X_v)
$$

where $X_{v, sat}$ is the mass fraction of saturated vapor, calculated as:

$$
X_{v, sat} = x_{v, sat} \frac{M_v}{x_{v, sat}M_v+(1-x_{v, sat})M_g}
$$

with molecular weight of vapor $M_v$, molecular weight of gas $M_g$. The mole fraction of saturated vapor is calculated as

$$
x_{v, sat} = \frac{p_{v, sat} (T)}{pstat}
$$

where $p_{v, sat}$ is the saturation pressure of vapor at each cell's temperature $T$, $p_{stat}$ is the static pressure of the fluid.

$$
\nabla X_v = \nabla (\alpha_1 X_{v, sat} + (1 - \alpha_1) X_v)
$$

## Energy Equation

A convective as well as conductive term is required. Moreover, a negative source term due to the enthalpy of vaporization needs to be added:

$$
\frac{\partial}{\partial t} (\rho c_p T) + \nabla \cdot (\rho c_p \mathbf{U} T) = \nabla \cdot (\lambda \nabla T) - \Delta h_v \dot m '''
$$

Local fluid properties such as density $\rho$, viscosity $\nu$ and heat conductivity $\lambda$ are calculated using the one-field approach (only one property for all phases) and volume averaged, e.g.:

$$
\rho = \sum \rho_i \alpha_i
$$

Only the specific heat transfer capacity $c_p$ is mass averaged:

$$
c_p = \frac{\sum_i \rho_i \alpha_i c_{p, i}}{\sum_i \rho_i \alpha_i}
$$

So far, all material properties are assumed to be constant. Temperature-dependent material properties could be added in the future for improved results.

## Modification of Governing Equations

### Continuity Equation

$$
\nabla \mathbf{U} = -\dot m ''' (\frac{1}{\rho_v} - \frac{1}{\rho_l})
$$

### Volume-of-Fluid Equations

$$
\frac{\partial \alpha_l}{\partial t} = \nabla \cdot (\mathbf {U} \alpha_l) = - \frac{\dot m '''}{\rho_l}
$$

$$
\frac{\partial \alpha_v}{\partial t} = \nabla \cdot (\mathbf {U} \alpha_v) = \frac{\dot m '''}{\rho_v}
$$

# compressibleInterFoam

## Structure

+ compressibleInterDyMFoam
+ compressibleInterFilmFoam
+ Make
+ surfaceTensionModels
+ twoPhaseMixtureThermo
+ VoFphaseCompressibleTurbulenceModels
+ Allwclean
+ Allwmake
+ alphaSuSp.H
+ compressibleAlphaEqnSubCycle..H
+ compressibleInterFoam.C
+ createFields.H
+ pEqn.H
+ rhofs.H
+ TEqn.H
+ UEqn.H

## twoPhaseMixtureThermo

### twoPhaseMixtureThermo.H

