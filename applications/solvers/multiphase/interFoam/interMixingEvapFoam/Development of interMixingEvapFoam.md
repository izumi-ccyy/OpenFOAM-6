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

