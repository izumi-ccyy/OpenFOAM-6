# TDAC

## chemistryReductionMethod

### chemistryReductionMethod.H

some important variables

```cpp
protected:

    const IOdictionary& dict_;

    //- Dictionary that store the algorithm data
    const dictionary coeffsDict_;

    //- Is mechanism reduction active?
    Switch active_;

    //- Switch to select performance logging
    Switch log_;

    TDACChemistryModel<CompType, ThermoType>& chemistry_;

    //- List of active species (active = true)
    List<bool> activeSpecies_;

    //- Number of active species
    label NsSimp_;

    //- Number of species
    const label nSpecie_;

    //- Tolerance for the mechanism reduction algorithm
    scalar tolerance_;
```

some important functions

```cpp
//- Is mechanism reduction active?
inline bool active() const;

//- Is performance data logging enabled?
inline bool log() const;

//- Return the active species
inline const List<bool>& activeSpecies() const;

//- Return the number of active species
inline label NsSimp();

//- Return the initial number of species
inline label nSpecie();

//- Return the tolerance
inline scalar tolerance() const;

//- Reduce the mechanism
virtual void reduceMechanism
(
    const scalarField &c,
    const scalar T,
    const scalar p
) = 0;
```

## DRG

### DRG.H

```cpp
//- List of label for the search initiating set
    labelList searchInitSet_;
```

### DRG.C

#### constructor

```cpp
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
```

#### DynamicList

`src/OpenFOAM/containers/Lists/DynamicList/DynamicList.H`

```cpp
 Description
     A 1D vector of objects of type <T> that resizes itself as necessary to
     accept the new objects.
 
     Internal storage is a compact array and the list can be shrunk to compact
     storage. The increase of list size is controlled by three template
     parameters, which allows the list storage to either increase by the given
     increment or by the given multiplier and divider (allowing non-integer
     multiples).
```