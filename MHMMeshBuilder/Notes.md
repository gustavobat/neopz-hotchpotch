# ```MHMixedMeshControl``` today

## Constructor
- Enters ```MHMeshControl``` constructor. Here, ```fGMesh``` is initialized with ```gmesh```
passed as a parameter. The variables ```fGeoToMHMDomain```, ```fMHMtoSubCMesh```,
 ```fGlobalSystemSize```, ```fGlobalSystemWithLocalCondensationSize``` and
```fNumeq``` are initialized with zero/default values.
Inside the function we have:
```
    fCMesh = new TPZCompMesh(fGMesh);
    fCMesh->SetDimModel(fGMesh->Dimension());
    fPressureFineMesh = fCMesh;
```
It goes back to the ```MHMixedMeshControl``` constructor, that fills the following variables:
```
    fFluxMesh = new TPZCompMesh(gmesh);
    fFluxMesh->SetDimModel(gmesh->Dimension());
    fCMesh = new TPZMultiphysicsCompMesh(gmesh);
    fCMesh->SetDimModel(gmesh->Dimension());
```

## Coarse indexes
The coarse indexes need to be computed outside of the class and then passed to it. 
In the example I've seen, the mesh is refined uniformly from 4 original squares elements. 
These 4 elements are the 'coarse indexes'. The vector that stores is has 4 elements: 0, 1, 2, 3
that refer to these original elements.
The mesh has a total of 140 elements.

## The ```DefinePartitionByCoarseIndices``` method
Mm... ok, now things get real.

The previously computed coarseindices vector is the only parameter that this functions takes.
A variable that contains the number of the coarse indices is created. In this case it is 4.
A variable that holds the number of elements in the ```fGmesh``` is also created.
```fGeoToMHMDomain``` is resized to the number of elements in ```fGmesh``` and filled with -1.
 The first 4 elements of ```fGeoToMHMDomain``` are filled with 0, 1, 2, 3. The variable
 ```fMHMToSubCMesh``` is a ```std::map<int64_t, int64_t>``` and is filled with [0, -1] ... [3, -1].
 Then for every other element, ```fGeoToMHMDomain[el]``` is assigned as the ```fGeoToMHMDomain``` 
 value of its eldest ancestor.
 
 The function calls the ```CreateSkeletonElements()``` method.
 
## The ```CreateSkeletonElements()``` method
It creates face elements on the interfaces between coarse element indexes.
The variable ```fInterfaces``` needs to be empty (it's checked at the beginning ot the method) and
```fSkeletonMatId``` must be greater or equal than 0 (I don't know why it is necessary).
A temporary computational mesh is created by the method ```CriaMalhaTemporaria()```, which does
the following:

- Resets ```fGMesh``` reference;
- Iterates through mesh elements, filtering 2D elements;
- Inserts element mat id into set variable;
- Filters (d-1)-dimensional sides
- Gets side neighbour, if ```neighbour.Element()``` is (d-1)-dimensional, inserts its material id
into a bcIds set, ends first for loop
- Creates a dummy compmesh, inserts dummy materials for the IDs contained in both sets,
sets all create functions discontinuous
- Creates CompEl for each element (0, 1, 2, 3) of ```fMHMtoSubCMesh```, checking its a volume 
element. The compel is created and if it has a boundary condition, a (d-1)-dimensional comp 
element is created
- Calls ```TPZCreateApproximationSpace::CreateInterfaces(*cmesh)``` and 
```cmesh->ExpandSolution()```, this creates interface GeoElements
- Returns cmesh to ```CreateSkeletonElements()``` method

An iteration occurs over mesh elements.



# The new MHMMeshBuilder class

