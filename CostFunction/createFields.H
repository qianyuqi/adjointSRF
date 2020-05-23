    Info<< "Reading field p\n" << endl;
    volScalarField p
    (
        IOobject
        (
            "p",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

    Info<< "Reading field U\n" << endl;
    volVectorField U
    (
        IOobject
        (
            "U",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

    #include "createPhi.H"


    label pRefCell = 0;
    scalar pRefValue = 0.0;
    setRefCell(p, mesh.solutionDict().subDict("SIMPLE"), pRefCell, pRefValue);

dictionary optFunc = mesh.solutionDict().subDict("objectiveFunctionDict");

int nObjPatch =
    optFunc.lookupOrDefault<scalar>("numberObjectivePatches",0);

int objFunction =
    optFunc.lookupOrDefault<scalar>("objectiveFunction", 0);

wordList objPatchNames = optFunc.lookup("objectivePatchesNames");

Info<<"Initializing objective function calculation:"<<endl;
Info<<"The objective function chosen is:"<<objFunction<<endl;
Info<<"Name of the patches for which the cost function will be calculated:"<<objPatchNames<<endl;
Info<<"Number of patches:"<<nObjPatch<<endl;
if (nObjPatch==0) {
Info<<"^Warning:: no objective patches defined, default 0 used."<<endl;
}
if (objFunction==0) {
     Info<<"^Warning:: no objective function defined, default 0 used."<<endl;
}
Info<<"Objective function power dissipated chosen."<<endl;
label objPatchList[nObjPatch];
int iLoop;
for (iLoop=0; iLoop<nObjPatch; iLoop++)
{
    objPatchList[iLoop] =
         mesh.boundaryMesh().findPatchID(objPatchNames[iLoop]);
}
vector Udlookup = optFunc.lookupOrDefault<vector>("Uduserdefnodim",vector(1,0,0));



