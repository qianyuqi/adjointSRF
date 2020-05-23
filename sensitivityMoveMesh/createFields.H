    Info<< "Reading field sensitivity\n" << endl;

    volScalarField sensitivity
    (
        IOobject
        (
            "sensitivity",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::NO_WRITE
        ),
        mesh
    );
    
    Info<< "Reading transportProperties\n" << endl;

    IOdictionary transportProperties
    (
        IOobject
        (
            "transportProperties",
            runTime.constant(),
            mesh,
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::NO_WRITE
        )
    );

    
    Info<< "Reading coefficient Alpha\n" << endl;

    dimensionedScalar alpha
    (
        transportProperties.lookup("alpha")
    );
