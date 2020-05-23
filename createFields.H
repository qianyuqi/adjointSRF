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

    Info<< "Reading field Urel\n" << endl;
    volVectorField Urel
    (
        IOobject
        (
            "Urel",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

    Info<< "Reading field Uabs\n" << endl;
    volVectorField Uabs
    (
        IOobject
        (
            "Uabs",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

    Info<< "Reading field omega\n" << endl;
    volScalarField omega
    (
        IOobject
        (
            "omega",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

    Info<< "Reading/calculating face flux field phi\n" << endl;
    surfaceScalarField phi
    (
        IOobject
        (
            "phi",
            runTime.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        linearInterpolate(Urel) & mesh.Sf()
    );


    label pRefCell = 0;
    scalar pRefValue = 0.0;
    setRefCell(p, mesh.solutionDict().subDict("SIMPLE"), pRefCell, pRefValue);


    Info<< "Reading field pa\n" << endl;
    volScalarField pa
    (
        IOobject
        (
            "pa",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

    Info<< "Reading field Ua\n" << endl;
    volVectorField Ua
    (
        IOobject
        (
            "Ua",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

    #include "createPhia.H"


    label paRefCell = 0;
    scalar paRefValue = 0.0;
    setRefCell
    (
        pa,
        mesh.solutionDict().subDict("SIMPLE"),
        paRefCell,
        paRefValue
    );


    singlePhaseTransportModel laminarTransport(Urel, phi);

    autoPtr<incompressible::RASModel> turbulence
    (
        incompressible::RASModel::New(Urel, phi, laminarTransport)
    );
    Info<< "Creating SRF model\n" << endl;
    autoPtr<SRF::SRFModel> SRF(SRF::SRFModel::New(Ua));


    dimensionedScalar zeroSensitivity("0", dimVelocity*dimVelocity, 0.0);
    dimensionedScalar zeroAlpha("0", dimless/dimTime, 0.0);

/*    dimensionedScalar lambda(laminarTransport.lookup("lambda"));
    dimensionedScalar alphaMax(laminarTransport.lookup("alphaMax"));

    const labelList& inletCells = mesh.boundary()["inlet"].faceCells();
    //const labelList& outletCells = mesh.boundary()["outlet"].faceCells();

    volScalarField alpha
    (
        IOobject
        (
            "alpha",
            runTime.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        lambda*max(Ua & U, zeroSensitivity)
    );
    zeroCells(alpha, inletCells);
    //zeroCells(alpha, outletCells);
*/
	dimensionedScalar nu(laminarTransport.lookup("nu"));
	dimensionedScalar gamma(laminarTransport.lookup("gamma"));

    const labelList& inletCells = mesh.boundary()["inlet"].faceCells();
    //const labelList& outletCells = mesh.boundary()["outlet"].faceCells();

	volScalarField sens
	(
		IOobject
		(
			"sensitivity",
			runTime.timeName(),
			mesh,
			IOobject::READ_IF_PRESENT,
			IOobject::AUTO_WRITE
		),
		(Ua&Urel)*0
	);

