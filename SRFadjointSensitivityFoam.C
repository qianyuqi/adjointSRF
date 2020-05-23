/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2012 OpenFOAM Foundation
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

Application
    ajointShapeOptimizationFoam

Description
    Steady-state solver for incompressible, turbulent flow of non-Newtonian
    fluids with optimisation of duct shape by applying "blockage" in regions
    causing pressure loss as estimated using an adjoint formulation.

    References:
    \verbatim
        "Implementation of a continuous adjoint for topology optimization of
         ducted flows"
        C. Othmer,
        E. de Villiers,
        H.G. Weller
        AIAA-2007-3947
        http://pdf.aiaa.org/preview/CDReadyMCFD07_1379/PV2007_3947.pdf
    \endverbatim

    Note that this solver optimises for total pressure loss whereas the
    above paper describes the method for optimising power-loss.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "singlePhaseTransportModel.H"
#include "RASModel.H"
#include "SRFModel.H"
#include "simpleControl.H"

template<class Type>
void zeroCells
(
    GeometricField<Type, fvPatchField, volMesh>& vf,
    const labelList& cells
)
{
    forAll(cells, i)
    {
        vf[cells[i]] = pTraits<Type>::zero;
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "setRootCase.H"

    #include "createTime.H"
    #include "createMesh.H"
    #include "createFields.H"
    #include "initContinuityErrs.H"
    #include "initAdjointContinuityErrs.H"

    simpleControl simple(mesh);

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nStarting time loop\n" << endl;

    while (simple.loop())
    {
        Info<< "Time = " << runTime.timeName() << nl << endl;

/*        laminarTransport.lookup("lambda") >> lambda;

            // Explicitly relax pressure for momentum corrector
            p.relax();

            // Momentum corrector
            Urel = Urel - rAU*fvc::grad(p);
            Urel.correctBoundaryConditions();
        }
*/
        // Adjoint Pressure-velocity SIMPLE corrector
        {
            // Adjoint Momentum predictor


            volVectorField adjointTransposeConvection((fvc::grad(Ua) & Urel));
            //volVectorField adjointTransposeConvection
            //(
            //    fvc::reconstruct
            //    (
            //        mesh.magSf()*(fvc::snGrad(Ua) & fvc::interpolate(U))
            //    )
            //);

            zeroCells(adjointTransposeConvection, inletCells);

            tmp<fvVectorMatrix> UaEqn
            (
                fvm::div(-phi, Ua)
              - adjointTransposeConvection
              + turbulence->divDevReff(Ua)
	      - SRF->Fcoriolis()
//	      -2.0*(omega ^ Ua) 
              
            );

            UaEqn().relax();

            solve(UaEqn() == -fvc::grad(pa));

            pa.boundaryField().updateCoeffs();
            volScalarField rAUa(1.0/UaEqn().A());
            Ua = rAUa*UaEqn().H();
            UaEqn.clear();
            phia = fvc::interpolate(Ua) & mesh.Sf();
            adjustPhi(phia, Ua, pa);

            // Non-orthogonal pressure corrector loop
            while (simple.correctNonOrthogonal())
            {
                fvScalarMatrix paEqn
                (
                    fvm::laplacian(rAUa, pa) == fvc::div(phia)
                );

                paEqn.setReference(paRefCell, paRefValue);
                paEqn.solve();

                if (simple.finalNonOrthogonalIter())
                {
                    phia -=   paEqn.flux();
                }
            }

            #include "adjointContinuityErrs.H"

            // Explicitly relax pressure for adjoint momentum corrector
            pa.relax();

            // Adjoint momentum corrector
            Ua -= rAUa*fvc::grad(pa);
            Ua.correctBoundaryConditions();
        }

//        turbulence->correct();
	word patchName="innerBox";
	label patchID=mesh.boundaryMesh().findPatchID(patchName);
	error err("Error!\n");
	if(-1==patchID)
	 err.exit();
	const fvPatch& cPatch=mesh.boundary()[patchID];
	const labelUList& cellsPatch=cPatch.faceCells();
	forAll(cellsPatch,cellI)
	 sens[cellsPatch[cellI]]=(Ua[cellsPatch[cellI]]&Urel[cellsPatch[cellI]])*mesh.V()[cellsPatch[cellI]];
	forAll(cPatch,faceI)
	{
	 label faceCellI=cPatch.faceCells()[faceI];
	 sens.boundaryField()[patchID][faceI]=sens[faceCellI];
	}
	scalar maxSens(gMax(sens));
	scalar minSens(gMin(sens));
	maxSens=(maxSens*maxSens>=minSens*minSens ? fabs(maxSens):fabs(minSens));
	sens=sens/(maxSens+SMALL);
        runTime.write();

        Info<< "ExecutionTime = "
            << runTime.elapsedCpuTime()
            << " s\n\n" << endl;
    }

    Info<< "End\n" << endl;

    return(0);
}


// ************************************************************************* //
