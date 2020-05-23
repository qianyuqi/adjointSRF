/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
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

#include "adjointOutletVelocityFvPatchVectorField.H"
#include "volFields.H"
#include "addToRunTimeSelectionTable.H"
#include "surfaceFields.H"
#include "fvPatchFieldMapper.H"
#include "RASModel.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::adjointOutletVelocityFvPatchVectorField::
adjointOutletVelocityFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchVectorField(p, iF)
{}


Foam::adjointOutletVelocityFvPatchVectorField::
adjointOutletVelocityFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchVectorField(p, iF)
{
    fvPatchVectorField::operator=(vectorField("value", dict, p.size()));
}


Foam::adjointOutletVelocityFvPatchVectorField::
adjointOutletVelocityFvPatchVectorField
(
    const adjointOutletVelocityFvPatchVectorField& ptf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchVectorField(ptf, p, iF, mapper)
{}


Foam::adjointOutletVelocityFvPatchVectorField::
adjointOutletVelocityFvPatchVectorField
(
    const adjointOutletVelocityFvPatchVectorField& pivpvf,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchVectorField(pivpvf, iF)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// Update the coefficients associated with the patch field
void Foam::adjointOutletVelocityFvPatchVectorField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    const fvsPatchField<scalar>& phiap =
        patch().lookupPatchField<surfaceScalarField, scalar>("phia");

    const fvPatchField<vector>& Up =
        patch().lookupPatchField<volVectorField, vector>("Urel");

    const fvPatchField<vector>& Upabs =
        patch().lookupPatchField<volVectorField, vector>("Uabs");

    const fvPatchField<vector>& Uap =
    	patch().lookupPatchField<volVectorField, vector>("Ua");

    const fvsPatchField<scalar>& phip =
    	patch().lookupPatchField<surfaceScalarField, scalar>("phi");

    const incompressible::RASModel& rasModel =
    	db().lookupObject<incompressible::RASModel>("RASProperties");

    scalarField nueff = rasModel.nuEff()().boundaryField()[patch().index()];

    const scalarField& deltainv = 
    	patch().deltaCoeffs(); // dist^(-1) 

//Primal velocity, mag of normal component and tangential component
    scalarField Up_ns = phip/patch().magSf();

    vectorField Upabs_t = Upabs - (phip * patch().Sf())/(patch().magSf()*patch().magSf());

//Tangential component of adjoint velocity in neighbouring node
    vectorField Uaneigh = Uap.patchInternalField();
    vectorField Uaneigh_n = (Uaneigh & patch().nf())*patch().nf();
//vectorField Uaneigh_n = (Uaneigh & normal)*normal;
//vectorField Uaneigh_n = (Uaneigh & patch().nf());
    vectorField Uaneigh_t = Uaneigh - Uaneigh_n;

    vectorField Uap_t =  ((Up_ns*Upabs_t) + nueff*deltainv*Uaneigh_t) / (Up_ns+nueff*deltainv) ;

    vectorField Uap_n = (phiap * patch().Sf())/(patch().magSf()*patch().magSf());

    operator==(Uap_t+Uap_n);

    fixedValueFvPatchVectorField::updateCoeffs();
}


void Foam::adjointOutletVelocityFvPatchVectorField::write(Ostream& os) const
{
    fvPatchVectorField::write(os);
    writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchVectorField,
        adjointOutletVelocityFvPatchVectorField
    );
}


// ************************************************************************* //
