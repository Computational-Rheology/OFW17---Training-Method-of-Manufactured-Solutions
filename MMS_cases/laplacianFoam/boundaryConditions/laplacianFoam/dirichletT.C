/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2020 OpenCFD Ltd.
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


    Note:
        File adapted by:
            Computational Rheology Group at the Institute for Polymers and Composites, University of Minho, Portugal, 2022
        If you find any problems/mistakes or have any improvement recommendations, please contact: info@crheo.org

\*---------------------------------------------------------------------------*/

#include "dirichletT.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::dirichletT::dirichletT
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(p, iF)
{}


Foam::dirichletT::dirichletT
(
    const dirichletT& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchScalarField(ptf, p, iF, mapper)
{}


Foam::dirichletT::dirichletT
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchScalarField(p, iF, dict, false)
{
    if (dict.found("value"))
    {
        fvPatchField<scalar>::operator=
        (
            scalarField("value", dict, p.size())
        );
    }
    else
    {
        //evaluate(Pstream::commsTypes::blocking);
    }
}


Foam::dirichletT::dirichletT
(
    const dirichletT& tppsf
)
:
    fixedValueFvPatchScalarField(tppsf)
{}


Foam::dirichletT::dirichletT
(
    const dirichletT& tppsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(tppsf, iF)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
/*
void Foam::dirichletT::autoMap
(
    const fvPatchFieldMapper& m
)
{

    FatalErrorInFunction
    << "autoMap not implemented"<< nl
    << exit(FatalError);

}


void Foam::dirichletT::rmap
(
    const fvPatchScalarField& ptf,
    const labelList& addr
)
{
    FatalErrorInFunction
    << "rmap not implemented"<< nl
    << exit(FatalError);
}*/


void Foam::dirichletT::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    // Get reference to patch face centers
    const vectorField& Cf = patch().Cf();

    scalarField result(Cf.size(), Zero);

    forAll(result, faceI)
    {
        const scalar x = Cf[faceI].x(); // Get the x-coordinate of the face
        const scalar y = Cf[faceI].y(); // Get the y-coordinate of the face
        
        result[faceI] = 225.0 + 150*Foam::cos(x*x + y*y); // Apply T(x,y) from MMS
    }

    operator==(result);

    fixedValueFvPatchScalarField::updateCoeffs();
}


void Foam::dirichletT::write(Ostream& os) const
{
    fvPatchScalarField::write(os);
    writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchScalarField,
        dirichletT
    );
}

// ************************************************************************* //
