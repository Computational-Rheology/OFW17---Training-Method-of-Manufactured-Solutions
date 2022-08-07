/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2019 OpenCFD Ltd.
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
    laplacianFoam

Group
    grpBasicSolvers

Description
    Laplace equation solver for a scalar quantity.

    \heading Solver details
    The solver is applicable to, e.g. for thermal diffusion in a solid.  The
    equation is given by:

    \f[
        \ddt{T}  = \div \left( D_T \grad T \right)
    \f]

    Where:
    \vartable
        T     | Scalar field which is solved for, e.g. temperature
        D_T   | Diffusion coefficient
    \endvartable

    \heading Required fields
    \plaintable
        T     | Scalar field which is solved for, e.g. temperature
    \endplaintable


    Notes:

    File adapted developed by:
          Computational Rheology Group at the Institute for Polymers and Composites, University of Minho, Portugal, 2022
    If you find any problems/mistakes or have any improvement recommendations, please contact: info@crheo.org

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "fvOptions.H"
#include "simpleControl.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::addNote
    (
        "Laplace equation solver for a scalar quantity."
    );

    #include "postProcess.H"

    #include "addCheckCaseOptions.H"
    #include "setRootCaseLists.H"
    #include "createTime.H"
    #include "createMesh.H"

    simpleControl simple(mesh);

    #include "createFields.H"

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nCalculating temperature distribution\n" << endl;

    while (simple.loop())
    {
        Info<< "Time = " << runTime.timeName() << nl << endl;

        while (simple.correctNonOrthogonal())
        {
            // fvScalarMatrix TEqn // Original Formulation
            // (
            //     fvm::ddt(T) - fvm::laplacian(DT, T) 
            //  ==
            //     fvOptions(T)
            // );

            fvScalarMatrix TEqn
            (
                fvm::ddt(T) - fvm::laplacian(DT, T) 
             ==
                fvOptions(T) + source
            );


            fvOptions.constrain(TEqn);
            TEqn.solve();
            fvOptions.correct(T);
        }



        // MMS


            // Calculate error norms
            // Get reference to cell volumes
        const scalarField& V = mesh.V();

        // Compute the difference between analytical and computed values
        err = mag(analyticalSolution - T);

        // Compute error norms
        Info << "L1 norm is: "      << gSum(err * V) / gSum(V)                     << endl;
        Info << "L2 norm is: "      << Foam::sqrt( gSum(err * err * V) / gSum(V) ) << endl;
        Info << "Linf norm is: "    << gMax(err)                                   << endl;
        
        #include "write.H"

        runTime.printExecutionTime(Info);
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
