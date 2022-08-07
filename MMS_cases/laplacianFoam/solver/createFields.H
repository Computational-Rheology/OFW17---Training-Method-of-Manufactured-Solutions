Info<< "Reading field T\n" << endl;

volScalarField T
(
    IOobject
    (
        "T",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);


Info<< "Reading diffusivity DT\n" << endl;

volScalarField DT
(
    IOobject
    (
        "DT",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar(dimViscosity, Zero)
);

if (!DT.headerOk())
{
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
    DT = dimensionedScalar("DT", dimViscosity, transportProperties);
}

#include "createFvOptions.H"




/* ******************************************************************** 
    
                                MMS

   ******************************************************************** */

// Reference to cell centers
const volVectorField& C = mesh.C();

// Reference to face centers
const surfaceVectorField& Cf = mesh.Cf();

// Field to store the analytical solution (MMS)
volScalarField analyticalSolution
(
    IOobject
    (
        "analyticalSolution",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("analyticalSolution_", dimTemperature, Zero)
);

// Field to store the source term from MMS (S_T)
volScalarField source
(
    IOobject
    (
        "source",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("source_", dimensionSet(0, 0, -1, 1, 0, 0, 0), Zero)
);

// Field to store the difference between the calculated and manufactured solution
volScalarField err
(
    IOobject
    (
        "err",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar(dimTemperature, Zero)
);


// Analytical solution

// Loop through internal cells
forAll(analyticalSolution, cellI)
{
    const scalar x = C[cellI].x(); // x- coordinate of cell center
    const scalar y = C[cellI].y(); // y- coordinate of cell center
    
    analyticalSolution[cellI] = 225.0 + 150*Foam::cos(x*x + y*y); // T(x,y) from MMS
}

// Loop through boundary faces (not acctually need, but stays for reference on how to loop through the boundaries)
forAll(analyticalSolution.boundaryField(), patchI)
{
    forAll(analyticalSolution.boundaryField()[patchI], faceI)
    {
        const scalar x = Cf.boundaryField()[patchI][faceI].x(); // x- coordinate of face center
        const scalar y = Cf.boundaryField()[patchI][faceI].y(); // y- coordinate of face center

        analyticalSolution.boundaryFieldRef()[patchI][faceI] = 225.0 + 150 * Foam::cos(x*x + y*y); // T(x,y) from MMS
    }
}


// Loop though the source term
forAll(source, cellI)
{
    const scalar x = C[cellI].x(); // x- coordinate of cell center
    const scalar y = C[cellI].y(); // y- coordinate of cell center

    // S_T from MMS
    source[cellI] =  (    0.0060000000000000001*Foam::sin(x*x + y*y) 
                        + 0.0060000000000000001*(x*x)*Foam::cos(x*x + y*y) 
                        + 0.0060000000000000001*(y*y)*Foam::cos(x*x + y*y)
                     );

}

