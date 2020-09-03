/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
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

#include "nutGammat.H"
#include "fvOptions.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace incompressible
{
namespace RASModels
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(nutGammat, 0);
addToRunTimeSelectionTable(RASModel, nutGammat, dictionary);

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void nutGammat::correctNut()
{
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

nutGammat::nutGammat
(
    const geometricOneField& alpha,
    const geometricOneField& rho,
    const volVectorField& U,
    const surfaceScalarField& alphaRhoPhi,
    const surfaceScalarField& phi,
    const transportModel& transport,
    const word& propertiesName,
    const word& type
)
:
    eddyViscosity<incompressible::RASModel>
    (
        type,
        alpha,
        rho,
        U,
        alphaRhoPhi,
        phi,
        transport,
        propertiesName
    ),

    //Dummy declarations required by OpenFOAM's baseline kEpsilon
    //turbulence model used as a template
    k_
    (
        IOobject
        (
            IOobject::groupName("k", alphaRhoPhi.group()),
            runTime_.timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
	    mesh_,
        dimensionedScalar("k",dimensionSet(0,2,-2,0,0,0,0),1.0)
    ),

    epsilon_
    (
        IOobject
        (
            IOobject::groupName("epsilon", alphaRhoPhi.group()),
            runTime_.timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
	    mesh_,
        dimensionedScalar("epsilon",dimensionSet(0,2,-3,0,0,0,0),1.0)
    ),

    couplingFactor_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "couplingFactor",
            this->coeffDict_,
            0.0
        )
    ),

    implicitGamma_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "implicitGamma",
            this->coeffDict_,
            0.0
        )
    ),

    nut_
    (
        IOobject
        (
            IOobject::groupName("nut", alphaRhoPhi.group()),
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh_
    ),

    t_
    (
        IOobject
        (
            IOobject::groupName("t", alphaRhoPhi.group()),
            runTime_.timeName(),
            mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_
    ),

    gamma_
    (
        IOobject
        (
            IOobject::groupName("gamma_t", alphaRhoPhi.group()),
            runTime_.timeName(),
            mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_
    )

{
    if (type == typeName)
    {
        
        if (couplingFactor_.value() < 0.0 || couplingFactor_.value() > 1.0)
        {
            FatalErrorInFunction
                << "couplingFactor = " << couplingFactor_
                << " is not in range 0 - 1" << nl
                << exit(FatalError);
        }

        if (implicitGamma_.value() != 0.0 && implicitGamma_.value() != 1.0)
        {
            FatalErrorInFunction
                << "implicitGamma = " << implicitGamma_
                << " is different than 0.0 or 1.0" << nl
                << exit(FatalError);
        }

        printCoeffs(type);
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool nutGammat::read()
{
    if (eddyViscosity<incompressible::RASModel>::read())
    {   
        couplingFactor_.readIfPresent(this->coeffDict());
        implicitGamma_.readIfPresent(this->coeffDict());
        return true;
    }
    else
    {
        return false;
    }
}

void nutGammat::correct()
{
    if (!turbulence_)
    {
        return;
    }

    // Local references
    const alphaField& alpha = this->alpha_;
    const rhoField& rho = this->rho_;
    const surfaceScalarField& alphaRhoPhi = this->alphaRhoPhi_;
    const volVectorField& gamma = this->gamma_;
    volVectorField& t = this->t_;
    fv::options& fvOptions(fv::options::New(this->mesh_));

    eddyViscosity<incompressible::RASModel>::correct();

    // treating part of gamma as a source term through fvm::SuSp
    // gamma is decomposed into a normal to t vector and a orthogonal to t
    // vector. The normal to t vector is decomposed into kappag*t.
    tmp<volScalarField> tmagsqrt(magSqr(t));
    const volScalarField& magsqrt = tmagsqrt();

    dimensionedScalar magmin("",magsqrt.dimensions(),SMALL);

    tmp<volScalarField> tkappa((gamma & t)/max(magsqrt,magmin));
    const volScalarField& kappa = tkappa();
    
    tmp<volVectorField> texp_gamma(gamma - implicitGamma_*kappa*t);
    const volVectorField& exp_gamma = texp_gamma();

    // Adapted Reynolds Force Vector t equation
    tmp<fvVectorMatrix> tEqn
    (
        fvm::ddt(alpha, rho, t)
      + fvm::div(alphaRhoPhi, t)
      - fvm::SuSp(implicitGamma_*alpha*rho*kappa,t)
      - fvm::laplacian(alpha*rho*nuEff(), t)
      ==
        alpha*rho*exp_gamma
      + fvOptions(alpha, rho, t)
    );

    tEqn.ref().relax();
    fvOptions.constrain(tEqn.ref());
    solve(tEqn);
    fvOptions.correct(t);
}

Foam::tmp<Foam::fvVectorMatrix>
nutGammat::divDevRhoReff
(
    volVectorField& U
) const
{

    if (couplingFactor_.value() > 0.0)
    {
        return
        (
            fvc::laplacian
            (
                (1.0 - couplingFactor_)*this->alpha_*rho_*this->nut(),
                U,
                "laplacian(nuEff,U)"
            )
            + fvc::div
            (
                couplingFactor_
                *this->alpha_*rho_*this->nut()*fvc::grad(U),
                "div(devRhoReff)"
            )
            + t_
            - fvc::div(this->alpha_*rho_*this->nu()*dev2(T(fvc::grad(U))))
            - fvm::laplacian(this->alpha_*rho_*this->nuEff(), U)
        );
    }
    else
    {
        return
        (
            fvc::laplacian
            (
                this->alpha_*rho_*this->nut(),
                U,
                "laplacian(nuEff,U)"
            )
            + t_
            - fvc::div(this->alpha_*rho_*this->nu()*dev2(T(fvc::grad(U))))
            - fvm::laplacian(this->alpha_*rho_*this->nuEff(), U)
        );
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RASModels
} // End namespace incompressible
} // End namespace Foam

// ************************************************************************* //
