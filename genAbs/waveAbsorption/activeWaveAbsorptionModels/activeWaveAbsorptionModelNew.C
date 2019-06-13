/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2016 OpenCFD Ltd.
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

#include "activeWaveAbsorptionModel.H"
#include "fvMesh.H"

Foam::autoPtr<Foam::activeWaveAbsorptionModel> Foam::activeWaveAbsorptionModel::New
(
    const fvMesh& mesh,
    const polyPatch& patch,
    const dictionary& dict
)
{
    word theory = (dict.lookupOrDefault<word>("theory", "" ));
    
    Info<< "Selecting active wave absorption model " << theory << endl;

    typename patchConstructorTable::iterator cstrIter = patchConstructorTablePtr_->find(theory);

    if (cstrIter == patchConstructorTablePtr_->end())
    {
        FatalError
            << "Unknown activeWaveAbsorptionModel type " << theory << nl
            << "Valid activeWaveAbsorptionModel types are:" << nl
            << patchConstructorTablePtr_->sortedToc()
            << exit(FatalError);
    }

    return autoPtr<activeWaveAbsorptionModel>(cstrIter()(mesh, patch, dict));
}


Foam::tmp<Foam::activeWaveAbsorptionModel> Foam::activeWaveAbsorptionModel::lookupOrCreate
(
    const polyPatch& patch,
    const fvMesh& mesh,
    const dictionary& dict
)
{
    word AWA_ID = modelID(patch.name());

    if (!mesh.foundObject<activeWaveAbsorptionModel>(AWA_ID))
    {
        autoPtr<activeWaveAbsorptionModel> 
            model(activeWaveAbsorptionModel::New(mesh, patch, dict));
        activeWaveAbsorptionModel* activeWaveAbsorptionModelPtr = model.ptr();
        activeWaveAbsorptionModelPtr->store();
        activeWaveAbsorptionModelPtr->info(Info);
    }

    return mesh.lookupObject<activeWaveAbsorptionModel>(AWA_ID);
}

// ************************************************************************* //
