/*---------------------------------------------------------------------------*\
License
    This file is part of olaFlow Project.

    olaFlow is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    olaFlow is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with olaFlow.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*\
| olaFlow Project                                       ll                    |
|                                                       l l                   |
|   Coder: Pablo Higuera Caubilla                 ooo   l l     aa            |
|   Bug reports: olaFlowCFD@gmail.com            o   o  l l    a  a           |
|                                                o   o  ll   l a  aa  aa      |
|                                                 ooo    llll   aa  aa        |
|                                                                             |
|                                                FFFFF L     OOOOO W   W      |
|                                                F     L     O   O W   W      |
|                                                FFFF  L     O   O W W W      |
|                                                F     L     O   O WW WW      |
|                                                F     LLLLL OOOOO W   W      |
|   -----------------------------------------------------------------------   |
| References:                                                                 |
|                                                                             |
| - Enhancing active wave absorption in RANS models                           |
|    Higuera, P. (2018)                                                       |
|    https://arxiv.org/abs/1810.03492                                         |
|                                                                             |
\*---------------------------------------------------------------------------*/

#include "waveAbsorptionVelocityFvPatchVectorField.H"
#include "volFields.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "surfaceFields.H"
#include "activeWaveAbsorptionModel.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::
waveAbsorptionVelocityFvPatchVectorField::
waveAbsorptionVelocityFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchField<vector>(p, iF),
    AWADictName_(""),
    AWADict_(),
    allCheck_(false)
{}


Foam::
waveAbsorptionVelocityFvPatchVectorField::
waveAbsorptionVelocityFvPatchVectorField
(
    const waveAbsorptionVelocityFvPatchVectorField& ptf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchField<vector>(ptf, p, iF, mapper),
    AWADictName_(ptf.AWADictName_),
    AWADict_(ptf.AWADict_),
    allCheck_(ptf.allCheck_)
{}


Foam::
waveAbsorptionVelocityFvPatchVectorField::
waveAbsorptionVelocityFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchField<vector>(p, iF, dict),
    AWADictName_(dict.lookupOrDefault<word>("AWADictName", "")),
    AWADict_(dict.subOrEmptyDict("AWADict")),
    allCheck_(dict.lookupOrDefault<bool>("allCheck", false ))
{
    if (!dict.isDict("AWADict"))
    {
        AWADict_ = dict;
    }

    fixedValueFvPatchField<vector>::operator==(this->patchInternalField());
    fixedValueFvPatchField<vector>::updateCoeffs();
}


Foam::
waveAbsorptionVelocityFvPatchVectorField::
waveAbsorptionVelocityFvPatchVectorField
(
    const waveAbsorptionVelocityFvPatchVectorField& ptf
)
:
    fixedValueFvPatchField<vector>(ptf),
    AWADictName_(ptf.AWADictName_),
    AWADict_(ptf.AWADict_),
    allCheck_(ptf.allCheck_)
{}


Foam::
waveAbsorptionVelocityFvPatchVectorField::
waveAbsorptionVelocityFvPatchVectorField
(
    const waveAbsorptionVelocityFvPatchVectorField& ptf,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchField<vector>(ptf, iF),
    AWADictName_(ptf.AWADictName_),
    AWADict_(ptf.AWADict_),
    allCheck_(ptf.allCheck_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::waveAbsorptionVelocityFvPatchVectorField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    Info << "Active wave absorption BC on patch " << this->patch().name() << endl;

    // Build the AWA dictionary
    if (!allCheck_)
    {
        #include "initialiseGeometry.H"
        #include "buildAWADict.H"
        allCheck_ = true;
    }

    // Build the AWA object
    tmp<activeWaveAbsorptionModel> AWAmodel
    (
        activeWaveAbsorptionModel::lookupOrCreate
        (
            patch().patch(),
            patch().boundaryMesh().mesh(),
            AWADict_
        )
    );

    activeWaveAbsorptionModel& model = 
        const_cast<activeWaveAbsorptionModel&>(AWAmodel());
    model.correct();
    scalarList refWaterDepths = AWADict_.lookupOrDefault<List<scalar> >
            ("initialWaterDepths", List<scalar> (1, -1.0));
    operator == (model.setVelocity(true, refWaterDepths));
    fixedValueFvPatchField<vector>::updateCoeffs();
}


void Foam::waveAbsorptionVelocityFvPatchVectorField::
write(Ostream& os) const
{
    fvPatchField<vector>::write(os);

    #if OFVERSION >= 1712
        os.writeEntryIfDifferent<word>("AWADictName", "", AWADictName_);
    #else
        writeEntryIfDifferent<word>(os, "AWADictName", "", AWADictName_);
    #endif

    os.writeKeyword("AWADict") << AWADict_ << token::END_STATEMENT << nl;
    os.writeKeyword("allCheck") << allCheck_ << token::END_STATEMENT << nl;

    writeEntry("value", os);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
   makePatchTypeField
   (
       fvPatchVectorField,
       waveAbsorptionVelocityFvPatchVectorField
   );
}


// ************************************************************************* //
