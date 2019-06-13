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

Class
    Foam::waveAbsorptionVelocityFvPatchVectorField

Description
    Wave absorption boundary condition based on shallow water theory and on a
    2D approach. Works both in 2D and 3D and for waves out of the shallow water
    regime.

    Example of the boundary condition specification:
    @verbatim
    inlet
    {
        type            waveAbsorption2DVelocity;
        nPaddles        1;
        nEdgeMin        0;
        nEdgeMax        0;
        absorptionDir   400; //I.e. calculate perpendicular to boundary
        value           uniform (0 0 0); // placeholder
    }
    @endverbatim

Note
    - The value is positive inwards
    - absorptionDir in degrees, and should point inside the domain

SourceFiles
    waveAbsorptionVelocityFvPatchVectorField.C

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

#ifndef waveAbsorptionVelocityFvPatchVectorField_H
#define waveAbsorptionVelocityFvPatchVectorField_H

#include "fixedValueFvPatchFields.H"
#include "mathematicalConstants.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
/*---------------------------------------------------------------------------*\
               Class waveAbsorptionVelocityFvPatch Declaration
\*---------------------------------------------------------------------------*/

class waveAbsorptionVelocityFvPatchVectorField
:
    public fixedValueFvPatchVectorField
{
    // Private data

        //- Name of the AWA dictionary to read the absorption parameters from
        word AWADictName_;

        //- AWA parameters to start or restart the AWA model
        dictionary AWADict_;

        //- BC parameters have been checked for errors
        bool allCheck_;

public:

   //- Runtime type information
   TypeName("waveAbsorptionVelocity");


   // Constructors

        //- Construct from patch and internal field
        waveAbsorptionVelocityFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        waveAbsorptionVelocityFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //  waveAbsorptionVelocityFvPatchVectorField
        //  onto a new patch
        waveAbsorptionVelocityFvPatchVectorField
        (
            const waveAbsorptionVelocityFvPatchVectorField&,
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        waveAbsorptionVelocityFvPatchVectorField
        (
            const waveAbsorptionVelocityFvPatchVectorField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchVectorField> clone() const
        {
            return tmp<fvPatchVectorField>
            (
                new waveAbsorptionVelocityFvPatchVectorField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        waveAbsorptionVelocityFvPatchVectorField
        (
            const waveAbsorptionVelocityFvPatchVectorField&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchVectorField> clone
        (
            const DimensionedField<vector, volMesh>& iF
        ) const
        {
            return tmp<fvPatchVectorField>
            (
                new waveAbsorptionVelocityFvPatchVectorField
                    (*this, iF)
            );
        }


    // Member functions

        //- Update the coefficients associated with the patch field
        virtual void updateCoeffs();

        //- Write
        virtual void write(Ostream&) const;

    // Other common member functions
        #include "myMemberFun.H"

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
