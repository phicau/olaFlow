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
    Foam::waveAbsorptionVelocity3DFvPatchVectorField

Description
    Wave absorption boundary condition based on shallow water theory and on a
    3D approach. Works both in 2D and 3D and for waves out of the shallow water
    regime.

    Example of the boundary condition specification:
    @verbatim
    inlet
    {
        type            waveAbsorption3DVelocity;
        nPaddles        1;
        nEdgeMin        0;
        nEdgeMax        0;
        value           uniform (0 0 0); // placeholder
    }
    @endverbatim

Note
    - The value is positive inwards

SourceFiles
    waveAbsorptionVelocity3DFvPatchVectorField.C

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
| - Realistic wave generation and active wave absorption for Navier-Stokes    |
|    models: Application to OpenFOAM.                                         |
|    Higuera, P., Lara, J.L. and Losada, I.J. (2013)                          |
|    Coastal Engineering, Vol. 71, 102-118.                                   |
|    http://dx.doi.org/10.1016/j.coastaleng.2012.07.002                       |
|                                                                             |
| - Simulating coastal engineering processes with OpenFOAM                    |
|    Higuera, P., Lara, J.L. and Losada, I.J. (2013)                          |
|    Coastal Engineering, Vol. 71, 119-134.                                   |
|    http://dx.doi.org/10.1016/j.coastaleng.2012.06.002                       |
|                                                                             |
\*---------------------------------------------------------------------------*/

#ifndef waveAbsorptionVelocity3DFvPatchVectorField_H
#define waveAbsorptionVelocity3DFvPatchVectorField_H

#include "fixedValueFvPatchFields.H"
#include "mathematicalConstants.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
/*---------------------------------------------------------------------------*\
               Class 3DWaveAbsorptionVelocityFvPatch Declaration
\*---------------------------------------------------------------------------*/

class waveAbsorptionVelocity3DFvPatchVectorField
:
    public fixedValueFvPatchVectorField
{
    // Private data

        //- Number of paddles
        label nPaddles_;

        //- Initial water depth (referece)
        scalarList initialWaterDepths_;

        //- Absorption angle for each paddle
        scalarList meanAngles_;

        //- Z-Span for each paddle
        scalarList zSpanL_;
        
        //- Number of paddles in the starting part of the domain with 2D
        // version of absorption and out flux only
        label nEdgeMin_;

        //- Number of paddles in the ending part of the domain with 2D
        // version of absorption and out flux only
        label nEdgeMax_;

        //- BC has been checked for errors
        bool allCheck_;

public:

   //- Runtime type information
   TypeName("waveAbsorption3DVelocity");


   // Constructors

        //- Construct from patch and internal field
        waveAbsorptionVelocity3DFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        waveAbsorptionVelocity3DFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //  waveAbsorptionVelocity3DFvPatchVectorField
        //  onto a new patch
        waveAbsorptionVelocity3DFvPatchVectorField
        (
            const waveAbsorptionVelocity3DFvPatchVectorField&,
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const fvPatchFieldMapper&
        );

        #if OFFLAVOUR == 3 && OFVERSION >= 900
        #else
        //- Construct as copy
        waveAbsorptionVelocity3DFvPatchVectorField
        (
            const waveAbsorptionVelocity3DFvPatchVectorField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchVectorField> clone() const
        {
            return tmp<fvPatchVectorField>
            (
                new waveAbsorptionVelocity3DFvPatchVectorField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        waveAbsorptionVelocity3DFvPatchVectorField
        (
            const waveAbsorptionVelocity3DFvPatchVectorField&,
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
                new waveAbsorptionVelocity3DFvPatchVectorField
                    (*this, iF)
            );
        }
        #endif


    // Member functions

        //- Update the coefficients associated with the patch field
        virtual void updateCoeffs();

        //- Write
        virtual void write(Ostream&) const;

    // Other common member functions
        #include "memberFun.H"

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
