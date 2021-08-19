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
    Foam::waveAbsorptionVelocity2DFvPatchVectorField

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
    waveAbsorptionVelocity2DFvPatchVectorField.C

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

#ifndef waveAbsorptionVelocity2DFvPatchVectorField_H
#define waveAbsorptionVelocity2DFvPatchVectorField_H

#include "fixedValueFvPatchFields.H"
#include "mathematicalConstants.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
/*---------------------------------------------------------------------------*\
               Class 2DWaveAbsorptionVelocityFvPatch Declaration
\*---------------------------------------------------------------------------*/

class waveAbsorptionVelocity2DFvPatchVectorField
:
    public fixedValueFvPatchVectorField
{
    // Private data

        //- Number of paddles
        label nPaddles_;

        //- Initial water depths (referece)
        scalarList initialWaterDepths_;

        //- Direction of absorption (perpendicular to the wall, inward)
        scalar absorptionDir_;

        //- Current velocity vector
        vector uCurrent_;
        
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
   TypeName("waveAbsorption2DVelocity");


   // Constructors

        //- Construct from patch and internal field
        waveAbsorptionVelocity2DFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        waveAbsorptionVelocity2DFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //  waveAbsorptionVelocity2DFvPatchVectorField
        //  onto a new patch
        waveAbsorptionVelocity2DFvPatchVectorField
        (
            const waveAbsorptionVelocity2DFvPatchVectorField&,
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const fvPatchFieldMapper&
        );

        #if OFFLAVOUR == 3 && OFVERSION >= 900
        #else
        //- Construct as copy
        waveAbsorptionVelocity2DFvPatchVectorField
        (
            const waveAbsorptionVelocity2DFvPatchVectorField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchVectorField> clone() const
        {
            return tmp<fvPatchVectorField>
            (
                new waveAbsorptionVelocity2DFvPatchVectorField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        waveAbsorptionVelocity2DFvPatchVectorField
        (
            const waveAbsorptionVelocity2DFvPatchVectorField&,
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
                new waveAbsorptionVelocity2DFvPatchVectorField
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
