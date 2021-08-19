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
    Foam::waveAlphaFvPatchScalarField

Description
    Describes a volumetric/mass flow normal Scalar boundary condition by its
    magnitude as an integral over its area.

    The basis of the patch (volumetric or mass) is determined by the
    dimensions of the flux, phi.
    The current density is used to correct the Alpha when applying the
    mass basis.

    Example of the boundary condition specification:
    @verbatim
    inlet
    {
        type            waveAlpha;
        value           uniform (0 0 0); // placeholder
    }
    @endverbatim

Note
    - The value is positive inwards
    - May not work correctly for transonic inlets
    - Strange behaviour with potentialFoam since the U equation is not solved

SourceFiles
    waveAlphaFvPatchScalarField.C

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

#ifndef waveAlphaFvPatchScalarField_H
#define waveAlphaFvPatchScalarField_H

#include "fixedValueFvPatchFields.H"
#include "mathematicalConstants.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
/*---------------------------------------------------------------------------*\
               Class waveAlphaFvPatch Declaration
\*---------------------------------------------------------------------------*/

class waveAlphaFvPatchScalarField
:
    public fixedValueFvPatchScalarField
{
    // Private data

        //- Wave period (seconds)
        scalar wavePeriod_;

        //- Wave periods for irregular case (seconds)
        scalarList wavePeriods_;

        //- Wave height (meters)
        scalar waveHeight_;

        //- Wave heights for irregular case (meters)
        scalarList waveHeights_;

        //- Wave length (meters)
        scalar waveLength_;

        //- Wave lengths for irregular case (meters)
        scalarList waveLengths_;

        //- Water depth (meters)
        scalar waterDepth_;

        //- Wave phase (radians)
        scalar wavePhase_;

        //- Wave phases for irregular case (radians)
        scalarList wavePhases_;
        
        //- Time lag to match the phase
        scalar timeLag_;
        
        //- Time lag to match the phase for irregular case
        scalarList timeLags_;

        //- Number of solitary waves
        label nSolitaryWaves_;

        //- Lambda parameter
        scalar lambdaStokesV_;

        //- m Cnoidal 
        scalar mCnoidal_;

        //- Parameter a - Stokes III
        scalar aE_SIII_;

        //- Parameter kl - Stokes III
        scalar klE_SIII_;

        //- Mean fluid speed in frame of reference (stream function)
        scalar uMean_;

        //- Stream function Bj coefficients
        scalarList Bjs_;

        //- Stream function Bj coefficients
        scalarList Ejs_;

        //- Generation + Absorption at the same time
        bool genAbs_;

        //- Second order irregular generation
        bool secondOrder_;

        //- Number of different paddles (for absorption)
        label nPaddles_;

        //- Fade in time (linear 0-1)
        scalar tSmooth_;

        //- Tuning factor
        scalar tuningFactor_;

        //- Number of different components (for irregular)
        label nComp_;

        //- Dictionary name
        word waveDictName_;

        //- Regular or Irregular
        word waveType_;

        //- Name of the theory
        word waveTheory_;

        //- Name of the theory
        word waveTheoryOrig_;

        //- BC has been checked for errors
        bool allCheck_;

        //- Direction of absorption (degrees, from X axis)
        scalar absDir_;

        //- Direction of propagation (degrees, from X axis)
        scalar waveDir_;

        //- Direction of propagation for irregular case (degrees, from X axis)
        scalarList waveDirs_;

        //- Time series for interpolation
        scalarList timeSeries_;

        //- Paddle position series for interpolation
        scalarList paddlePosition_;

        //- Horizontal velocity series for interpolation
        scalarList paddleVelocityU_;

        //- Vertical velocity series for interpolation
        scalarList paddleVelocityW_;

        //- Free surface series at the paddle for interpolation
        scalarList paddleEta_;


public:

   //- Runtime type information
   TypeName("waveAlpha");


   // Constructors

        //- Construct from patch and internal field
        waveAlphaFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        waveAlphaFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //  waveAlphaFvPatchScalarField
        //  onto a new patch
        waveAlphaFvPatchScalarField
        (
            const waveAlphaFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        #if OFFLAVOUR == 3 && OFVERSION >= 900
        #else
        //- Construct as copy
        waveAlphaFvPatchScalarField
        (
            const waveAlphaFvPatchScalarField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchScalarField> clone() const
        {
            return tmp<fvPatchScalarField>
            (
                new waveAlphaFvPatchScalarField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        waveAlphaFvPatchScalarField
        (
            const waveAlphaFvPatchScalarField&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchScalarField> clone
        (
            const DimensionedField<scalar, volMesh>& iF
        ) const
        {
            return tmp<fvPatchScalarField>
            (
                new waveAlphaFvPatchScalarField(*this, iF)
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
