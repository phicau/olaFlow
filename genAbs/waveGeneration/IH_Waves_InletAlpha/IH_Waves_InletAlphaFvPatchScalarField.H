/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2006-2010 OpenCFD Ltd.
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

Class
    Foam::IH_Waves_InletAlphaFvPatchScalarField

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
        type            IH_Waves_InletAlpha;
        value           uniform (0 0 0); // placeholder
    }
    @endverbatim

Note
    - The value is positive inwards
    - May not work correctly for transonic inlets
    - Strange behaviour with potentialFoam since the U equation is not solved

SourceFiles
    IH_Waves_InletAlphaFvPatchScalarField.C

\*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*\
| Work developed at IH Cantabria       IIIII H   H FFFFF OOOOO AAAAA M   M    |
|                                        I   H   H F     O   O A   A MM MM    |
|   Coder: Pablo Higuera Caubilla        I   HHHHH FFFF  O   O AAAAA M M M    |
|   Bug reports: higuerap@unican.es      I   H   H F     O   O A   A M   M    |
|                                      IIIII H   H F     OOOOO A   A M   M    |
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

#ifndef IH_Waves_InletAlphaFvPatchScalarField_H
#define IH_Waves_InletAlphaFvPatchScalarField_H

#include "fixedValueFvPatchFields.H"
#include "mathematicalConstants.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
/*---------------------------------------------------------------------------*\
               Class IH_Waves_InletAlphaFvPatch Declaration
\*---------------------------------------------------------------------------*/

class IH_Waves_InletAlphaFvPatchScalarField
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

        //- Lambda parameter
        scalar lambdaStokesV_;

        //- m Cnoidal 
        scalar mCnoidal_;

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

        //- Paddle velocity series for interpolation
        scalarList paddleVelocity_;

        //- Free surface series at the paddle for interpolation
        scalarList paddleEta_;


public:

   //- Runtime type information
   TypeName("IH_Waves_InletAlpha");


   // Constructors

        //- Construct from patch and internal field
        IH_Waves_InletAlphaFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        IH_Waves_InletAlphaFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //  IH_Waves_InletAlphaFvPatchScalarField
        //  onto a new patch
        IH_Waves_InletAlphaFvPatchScalarField
        (
            const IH_Waves_InletAlphaFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        IH_Waves_InletAlphaFvPatchScalarField
        (
            const IH_Waves_InletAlphaFvPatchScalarField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchScalarField> clone() const
        {
            return tmp<fvPatchScalarField>
            (
                new IH_Waves_InletAlphaFvPatchScalarField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        IH_Waves_InletAlphaFvPatchScalarField
        (
            const IH_Waves_InletAlphaFvPatchScalarField&,
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
                new IH_Waves_InletAlphaFvPatchScalarField(*this, iF)
            );
        }


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
