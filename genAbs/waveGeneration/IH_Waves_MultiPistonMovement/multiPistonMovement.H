/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2012 OpenFOAM Foundation
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

Description
    fixedValue pointPatchField.

    - without state

SourceFiles
    multiPistonMovement.C

\*---------------------------------------------------------------------------*/

#ifndef fixedValuePointPatchTemplateVectorField_H
#define fixedValuePointPatchTemplateVectorField_H

#include "fixedValuePointPatchFields.H"

#include "mathematicalConstants.H"
#include "fixedValueFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class IH_Waves_MultiPistonMovement
:
    public fixedValuePointPatchField<vector>
{
    // Private data

        //- Dictionary name
        word multiPistonDictName_;

        //- Time series for interpolation
        scalarList timeSeries_;

        //- Number of different paddles (for absorption)
        label nPaddles_;

        //- Paddle position series for interpolation
        List<List<scalar> > paddlePosition_;

        //- Free surface series at the paddle for interpolation
        List<List<scalar> > paddleEta_;

        //- Initial water depth (meters)
        scalar initialWaterDepth_;

        //- Perpendicular direction to the paddle
        vector meanAngle_;

        //- Generation + Absorption at the same time
        bool genAbs_;

        //- Drift Prevention System acting
        List<bool> DPS_;

        //- Maximum allowed stroke
        scalar maxStroke_;

        //- DPS period
        scalar DPST_;

        //- Drift Prevention System sign (+/-1)
        scalarList DPSsign_;

        //- DPS starting time
        scalarList DPStIni_;

        //- Instantaneous paddle position correction 
        scalarList instDPSCorrection_;

        //- Cumulative paddle position correction 
        scalarList cumDPSCorrection_;

        //- Cumulative paddle position correction
        scalarList cumAbsCorrection_;

        //- Fade in time (linear 0-1)
        scalar tSmooth_;

        //- Tuning factor
        scalar tuningFactor_;

        //- BC has been checked for first time
        bool allCheck_;

public:

    //- Runtime type information
    TypeName("IH_Waves_MultiPistonMovement");


    // Constructors

        //- Construct from patch and internal field
        IH_Waves_MultiPistonMovement
        (
            const pointPatch&,
            const DimensionedField<vector, pointMesh>&
        );

        //- Construct from patch, internal field and dictionary
        IH_Waves_MultiPistonMovement
        (
            const pointPatch&,
            const DimensionedField<vector, pointMesh>&,
            const dictionary&,
            const bool valueRequired=true
        );

        //- Construct by mapping a copy onto a new patch
        IH_Waves_MultiPistonMovement
        (
            const IH_Waves_MultiPistonMovement&,
            const pointPatch&,
            const DimensionedField<vector, pointMesh>&,
            const pointPatchFieldMapper&
        );

        //- Construct as copy
        IH_Waves_MultiPistonMovement
        (
            const IH_Waves_MultiPistonMovement&
        );

        //- Construct and return a clone
        virtual autoPtr< pointPatchField<vector> > clone() const
        {
            return autoPtr< pointPatchField<vector> >
            (
                new IH_Waves_MultiPistonMovement(*this)
            );
        }

        //- Construct as copy setting internal field reference
        IH_Waves_MultiPistonMovement
        (
            const IH_Waves_MultiPistonMovement&,
            const DimensionedField<vector, pointMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual autoPtr< pointPatchField<vector> > clone
        (
            const DimensionedField<vector, pointMesh>& iF
        ) const
        {
            return autoPtr< pointPatchField<vector> >
            (
                new IH_Waves_MultiPistonMovement(*this, iF)
            );
        }


    //- Destructor
    virtual ~IH_Waves_MultiPistonMovement();


    // Member functions

        //- Update the coefficients associated with the patch field
        virtual void updateCoeffs();

        //- Write
        virtual void write(Ostream&) const;

        #include "pistonMemberFun.H"

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

