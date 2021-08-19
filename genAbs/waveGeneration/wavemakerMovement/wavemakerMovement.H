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

Description
    fixedValue pointPatchField.

    - without state

SourceFiles
    wavemakerMovement.C

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
| - Three-dimensional numerical wave generation with moving boundaries        |
|    Higuera, P., Losada, I.J. and Lara, J.L. (2015)                          |
|    Coastal Engineering, Vol. 101, 35-47.                                    |
|    http://dx.doi.org/10.1016/j.coastaleng.2015.04.003                       |
|                                                                             |
\*---------------------------------------------------------------------------*/

#ifndef fixedValuePointPatchTemplateVectorField_H
#define fixedValuePointPatchTemplateVectorField_H

#include "fixedValuePointPatchFields.H"

#include "mathematicalConstants.H"
#include "fixedValueFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class wavemakerMovement
:
    #if OFFLAVOUR == 1
        public FixedValuePointPatchField
        <
            pointPatchField,
            pointMesh,
            pointPatch,
            DummyMatrix,
            vector
        >
    #else
        public fixedValuePointPatchField<vector>
    #endif
{
    // Private data

        //- Dictionary name
        word wavemakerDictName_;

        //- Piston, Flap, Mixed
        word wavemakerType_;

        //- Time series for interpolation
        scalarList timeSeries_;

        //- Number of different paddles (for absorption)
        label nPaddles_;

        //- Paddle position series for interpolation
        List<List<scalar> > paddlePosition_;

        //- Paddle tilt angle series for interpolation
        List<List<scalar> > paddleTilt_;

        //- Free surface series at the paddle for interpolation
        List<List<scalar> > paddleEta_;

        //- Initial water depth (meters)
        scalar initialWaterDepth_;

        //- Perpendicular direction to the paddle
        vector meanAngle_;

        //- Z location of the flap paddle hinge (meters)
        scalar hingeHeight_;

        //- X location of the flap paddle hinge (meters)
        scalar hingeLocation_;

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

        //- Paddle tilt angle (previous time step)
        scalarList tiltOld_;

        //- Fade in time (linear 0-1)
        scalar tSmooth_;

        //- Tuning factor
        scalar tuningFactor_;

        //- BC has been checked for first time
        bool allCheck_;

public:

    //- Runtime type information
    TypeName("wavemakerMovement");


    // Constructors

        //- Construct from patch and internal field
        wavemakerMovement
        (
            const pointPatch&,
            const DimensionedField<vector, pointMesh>&
        );

        //- Construct from patch, internal field and dictionary
        wavemakerMovement
        (
            const pointPatch&,
            const DimensionedField<vector, pointMesh>&,
            const dictionary&,
            const bool valueRequired=true
        );

        //- Construct by mapping a copy onto a new patch
        wavemakerMovement
        (
            const wavemakerMovement&,
            const pointPatch&,
            const DimensionedField<vector, pointMesh>&,
            #if OFFLAVOUR == 1
                const PointPatchFieldMapper&
            #else
                const pointPatchFieldMapper&
            #endif
        );

        #if OFFLAVOUR == 3 && OFVERSION >= 900
        #else
        //- Construct as copy
        wavemakerMovement
        (
            const wavemakerMovement&
        );

        //- Construct and return a clone
        virtual autoPtr< pointPatchField<vector> > clone() const
        {
            return autoPtr< pointPatchField<vector> >
            (
                new wavemakerMovement(*this)
            );
        }

        //- Construct as copy setting internal field reference
        wavemakerMovement
        (
            const wavemakerMovement&,
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
                new wavemakerMovement(*this, iF)
            );
        }
        #endif


    //- Destructor
    virtual ~wavemakerMovement();


    // Member functions

        //- Update the coefficients associated with the patch field
        virtual void updateCoeffs();

        //- Write
        virtual void write(Ostream&) const;

        #include "wavemakerMemberFun.H"

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

