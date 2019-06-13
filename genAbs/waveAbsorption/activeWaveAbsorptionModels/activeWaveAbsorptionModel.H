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

Class
    Foam::activeWaveAbsorptionModel

Description
    Base class for active wave absorption models

\*---------------------------------------------------------------------------*/

#ifndef activeWaveAbsorptionModel_H
#define activeWaveAbsorptionModel_H

#include "autoPtr.H"
#include "runTimeSelectionTables.H"
#include "mathematicalConstants.H"

#include "IOdictionary.H"
#include "dictionary.H"
#include "pointField.H"
#include "scalarField.H"
#include "vectorField.H"
#include "polyBoundaryMesh.H"
#include "SubField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class fvMesh;
class polyPatch;

/*---------------------------------------------------------------------------*\
                Class activeWaveAbsorptionModel Declaration
\*---------------------------------------------------------------------------*/

class activeWaveAbsorptionModel
:
    public refCount,
    public IOdictionary
{
protected:

    // Protected data

        //- Reference to the mesh database
        const fvMesh& mesh_;

        //- Reference to the polyPatch
        const polyPatch& patch_;

        //- Reference to the input dictionary
        const dictionary& dict_;

        //- Reference to the AWA type
        word theory_;

        //- Number of paddles
        label nPaddles_;

        //- Number of paddles in the starting part of the domain with 2D
        // version of absorption and out flux only
        label nEdgeMin_;

        //- Number of paddles in the ending part of the domain with 2D
        // version of absorption and out flux only
        label nEdgeMax_;
        
        //- X coordinate of paddle group centres
        scalarList xPaddleGroup_;
        
        //- Y coordinate of paddle group centres
        scalarList yPaddleGroup_;

        //- Direction of absorption (perpendicular to the wall, inward)
        scalar absDir_;

        //- Current velocity vector
        vector uCurrent_;
        
        //- Addressing from face to paddle group index
        labelList faceGroup_;
        
        //- Absorption angle for each paddle
        scalarList meanAngles_;

        //- Z-Span for each paddle
        scalarList zSpanL_;

        //- Initial water depths (reference)
        scalarList initialWaterDepths_;

        //- Measured water depths
        scalarList waterDepths_;


    // Protected Member Functions

        //- Initialise
        virtual void initialise();

        //- Water level
        virtual scalarList calcWaterLevels();


public:

    //- Runtime type information
    TypeName("activeWaveAbsorptionModel");

    // Declare runtime constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            activeWaveAbsorptionModel,
            patch,
            (
                const fvMesh& mesh,
                const polyPatch& patch,
                const dictionary& dict
            ),
            (mesh, patch, dict)
        );


    //- Constructor
    activeWaveAbsorptionModel
    (
        const fvMesh& mesh,
        const polyPatch& patch,
        const dictionary& dict
    );


    // Selectors

        //- Return a reference to the selected wave model
        static autoPtr<activeWaveAbsorptionModel> New
        (
            const fvMesh& mesh,
            const polyPatch& patch,
            const dictionary& dict
        );

        //- Lookup activeWaveAbsorptionModel from database, or create new
        static tmp<activeWaveAbsorptionModel> lookupOrCreate
        (
            const polyPatch& patch,
            const fvMesh& mesh,
            const dictionary& dict
        );


    //- Destructor
    virtual ~activeWaveAbsorptionModel();


    // Public Member Functions

        //- Return the unique model ID
        static word modelID(const word& patchName);

        //- Read from dictionary
        virtual bool readDict(const dictionary& overrideDict);

        //- Correct the model
        virtual void correct();

        //- Calculate the absorption model alpha
        virtual scalarField setAlpha();

        //- Calculate the absorption model velocity
        virtual vectorField setVelocity
            (bool isPureAWA, scalarList refWaterDepths);

        //- Return the reference to X coordinate of paddle group centres
        virtual const scalarList& xPaddleCentres() const;

        //- Return the reference to Y coordinate of paddle group centres
        virtual const scalarList& yPaddleCentres() const;

        //- Info
        virtual void info(Ostream& os) const;

    // Other common member functions
        #include "newMemberFun.H"
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
