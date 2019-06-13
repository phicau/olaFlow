// Geometric variables & constants
const volScalarField& alpha = 
    db().lookupObject<volScalarField>(alphaName());

const fvMesh& mesh = alpha.mesh();
const word& patchName = this->patch().name();
const label patchID = mesh.boundaryMesh().findPatchID(patchName);
const label nF = patch().faceCells().size();

const scalarField alphaCell = 
    alpha.boundaryField()[patchID].patchInternalField();

// 3D Variables
const vector cMin = gMin(patch().patch().localPoints());
const vector cMax = gMax(patch().patch().localPoints());
const vector cSpan = cMax - cMin;

scalar dMin = 0.0;
scalar dSpan = 0.0;
const scalarField patchD = patchDirection( cSpan, &dMin, &dSpan );

// Default values
word theory = "";
label nPaddles = 1;
scalarList initialWaterDepths = List<scalar> (1, -1.0);
scalar absDir = 666.;
scalar wavePeriod = -1.;
vector uCurrent = vector::zero;
scalarList meanAngles = List<scalar> (1, -1.0);
scalarList zSpanL = List<scalar> (1, -1.0);
label nEdgeMin = 0;
label nEdgeMax = 0;

// Check if data is provided directly in the BC dict
if ( AWADictName_ == "" )
{
    theory = AWADict_.lookupOrDefault<word>("theory", "");
    nPaddles = AWADict_.lookupOrDefault<label>("nPaddles", 1);
    initialWaterDepths = 
        AWADict_.lookupOrDefault<List<scalar> >
            ("initialWaterDepths", List<scalar> (1, -1.0));
    absDir = AWADict_.lookupOrDefault<scalar>("absDir", 666.);
    wavePeriod = AWADict_.lookupOrDefault<scalar>("wavePeriod", -1.);
    uCurrent = AWADict_.lookupOrDefault<vector>("uCurrent", vector::zero);
    nEdgeMin = AWADict_.lookupOrDefault<label>("nEdgeMin", 0);
    nEdgeMax = AWADict_.lookupOrDefault<label>("nEdgeMax", 0);
}
else
{
    // Define dictionary
    IOdictionary fileAWADict
    (
        IOobject
        (
            AWADictName_,
            this->db().time().constant(),
            this->db(),
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        )
    );

    theory = fileAWADict.lookupOrDefault<word>("theory", "");
    nPaddles = fileAWADict.lookupOrDefault<label>("nPaddles", 1);
    initialWaterDepths = 
        fileAWADict.lookupOrDefault<List<scalar> >
            ("initialWaterDepths", List<scalar> (1, -1.0));
    absDir = fileAWADict.lookupOrDefault<scalar>("absDir", 666.);
    wavePeriod = fileAWADict.lookupOrDefault<scalar>("wavePeriod", -1.);
    uCurrent = fileAWADict.lookupOrDefault<vector>("uCurrent", vector::zero);
    nEdgeMin = fileAWADict.lookupOrDefault<label>("nEdgeMin", 0);
    nEdgeMax = fileAWADict.lookupOrDefault<label>("nEdgeMax", 0);
}

// Check for common errors
if (theory == "")
{
    Warning << 
    "Active wave absorption theory not defined in boundary condition.\n" << 
    "The expected keyword is: theory\n" << 
    "Using theory = shallowWaterAWA by default.\n" << endl;

    theory = "shallowWaterAWA";
}

// Check if the value of nPaddles is correct for the number of columns
if (nPaddles < 1)
{
    FatalError
        << "Check nPaddles value."
        << exit(FatalError);
}

if ( nPaddles > 1 )
{
    nPaddles = decreaseNPaddles( nPaddles, patchD, dMin, dSpan );
    reduce(nPaddles, minOp<label>());
}

// Check nEdges
if ( nEdgeMin < 0 || nEdgeMax < 0 )
{
    FatalError
        << "Check nEdgeMin/Max value."
        << exit(FatalError); 
}

if ( nEdgeMin + nEdgeMax > nPaddles )
{
    FatalError
        << "Check: nEdges > nPaddles."
        << exit(FatalError); 
}

// Grouping part
labelList faceGroup = Foam::labelList(nF, 0);
scalarList dBreakPoints = Foam::scalarList(nPaddles+1, dMin); 
scalarList xGroup = Foam::scalarList(nPaddles, 0.0);
scalarList yGroup = Foam::scalarList(nPaddles, 0.0);

for (label i=0; i<nPaddles; i++)
{
    // Breakpoints, X & Y centre of the paddles
    dBreakPoints[i+1] = dMin + dSpan/(nPaddles)*(i+1);
    xGroup[i] = cMin[0] + cSpan[0]/(2.0*nPaddles) + cSpan[0]/(nPaddles)*i;
    yGroup[i] = cMin[1] + cSpan[1]/(2.0*nPaddles) + cSpan[1]/(nPaddles)*i;
}

forAll(patchD, patchCells) 
{
    for (label i=0; i<nPaddles; i++)
    {
        if ( (patchD[patchCells]>=dBreakPoints[i]) && 
            (patchD[patchCells]<dBreakPoints[i + 1]) )
        {
            faceGroup[patchCells] = i + 1; // Group of each face
            continue;
        }
    }      
}

// Calculate Z bounds of the faces
scalarField zSup, zInf;
faceBoundsZ( &zSup, &zInf );
// Z-span paddlewise
zSpanL = zSpanList( faceGroup, zInf, zSup, nPaddles );
//inlinePrint( "Z-span paddlewise", zSpanL );

// Calculate initial water depth
if ( gMin(initialWaterDepths) < 0.0 )
{
    initialWaterDepths = calcWL(alphaCell, faceGroup, zSpanL, nPaddles);
    inlinePrint( "Initial water depths for absorption", 
        initialWaterDepths );
}

// Absorption direction part
meanAngles = initialWaterDepths;

// Check if absorption is directional
if ( absDir > 360.0 ) // Automatic
{
    meanAngles = meanPatchDirs( faceGroup, nPaddles );
}
else // Fixed
{
    meanAngles = absDir*PI()/180.0;
}
//inlinePrint( "Paddle angles", meanAngles );
