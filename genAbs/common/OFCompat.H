#pragma once

#include "scalar.H"
#include "word.H"
#include "Ostream.H"

// Centralised OpenFOAM version/flavour compatibility adapters.
// All #if OFVERSION / OFFLAVOUR guards in the codebase are defined here only.
// Source files and headers call through OFCompat:: instead of branching directly.

namespace OFCompat {

// --- PI constant --------------------------------------------------------------
inline Foam::scalar pi() noexcept
{
#   if OFVERSION >= 200
        return Foam::constant::mathematical::pi;
#   else
        return Foam::mathematicalConstant::pi;
#   endif
}

// --- Alpha field name ---------------------------------------------------------
inline Foam::word alphaName() noexcept
{
#   if OFVERSION >= 230
        return "alpha.water";
#   else
        return "alpha1";
#   endif
}

// --- writeEntryIfDifferent ----------------------------------------------------
// Wraps the member-function form (>= 1712) vs the free-function form (< 1712).
template<class T>
void writeEntryIfDifferent(
    Foam::Ostream& os,
    const Foam::word& key,
    const T& defaultVal,
    const T& val)
{
#   if OFVERSION >= 1712
        os.writeEntryIfDifferent<T>(key, defaultVal, val);
#   else
        Foam::writeEntryIfDifferent<T>(os, key, defaultVal, val);
#   endif
}

// --- writeListEntry -----------------------------------------------------------
// Wraps the list/field write-entry API change introduced in OpenFOAM.org v7.
// Guard uses both OFFLAVOUR and OFVERSION: OpenFOAM.com (OFFLAVOUR=2) uses 4-digit
// version numbers that are all >= 700, but still uses the old member form.
template<class T>
void writeListEntry(
    Foam::Ostream& os,
    const Foam::word& key,
    const T& val)
{
#   if OFFLAVOUR == 3 && OFVERSION >= 700
        writeEntry(os, key, val);
#   else
        val.writeEntry(key, os);
#   endif
}

// --- Component ref accessor ---------------------------------------------------
// .ref() was added in OpenFOAM 4.0 for non-const access to geometric field
// components.  Used in Dev AWA models (integrated in Phase 6).
// Usage: OFCompat::componentRef(UPatch, 2) = scalar(0);
template<class Field>
decltype(auto) componentRef(Field& f, Foam::label cmpt)
{
#   if OFVERSION >= 400
        return f.component(cmpt).ref();
#   else
        return f.component(cmpt);
#   endif
}

} // namespace OFCompat


// --- Point patch base class aliases -------------------------------------------
// Only compiled when fixedValuePointPatchFields.H has already been included
// (wavemakerMovement.H includes it before OFCompat.H so the aliases are
// available in that class declaration).
#ifdef fixedValuePointPatchFields_H
#   if OFFLAVOUR == 1
        namespace OFCompat {
            // FOAM-Extend uses a 5-parameter template form.
            template<class T>
            using FixedValuePointPatch =
                Foam::FixedValuePointPatchField<
                    Foam::pointPatchField,
                    Foam::pointMesh,
                    Foam::pointPatch,
                    Foam::DummyMatrix,
                    T>;
            using PointPatchMapper = Foam::PointPatchFieldMapper;
        }
#   else
        namespace OFCompat {
            template<class T>
            using FixedValuePointPatch = Foam::fixedValuePointPatchField<T>;
            using PointPatchMapper = Foam::pointPatchFieldMapper;
        }
#   endif
#endif
