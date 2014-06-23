// $Id$
//
// Copyright (c) 2013,
// Andrew Larkoski, Simone Marzani, Gregory Soyez, and Jesse Thaler
//
//----------------------------------------------------------------------
// This file is part of FastJet contrib.
//
// It is free software; you can redistribute it and/or modify it under
// the terms of the GNU General Public License as published by the
// Free Software Foundation; either version 2 of the License, or (at
// your option) any later version.
//
// It is distributed in the hope that it will be useful, but WITHOUT
// ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
// or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public
// License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this code. If not, see <http://www.gnu.org/licenses/>.
//----------------------------------------------------------------------

#ifndef __FASTJET_CONTRIB_SOFTDROP_HH__
#define __FASTJET_CONTRIB_SOFTDROP_HH__

#include <fastjet/internal/base.hh>
#include <fastjet/tools/Transformer.hh>
#include <fastjet/WrappedStructure.hh>
#include <fastjet/LimitedWarning.hh>

FASTJET_BEGIN_NAMESPACE      // defined in fastjet/internal/base.hh

namespace contrib{

// a few forward declarations
class SoftDropTagger;
class SoftDropTaggerStructure;

//------------------------------------------------------------------------
/// \class SoftDropTagger
/// Class that helps perform the SoftDrop declustering introduced
/// by Andrew Larkoski, Simone Marzani, Gregory Soyez, and Jesse
/// Thaler in arXiv:YYMM.NNNN.
///
/// The tagger proceeds as follows:
///
///  0. recluster the jet with the Cambridge/Aachen algorithm
///     (or any of the generalized kT algorithms.)
///
///  1. undo the last step of the clustering step j -> j1 + j2.
///  
///  2. if the two following conditions are satisfied 
///       (i) there is a soft frop, i.e. \f${\rm
///           min}(p_{tj1},p_{tj2})/(p_{tj1},p_{tj2}) > z_{/rm cut}
///           (\Delta R_{j1,j2}/R_0)^\beta\f$,
///      (ii) there is a mass drop, i.e. \f$m_{j1}/m_{j} < \mu\f$
///     keep j as the result of the tagger (with j1 and j2 its 2
///     subjets)
///
///     In the above conditions, zcut, beta and mu are parameters of
///     the tagger and R0 is the radius parameter of the original jet.
///
///  3. otherwise, redefine j to be equal to the largest pt of j1 and
///     j2 and return to step 1.
///
/// \section desc Options
/// 
/// The constructor has the following arguments:
///   - The first  argument is the angular exponent beta 
///   - The second argument is the soft-drop cut (zcut) [0.1 by default]
///   - The third  argument is the minimal mass drop that is required
///      (mu) [1 by default]
///   - The fourth argument is the exponent for the generalized
///      kt algorithm [0 = cambridge_aachen by default]
///
/// Note that the arguments are in an opposite order as the
/// MassDropTagger. This allows to keep mu=1 as a default.
///
/// \section input Input conditions
/// 
///  - the original jet must have constituents
///
///  - the current version requires the jet to have an associated
///    ClusterSequence.
///
///    Note that for future versions we intend to also support more
///    generic jets [e.g. CompositeJets]. [Note for developpers: in
///    practice, this means being a bit more careful with the
///    reculstering + allowing the user to specify R0 explicitly]
///
///
/// \section output Output/structure
/// 
///  - the 2 subjets are kept as pieces if some substructure is found,
///    otherwise a single 0-momentum piece is returned
///  - the 'mu' and 'z/(\Delta R_{j1,j2}/R_0)^\beta' values
///    corresponding to the unclustering step that passed the tagger's
///    cuts
///
/// See also example.cc for a usage example.
///
///----------------------------------------------------------------------
/// Open question(s) [for developpers]: 
///  - Q: R0 can be undefined in some cases [e.g. if the input jet is a
///       composite jet]?
///    A: In a future version, support more generic inputs. For these
///       cases, add a 2nd ctor where one specifies R0.
///
///  - Q: do we have a set_subtractor for interaction with PU subtraction?
///    A: Yes [to be implemented!]
///----------------------------------------------------------------------
///
class SoftDropTagger : public Transformer{
public:
  /// default ctor
  ///
  /// The parameters are the following:
  ///   \param beta  the angular exponent to use for the soft drop criterion
  ///   \param zcut  the energy cut of the soft drop criterion [0.1 by default]
  ///   \param mu    the minimal mass-drop required [1 by default]
  ///   \param p     the exponent in generalized kt [0 = CA by default]
  SoftDropTagger(const double beta,
		 const double zcut=0.1, 
		 const double mu=1.0,
       const double p=0.0)
    : _beta(beta), _zcut(zcut), _mu(mu), _p(p){};

  /// returns a textual description of the tagger
  virtual std::string description() const;

  /// runs the tagger on the given jet and
  /// returns the tagged PseudoJet if successful, a PseudoJet==0 otherwise
  ///   \param jet   the PseudoJet to tag
  /// Note that standard usage is through the () operator rather than
  /// directly calling "result"
  virtual PseudoJet result(const PseudoJet & jet) const;

  /// the type of the associated structure
  typedef SoftDropTaggerStructure StructureType;

private:
  double _beta; ///< the angular exponent for the soft-drop condition
  double _zcut; ///< the energy parameter for the soft-drop condition
  double _mu;   ///< the parameter of the mass-drop condition
  double _p;    ///< the value of the generalized kt exponent

  /// recluster the input jet with the Cambridge-Aachen algorithm
  PseudoJet _recluster(const PseudoJet &jet) const;

  static LimitedWarning _warnings_recluster;
};

//------------------------------------------------------------------------
/// \class SoftDropTaggerStructure
/// the structure returned by the BetaMassDropTagger transformer.
///
/// See the BetaMassDropTagger class description for the details of what
/// is inside this structure
///
class SoftDropTaggerStructure : public WrappedStructure{
public:
  /// ctor with initialisation
  ///  \param pieces  the pieces of the created jet
  ///  \param rec     the recombiner from the underlying cluster sequence
  SoftDropTaggerStructure(const PseudoJet & result_jet) :
  WrappedStructure(result_jet.structure_shared_ptr()), 
  _soft_drop(0.0), _mu(0.0), _z(0.0), _Rg(0.0) {}

  /// returns the value of z/(\Delta R_{j1,j2}/R_0)^\beta for the
  /// splitting that triggered the soft-drop condition
  inline double soft_drop() const {return _soft_drop;}

  /// returns the mass-drop ratio, pieces[0].m()/jet.m(), for the splitting
  /// that triggered the soft-drop condition
  inline double mu() const{return _mu;}
  
  /// returns the energy fraction z for the splitting
  /// where the soft-drop condition was satisfied
  inline double z() const{return _z;}

  /// returns the angle between the two groomed subjets
  /// at the end of the soft-drop procedure
  inline double Rg() const{return _Rg;}
  

protected:
  double _soft_drop; ///< the value of the soft-drop parameter
  double _mu;        ///< the value of the mass-drop parameter
  double _z;         ///< the value of the energy fraction
  double _Rg;        ///< the value of the groomed subjet separation

  // allow the tagger to set these
  friend class SoftDropTagger;
};

} // namespace contrib

FASTJET_END_NAMESPACE

#endif  // __FASTJET_CONTRIB_SOFTDROP_HH__
