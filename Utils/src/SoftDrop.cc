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

#include "BaconProd/Utils/interface/SoftDrop.hh"
#include <fastjet/ClusterSequenceAreaBase.hh>
#include <fastjet/ClusterSequenceActiveAreaExplicitGhosts.hh>
#include <sstream>

FASTJET_BEGIN_NAMESPACE      // defined in fastjet/internal/base.hh

namespace contrib{

using namespace std;

LimitedWarning SoftDropTagger::_warnings_recluster;

//----------------------------------------------------------------------
// SoftDropTagger class implementation
//----------------------------------------------------------------------

//------------------------------------------------------------------------
// description of the tagger
string SoftDropTagger::description() const{ 
  ostringstream oss;
  oss << "SoftDropTagger with beta=" << _beta
      << ", zcut=" << _zcut
      << "and mu=" << _mu;
  return oss.str();
}

//------------------------------------------------------------------------
// returns the tagged PseudoJet if successful, 0 otherwise
//  - jet   the PseudoJet to tag
PseudoJet SoftDropTagger::result(const PseudoJet & jet) const{
  //--------------------------------------------------------------
  // check that the jet has constituents
  if (!jet.has_constituents()){
    throw Error("SoftDropTagger can only be applied to jet with constituents");
    return PseudoJet();
  }

  //--------------------------------------------------------------
  // get the original cliustering radius R0
  if (!jet.has_associated_cluster_sequence()){
    throw Error("SoftDropTagger can only be applied on jets resulting from a previous clustering [i.e. jets with an associated cluster sequence].");
    return PseudoJet();
  }
  double R0 = jet.associated_cs()->jet_def().R();

  //--------------------------------------------------------------
  // Recluster the jet with the requested jet alg (CA by default)
  PseudoJet j = _recluster(jet);

  //--------------------------------------------------------------
  // now do the declustering until the 2 conditions (soft-drop and
  // mass-drop) are met.
  PseudoJet j1, j2;
  bool had_parents;

  double pt1=0.0, pt2=0.0, theta12_squared=0.0;

  while ((had_parents = j.has_parents(j1,j2))) {
    // pre-compute the transverse momenta
    pt1 = j1.pt();
    pt2 = j2.pt();
    // and the angle (including the normalisation)
    theta12_squared = j1.squared_distance(j2)/(R0*R0); 

    // Check mass if desired
    bool passMassDrop = (_mu == 1.0 || (max(j1.m2(), j2.m2()) < _mu*_mu*j.m2()));
     
    // check the 2 conditions
    if ( (min(pt1,pt2) > (pt1+pt2)*_zcut*pow(theta12_squared,0.5*_beta)) && passMassDrop)
      // conditions satisfied: exit the loop
      break;
    else
      // recurse into the subjet with the largest transverse momentum
      j = (j1.pt2() < j2.pt2()) ? j2 : j1;
  }
  
  // create the result and its structure
  PseudoJet result_local = j;
  SoftDropTaggerStructure * s = new SoftDropTaggerStructure(result_local);

  if (!had_parents){
    // no substructure found, return the "leading" particle
    // [Note that this behaviours differs from the MassDropTagger. In
    // this regard, the SoftDropTagger behaves more like a groomer
    // than a tagger]
    s->_mu = 0.0;
    s->_soft_drop = 0.0;
    s->_z = 0.0;
    s->_Rg = 0.0;
  } else {
    s->_mu = (j.m2()!=0.0) ? sqrt(j1.m2()/j.m2()) : 0.0;
    s->_soft_drop = min(pt1,pt2)/(pt1+pt2)*pow(theta12_squared,-0.5*_beta);
    s->_z = min(pt1,pt2)/(pt1+pt2);
    s->_Rg = j1.delta_R(j2);
  }

  result_local.set_structure_shared_ptr(SharedPtr<PseudoJetStructureBase>(s));

  return result_local;
}

//------------------------------------------------------------------------
// rcluster the jet with the Cambridge/Aachen algorithm (or any
// generalized kt algorithm)
// Note that we assume that the jet has constituents and an associated
// cluster sequence
PseudoJet SoftDropTagger::_recluster(const PseudoJet & jet) const{
  // Before everything else, if the jet results from a C/A clustering,
  // (and we have the exponent p = 0.0) we do not need to do anything!
  const JetDefinition & original_jet_def = jet.associated_cs()->jet_def();
  if (original_jet_def.jet_algorithm() == cambridge_algorithm && _p == 0.0)
    return jet;

  // in all other cases, recluster the constituents.
  // 
  // Prepare the cluster sequence and the jet definition.
  // If the input jet has been obtained with a user-defined
  // recombiner, use it also for the re-clustering.
  ClusterSequence *cs;
  JetDefinition jet_def(genkt_algorithm,
                        JetDefinition::max_allowable_R,
                        _p);
  jet_def.set_recombiner(original_jet_def.recombiner());

  // get the jet constituents
  vector<PseudoJet> constituents = jet.constituents();
  if (constituents.size() == 0) return PseudoJet();

  // if we have explicit ghosts, do the clustering keeping area info
  if ((jet.has_area()) &&
      (jet.validated_csab()->has_explicit_ghosts())){
    // we have area readily available from explicit ghosts; split the
    // constituents into ghosts and regular particles and user the
    // former to keep the area information in the new clustering
    std::vector<PseudoJet> particles, ghosts;
    SelectorIsPureGhost().sift(constituents, ghosts, particles);

    // get the ghost area (if there is no ghost in the jet, any number
    // would do here!)
    double ghost_area = ghosts.size() ? ghosts[0].area() : 0.01;

    // do the clustering using the particles and ghosts
    cs = new ClusterSequenceActiveAreaExplicitGhosts(particles, jet_def, 
						     ghosts, ghost_area);
  } else {  // no area info or area without explicit ghosts
    // use a regular clustering without area info
    cs = new ClusterSequence(constituents, jet_def);
  }

  // keep clustering info available
  vector<PseudoJet> new_jets = cs->inclusive_jets();
  
  // print a warning if we get more than one jet
  if (new_jets.size()!=1){
    _warnings_recluster.warn("SoftDropTagger::_recluster: more than one jet obtained after re-clustering. Keeping the hardest");
    new_jets = sorted_by_pt(new_jets);
  }

  // make sure we keep the re-clustering sequence alive
  cs->delete_self_when_unused();
  
  return new_jets[0];
}

} // namespace contrib

FASTJET_END_NAMESPACE
