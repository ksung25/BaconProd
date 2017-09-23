#include "BaconProd/Utils/interface/RecursiveSoftDrop.hh"

using namespace std;

FASTJET_BEGIN_NAMESPACE      // defined in fastjet/internal/base.hh

// return the symmetry cut = zcut (R_ij/R0)^beta for specified R0 value
double RecursiveSoftDrop::symmetry_cut_fn(double square_distance_ij, double R0sqr) const {
  return _zcut * pow(square_distance_ij/R0sqr, 0.5*_beta);
}

// return the symmetry cut = zcut (R_ij/R0)^beta
double RecursiveSoftDrop::symmetry_cut_fn(double square_distance_ij) const {
  return symmetry_cut_fn(square_distance_ij, _R0sqr);
}

// this routine applies the Soft Drop criterion recursively on the
// CA tree until we find n subjets (or until it converges), and
// adds them together into a groomed PseudoJet
PseudoJet RecursiveSoftDrop::result(const PseudoJet &jet) const {
  // start by reclustering jet with C/A algorithm
  PseudoJet ca_jet =
    contrib::Recluster(cambridge_algorithm, JetDefinition::max_allowable_R)(jet);

  // create a priority_queue containing the d_ij values and their corresponding subjet
  // std::priority_queue< std::pair<double, PseudoJet>,
  // 		       std::vector< pair<double, PseudoJet> >,
  // 		       CompareJetsWithDeltaRsqr >  subjets;
  // create a priority queue containing the subjets and a comparison definition
  priority_queue<PseudoJet, vector<PseudoJet>, CompareJetsWithDeltaRsqr>  subjets;
  
  // alternative options for comparison with jet pt
  // std::priority_queue<PseudoJet, std::vector<PseudoJet>,
  // 		      CompareJetsWithPtLargest>  subjets;
  // std::priority_queue<PseudoJet, std::vector<PseudoJet>,
  // 		      CompareJetsWithPtSmallest>  subjets;
  // vector<PseudoJet> solitons;
  
  // initialize counter to 1 subjet (i.e. the full ca_jet)
  subjets.push(ca_jet);
  int counter = 1;
  // initialize max value for subjets
  int max_njet = ca_jet.constituents().size();
  
  // initialize these vectors for use when _keep_structure is enabled
  vector<double> Rg;
  vector<double> zg;
  
  // loop over C/A tree until we reach the appropriate number of subjets
  double R0sqr = _R0sqr;
  PseudoJet piece1, piece2;
  while (continue_grooming(counter-1)) {
    if (subjets.top().has_parents(piece1,piece2)) {
      double sym = piece1.pt() + piece2.pt();
      if (sym == 0) break; 
      sym = min(piece1.pt(), piece2.pt()) / sym;
      double sqr_dist = piece1.squared_distance(piece2);
      double symmetry_cut = symmetry_cut_fn(sqr_dist, R0sqr);
      // double symmetry_cut = _zcut * pow(sqr_dist/R0sqr, 0.5*_beta);
      if (sym > symmetry_cut) {
	subjets.pop();
	subjets.push(piece1);
	subjets.push(piece2);
	if (_dynamic_R0)
	  R0sqr = sqr_dist;
	if (_keep_structure) {
	  Rg.push_back(sqrt(sqr_dist/_R0sqr));
	  zg.push_back(sym);
	}
	++counter;
      } else {
	subjets.pop();
	int choice = piece1.pt2() > piece2.pt2() ? 1 : 2;
	PseudoJet subjet = (choice == 1) ? piece1 : piece2;
	PseudoJet discarded = (choice == 1) ? piece2 : piece1;
	subjets.push(subjet);
	max_njet -= discarded.constituents().size();
      }
    }
    // else if (counter < max_njet) {
    // 	solitons.push_back(subjets.top());
    // 	subjets.pop();
    // }
    if (counter == max_njet) break;
  }
  // now create and fill the pieces vector with all the subjets
  vector<PseudoJet> pieces(subjets.size());
  copy(&(subjets.top()), &(subjets.top()) + subjets.size(), &pieces[0]);
  // vector<PseudoJet> pieces(subjets.size() + solitons.size());
  // copy(&(subjets.top()), &(subjets.top()) + subjets.size(), &pieces[0]);
  // result.insert(pieces.end(), solitons.begin(), solitons.end() );
  
  // create the final result jet
  PseudoJet result = join(pieces);

  // get the delta R value of the final step
  double lastDeltaR = 0.0;
  if (subjets.top().has_parents(piece1,piece2))
    lastDeltaR = piece1.delta_R(piece2);
  
  // and lastly, set up the associated StructureType
  StructureType * structure = new StructureType(result);
  structure->_delta_R = lastDeltaR;
  if (_keep_structure) {
    structure->_has_structure = true;
    structure->_Rg = Rg;
    structure->_zg = zg;
  }
  result.set_structure_shared_ptr(SharedPtr<PseudoJetStructureBase>(structure));

  return result;
}

//----------------------------------------------------------------------
// this routine applies the Soft Drop criterion recursively with a "depth" N,
// i.e. applying it N times on each branch that passes (or until it
// converges), and returns a groomed PseudoJet
PseudoJet RecursiveSoftDropSameDepth::result(const PseudoJet &jet) const {
  // start by reclustering jet with C/A algorithm
  PseudoJet ca_jet =
    contrib::Recluster(cambridge_algorithm, JetDefinition::max_allowable_R)(jet);
  vector<PseudoJet> pieces;
  recursive_steps(ca_jet, pieces, 0);
  return join(pieces);
}


// this recursive method goes through the jet's parent history
// until we have applied Soft Drop down to depth _n.
void RecursiveSoftDropSameDepth::recursive_steps(const PseudoJet& jet,
						 vector<PseudoJet> & pieces,
						 int current_N) const {
  if (continue_grooming(current_N)) {
    PseudoJet piece1, piece2;
    if (jet.has_parents(piece1,piece2)) {
      double sym = piece1.pt() + piece2.pt();
      if (sym == 0) return; 
      sym = min(piece1.pt(), piece2.pt()) / sym;
      double sqr_dist = piece1.squared_distance(piece2);
      double symmetry_cut = symmetry_cut_fn(sqr_dist);
      //double symmetry_cut = _zcut * pow(sqr_dist/R0sqr, 0.5*_beta);
      if (sym > symmetry_cut) {
	// found two prongs, continue applying Soft Drop to both of
	// them, increasing current_N value by 1
	recursive_steps(piece1, pieces, current_N + 1);
	recursive_steps(piece2, pieces, current_N + 1);
      } else {
	// didn't find two prongs, keep going but don't increase current_N
	int choice = piece1.pt2() > piece2.pt2() ? 1 : 2;
	PseudoJet subjet = (choice == 1) ? piece1 : piece2;
	recursive_steps(subjet, pieces, current_N);
      }
    } else {
      // if we arrive here, jet has no parents, so we can stop recursing
      pieces.push_back(jet);
    }
  } else {
    // if we arrive here, we have arrived at "grooming depth" N, so we stop
    pieces.push_back(jet);
  }
}
FASTJET_END_NAMESPACE
