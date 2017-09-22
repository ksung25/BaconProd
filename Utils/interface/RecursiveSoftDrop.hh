#ifndef __RECURSIVESOFTDROP_HH__
#define __RECURSIVESOFTDROP_HH__
#include "fastjet/tools/Transformer.hh"
#include "fastjet/WrappedStructure.hh"
#include "fastjet/contrib/Recluster.hh"

#include <iostream>
#include <queue>
#include <vector>

FASTJET_BEGIN_NAMESPACE

//------------------------------------------------------------------------
/// \class RecursiveSoftDrop
/// An implementation of the RecursiveSoftDrop.
///
/// Recursive Soft Drop will recursively groom a jet, removing
/// particles that fail the criterion
/// \f[
///     z > z_{\rm cut} (\theta/R0)^\beta
/// \f]
/// until n subjets have been found.
//----------------------------------------------------------------------
class RecursiveSoftDrop : public Transformer {
public:
  /// minimal constructor, which takes a parameter beta and zcut,
  /// as well as optional parameters R0 (default 1.0), n (default -1,
  /// i.e. infinity), dynamic_R0 and keep_struct)
  /// to JetDefinition::max_allowable_R (practically equivalent to
  /// infinity) and also tries to use a recombiner based on the one in
  /// the jet definition of the particular jet being Soft Dropped.
  ///
  ///  \param jet_alg     the jet algorithm for the internal clustering
  ///  \param zcut        the value for zcut
  ///  \param beta        the value for beta
  ///  \param R0          the value for R0
  ///  \param n           the value for n (n = 0 means infinite)
  ///  \param dynamic_R0  the flag for dynamic choice of R0  
  ///  \param keep_struct the flag for the level of detail of StructureType
  RecursiveSoftDrop(double beta, double zcut, double R0 = 1.0, int n = -1,
		    bool dynamic_R0 = false, bool keep_struct = false) 
    : _beta(beta), _zcut(zcut), _R0sqr(R0*R0), _n(n),
      _dynamic_R0(dynamic_R0), _keep_structure(keep_struct) {}

  // associated structure for groomed jets
  class StructureType;
  
  // definitions needed for comparison of subjets
  struct CompareJetsWithDeltaRsqr {
    // return the squared Delta R value between the two subjets
    // is there a better way of doing this by saving the cluster sequence in Recluster??
    double jet_deltaRsqr(const PseudoJet& jet) const {
      PseudoJet piece1, piece2;
      if (jet.has_parents(piece1,piece2))
	return piece1.squared_distance(piece2);
      return 0.0;
    }
    
    bool operator ()(const PseudoJet& j1, const PseudoJet& j2) const {
      return jet_deltaRsqr(j1) < jet_deltaRsqr(j2);
    }

    // alternative definition using a pair with dij and PseudoJet
    // bool operator ()(const std::pair<double, PseudoJet>& p1, const std::pair<double, PseudoJet>& p2) const { 
    //   if (p1.first != p2.first) return p1.first < p2.first;
    //   return p1.second.pt() < p2.second.pt();
    // }
  };

  // // alternative definition for comparison of subjets in terms of kt distance
  // struct CompareJetsWithKtDist {
  //   // return the squared Delta R value between the two subjets
  //   // is there a better way of doing this by saving the cluster sequence in Recluster??
  //   double jet_ktdist(const PseudoJet& jet) const {
  //     PseudoJet piece1, piece2;
  //     if (jet.has_parents(piece1,piece2))
  // 	return piece1.kt_distance(piece2);
  //     return 0.0;
  //   }
    
  //   bool operator ()(const PseudoJet& j1, const PseudoJet& j2) const {
  //     return jet_ktdist(j1) < jet_ktdist(j2);
  //   }

  //   // alternative definition using a pair with dij and PseudoJet
  //   // bool operator ()(const std::pair<double, PseudoJet>& p1, const std::pair<double, PseudoJet>& p2) const { 
  //   //   if (p1.first != p2.first) return p1.first < p2.first;
  //   //   return p1.second.pt() < p2.second.pt();
  //   // }
  // };

  // // alternative definition for comparison in terms of jet pt
  // struct CompareJetsWithPtLargest {    
  //   bool operator ()(const PseudoJet& j1, const PseudoJet& j2) const {
  //     //cout << j1.pt() << ", " << j2.pt() << " / " << (j1.pt() < j2.pt()) << endl;
  //     return j1.pt() < j2.pt();
  //   }
  // };
  // // alternative definition for comparison in terms of smallest jet pt
  // struct CompareJetsWithPtSmallest {    
  //   bool operator ()(const PseudoJet& j1, const PseudoJet& j2) const {
  //     return j1.pt() > j2.pt();
  //   }
  // };

  /// action on a single jet:
  /// this routine applies the Soft Drop criterion recursively on the
  /// CA tree until we find n subjets (or until it converges), and
  /// adds them together into a groomed PseudoJet
  virtual PseudoJet result(const PseudoJet &jet) const;
  
  /// description
  virtual std::string description() const {
    std::ostringstream oss;
    oss << "Recursive Soft Drop with: " << param_description() 
	<< ", " << ( (_dynamic_R0) ? "dynamic R0" : "fixed R0" );
    return oss.str();
  }
  
protected:  
  /// return the symmetry cut zcut(R_ij/R0)^beta for specified R0 value
  double symmetry_cut_fn(double square_distance_ij, double R0sqr) const;
  
  /// return the symmetry cut zcut(R_ij/R0)^beta
  double symmetry_cut_fn(double square_distance_ij) const;

  std::string param_description() const {
    std::ostringstream oss;
    oss << " N = " << _n << ", zcut = " << _zcut
	<< ", beta = " << _beta << ", R0 = " << sqrt(_R0sqr);
    return oss.str();
  }

  /// return false if we reached desired layer of grooming _n
  bool continue_grooming(int current_n) const {
    return ((_n < 0) or (current_n < _n));
  }
  
private:
  double _beta;         ///< the power of the angular distance to be used
                        ///< in the symmetry condition
  double _zcut;         ///< the value of zcut (the prefactor in the asymmetry cut)
  double _R0sqr;        ///< normalisation of the angular distance
                        ///< (typically set to the jet radius, 1 by default)
  int    _n;            ///< the value of n
  bool _dynamic_R0;     ///< flag for dynamic choice of R0
  bool _keep_structure; ///< flag for how much information to store in StructureType object
};
  
//----------------------------------------------------------------------
/// class to hold the structure of a jet tagged by RecursiveSymmetryCutBase.
class RecursiveSoftDrop::StructureType : public WrappedStructure {
public:
  StructureType(const PseudoJet & j) :
    WrappedStructure(j.structure_shared_ptr()),
    _has_structure(false) // by default, do not store full structure
  {}

  // last delta R value
  double deltaR() const {return _delta_R;}
  
  // information branching history
  std::vector<double> Rg() const {
    if (!_has_structure)
      throw Error("RecursiveSoftDrop::StructureType: Verbose structure must be turned on to get Rg() values.");
    return _Rg;
  };
  
  std::vector<double> zg() const {
    if (!_has_structure)
      throw Error("RecursiveSoftDrop::StructureType: Verbose structure must be turned on to get zg() values.");
    return _zg;
  };
  
  
private:
  double _delta_R;
  friend class RecursiveSoftDrop;

  // additional information on branches
  bool _has_structure;
  
  // information about dropped values
  std::vector<double> _Rg; ///< the delta R of the subjets passing the SD criterion
  std::vector<double> _zg; ///< the energy fraction of the subjets passing the SD criterion
  
};

//----------------------------------------------------------------------
/// class to provide alternative Recursive Soft Drop algorithm, to
/// remove non-global contributions in N-prong observables.
/// It is equivalent to RSD in the N->infty limit, but for finite N it
/// applies equal numbers of Soft Drop grooming to each branch.
class RecursiveSoftDropSameDepth : public RecursiveSoftDrop {
public:
  RecursiveSoftDropSameDepth(double beta, double zcut, double R0 = 1.0, int n = -1) 
    : RecursiveSoftDrop(beta, zcut, R0, n, false, false) {}
  
  /// action on a single jet:
  /// this routine applies the Soft Drop criterion recursively with a "depth" N,
  /// i.e. applying it N times on each branch that passes (or until it
  /// converges), and returns a groomed PseudoJet
  virtual PseudoJet result(const PseudoJet &jet) const;

  /// description
  virtual std::string description() const {
    std::ostringstream oss;
    oss << "Same Depth Recursive Soft Drop with: " << param_description() << "";
    return oss.str();
  }

private:
  /// this recursive method goes through the jet's parent history
  /// until we have applied Soft Drop down to depth _n.
  void recursive_steps(const PseudoJet & jet, std::vector<PseudoJet>& pieces,
		       int current_N) const;
};

FASTJET_END_NAMESPACE      // defined in fastjet/internal/base.hh
#endif // __RECURSIVESOFTDROP_HH__
