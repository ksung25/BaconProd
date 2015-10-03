#ifndef BACONPROD_UTILS_TRIGGERTOOLS_HH
#define BACONPROD_UTILS_TRIGGERTOOLS_HH

#include "BaconAna/DataFormats/interface/BaconAnaDefs.hh"
#include "BaconAna/Utils/interface/TriggerRecord.hh"
#include "DataFormats/HLTReco/interface/TriggerEvent.h"
#include <vector>

namespace baconhep {

class TriggerTools
{
public:
  static TriggerObjects  matchHLT(const double eta, const double phi, 
				  const std::vector<TriggerRecord> &triggerRecords,
				  const trigger::TriggerEvent &triggerEvent);
};

}
#endif
