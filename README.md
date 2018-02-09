BaconProd
=========

Package for producing bacon files

 * Depends on BaconAna package
 * Place package in `$CMSSW_BASE/src` area

All objects are declared in BaconAna and filled in BaconProd/Ntupler, see e.g.:

[TJet for Jets](https://github.com/ksung25/BaconAna/blob/master/DataFormats/interface/TJet.hh)

[FillerJet for Jets](https://github.com/ksung25/BaconProd/blob/master/Ntupler/src/FillerJet.cc)

Setup
----------

Run setup script after setting CMSSW environment:

```Shell
source BaconProd/scripts/setup_prod.sh
```

```
$ git clone https://github.com/ksung/BaconProd.git
$ git clone https://github.com/ksung/BaconAna.git
$ scram b -j 10
```

Production
----------

Is done through the following config scripts:

```
$ makingBacon_MC_25ns_MINIAOD.py
$ makingBacon_Data_25ns_MINIAOD.py
```

For crab configuration see e.g.:
----------

```
lxplus: /afs/cern.ch/work/p/pharris/public/bacon/prod/CMSSW_8_0_20/src/BaconProd/Ntupler/crab
lpc: /uscms_data/d3/cmantill/bacon/CMSSW_9_4_0_patch1/src/BaconProd/Ntupler/crab
```

For running a list of samples e.g. mc.txt:
```
./runList.sh mc.txt
```

Modify `run.sh` and `crab_template*.py` according to your preferences.
