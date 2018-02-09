BaconProd
=========

Package for producing bacon files

 * Depends on BaconAna package
 * Place package in `$CMSSW_BASE/src` area

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

Production is done through the following config scripts:
----------

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

