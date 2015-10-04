#!/bin/bash

if test -z $CMSSW_VERSION; then
  echo "[BaconProd] Need CMSSW project area setup!";
  echo
  return 0;
fi

CURRDIR=$PWD
PATCHDIR=/afs/cern.ch/work/k/ksung/public/Development/Run2Packages
cd $CMSSW_BASE/src

cp -r ${PATCHDIR}/* ./

echo
echo "[BaconProd] Setup complete!"

cd $CURRDIR

return 1;
