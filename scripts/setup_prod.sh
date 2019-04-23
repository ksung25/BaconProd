#!/bin/bash

if test -z $CMSSW_VERSION; then
  echo "[BaconProd] Need CMSSW project area setup!";
  echo
  return 0;
fi

git cms-merge-topic cmantill:baconprod-10213-v15

echo
echo "[BaconProd] Setup complete!"

cd $CURRDIR

return 1;
