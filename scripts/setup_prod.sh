#!/bin/bash

if test -z $CMSSW_VERSION; then
  echo "[BaconProd] Need CMSSW project area setup!";
  echo
  return 0;
fi

if test -z "$CVSROOT"; then
  echo "[BaconProd] Need to set CVSROOT!"
  echo "[BaconProd] Something like:"
  echo "            export CVSROOT=:ext:<cern-user-account>@lxplus.cern.ch:/afs/cern.ch/user/c/cvscmssw/public/CMSSW"
  echo
  return 0;
fi


CURRDIR=$PWD
PATCHDIR=$CMSSW_BASE/src/BaconProd/patch
cd $CMSSW_BASE/src


### MuScleFit corrections for muons
echo
echo "[BaconProd] Checking out MuScleFit package..."
cvs co -r muscle_v4_2_0 UserCode/scasasso/MuScleFit/Calibration
mv UserCode/scasasso/MuScleFit ./
rm -rf UserCode


### Electron MVA ID
echo
echo "[BaconProd] Checking out Electron MVA ID package..."
cvs co -r V00-00-30-02 UserCode/EGamma/EGammaAnalysisTools
mv UserCode/EGamma ./
rm -rf UserCode
mv EGamma/EGammaAnalysisTools/test/BuildFile.xml EGamma/EGammaAnalysisTools/test/BuildFile.xmlSilent
mv EGamma/EGammaAnalysisTools/plugins/BuildFile.xml EGamma/EGammaAnalysisTools/plugins/BuildFile.xmlSilent
cp $PATCHDIR/EGamma/EGammaAnalysisTools/interface/EGammaMvaEleEstimator.h.53Xpatch EGamma/EGammaAnalysisTools/interface/EGammaMvaEleEstimator.h


### MET filters
echo
echo "[BaconProd] Checking out packages for MET filters..."
cvs co -r V00-00-13-01 RecoMET/METFilters
cvs co -r V00-00-08    RecoMET/METAnalyzers


### Jet/MET packages
echo
echo "[BaconProd] Checking out jet/MET packages..."
git clone https://github.com/nhanvtran/JetTools.git
cp $PATCHDIR/JetTools/AnalyzerToolbox/python/njettinessadder_cfi.py JetTools/AnalyzerToolbox/python/

cp -r /afs/cern.ch/work/p/pharris/public/tmp/CMSSW_5_3_13/src/RecoMET/METPUSubtraction RecoMET
cp -r /afs/cern.ch/work/p/pharris/public/tmp/CMSSW_5_3_13/src/DataFormats .
cp -r /afs/cern.ch/work/p/pharris/public/tmp/CMSSW_5_3_13/src/RecoJets .
cp -r /afs/cern.ch/work/p/pharris/public/tmp/CMSSW_5_3_13/src/CondFormats .

git clone https://github.com/violatingcp/Jets_Short.git
mv Jets_Short/RecoJets/JetProducers/data/*.xml RecoJets/JetProducers/data/
rm -rf Jets_Short

cp $PATCHDIR/RecoMET/METPUSubtraction/python/mvaPFMET_leptons_cff.py RecoMET/METPUSubtraction/python/
cp $PATCHDIR/RecoMET/METPUSubtraction/data/*Sep*.root                RecoMET/METPUSubtraction/data/

### clean up large unncessary files to fit into CRAB input sandbox (100MB)
rm -f RecoJets/JetProducers/data/TMVAClassificationCategory_JetID_53X_chs_Dec2012.weights.xml.~1.1.~
rm -f JetTools/AnalyzerToolbox/data/TMVAClassificationCategory_JetID_53X_Dec2012.weights.xml
rm -f JetTools/AnalyzerToolbox/data/TMVAClassification_5x_BDT_chsSimpleNoVtxCat.weights.xml
rm -f JetTools/AnalyzerToolbox/data/TMVAClassification_5x_BDT_simpleNoVtxCat.weights.xml
rm -f JetTools/AnalyzerToolbox/data/TMVAClassification_PuJetIdOptMVA.weights.xml

### remove any compiled python scrips from copied directories...
rm -f */*/python/*.pyc

echo
echo "[BaconProd] Setup complete!"

cd $CURRDIR

return 1;
