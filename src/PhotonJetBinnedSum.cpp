#include <PhotonJetBinnedSum.hpp>
#include <Rebin.hpp>

#include <cmath>
#include <memory>
#include <sstream>
#include <TVectorD.h>

#include <TFile.h>


PhotonJetBinnedSum::PhotonJetBinnedSum(std::string const &fileName,
  PhotonJetBinnedSum::Method method_):
    method(method_)
{
    std::string methodLabel;
    
    if (method == Method::PtBal)
        methodLabel = "Bal";
    else if (method == Method::MPF)
        methodLabel = "MPF";
    
    
    std::unique_ptr<TFile> inputFile(TFile::Open(fileName.c_str()));
    
    if (not inputFile or inputFile->IsZombie())
    {
        std::ostringstream message;
        message << "PhotonJetBinnedSum::PhotonJetBinnedSum: Failed to open file \"" <<
          fileName << "\".";
        throw std::runtime_error(message.str());
    }
    
    auto ptThreshold = dynamic_cast<TVectorD *>(inputFile->Get(("MC_MinPt" + methodLabel).c_str()));

    if (not ptThreshold or ptThreshold->GetNoElements() != 1)
     {
        std::ostringstream message;
        message << "PhotonJetBinnedSum::PhotonJetBinnedSum: Failed to read jet pt threshold " <<
           "from file \"" << fileName << "\".";
        throw std::runtime_error(message.str());
      }

    jetPtMin = (*ptThreshold)[1];
  
    simBalProfile.reset(dynamic_cast<TProfile *>(inputFile->Get(
      ("MC_new" + methodLabel + "_vs_ptphoton").c_str())));
    balProfile.reset(dynamic_cast<TProfile *>(inputFile->Get(
      ("DATA_new" + methodLabel + "_vs_ptphoton").c_str())));
    ptPhoton.reset(dynamic_cast<TH1 *>(inputFile->Get("DATA_phopt_for_nevts")));
    ptPhotonProfile.reset(dynamic_cast<TProfile *>(inputFile->Get("DATA_ptphoton_vs_ptphoton")));
    ptJetSumProj.reset(dynamic_cast<TH2 *>(inputFile->Get("DATA_Skl_phopt_vs_jetpt")));
    ptJet2DProfile.reset(dynamic_cast<TProfile2D *>(inputFile->Get("DATA_jetpt_phopt_vs_jetpt")));
    
    
    simBalProfile->SetDirectory(nullptr);
    balProfile->SetDirectory(nullptr);
    ptPhoton->SetDirectory(nullptr);
    ptPhotonProfile->SetDirectory(nullptr);
    ptJetSumProj->SetDirectory(nullptr);
    ptJet2DProfile->SetDirectory(nullptr);
    
    inputFile->Close();
    
    
    // Compute combined (squared) uncertainty on the balance observable in data and simulation.
    //The data profile is rebinned with the binning used for simulation. This is done assuming that
    //bin edges of the two binnings are aligned, which should normally be the case.
    std::unique_ptr<TH1> balRebinned(balProfile->Rebin(simBalProfile->GetNbinsX(), "",
      simBalProfile->GetXaxis()->GetXbins()->GetArray()));
    
    for (int i = 1; i <= simBalProfile->GetNbinsX(); ++i)
    {
        double const unc2 = std::pow(simBalProfile->GetBinError(i), 2) +
          std::pow(balRebinned->GetBinError(i), 2);
        totalUnc2.emplace_back(unc2);
    }
    
    
    recompBal.resize(simBalProfile->GetNbinsX());
}


unsigned PhotonJetBinnedSum::GetDim() const
{
    return simBalProfile->GetNbinsX();
}


double PhotonJetBinnedSum::Eval(JetCorrBase const &corrector, Nuisances const &nuisances) const
{
    UpdateBalance(corrector, nuisances);
    double chi2 = 0.;
    
    for (int photonBinIndex = 1; photonBinIndex <= simBalProfile->GetNbinsX(); ++photonBinIndex)
    {
        double const meanBal = recompBal[photonBinIndex - 1];
        double const simMeanBal = simBalProfile->GetBinContent(photonBinIndex);
        chi2 += std::pow(meanBal - simMeanBal, 2) / totalUnc2[photonBinIndex - 1];
    }
    
    return chi2;
}


double PhotonJetBinnedSum::ComputeMPF(FracBin const &ptPhotonStart, FracBin const &ptPhotonEnd,
  JetCorrBase const &corrector, Nuisances const &nuisances) const
{
    
    // Find the bin in jet pt that includes the value of pt that, after the current correction,
    // would give the nominal minimal pt threshold. Compute also the fraction of this bin that
    // should included in the sum.
    TAxis const *ptJetAxis = ptJetSumProj->GetYaxis();
    double const uncorrJetPtMin = corrector.UndoCorr(jetPtMin);
    int const startBin = ptJetAxis->FindBin(uncorrJetPtMin);
    double const fracStartBin = 1. - (uncorrJetPtMin - ptJetAxis->GetBinLowEdge(startBin)) /
      ptJetAxis->GetBinWidth(startBin);
    
    
    double sumBal = 0., sumWeight = 0.,  sumJets = 0.;
    
    for (unsigned photonBinIndex = ptPhotonStart.index; photonBinIndex <= ptPhotonEnd.index;
      ++photonBinIndex)
    {
        double const numEvents = ptPhoton->GetBinContent(photonBinIndex);
        
        if (numEvents == 0)
            continue;
        
        double const meanPhotonPt = ptPhotonProfile->GetBinContent(photonBinIndex) *
          (1 + nuisances.photonScale);
        
        
        // Recompute mean value for the MPF observable in data by summing over all jet pt bins
        sumBal += balProfile->GetBinContent(photonBinIndex) * numEvents ;
        sumWeight += numEvents ;
        
        sumJets = 0.;
        
        for (int jetBinIndex = startBin; jetBinIndex <= ptJetSumProj->GetNbinsY(); ++jetBinIndex)
        {
            double const s = ptJetSumProj->GetBinContent(photonBinIndex, jetBinIndex);
            
            if (s == 0.)
                continue;
            
            double const meanJetPt = ptJet2DProfile->GetBinContent(photonBinIndex, jetBinIndex);
            
            if (jetBinIndex == startBin)
                sumJets -= s * (1. - corrector.Eval(meanJetPt)) * fracStartBin;
            else
                sumJets -= s * (1. - corrector.Eval(meanJetPt));
        }
        
        sumJets /= meanPhotonPt;
        sumBal += sumJets;
    }
    
    return sumBal / sumWeight;
}


double PhotonJetBinnedSum::ComputePtBal(FracBin const &ptPhotonStart, FracBin const &ptPhotonEnd,
 JetCorrBase const &corrector, Nuisances const &nuisances) const
{
    
    // Find the bin in jet pt that includes the value of pt that, after the current correction,
    // would give the nominal minimal pt threshold. Compute also the fraction of this bin that
    // should included in the sum.
    TAxis const *ptJetAxis = ptJetSumProj->GetYaxis();
    double const uncorrJetPtMin = corrector.UndoCorr(jetPtMin);
    int const startBin = ptJetAxis->FindBin(uncorrJetPtMin);
    double const fracStartBin = 1. - (uncorrJetPtMin - ptJetAxis->GetBinLowEdge(startBin)) /
      ptJetAxis->GetBinWidth(startBin);
    
    
    // Recompute mean value for the balance observable in data by summing over all jet pt bins
    double meanBal = 0.;
    double sumWeight = 0.;
    
    for (unsigned photonBinIndex = ptPhotonStart.index; photonBinIndex <= ptPhotonEnd.index;
      ++photonBinIndex)
    {
        double const numEvents = ptPhoton->GetBinContent(photonBinIndex);
        
        if(numEvents == 0)
            continue;
        
        sumWeight += numEvents;
        double meanBalInBin = 0.;
        
        double const meanPhotonPt = ptPhotonProfile->GetBinContent(photonBinIndex) *
          (1 + nuisances.photonScale);
        
        for (int jetBinIndex = startBin; jetBinIndex <= ptJetSumProj->GetNbinsY(); ++jetBinIndex)
        {
            double const s = ptJetSumProj->GetBinContent(photonBinIndex, jetBinIndex);
            
            if (s == 0.)
                continue;
            
            double const meanJetPt = ptJet2DProfile->GetBinContent(photonBinIndex, jetBinIndex);
            
            if(jetBinIndex == startBin)
                meanBalInBin += s * corrector.Eval(meanJetPt) * fracStartBin;
            else
                meanBalInBin += s * corrector.Eval(meanJetPt);
        }
    
        meanBalInBin /= meanPhotonPt;
        meanBal += meanBalInBin;
    }
    
    meanBal /= sumWeight;
    return meanBal;
}

void PhotonJetBinnedSum::UpdateBalance(JetCorrBase const &corrector, Nuisances const &nuisances)
  const
{
    std::vector<double> simPtBinning;
    std::vector<double> dataPtBinning;
    
    for (int i = 1; i <= simBalProfile->GetNbinsX() + 1; ++i)
    {
        double const pt = simBalProfile->GetBinLowEdge(i);
        simPtBinning.emplace_back(pt);
    }
    
    for (int i = 1; i <= balProfile->GetNbinsX() + 1; ++i)
    {
        double const pt = balProfile->GetBinLowEdge(i);
        dataPtBinning.emplace_back(pt);
    }
    
    
    // Build a map from the simulation (wide) binning to the fine binning used in data
    auto binMap = mapBinning(dataPtBinning, simPtBinning);
    binMap.erase(0);
    binMap.erase(simBalProfile->GetNbinsX() + 1);
    
    for (auto const &binMapPair: binMap)
    {
        auto const &binIndex = binMapPair.first;
        auto const &binRange = binMapPair.second;
        
        double meanBal;
        
        if (method == Method::PtBal)
            meanBal = ComputePtBal(binRange[0], binRange[1], corrector, nuisances);
        else
            meanBal = ComputeMPF(binRange[0], binRange[1], corrector, nuisances);
        
        recompBal[binIndex - 1] = meanBal;
    }
}
