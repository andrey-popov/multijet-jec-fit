#include <NewPhotonJet.hpp>

#include <cmath>
#include <memory>

#include <TFile.h>

using namespace std::string_literals;


NewPhotonJet::NewPhotonJet(std::string const &fileName, Method method)
{
    if (method == Method::MPF)
    {
        std::ostringstream message;
        message << "NewPhotonJet:: NewPhotonJet: MPF method is not implemented yet.";
        throw std::runtime_error(message.str());
    }
    
    
    std::string methodLabel;

    if (method == Method::PtBal)
        methodLabel = "Bal";
    else if (method == Method::MPF)
        methodLabel = "MPF";
    
    TFile inputFile(fileName.c_str());
    
    if (inputFile.IsZombie())
    {
        std::ostringstream message;
        message << "NewPhotonJet:: NewPhotonJet: Failed to open file \"" << fileName << "\".";
        throw std::runtime_error(message.str());
    }
    
    
    simBalProfile.reset(dynamic_cast<TProfile *>(inputFile.Get(
      ("MC_new"s + methodLabel + "_vs_ptphoton").c_str())));
    balProfile.reset(dynamic_cast<TProfile *>(inputFile.Get(
      ("DATA_new"s + methodLabel + "_vs_ptphoton").c_str())));
    ptPhoton.reset(dynamic_cast<TH1 *>(inputFile.Get("DATA_phopt_for_nevts")));
    ptPhotonProfile.reset(dynamic_cast<TProfile *>(inputFile.Get("DATA_ptphoton_vs_ptphoton")));
    ptJetSumProj.reset(dynamic_cast<TH2 *>(inputFile.Get("DATA_Skl_phopt_vs_jetpt")));
    ptJet2DProfile.reset(dynamic_cast<TProfile2D *>(inputFile.Get("DATA_jetpt_phopt_vs_jetpt")));

    simBalProfile->SetDirectory(nullptr);
    balProfile->SetDirectory(nullptr);
    ptPhoton->SetDirectory(nullptr);
    ptPhotonProfile->SetDirectory(nullptr);
    ptJetSumProj->SetDirectory(nullptr);
    ptJet2DProfile->SetDirectory(nullptr);
    
    inputFile.Close();
    
    
    for (int i = 1; i <= simBalProfile->GetNbinsX(); ++i)
    {
        double const unc2 = std::pow(simBalProfile->GetBinError(i), 2) +
            std::pow(balProfile->GetBinError(i), 2);
        totalUnc2.emplace_back(unc2);
    }
}


unsigned NewPhotonJet::GetDim() const
{
    return simBalProfile->GetNbinsX();
}


double NewPhotonJet::Eval(JetCorrBase const &corrector, Nuisances const &nuisances) const
{
    double chi2 = 0.;
    double const jetPtMin = 30.;
    
    for (int photonBinIndex = 1; photonBinIndex <= simBalProfile->GetNbinsX(); ++photonBinIndex)
    {
        double const numEvents = ptPhoton->GetBinContent(photonBinIndex);
        
        if (numEvents == 0.)
            continue;
        
        double const meanPhotonPt =
          ptPhotonProfile->GetBinContent(photonBinIndex) * (1 + nuisances.photonScale);
        
        
        // Find the bin in jet pt that includes the value of pt that, after the current correction,
        //would give the nominal minimal pt threshold. Compute also the fraction of this bin that
        //should included in the sum.
        TAxis const *ptJetAxis = ptJetSumProj->GetYaxis();
        double const uncorrJetPtMin = corrector.UndoCorr(jetPtMin);
        int const startBin = ptJetAxis->FindBin(uncorrJetPtMin);
        double const fracStartBin = 1. - (uncorrJetPtMin - ptJetAxis->GetBinLowEdge(startBin)) /
          ptJetAxis->GetBinWidth(startBin);
        
        
        // Recompute mean value for the balance observable in data by summing over all jet pt bins
        double meanBal = 0.;
        
        for (int jetBinIndex = startBin; jetBinIndex <= ptJetSumProj->GetNbinsY(); ++jetBinIndex)
        {
            double const s = ptJetSumProj->GetBinContent(photonBinIndex, jetBinIndex);
            
            if (s == 0.)
                continue;
            
            double const meanJetPt = ptJet2DProfile->GetBinContent(photonBinIndex, jetBinIndex);
            
            if(jetBinIndex == startBin)
                meanBal += s * corrector.Eval(meanJetPt) * fracStartBin;
            else
                meanBal += s * corrector.Eval(meanJetPt);
        }
        
        meanBal /= meanPhotonPt * numEvents;
        
        
        double const simMeanBal = simBalProfile->GetBinContent(photonBinIndex);
        chi2 += std::pow(meanBal - simMeanBal, 2) / totalUnc2[photonBinIndex - 1];
    }
    
    return chi2;
}
