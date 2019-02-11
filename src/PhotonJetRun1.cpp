#include <PhotonJetRun1.hpp>

#include <cmath>
#include <memory>
#include <sstream>

#include <TFile.h>
#include <TGraphErrors.h>

using namespace std::string_literals;


PhotonJetRun1::PhotonJetRun1(std::string const &fileName, Method method,
  NuisanceDefinitions &nuisanceDefs):
    photonScaleVar(0.01)  // a dummy value
{
    std::string methodLabel;
    
    if (method == Method::PtBal)
        methodLabel = "PtBal";
    else if (method == Method::MPF)
        methodLabel = "MPF";
    
    std::unique_ptr<TFile> inputFile(TFile::Open(fileName.c_str()));
    
    if (not inputFile or inputFile->IsZombie())
    {
        std::ostringstream message;
        message << "PhotonJetRun1::PhotonJetRun1: Failed to open file \"" << fileName << "\".";
        throw std::runtime_error(message.str());
    }
    
    std::unique_ptr<TGraphErrors> extrapRatio(dynamic_cast<TGraphErrors *>(
      inputFile->Get(("resp_"s + methodLabel + "chs_extrap_a30_eta00_13").c_str())));
    
    inputFile->Close();
    
    
    bins.reserve(extrapRatio->GetN());
    
    for (int i = 0; i < extrapRatio->GetN(); ++i)
    {
        double x, y;
        extrapRatio->GetPoint(i, x, y);
        
        PtBin bin;
        bin.ptPhoton = x;
        bin.balanceRatio = y;
        bin.unc2 = pow(extrapRatio->GetErrorY(i), 2);
        
        bins.push_back(bin);
    }


    nuisanceDefs.Register("PhotonScale");
}


unsigned PhotonJetRun1::GetDim() const
{
    return bins.size();
}


double PhotonJetRun1::Eval(JetCorrBase const &corrector, Nuisances const &nuisances) const
{
    double chi2 = 0.;
    
    for (auto const &bin: bins)
    {
        // Correct the balance ratio and photon pt for the potential offset in the photon pt scale
        double const photonScaleFactor = 1 + photonScaleVar * nuisances["PhotonScale"];
        double const balanceRatioCorr = bin.balanceRatio / photonScaleFactor;
        double const ptPhoton = bin.ptPhoton * photonScaleFactor;
        
        // Assume that pt of the jet is the same as pt of the photon
        chi2 += std::pow(balanceRatioCorr - 1 / corrector.Eval(ptPhoton), 2) / bin.unc2;
    }
    
    return chi2;
}

