#include <PhotonJet.hpp>

#include <cmath>

#include <TFile.h>

#include <iostream>

using namespace std::string_literals;


PhotonJet::PhotonJet(std::string const &fileName, Method method)
{
    // TODO: Read input file and fill the vector bins with information from it. Whether inputs for
    // pt balancing or MPF are read, is controlled by the second argument. If inputs are provided
    // for different values of alpha instead of ones extrapolated to alpha -> 0, the extrapolation
    // should be performed here.


    std::string methodLabel;
    
    if (method == Method::PtBal)
        methodLabel = "PtBal";
    else if (method == Method::MPF)
        methodLabel = "MPF";

    TFile inputFile(fileName.c_str());
    
    if (inputFile.IsZombie())
    {
        std::ostringstream message;
        message << "Failed to open file \"" << fileName << "\".";
        throw std::runtime_error(message.str());
    }

ExtrapRatio.reset(dynamic_cast<TGraphErrors *>(inputFile.Get(("resp_"s + methodLabel + "chs_extrap_a30_eta00_13").c_str())));

double xtmp, ytmp;

PtBin bin;

for (int i=1; i<(ExtrapRatio.get())->GetN(); i++) {
//get values from the graph
(ExtrapRatio.get())->GetPoint(i, xtmp, ytmp);
//fill PtBin structure with values
bin.ptPhoton = xtmp;
bin.balanceRatio = ytmp;
bin.unc2 = (pow((ExtrapRatio.get())->GetErrorY(i),2));
//push back in the vector of PtBins
if((ExtrapRatio.get())->GetErrorY(i)>0.) bins.push_back(bin);
}


}


unsigned PhotonJet::GetDim() const
{
    return bins.size();
}


double PhotonJet::Eval(JetCorrBase const &corrector, Nuisances const &nuisances) const
{
    double chi2 = 0.;
    
    for (auto const &bin: bins)
    {
        // Correct the balance ratio for the potential offset in the photon pt scale
        double const balanceRatioCorr = bin.balanceRatio / (1 + nuisances.photonScale);
        
        // Assume that pt of the jet is the same as pt of the photon
        chi2 += std::pow(balanceRatioCorr - 1/(corrector.Eval(bin.ptPhoton*(1 + nuisances.photonScale))), 2) / bin.unc2;
    }
    
    return chi2;
}
