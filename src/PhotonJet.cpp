#include <PhotonJet.hpp>

#include <cmath>


PhotonJet::PhotonJet(std::string const &fileName, Method method)
{
    // TODO: Read input file and fill the vector bins with information from it. Whether inputs for
    // pt balancing or MPF are read, is controlled by the second argument. If inputs are provided
    // for different values of alpha instead of ones extrapolated to alpha -> 0, the extrapolation
    // should be performed here.
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
        chi2 += std::pow(balanceRatioCorr - corrector.Eval(bin.ptPhoton), 2) / bin.unc2;
    }
    
    return chi2;
}
