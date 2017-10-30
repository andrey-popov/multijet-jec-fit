#pragma once

#include <FitBase.hpp>

#include <TH1.h>
#include <TH2.h>
#include <TProfile.h>
#include <TProfile2D.h>

#include <vector>


/**
 * \class NewPhotonJet
 * \brief Implements computation of deviation of data from expectation in the photon + jet analysis
 * 
 * The deviation is computed as a chi^2 distance,
 *   chi^2 = sum_i (B^{Data}_i - B^{Sim}_i)^2 / (sigma^{Data}_i^2 + sigma^{Sim}_i^2),
 * where B_i is the mean balance observable in bin i in pt of the photon and sigma_i is its
 * statistical uncertainty. In data the mean balance observable is recomputed for the given jet
 * correction following an approach similar to the multijet analysis.
 * 
 * Changes of photon pt scale in data are propagated into the pt of the photon.
 */
class NewPhotonJet: public DeviationBase
{
public:
    /// Supported methods of computation
    enum class Method
    {
        PtBal,
        MPF
    };
    
public:
    /// Constructor
    NewPhotonJet(std::string const &fileName, Method method);
    
public:
    /**
     * \brief Returns dimensionality of the deviation
     * 
     * Implemented from DeviationBase.
     */
    virtual unsigned GetDim() const override;
    
    /**
     * \brief Evaluates the deviation with the given jet corrector and set of nuisances
     * 
     * Implemented from DeviationBase.
     */
    virtual double Eval(JetCorrBase const &corrector, Nuisances const &nuisances) const override;
    
private:
    /// Profiles of the balance observable in data and simulation
    std::unique_ptr<TProfile> balProfile, simBalProfile;
    
    /// Distribution of the pt of the photon in data
    std::unique_ptr<TH1> ptPhoton;
    
    /// Profile of the pt of the photon in data
    std::unique_ptr<TProfile> ptPhotonProfile;
    
    /// Sum of projections of pt of jets in bins of pt of the photon and jets
    std::unique_ptr<TH2> ptJetSumProj;
    
    /// 2D profile of pt of jets
    std::unique_ptr<TProfile2D> ptJet2DProfile; 
    
    /**
     * \brief Squared uncertainty on the difference between mean balance observables in data
     * and simulation
     */
    std::vector<double> totalUnc2;
};
