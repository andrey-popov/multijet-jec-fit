#pragma once

#include <FitBase.hpp>

#include <Nuisances.hpp>

#include <set>
#include <vector>


/**
 * \class PhotonJetRun1
 * \brief Implements computation of deviation of data from expectation in the photon + jet analysis
 * 
 * The deviation is computed as a chi^2 distance with
 *   chi^2 = sum_i (B^{Data}_i / B^{Sim}_i - 1 / corr(p_i))^2 / sigma_i^2,
 * where B_i is the mean balance observable in bin i, p_i is the mean pt of the photon, and sigma_i
 * is the statistical uncertainty on the ratio. The mean balance observables are extrapolated to
 * alpha = 0.
 * 
 * Changes of photon pt scale in data are propagated into the ratio of balance observables and the
 * pt of the photon.
 */
class PhotonJetRun1: public MeasurementBase
{
public:
    /// Supported methods of computation
    enum class Method
    {
        PtBal,
        MPF
    };
    
private:
    /// Auxiliary structure to aggregate information about a single bin in photon pt
    struct PtBin
    {
        /// Mean photon pt
        double ptPhoton;
        
        /// Ratio between balance observables in data and simulation
        double balanceRatio;
        
        /// (Squared) statistical uncertainty on balanceRatio
        double unc2;
    };
    
public:
    /// Constructor
    PhotonJetRun1(std::string const &fileName, Method method, NuisanceDefinitions &nuisanceDefs);
    
public:
    /**
     * \brief Returns dimensionality of the deviation
     * 
     * Implemented from MeasurementBase.
     */
    virtual unsigned GetDim() const override;
    
    /**
     * \brief Evaluates the deviation with the given jet corrector and set of nuisances
     * 
     * Implemented from MeasurementBase.
     */
    virtual double Eval(JetCorrBase const &corrector, Nuisances const &nuisances) const override;
    
private:
    /// Input data in bins of photon pt
    std::vector<PtBin> bins;
    
    /**
     * \brief Size of the variation in the photon pt scale
     * 
     * This is the relative change in the pt scale when the value of the corresponding nuisance
     * parameter is +1.
     */
    double photonScaleVar;
};
