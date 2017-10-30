#pragma once

#include <FitBase.hpp>

#include <vector>


/**
 * \class ZJetRun1
 * \brief Implements computation of deviation of data from expectation in the Z + jet analysis
 * 
 * The deviation is computed as a chi^2 distance with
 *   chi^2 = sum_i (B^{Data}_i / B^{Sim}_i - 1 / corr(p_i))^2 / sigma_i^2,
 * where B_i is the mean balance observable in bin i, p_i is the mean pt of the Z boson, and
 * sigma_i is the statistical uncertainty on the ratio. The mean balance observables are
 * extrapolated to alpha = 0.
 */
class ZJetRun1: public DeviationBase
{
public:
    /// Supported methods of computation
    enum class Method
    {
        PtBal,
        MPF
    };
    
private:
    /// Auxiliary structure to aggregate information about a single bin in pt of Z boson
    struct PtBin
    {
        /// Mean pt of Z boson
        double ptZ;
        
        /// Ratio between balance observables in data and simulation
        double balanceRatio;
        
        /// (Squared) statistical uncertainty on balanceRatio
        double unc2;
    };
    
public:
    /// Constructor
    ZJetRun1(std::string const &fileName, Method method);
    
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
    /// Input data in bins of pt of Z boson
    std::vector<PtBin> bins;
};
