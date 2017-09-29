#pragma once

#include <FitBase.hpp>

#include <vector>


/**
 * \class PhotonJet
 * \brief Implements computation of deviation of data from expectation in the photon + jet analysis
 * 
 * At the moment this is only a template for the actual implementation.
 */
class PhotonJet: public DeviationBase
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
    PhotonJet(std::string const &fileName, Method method);
    
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
    /// Input data in bins of photon pt
    std::vector<PtBin> bins;
};
