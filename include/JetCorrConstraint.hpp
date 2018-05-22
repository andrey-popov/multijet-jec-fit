#pragma once

#include <FitBase.hpp>


/**
 * \class JetCorrConstraint
 * \brief Fake measurement to constrain jet correction at given pt
 * 
 * This fake measurement sets, within an uncertainty, the value for the correction at given
 * pt scale. It can be useful to fit correction in the multijet analysis alone as it is not very
 * sensitive to the absolute scale of the correction.
 */
class JetCorrConstraint: public MeasurementBase
{
public:
    /// Constructor
    JetCorrConstraint();
    
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
    virtual double Eval(JetCorrBase const &corrector, Nuisances const &) const override;
    
    /// Sets new parameters for the constraint
    void SetParameters(double ptRef, double targetCorrection, double relUncertainty);
    
private:
    /// Transverse momentum at which the constraint is set, GeV
    double ptRef;
    
    /// Target correction at given pt scale
    double targetCorrection;
    
    /// Relative uncertainty on the correction
    double relUncertainty;
};
