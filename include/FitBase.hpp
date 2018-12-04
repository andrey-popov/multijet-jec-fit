#pragma once

#include <Nuisances.hpp>

#include <memory>
#include <ostream>
#include <set>
#include <vector>


/**
 * \class JetCorrBase
 * \brief Base class to for a jet correction
 * 
 * The correction is multiplicative. Its parameters are stored as data members.
 */
class JetCorrBase
{
public:
    /// Constructor from the number of parameters
    JetCorrBase(unsigned numParams);
    
    /// Trivial virtual destructor
    virtual ~JetCorrBase() noexcept;
    
public:
    /// Applies the correction to the given jet pt, returning the corrected pt
    double Apply(double pt) const;
    
    /// Returns number of parameters of the correction
    unsigned GetNumParams() const;
    
    /// Returns set values for parameters of the correction
    std::vector<double> const &GetParams() const;
    
    /**
     * \brief Evaluates the correction for the given jet pt
     * 
     * To be implemented in a derived class.
     */
    virtual double Eval(double pt) const = 0;
    
    /**
     * \brief Updates parameters of the correction
     * 
     * Throws an exception if the given number of parameters does not match the expected one.
     */
    void SetParams(std::vector<double> const &newParams);
    
    /**
     * \brief Updates parameters of the correction
     * 
     * Reads a number GetNumParams() parameters from the buffer. The version that accepts an
     * std::vector should usually be preferred.
     */
    void SetParams(double const *newParams);
    
    /**
     * \brief Inverts jet correction
     * 
     * Returns uncorrected pt such that ptUncorr * corr(ptUncorr) recovers the given corrected pt.
     * The computation is done iteratively and stops when the corrected pt is reproduced with the
     * specified relative tolerance.
     */
    virtual double UndoCorr(double pt, double tolerance = 1e-10) const;
    
protected:
    /// Current parameters of the correction
    std::vector<double> parameters;
};


/// An utility function to print out parameters of a jet correction
std::ostream &operator<<(std::ostream &os, JetCorrBase const &corrector);


/**
 * \class MeasurementBase
 * \brief Base class to describe an analysis
 * 
 * For the purpose of determination of the jet correction, an analysis is represented by a
 * deviation of data (with a hypothesized jet correction applied) from expectation, which is
 * intended to be used to construct the loss function for the fit. A derived class should implement
 * computation of the deviation in a single analysis.
 */
class MeasurementBase
{
public:
    virtual ~MeasurementBase() = default;
    
public:
    /**
     * \brief Returns dimensionality of the deviation
     * 
     * For a chi^2 measure the dimensionality is the number of individual chi^2 terms in the sum.
     * 
     * To be implemented in a derived class.
     */
    virtual unsigned GetDim() const = 0;
    
    /**
     * \brief Returns names of nuisance parameters
     * 
     * In the default implementation, an empty set is returned.
     */
    virtual std::set<std::string> GetNuisances() const;
    
    /**
     * \brief Evaluates the deviation with the given jet corrector and set of nuisances
     * 
     * To be implemented in a derived class.
     */
    virtual double Eval(JetCorrBase const &corrector, Nuisances const &nuisances) const = 0;
};


/**
 * \class CombLossFunction
 * \brief Loss function computed from multiple measurements
 * 
 * Deviations computed for different measurements are summed up resulting in a combined loss
 * function, which can then be minimized to fit parameters of the jet correction. The interface for
 * a minimization package is implemented with method EvalRawInput.
 * 
 * In the context of minimization of the loss function, nuisance parameters are split into two
 * groups. Marginalized nuisances are allowed to float in the fit, and thus they are treated in the
 * same way as the sought-for parameters of the jet correction. Remaining nuisances are fixed as
 * external parameters. They are not modified with EvalRawInput and should instead be set using
 * method SetExternalNuisances. In this implementation no marginalized nuisances are included; they
 * can be added in a derived class.
 */
class CombLossFunction
{
public:
    /**
     * \brief Constructor
     * 
     * JetCorrBase object is owned by this.
     */
    CombLossFunction(std::unique_ptr<JetCorrBase> &&corrector);
    
    /**
     * \brief Constructor
     * 
     * JetCorrBase object is owned by this.
     */
    CombLossFunction(JetCorrBase *corrector);
    
    virtual ~CombLossFunction() = default;
    
public:
    /**
     * \brief Adds a new measurement that will contribute to the loss function
     * 
     * Provided object is not owned by this.
     */
    void AddMeasurement(MeasurementBase const *measurement);
    
    /**
     * \brief Returns the number of degrees of freedom
     * 
     * The number of degrees of freedom is computed as the sum of dimensionality of all included
     * deviations minus the number of parameters to be fitted.
     */
    unsigned GetNDF() const;
    
    /**
     * \brief Returns the number of parameters to be fitted
     * 
     * Computed as the number of parameters of the jet correction plus the number of nuisances to
     * be marginalized in the fit.
     */
    virtual unsigned GetNumParams() const;
    
    /// Returns internal copy of nuisance parameters
    Nuisances const &GetNuisances() const;
    
    /// Returns internal copy of nuisance parameters
    Nuisances &GetNuisances();
    
    /// Wrapper for EvalRawInput that checks the size of the given vector
    double Eval(std::vector<double> const &x) const;
    
    /**
     * \brief Computes the combined loss function for given correction parameters and nuisances
     * 
     * This method is intended to be used outside of the context of the fit, for instance, to
     * perform a scan over some parameters.
     * 
     * Provided values of nuisances are saved by calling method SetExternalNuisances.
     */
    double Eval(std::vector<double> const &corrParams, Nuisances const &nuisances) const;
    
    /**
     * \brief Interface to evaluate the combined loss function in the fit
     * 
     * Evaluates the combined loss function for the given point. The argument is a pointer to an
     * array, which contains values of the parameters of the jet correction followed by
     * marginalized nuisances. Values for other nuisance parameters, which are externalized in the
     * fit, are taken from the internal state set by method SetExternalNuisances.
     * 
     * In this implementation no marginalized nuisances are included, which can be changed in a
     * derived class.
     */
    virtual double EvalRawInput(double const *x) const;
    
protected:
    /// Jet corrector object
    std::unique_ptr<JetCorrBase> corrector;
    
    /// Nuisance parameters
    mutable Nuisances nuisances;
    
    /// Non-owning pointers to individual contributing measurements
    std::vector<MeasurementBase const *> measurements;
};
