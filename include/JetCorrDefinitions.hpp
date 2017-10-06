/**
 * \file JetCorrDefinitions.hpp
 * 
 * This file defines standard forms for jet corrections.
 */


#pragma once

#include <FitBase.hpp>

#include <array>


/**
 * \class JetCorrStableLogLin
 * \brief Log-linear correction that smoothly approaches unity given minimal pt threshold
 * 
 * The correction is unity at and below a given pt threshold and proportional to log(pt) for large
 * values of pt. It transitions smoothely between the two regimes. A nice property of this
 * correction is that jets cannot migrate through the threshold in any direction. The only free
 * parameter determines the slope of the log-linear function for large values of pt. This
 * correction is useful for tests of the multijet method. Introduced in [1].
 * [1] https://indico.cern.ch/event/646599/contributions/2627202/attachments/1476692/2287829/2017.06.14_JERC_Balance-recomputation_v1.pdf#9
 */
class JetCorrStableLogLin: public JetCorrBase
{
public:
    /// Construct from a pt threshold
    JetCorrStableLogLin(double ptMin = 15.);
    
public:
    /**
     * \brief Computes correction for a jet with given pt
     * 
     * Implemented from JetCorrBase.
     */
    virtual double Eval(double pt) const override;
    
private:
    /// Threshold below which the correction is unity
    double ptMin;
};


/**
 * \class JetCorrStd2P
 * \brief Standard two-parameter correction.
 * 
 * Standard two-parameter correction used with 2016 data. Functional form is taken from [1], with
 * the exception that a unity is added to the correction, which effectively redefines the first
 * parameter. Because of this, setting all free parameters to zero results in a unity correction.
 * 
 * In addition to the two free parameters, three additional parameters configure the single-pion
 * response function, which is used as a building block for the correction. These additional
 * parameters are set by dedicated method and cannot float in the fit.
 * 
 * [1] https://github.com/miquork/jecsys/blob/Summer16_23Sep2016V4/globalFitL3Res.C#L128
 */
class JetCorrStd2P: public JetCorrBase
{
public:
    /**
     * \brief Constructor
     * 
     * Fixed parameters for single-pion response are initialized to values used in JEC on version
     * Summer16_03Feb2017H_V3.
     */
    JetCorrStd2P();
    
public:
    /**
     * \brief Computes correction for a jet with given pt
     * 
     * Implemented from JetCorrBase.
     */
    virtual double Eval(double pt) const override;
    
    /// Sets parameters of the single-pion response
    void SetParamsSPR(std::initializer_list<double> paramsSPR);
    
protected:
    /**
     * \brief Single-pion response in HCAL [1]
     * 
     * A building block of the full correction.
     * [1] https://github.com/miquork/jecsys/blob/Summer16_23Sep2016V4/globalFitL3Res.C#L85
     */
    double fSPR(double pt) const;
    
protected:
    /// Reference pt scale
    double ptRef;
    
    /// (Fixed) parameters for single-pion response
    std::array<double, 3> paramsSPR;
};


/**
 * \class JetCorrStd3P
 * \brief Standard three-parameter correction.
 * 
 * Standard three-parameter correction used with 2016 data. Functional form is taken from [1], with
 * the exception that a unity is added to the correction, which effectively redefines the first
 * parameter.
 * 
 * Implemented by extending class JetCorrStd2P. The correction is evaluated using the single-pion
 * response function defined in that class and also adding a contribution that reflects the
 * difference between L1FastJet and L1RC. It is described with two additional parameters, which
 * can be set using dedicated method.
 * 
 * [1] https://github.com/miquork/jecsys/blob/Summer16_23Sep2016V4/globalFitL3Res.C#L137-L138
 */
class JetCorrStd3P: public JetCorrStd2P
{
public:
    /**
     * \brief Constructor
     * 
     * Fixed parameters for building blocks of the correction are initialized to values used in
     * JEC of version Summer16_03Feb2017H_V3.
     */
    JetCorrStd3P();
    
public:
    /**
     * \brief Computes correction for a jet with given pt
     * 
     * Reimplemented from JetCorrStd2P.
     */
    virtual double Eval(double pt) const override;
    
    /// Sets parameters related to L1 corrections
    void SetParamsL1(std::initializer_list<double> paramsL1);
    
private:
    /**
     * \brief Difference between L1FastJet and L1RC (whatever it means).
     * 
     * A building block of the full correction. Functional form is taken from [1].
     * [1] https://github.com/miquork/jecsys/blob/Summer16_23Sep2016V4/globalFitL3Res.C#L88
     */
    double fL1(double pt) const;
    
private:
    /// (Fixed) parameters describing L1 corrections
    std::array<double, 2> paramsL1;
};
