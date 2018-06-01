#include <JetCorrDefinitions.hpp>

#include <algorithm>
#include <cmath>
#include <sstream>
#include <stdexcept>


JetCorrStableLogLin::JetCorrStableLogLin(double ptMin_):
    JetCorrBase(1),
    ptMin(ptMin_)
{}


double JetCorrStableLogLin::Eval(double pt) const
{
    double const b = 1.;
    return 1. + parameters[0] * std::log(pt / ptMin) +
      parameters[0] / b * (std::pow(pt / ptMin, -b) - 1);
}


JetCorrStd2P::JetCorrStd2P():
    JetCorrBase(2),
    ptRef(208.),
    paramsSPR({{1.03091, -0.051154, -0.154227}})  // Summer16_03Feb2017H_V3
{}


double JetCorrStd2P::Eval(double pt) const
{
    double response = 1. + parameters[0] + parameters[1] / 0.03 * (fSPR(pt) - fSPR(ptRef));
    return 1 / response;
}


void JetCorrStd2P::SetParamsSPR(std::initializer_list<double> paramsSPR_)
{
    if (paramsSPR_.size() != paramsSPR.size())
    {
        std::ostringstream message;
        message << "JetCorrStd2P::SetParamsSPR: Given " << paramsSPR_.size() <<
          " parameters while " << paramsSPR.size() << " expected.";
        throw std::runtime_error(message.str());
    }
    
    std::copy(paramsSPR_.begin(), paramsSPR_.end(), paramsSPR.begin());
}


double JetCorrStd2P::fSPR(double pt) const
{
    return std::max(0., paramsSPR[0] + paramsSPR[1] * std::pow(pt, paramsSPR[2]));
}


JetCorrStd3P::JetCorrStd3P():
    JetCorrStd2P(),
    paramsL1({{-0.196332, 0.307378}})  // Summer16_07AugBCDEFGH_V6
{
    // Change the number of free parameters from two to three
    parameters.resize(3);
}


double JetCorrStd3P::Eval(double pt) const
{
    double response = 1. + parameters[0] + \
        parameters[1] / 0.03 * (fSPR(pt) - fSPR(ptRef)) + \
        parameters[2] * (fL1(pt) - fL1(ptRef));
    return 1 / response;
}


void JetCorrStd3P::SetParamsL1(std::initializer_list<double> paramsL1_)
{
    if (paramsL1_.size() != paramsL1.size())
    {
        std::ostringstream message;
        message << "JetCorrStd2P::SetParamsL1: Given " << paramsL1_.size() <<
          " parameters while " << paramsL1.size() << " expected.";
        throw std::runtime_error(message.str());
    }
    
    std::copy(paramsL1_.begin(), paramsL1_.end(), paramsL1.begin());
}


double JetCorrStd3P::fL1(double pt) const
{
    return 1. - (paramsL1[0] + paramsL1[1] * std::log(pt)) / pt;
}
