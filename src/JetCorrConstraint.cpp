#include <JetCorrConstraint.hpp>

#include <cmath>


JetCorrConstraint::JetCorrConstraint():
    ptRef(208.),
    // Target correction and its uncertainty are taken from Summer16_07Aug2017GH_V10 for
    //pt = 208 GeV, eta = 0
    targetCorrection(1.0107136), relUncertainty(0.0078)
{}


unsigned JetCorrConstraint::GetDim() const
{
    return 1;
}


double JetCorrConstraint::Eval(JetCorrBase const &corrector, Nuisances const &) const
{
    double const correction = corrector.Eval(ptRef);
    return std::pow((1 - correction / targetCorrection) / relUncertainty, 2);
}


void JetCorrConstraint::SetParameters(double ptRef_, double targetCorrection_,
  double relUncertainty_)
{
    ptRef = ptRef_;
    targetCorrection = targetCorrection_;
    relUncertainty = relUncertainty_;
}
