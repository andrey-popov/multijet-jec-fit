#include <JetCorrConstraint.hpp>

#include <cmath>


JetCorrConstraint::JetCorrConstraint(double ptRef_, double targetCorrection_,
  double relUncertainty_):
    ptRef(ptRef_), targetCorrection(targetCorrection_), relUncertainty(relUncertainty_)
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
