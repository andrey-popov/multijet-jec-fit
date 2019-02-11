#include <Morphing.hpp>

#include <sstream>
#include <stdexcept>



PointMorph::PointMorph(double central_, double up_, double down_):
    central(central_), up(up_), down(down_)
{}


double PointMorph::Eval(double x) const
{
    return Morph(central, up, down, x);
}


double PointMorph::Morph(double central, double up, double down, double x)
{
    double const deltaUp = up - central;
    double const deltaDown = down - central;
    
    return central + ((deltaUp - deltaDown) / 2. + (deltaUp + deltaDown) / 2. * SmoothStep(x)) * x;
}


double PointMorph::operator()(double x) const
{
    return Eval(x);
}


double PointMorph::SmoothStep(double x)
{
    if (x >= 1.)
        return 1.;
    else if (x <= -1.)
        return -1.;
    
    double const x2 = x * x;
    return x * (x2 * (3 * x2 - 10) + 15) / 8.;
}



HistMorph::HistMorph(std::vector<double> const &central, std::vector<double> const &up,
  std::vector<double> const &down)
{
    if (central.size() != up.size() or central.size() != down.size())
    {
        std::ostringstream message;
        message << "HistMorph::HistMorph: Lengths of given vectors (" << central.size() << ", " <<
          up.size() << ", " << down.size() << ") do not match.";
        throw std::logic_error(message.str());
    }

    bins.reserve(central.size());

    for (unsigned i = 0; i < central.size(); ++i)
        bins.emplace_back(central[i], up[i], down[i]);
}


HistMorph::HistMorph(std::vector<double> const &up, std::vector<double> const &down)
{
    if (up.size() != down.size())
    {
        std::ostringstream message;
        message << "HistMorph::HistMorph: Lengths of given vectors (" << up.size() << ", " <<
          down.size() << ") do not match.";
        throw std::logic_error(message.str());
    }
    
    bins.reserve(up.size());

    for (unsigned i = 0; i < up.size(); ++i)
        bins.emplace_back(0., up[i], down[i]);
}


HistMorph::HistMorph(TH1 const &histCentral, TH1 const &histUp, TH1 const &histDown)
{
    unsigned const numBins = histCentral.GetNbinsX();
    
    if (unsigned(histUp.GetNbinsX()) != numBins or unsigned(histDown.GetNbinsX()) != numBins)
    {
        std::ostringstream message;
        message << "HistMorph::HistMorph: Numbers of bins in given histograms (" << numBins <<
          ", " << histUp.GetNbinsX() << ", " << histDown.GetNbinsX() << ") do not match.";
        throw std::logic_error(message.str());
    }
    
    bins.reserve(numBins);
    
    for (unsigned bin = 1; bin <= numBins; ++bin)
        bins.emplace_back(histCentral.GetBinContent(bin), histUp.GetBinContent(bin),
          histDown.GetBinContent(bin));
}


HistMorph::HistMorph(TH1 const &histUp, TH1 const &histDown)
{
    unsigned const numBins = histUp.GetNbinsX();
    
    if (unsigned(histDown.GetNbinsX()) != numBins)
    {
        std::ostringstream message;
        message << "HistMorph::HistMorph: Numbers of bins in given histograms (" << numBins <<
          ", " << histDown.GetNbinsX() << ") do not match.";
        throw std::logic_error(message.str());
    }
    
    bins.reserve(numBins);
    
    for (unsigned bin = 1; bin <= numBins; ++bin)
        bins.emplace_back(0., histUp.GetBinContent(bin), histDown.GetBinContent(bin));
}


double HistMorph::Eval(unsigned bin, double x) const
{
    if (bin >= bins.size())
    {
        std::ostringstream message;
        message << "HistMorph::Eval: Bin with index " << bin << " requested, but only " <<
          bins.size() << " bins are available.";
        throw std::out_of_range(message.str());
    }

    return bins[bin](x);
}

