#include <HistMorph.hpp>

#include <sstream>
#include <stdexcept>


HistMorph::HistMorph(std::vector<double> const &central_, std::vector<double> const &up_,
  std::vector<double> const &down_):
    central(central_), up(up_), down(down_)
{
    if (central.size() != up.size() or central.size() != down.size())
    {
        std::ostringstream message;
        message << "HistMorph::HistMorph: Lengths of given vectors (" << central.size() << ", " <<
          up.size() << ", " << down.size() << ") do not match.";
        throw std::logic_error(message.str());
    }
}


HistMorph::HistMorph(std::vector<double> const &up_, std::vector<double> const &down_):
    up(up_), down(down_)
{
    if (up.size() != down.size())
    {
        std::ostringstream message;
        message << "HistMorph::HistMorph: Lengths of given vectors (" << up.size() << ", " <<
          down.size() << ") do not match.";
        throw std::logic_error(message.str());
    }
    
    central.resize(up.size(), 0.);
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
    
    central.resize(numBins);
    up.resize(numBins);
    down.resize(numBins);
    
    for (unsigned bin = 1; bin <= numBins; ++bin)
    {
        central[bin - 1] = histCentral.GetBinContent(bin);
        up[bin - 1] = histUp.GetBinContent(bin);
        down[bin - 1] = histDown.GetBinContent(bin);
    }
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
    
    central.resize(numBins, 0.);
    up.resize(numBins);
    down.resize(numBins);
    
    for (unsigned bin = 1; bin <= numBins; ++bin)
    {
        up[bin - 1] = histUp.GetBinContent(bin);
        down[bin - 1] = histDown.GetBinContent(bin);
    }
}


double HistMorph::Eval(unsigned bin, double x) const
{
    if (bin >= central.size())
    {
        std::ostringstream message;
        message << "HistMorph::Eval: Bin with index " << bin << " requested, but only " <<
          central.size() << " bins are available.";
        throw std::out_of_range(message.str());
    }
    
    double const deltaUp = up[bin] - central[bin];
    double const deltaDown = down[bin] - central[bin];
    
    return central[bin] +
      ((deltaUp - deltaDown) / 2. + (deltaUp + deltaDown) / 2. * SmoothStep(x)) * x;
}


double HistMorph::SmoothStep(double x)
{
    if (x >= 1.)
        return 1.;
    else if (x <= -1.)
        return -1.;
    
    double const x2 = x * x;
    return x * (x2 * (3 * x2 - 10) + 15) / 8.;
}
