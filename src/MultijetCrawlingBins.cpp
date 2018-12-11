#include <MultijetCrawlingBins.hpp>

#include <TFile.h>
#include <TKey.h>
#include <TVectorD.h>

#include <algorithm>
#include <cmath>
#include <limits>
#include <sstream>
#include <stdexcept>


using namespace std::string_literals;


MultijetCrawlingBins::JetCache::JetCache(std::vector<double> const &meanPtLead_,
  std::vector<double> const &meanPtJet_, double thresholdStart_, double thresholdEnd_):
    meanPtLead(meanPtLead_), meanPtJet(meanPtJet_),
    thresholdStart(thresholdStart_), thresholdEnd(thresholdEnd_),
    ptLeadCorrections(meanPtLead.size(), 0.), ptJetCorrections(meanPtJet.size(), 0.),
    jetWeights(meanPtJet.size(), 0.), firstPtJetBin(1), lastPtJetBin(jetWeights.size())
{}


double MultijetCrawlingBins::JetCache::CorrectedMeanPtLead(unsigned bin) const
{
    return meanPtLead[bin - 1] * ptLeadCorrections[bin - 1];
}


double MultijetCrawlingBins::JetCache::CorrectionPtLead(unsigned bin) const
{
    return ptLeadCorrections[bin - 1];
}


double MultijetCrawlingBins::JetCache::CorrectionPtJet(unsigned bin) const
{
    return ptJetCorrections[bin - 1];
}


std::pair<unsigned, unsigned> MultijetCrawlingBins::JetCache::PtJetBinRange() const
{
    return {firstPtJetBin, lastPtJetBin};
}


void MultijetCrawlingBins::JetCache::Update(JetCorrBase const &corrector)
{
    for (unsigned i = 0; i < meanPtLead.size(); ++i)
        ptLeadCorrections[i] = corrector.Eval(meanPtLead[i]);
    
    for (unsigned i = 0; i < meanPtJet.size(); ++i)
    {
        ptJetCorrections[i] = corrector.Eval(meanPtJet[i]);
        jetWeights[i] = JetWeight(meanPtJet[i] * ptJetCorrections[i]);
    }
    
    // Find first bin for which the weight is not zero
    firstPtJetBin = 0;
    
    while (jetWeights[firstPtJetBin] == 0.)
        ++firstPtJetBin;
    
    // Convert to ROOT indexing convention
    ++firstPtJetBin;
    
    // The last bin in ROOT indexing convention
    lastPtJetBin = jetWeights.size();
}


double MultijetCrawlingBins::JetCache::Weight(unsigned bin) const
{
    return jetWeights[bin - 1];
}


double MultijetCrawlingBins::JetCache::JetWeight(double pt) const
{
    // Special treatment for a sharp threshold
    if (thresholdStart == thresholdEnd)
    {
        if (pt >= thresholdStart)
            return 1.;
        else
            return 0.;
    }
    
    double const x = (pt - thresholdStart) / (thresholdEnd - thresholdStart);
    
    if (x < 0.)
        return 0.;
    else if (x > 1.)
        return 1.;
    else
        return -2 * std::pow(x, 3) + 3 * std::pow(x, 2);
}



MultijetCrawlingBins::Chi2Bin::Chi2Bin(MultijetCrawlingBins::Method method, unsigned firstBin_,
  unsigned lastBin_, std::shared_ptr<TH1> ptLeadHist_, std::shared_ptr<TProfile> mpfProfile_,
  std::shared_ptr<TH2> sumProj_, std::shared_ptr<TSpline3> simBalSpline_, double unc2_):
    firstBin(firstBin_), lastBin(lastBin_),
    ptLeadHist(ptLeadHist_), mpfProfile(mpfProfile_), sumProj(sumProj_),
    simBalSpline(simBalSpline_), unc2(unc2_),
    jetCache(nullptr)
{
    if (method == MultijetCrawlingBins::Method::PtBal)
        meanBalanceCalc = &Chi2Bin::MeanPtBal;
    else
        meanBalanceCalc = &Chi2Bin::MeanMPF;
}


double MultijetCrawlingBins::Chi2Bin::Chi2() const
{
    return std::pow(MeanBalance() - MeanSimBalance(), 2) / unc2;
}


double MultijetCrawlingBins::Chi2Bin::MeanBalance() const
{
    return (this->*meanBalanceCalc)();
}


double MultijetCrawlingBins::Chi2Bin::MeanSimBalance() const
{
    double sumBal = 0., numEvents = 0.;
    
    for (unsigned binPtLead = firstBin; binPtLead <= lastBin; ++binPtLead)
    {
        double const n = ptLeadHist->GetBinContent(binPtLead);
        sumBal += SimBalance(jetCache->CorrectedMeanPtLead(binPtLead)) * n;
        numEvents += n;
    }
    
    return sumBal / numEvents;
}


std::pair<double, double> MultijetCrawlingBins::Chi2Bin::PtRange() const
{
    return {ptLeadHist->GetBinLowEdge(firstBin), ptLeadHist->GetBinLowEdge(lastBin + 1)};
}


double MultijetCrawlingBins::Chi2Bin::SimBalance(double ptLead) const
{
    return simBalSpline->Eval(std::log(ptLead));
}


void MultijetCrawlingBins::Chi2Bin::SetJetCache(JetCache const *jetCache_)
{
    jetCache = jetCache_;
}


double MultijetCrawlingBins::Chi2Bin::MeanMPF() const
{
    double sumBal = 0.;
    double numEvents = 0.;
    
    auto const ptJetBinRange = jetCache->PtJetBinRange();
    
    for (unsigned binPtLead = firstBin; binPtLead <= lastBin; ++binPtLead)
    {
        sumBal += mpfProfile->GetBinContent(binPtLead) * ptLeadHist->GetBinContent(binPtLead) / \
          jetCache->CorrectionPtLead(binPtLead);
        
        double sumJets = 0.;
        
        for (unsigned binPtJet = ptJetBinRange.first; binPtJet <= ptJetBinRange.second; ++binPtJet)
            sumJets -= sumProj->GetBinContent(binPtLead, binPtJet) * \
              (1 - jetCache->CorrectionPtJet(binPtJet)) * jetCache->Weight(binPtJet);
        
        sumBal += sumJets / jetCache->CorrectionPtLead(binPtLead);
        numEvents += ptLeadHist->GetBinContent(binPtLead);
    }
    
    return sumBal / numEvents;
}


double MultijetCrawlingBins::Chi2Bin::MeanPtBal() const
{
    double sumBal = 0.;
    double numEvents = 0.;
    
    auto const ptJetBinRange = jetCache->PtJetBinRange();
    
    for (unsigned binPtLead = firstBin; binPtLead <= lastBin; ++binPtLead)
    {
        double sumJets = 0.;
        
        for (unsigned binPtJet = ptJetBinRange.first; binPtJet <= ptJetBinRange.second; ++binPtJet)
            sumJets += sumProj->GetBinContent(binPtLead, binPtJet) * \
              jetCache->CorrectionPtJet(binPtJet) * jetCache->Weight(binPtJet);
        
        sumBal += sumJets / jetCache->CorrectionPtLead(binPtLead);
        numEvents += ptLeadHist->GetBinContent(binPtLead);
    }
    
    return sumBal / numEvents;
}



MultijetCrawlingBins::MultijetCrawlingBins(std::string const &fileName,
  MultijetCrawlingBins::Method method_):
    method(method_)
{
    std::string methodLabel;
    
    if (method == Method::PtBal)
        methodLabel = "PtBal";
    else if (method == Method::MPF)
        methodLabel = "MPF";
    
    
    std::unique_ptr<TFile> inputFile(TFile::Open(fileName.c_str()));
    
    if (not inputFile or inputFile->IsZombie())
    {
        std::ostringstream message;
        message << "MultijetCrawlingBins::MultijetCrawlingBins: Failed to open file \"" <<
          fileName << "\".";
        throw std::runtime_error(message.str());
    }
    
    
    // Read the jet pt threshold. It is not a free parameter and must be set to the same value as
    // used to construct the inputs. For the pt balance method, it affects the definition of the
    // balance observable in simulation (while in data it can be recomputed for any not too low
    // threshold). In the case of the MPF method, the definition of the balance observable in both
    // data and simulation is affected.
    std::unique_ptr<TVectorD> ptThreshold(dynamic_cast<TVectorD *>(
      inputFile->Get((methodLabel + "Threshold").c_str())));
    
    if (not ptThreshold)
    {
        std::ostringstream message;
        message << "MultijetCrawlingBins::MultijetCrawlingBins: Failed to read jet " <<
          "pt threshold from file \"" << fileName << "\".";
        throw std::runtime_error(message.str());
    }
    
    if (ptThreshold->GetNoElements() != 2)
    {
        std::ostringstream message;
        message << "MultijetCrawlingBins::MultijetCrawlingBins: Unexpected number of elements "
          "read for jet pt threshold.";
        throw std::runtime_error(message.str());
    }
    
    
    // Read target binning for computation of chi^2 and data histograms
    for (auto const &name: std::initializer_list<std::string>{"Binning", "PtLead", "PtLeadProfile",
      methodLabel + "Profile", "RelPtJetSumProj"})
    {
        if (not inputFile->Get(name.c_str()))
        {
            std::ostringstream message;
            message << "MultijetCrawlingBins::MultijetCrawlingBins: File \"" << fileName <<
              "\" does not contain required key \"" << name << "\".";
            throw std::runtime_error(message.str());
        }
    }
    
    std::unique_ptr<TVectorD> binning(dynamic_cast<TVectorD *>(inputFile->Get("Binning")));
    
    std::shared_ptr<TH1> ptLeadHist(dynamic_cast<TH1 *>(inputFile->Get("PtLead")));
    std::unique_ptr<TProfile> ptLeadProfile(dynamic_cast<TProfile *>(
      inputFile->Get("PtLeadProfile")));
    std::shared_ptr<TProfile> balProfile(dynamic_cast<TProfile *>(
      inputFile->Get((methodLabel + "Profile").c_str())));
    std::shared_ptr<TH2> sumProj(dynamic_cast<TH2 *>(inputFile->Get("RelPtJetSumProj")));
    
    ptLeadHist->SetDirectory(nullptr);
    ptLeadProfile->SetDirectory(nullptr);
    balProfile->SetDirectory(nullptr);
    sumProj->SetDirectory(nullptr);
    
    
    // Rebin TProfile with mean balance observable in data to the target binning.  It will be used
    // to obtain per-bin uncertainties.
    std::unique_ptr<TH1> balProfileRebinned(balProfile->Rebin(binning->GetNoElements() - 1,
      (balProfile->GetName() + "Rebinned"s).c_str(), binning->GetMatrixArray()));
    balProfileRebinned->SetDirectory(nullptr);
    
    
    // Read splines to compute mean value of the balance observable in simulation. They are
    // provided separately for different trigger bins.
    std::vector<std::pair<double, std::shared_ptr<TSpline3>>> simBalSplines;
    
    TIter fileIter(inputFile->GetListOfKeys());
    TKey *key;
    
    while ((key = dynamic_cast<TKey *>(fileIter())))
    {
        if (key->GetClassName() != "TDirectoryFile"s)
            continue;
        
        TDirectoryFile *directory = dynamic_cast<TDirectoryFile *>(key->ReadObj());
        
        for (auto const &name: std::initializer_list<std::string>{"Range", "Sim" + methodLabel})
        {
            if (not directory->Get(name.c_str()))
            {
                std::ostringstream message;
                message << "MultijetCrawlingBins::MultijetCrawlingBins: Directory \"" <<
                  directory->GetName() << "\" in file \"" << fileName <<
                  "\" does not contain required key \"" << name << "\".";
                throw std::runtime_error(message.str());
            }
        }
        
        std::unique_ptr<TVectorD> range(dynamic_cast<TVectorD *>(directory->Get("Range")));
        simBalSplines.emplace_back(std::make_pair((*range)[0], dynamic_cast<TSpline3 *>(
          directory->Get(("Sim" + methodLabel).c_str()))));
    }
    
    
    std::sort(simBalSplines.begin(), simBalSplines.end(),
      [](auto const &lhs, auto const &rhs){return (lhs.first < rhs.first);});
    
    inputFile->Close();
    
    
    // Find a number that is smaller than the width of any bin in pt of the leading jet in the
    // underlying histograms. It is used for the matching between the underlying binning and
    // the target chi^2 binning.
    double eps = std::numeric_limits<double>::infinity();
    
    for (int bin = 1; bin <= ptLeadHist->GetNbinsX(); ++bin)
    {
        double const binWidth = ptLeadHist->GetBinWidth(bin);
        
        if (eps > binWidth)
            eps = binWidth;
    }
    
    eps /= 2;
    
    if (eps <= 0.)
    {
        std::ostringstream message;
        message << "MultijetCrawlingBins::MultijetCrawlingBins: Found bins of zero width.";
        throw std::runtime_error(message.str());
    }
    
    
    // Construct chi^2 bins. Each one consists of one or (typically) more bins in pt of the leading
    // jet that are included in the range of a single bin in variable `binning`.
    for (int binChi2 = 1; binChi2 < binning->GetNoElements(); ++binChi2)
    {
        unsigned firstBin = ptLeadHist->FindFixBin((*binning)[binChi2 - 1] + eps);
        unsigned lastBin = ptLeadHist->FindFixBin((*binning)[binChi2] - eps);
        
        auto simBalSplineIt = std::lower_bound(simBalSplines.begin(), simBalSplines.end(),
          (*binning)[binChi2 - 1] + eps,
          [](auto const &lhs, double const &rhs){return (lhs.first < rhs);});
        --simBalSplineIt;
        
        chi2Bins.emplace_back(method, firstBin, lastBin, ptLeadHist,
          (method == MultijetCrawlingBins::Method::MPF) ? balProfile : nullptr,
          sumProj, simBalSplineIt->second, std::pow(balProfileRebinned->GetBinError(binChi2), 2));
    }
    
    if (chi2Bins.empty())
    {
        std::ostringstream message;
        message << "MultijetCrawlingBins::MultijetCrawlingBins: No data read from file \"" <<
          fileName << "\".";
        throw std::runtime_error(message.str());
    }
    
    chi2BinMask.assign(chi2Bins.size(), true);
    
    
    // Initialize the object to cache values of jet corrections
    std::vector<double> meanPtLead;
    meanPtLead.reserve(ptLeadHist->GetNbinsX());
    
    for (int bin = 1; bin <= ptLeadHist->GetNbinsX(); ++bin)
    {
        if (ptLeadHist->GetBinContent(bin) > 0.)
            meanPtLead.emplace_back(ptLeadProfile->GetBinContent(bin));
        else
            meanPtLead.emplace_back(ptLeadProfile->GetBinCenter(bin));
    }
    
    std::vector<double> meanPtJet;
    meanPtJet.reserve(sumProj->GetNbinsY());
    
    for (int bin = 1; bin <= sumProj->GetNbinsY(); ++bin)
        meanPtJet.emplace_back(sumProj->GetYaxis()->GetBinCenter(bin));
    
    jetCache.reset(new JetCache(meanPtLead, meanPtJet, (*ptThreshold)[0], (*ptThreshold)[1]));
    
    for (auto &chi2Bin: chi2Bins)
        chi2Bin.SetJetCache(jetCache.get());
}


unsigned MultijetCrawlingBins::GetDim() const
{
    return std::count(chi2BinMask.begin(), chi2BinMask.end(), true);
}


double MultijetCrawlingBins::Eval(JetCorrBase const &corrector, Nuisances const &) const
{
    jetCache->Update(corrector);
    double chi2 = 0.;
    
    for (unsigned i = 0; i < chi2Bins.size(); ++i)
    {
        if (chi2BinMask[i])
            chi2 += chi2Bins[i].Chi2();
    }
    
    return chi2;
}


std::pair<double, double> MultijetCrawlingBins::SetPtLeadRange(double minPt, double maxPt)
{
    // Construct an auxiliary vector of all boundaries between chi^2 bins. Assume that all bins are
    // adjacent.
    std::vector<double> edges;
    edges.reserve(chi2Bins.size() + 1);
    
    for (auto const &chi2Bin: chi2Bins)
        edges.emplace_back(chi2Bin.PtRange().first);
    
    edges.emplace_back(chi2Bins.back().PtRange().second);
    
    
    // Find closest edges
    unsigned iEdgeMin = std::lower_bound(edges.begin(), edges.end(), minPt) - edges.begin();
    
    if (iEdgeMin == edges.size())
        --iEdgeMin;
    else if (iEdgeMin > 0)
    {
        if (edges[iEdgeMin] - minPt > minPt - edges[iEdgeMin - 1])
            --iEdgeMin;
    }
    
    unsigned iEdgeMax = std::lower_bound(edges.begin(), edges.end(), maxPt) - edges.begin();
    
    if (iEdgeMax == edges.size())
        --iEdgeMax;
    else if (iEdgeMax > 0)
    {
        if (edges[iEdgeMax] - maxPt > maxPt - edges[iEdgeMax - 1])
            --iEdgeMax;
    }
    
    if (iEdgeMax <= iEdgeMin)
    {
        std::ostringstream message;
        message << "MultijetCrawlingBins::SetPtLeadRange: Requested range is too narrow.";
        throw std::runtime_error(message.str());
    }
    
    
    // Mask chi^2 bins outsize of the range
    for (unsigned i = 0; i < iEdgeMin; ++i)
        chi2BinMask[i] = false;
    
    for (unsigned i = iEdgeMin; i < iEdgeMax; ++i)
        chi2BinMask[i] = true;
    
    for (unsigned i = iEdgeMax; i < chi2BinMask.size(); ++i)
        chi2BinMask[i] = false;
    
    
    return {edges[iEdgeMin], edges[iEdgeMax]};
}
