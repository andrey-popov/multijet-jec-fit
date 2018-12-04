#include <MultijetCrawlingBins.hpp>

#include <TFile.h>
#include <TKey.h>
#include <TVectorD.h>

#include <algorithm>
#include <cmath>
#include <cstring>
#include <limits>
#include <sstream>
#include <stdexcept>


using namespace std::string_literals;


MultijetCrawlingBins::JetCache::JetCache(std::vector<double> const &meanPtLead_,
  std::vector<double> const &meanPtJet_, double thresholdStart_, double thresholdEnd_):
    meanPtLead(meanPtLead_), meanPtJet(meanPtJet_),
    thresholdStart(thresholdStart_), thresholdEnd(thresholdEnd_),
    ptLeadCorrections(meanPtLead.size(), 0.), ptJetCorrections(meanPtJet.size(), 0.),
    ptBalWeights(meanPtJet.size(), 0.), firstPtJetBin(1), lastPtJetBin(ptBalWeights.size())
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
        ptBalWeights[i] = JetWeight(meanPtJet[i] * ptJetCorrections[i]);
    }
    
    // Find first bin for which the weight is not zero
    firstPtJetBin = 0;
    
    while (ptBalWeights[firstPtJetBin] == 0.)
        ++firstPtJetBin;
    
    // Convert to ROOT indexing convention
    ++firstPtJetBin;
    
    // The last bin in ROOT indexing convention
    lastPtJetBin = ptBalWeights.size();
}


double MultijetCrawlingBins::JetCache::Weight(unsigned bin) const
{
    return ptBalWeights[bin - 1];
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



MultijetCrawlingBins::Chi2Bin::Chi2Bin(unsigned firstBin_, unsigned lastBin_,
  std::shared_ptr<TH1> ptLeadHist_, std::shared_ptr<TH2> sumProj_,
  std::shared_ptr<TSpline3> simBalSpline_, double unc2_):
    firstBin(firstBin_), lastBin(lastBin_),
    ptLeadHist(ptLeadHist_), sumProj(sumProj_),
    simBalSpline(simBalSpline_), unc2(unc2_),
    jetCache(nullptr)
{}


double MultijetCrawlingBins::Chi2Bin::Chi2() const
{
    return std::pow(MeanBalance() - MeanSimBalance(), 2) / unc2;
}


double MultijetCrawlingBins::Chi2Bin::MeanBalance() const
{
    // Compute mean balance observable in data according to selected algorithm.  At the moment only
    // pt balance is supported.
    return MeanPtBal();
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
    if (method == Method::MPF)
        throw std::runtime_error("MPF method is not yet implemented in MultijetCrawlingBins.");
    
    
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
    auto ptThreshold = dynamic_cast<TVectorD *>(
      inputFile->Get((methodLabel + "Threshold").c_str()));
    
    if (not ptThreshold)
    {
        std::ostringstream message;
        message << "MultijetCrawlingBins::MultijetCrawlingBins: Failed to read jet " <<
          "pt threshold from file \"" << fileName << "\".";
        throw std::runtime_error(message.str());
    }
    
    
    // Loop over directories in the input file. Each directory describes a single trigger bin.
    TIter fileIter(inputFile->GetListOfKeys());
    TKey *key;
    
    std::vector<double> meanPtLead, meanPtJet;
    
    while ((key = dynamic_cast<TKey *>(fileIter())))
    {
        if (strcmp(key->GetClassName(), "TDirectoryFile") != 0)
            continue;
        
        TDirectoryFile *directory = dynamic_cast<TDirectoryFile *>(key->ReadObj());
        
        for (auto const &name: std::initializer_list<std::string>{"Binning", "Sim" + methodLabel,
          "PtLead", "PtLeadProfile", methodLabel + "Profile", "RelPtJetSumProj"})
        {
            if (not directory->Get(name.c_str()))
            {
                std::ostringstream message;
                message << "MultijetCrawlingBins::MultijetCrawlingBins: Directory \"" <<
                  key->GetName() << "\" in file \"" << fileName <<
                  "\" does not contain required key \"" << name << "\".";
                throw std::runtime_error(message.str());
            }
        }
        
        
        // Read information about the current trigger bin
        auto binning = *dynamic_cast<TVectorD *>(directory->Get("Binning"));
        auto ptLeadProfile = dynamic_cast<TProfile *>(directory->Get("PtLeadProfile"));
        auto balProfile = dynamic_cast<TProfile *>(
          directory->Get((methodLabel + "Profile").c_str()));
        
        std::shared_ptr<TH1> ptLeadHist(dynamic_cast<TH1 *>(directory->Get("PtLead")));
        std::shared_ptr<TH2> sumProj(dynamic_cast<TH2 *>(directory->Get("RelPtJetSumProj")));
        std::shared_ptr<TSpline3> simBalSpline(
          dynamic_cast<TSpline3 *>(directory->Get(("Sim" + methodLabel).c_str())));
        
        ptLeadHist->SetDirectory(nullptr);
        sumProj->SetDirectory(nullptr);
        
        
        // Rebin TProfile with mean balance observable in data to the binning used in the
        // computation of the chi^2 in order to obtain per-bin uncertainties
        TH1 *balProfileRebinned = balProfile->Rebin(binning.GetNoElements() - 1,
          (balProfile->GetName() + "Rebinned"s).c_str(), binning.GetMatrixArray());
        
        
        // Find a number that is smaller than the width of any bin in pt of the leading jet in the
        // underlying histograms. It is used for the matching between the underlying binning and
        // the binning used in the computation of the chi^2.
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
        
        
        // Construct chi^2 bins. Each one consists of one or (typically) more bins in pt of the
        // leading jet that are included in the range of a single bin in variable `binning`.
        for (int binChi2 = 1; binChi2 < binning.GetNoElements(); ++binChi2)
        {
            unsigned firstBin = ptLeadHist->FindFixBin(binning[binChi2 - 1] + eps);
            unsigned lastBin = ptLeadHist->FindFixBin(binning[binChi2] - eps);
            
            chi2Bins.emplace_back(firstBin, lastBin, ptLeadHist, sumProj, simBalSpline,
              std::pow(balProfileRebinned->GetBinError(binChi2), 2));
        }
        
        
        // Save mean values of pt of the leading jet in each bin of the underlying binning. If the
        // vector is empty, first prefill it with positions of bin centres. Then update the
        // estimates for non-empty using more precise results from `ptLeadProfile`. Note that the
        // binning is shared among all trigger bins, but there is no overlap in filled bins.
        if (meanPtLead.empty())
        {
            meanPtLead.reserve(ptLeadProfile->GetNbinsX());
            
            for (int bin = 1; bin <= ptLeadProfile->GetNbinsX(); ++bin)
                meanPtLead.emplace_back(ptLeadProfile->GetBinCenter(bin));
        }
        
        for (int bin = 1; bin <= ptLeadProfile->GetNbinsX(); ++bin)
        {
            if (ptLeadProfile->GetBinContent(bin) > 0)
                meanPtLead[bin - 1] = ptLeadProfile->GetBinContent(bin);
        }
        
        
        
        // Save mean jet pt along the second axis of sumProj, which is needed for the JetCache
        // object. Assume the binning is the same for all triggers. Take positions of bin centres
        // as approximate mean values.
        if (meanPtJet.empty())
        {
            meanPtJet.reserve(sumProj->GetNbinsY());
            
            for (int bin = 1; bin <= sumProj->GetNbinsY(); ++bin)
                meanPtJet.emplace_back(sumProj->GetYaxis()->GetBinCenter(bin));
        }
    }
    
    inputFile->Close();
    
    if (chi2Bins.empty())
    {
        std::ostringstream message;
        message << "MultijetCrawlingBins::MultijetCrawlingBins: No data read from file \"" <<
          fileName << "\".";
        throw std::runtime_error(message.str());
    }
    
    
    // Make sure trigger bins and `meanPtLead` are sorted in pt
    std::sort(chi2Bins.begin(), chi2Bins.end(),
      [](auto const &left, auto const &right)
      {return (left.PtRange().first < right.PtRange().first);});
    std::sort(meanPtLead.begin(), meanPtLead.end());
    
    for (unsigned i = 0; i < chi2Bins.size() - 1; ++i)
        if (chi2Bins[i].PtRange().second > chi2Bins[i + 1].PtRange().first)
        {
            std::ostringstream message;
            message << "MultijetCrawlingBins::MultijetCrawlingBins: Overlapping chi^2 bins found.";
            throw std::runtime_error(message.str());
        }
    
    for (unsigned i = 0; i < meanPtLead.size() - 1; ++i)
        if (meanPtLead[i] == meanPtLead[i + 1])
        {
            std::ostringstream message;
            message << "MultijetCrawlingBins::MultijetCrawlingBins: Duplicate bins in pt of the "
              "leading jet found.";
            throw std::runtime_error(message.str());
        }
    
    chi2BinMask.assign(chi2Bins.size(), true);
    
    
    // Initialize the object to cache values of jet corrections
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
