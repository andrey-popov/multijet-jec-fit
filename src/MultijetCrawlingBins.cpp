#include <MultijetCrawlingBins.hpp>

#include <TFile.h>
#include <TKey.h>
#include <TVectorD.h>

#include <algorithm>
#include <cmath>
#include <limits>
#include <regex>
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
  std::shared_ptr<TH2> sumProj_, std::shared_ptr<Spline> simBalSpline_, double unc2_):
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


void MultijetCrawlingBins::Chi2Bin::AddDataSyst(unsigned nuisanceIndex, double up, double down)
{
    if (dataVariations.count(nuisanceIndex) > 0)
    {
        std::ostringstream message;
        message << "MultijetCrawlingBins::Chi2Bin::AddDataSyst: Systematic variation in data for "
          "the nuisance parameter with index " << nuisanceIndex << " has already been registered.";
        throw std::runtime_error(message.str());
    }

    dataVariations[nuisanceIndex] = PointMorph(0., up, down);
}


void MultijetCrawlingBins::Chi2Bin::AddSimSyst(unsigned nuisanceIndex, std::shared_ptr<Spline> up,
  std::shared_ptr<Spline> down)
{
    if (simVariations.count(nuisanceIndex) > 0)
    {
        std::ostringstream message;
        message << "MultijetCrawlingBins::Chi2Bin::AddSimSyst: Systematic variation in simulation "
          "for the nuisance parameter with index " << nuisanceIndex << " has already been "
          "registered.";
        throw std::runtime_error(message.str());
    }

    simVariations[nuisanceIndex] = std::array<std::shared_ptr<Spline>, 2>{up, down};
}


double MultijetCrawlingBins::Chi2Bin::Chi2(Nuisances const &nuisances) const
{
    return std::pow(MeanBalance(nuisances) - MeanSimBalance(nuisances), 2) / unc2;
}


double MultijetCrawlingBins::Chi2Bin::MeanBalance(Nuisances const &nuisances) const
{
    return (this->*meanBalanceCalc)(nuisances);
}


double MultijetCrawlingBins::Chi2Bin::MeanPt() const
{
    double sumPt = 0., numEvents = 0.;

    for (unsigned binPtLead = firstBin; binPtLead <= lastBin; ++binPtLead)
    {
        double const n = ptLeadHist->GetBinContent(binPtLead);
        sumPt += jetCache->CorrectedMeanPtLead(binPtLead) * n;
        numEvents += n;
    }

    return sumPt / numEvents;
}


double MultijetCrawlingBins::Chi2Bin::MeanSimBalance(Nuisances const &nuisances) const
{
    double sumBal = 0., numEvents = 0.;
    
    for (unsigned binPtLead = firstBin; binPtLead <= lastBin; ++binPtLead)
    {
        double const n = ptLeadHist->GetBinContent(binPtLead);
        sumBal += SimBalance(jetCache->CorrectedMeanPtLead(binPtLead), nuisances) * n;
        numEvents += n;
    }
    
    return sumBal / numEvents;
}


std::pair<double, double> MultijetCrawlingBins::Chi2Bin::PtRange() const
{
    return {ptLeadHist->GetBinLowEdge(firstBin), ptLeadHist->GetBinLowEdge(lastBin + 1)};
}


double MultijetCrawlingBins::Chi2Bin::SimBalance(double ptLead, Nuisances const &nuisances) const
{
    double const logPt = std::log(ptLead);
    double meanBalance = simBalSpline->Eval(logPt);


    // Apply systematic variations
    for (auto const &syst: simVariations)
    {
        // Reference up and down relative deviations for the current uncertainty
        double const up = syst.second[0]->Eval(logPt);
        double const down = syst.second[1]->Eval(logPt);

        // Interpolate between them and apply the resulting deviation
        meanBalance *= 1 + PointMorph::Morph(0, up, down, nuisances[syst.first]);
    }

    return meanBalance;
}


void MultijetCrawlingBins::Chi2Bin::SetJetCache(JetCache const *jetCache_)
{
    jetCache = jetCache_;
}


double MultijetCrawlingBins::Chi2Bin::Uncertainty() const
{
    return std::sqrt(unc2);
}


double MultijetCrawlingBins::Chi2Bin::MeanMPF(Nuisances const &nuisances) const
{
    // Compute nominal mean MPF balance
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
    
    double meanBalance = sumBal / numEvents;


    // Apply systematic variations
    for (auto const &syst: dataVariations)
    {
        double const deviation = syst.second(nuisances[syst.first]);
        meanBalance *= 1 + deviation;
    }

    return meanBalance;
}


double MultijetCrawlingBins::Chi2Bin::MeanPtBal(Nuisances const &nuisances) const
{
    // Compute nominal pt balance
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
    
    double meanBalance = sumBal / numEvents;


    // Apply systematic variations
    for (auto const &syst: dataVariations)
    {
        double const deviation = syst.second(nuisances[syst.first]);
        meanBalance *= 1 + deviation;
    }

    return meanBalance;
}



MultijetCrawlingBins::MultijetCrawlingBins(std::string const &fileName,
  MultijetCrawlingBins::Method method_, NuisanceDefinitions &nuisanceDefs,
  std::set<std::string> systToExclude):
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
    std::vector<std::pair<double, std::shared_ptr<Spline>>> simBalSplines;
    
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
        simBalSplines.emplace_back(std::make_pair((*range)[0], dynamic_cast<Spline *>(
          directory->Get(("Sim" + methodLabel).c_str()))));
    }
    
    
    std::sort(simBalSplines.begin(), simBalSplines.end(),
      [](auto const &lhs, auto const &rhs){return (lhs.first < rhs.first);});
    
    
    // Read systematic variations in data
    std::map<std::string, std::array<std::unique_ptr<TH1>, 2>> dataVariations;
    
    std::regex dataSystRegex("RelVar_" + methodLabel + "_(.+)Up", std::regex::extended);
    std::cmatch matchResult;
    fileIter = inputFile->GetListOfKeys();

    while ((key = dynamic_cast<TKey *>(fileIter())))
    {
        if (not std::regex_match(key->GetName(), matchResult, dataSystRegex))
            continue;

        std::string const systLabel(matchResult[1]);

        if (systToExclude.count(systLabel) > 0)
            continue;

        std::unique_ptr<TH1> histUp(dynamic_cast<TH1 *>(key->ReadObj()));
        std::unique_ptr<TH1> histDown(dynamic_cast<TH1 *>(inputFile->Get(
          ("RelVar_" + methodLabel + "_" + systLabel + "Down").c_str())));

        if (not histUp or not histDown)
        {
            std::ostringstream message;
            message << "MultijetCrawlingBins::MultijetCrawlingBins: Failed to read systematic "
              "variation \"" << systLabel << "\" for data.";
            throw std::runtime_error(message.str());
        }

        int const numBins = binning->GetNoElements() - 1;

        if (histUp->GetNbinsX() != numBins or histDown->GetNbinsX() != numBins)
        {
            std::ostringstream message;
            message << "MultijetCrawlingBins::MultijetCrawlingBins: Number of bins in histograms "
              "that define systematic variation \"" << systLabel << "\" in data, does not agree "
              "with the given chi^2 binning.";
            throw std::runtime_error(message.str());
        }

        histUp->SetDirectory(nullptr);
        histDown->SetDirectory(nullptr);

        dataVariations[systLabel] = std::array<std::unique_ptr<TH1>, 2>{
          std::move(histUp), std::move(histDown)};
    }


    // Read systematic variations in simulation. A single variation is descibed by an array of two
    // splines (which are wrapped into shared_ptr). The variations are associated with the lower
    // boundaries of the corresponding trigger bins, using an std::pair, and put into an ordered
    // vector. The vectors for different systematic uncertainties are aggregated in a map, whose
    // keys are the labels of the uncertainties.
    std::map<std::string, std::vector<std::pair<double, std::array<std::shared_ptr<Spline>, 2>>>>
      simVariations;

    std::regex simSystRegex("RelVar_Sim" + methodLabel + "_(.+)Up", std::regex::extended);
    fileIter = inputFile->GetListOfKeys();

    while ((key = dynamic_cast<TKey *>(fileIter())))
    {
        if (key->GetClassName() != "TDirectoryFile"s)
            continue;
        
        TDirectoryFile *directory = dynamic_cast<TDirectoryFile *>(key->ReadObj());
        std::unique_ptr<TVectorD> range(dynamic_cast<TVectorD *>(directory->Get("Range")));
        
        TIter dirIter(directory->GetListOfKeys());
        TKey *subKey;

        while ((subKey = dynamic_cast<TKey *>(dirIter())))
        {
            if (not std::regex_match(subKey->GetName(), matchResult, simSystRegex))
                continue;

            std::string const systLabel(matchResult[1]);

            if (systToExclude.count(systLabel) > 0)
                continue;

            std::shared_ptr<Spline> splineUp(dynamic_cast<Spline *>(subKey->ReadObj()));
            std::shared_ptr<Spline> splineDown(dynamic_cast<Spline *>(
              directory->Get(("RelVar_Sim" + methodLabel + "_" + systLabel + "Down").c_str())));

            if (not splineUp or not splineDown)
            {
                std::ostringstream message;
                message << "MultijetCrawlingBins::MultijetCrawlingBins: Failed to read systematic "
                  "variation \"" << systLabel << "\" for simulation.";
                throw std::runtime_error(message.str());
            }

            std::array<std::shared_ptr<Spline>, 2> splinePair{splineUp, splineDown};

            if (simVariations.find(systLabel) == simVariations.end())
                simVariations.emplace(std::piecewise_construct, std::forward_as_tuple(systLabel),
                  std::forward_as_tuple());

            simVariations[systLabel].emplace_back((*range)[0], splinePair);
        }
    }

    // Sort vectors for all systematic uncertainties in simulation according to the lower bounds of
    // the pt ranges of the corresponding trigger bins
    for (auto &syst: simVariations)
    {
        std::sort(syst.second.begin(), syst.second.end(),
          [](auto const &lhs, auto const &rhs){return (lhs.first < rhs.first);});
    }

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
        unsigned const splineIndex = std::distance(simBalSplines.begin(), simBalSplineIt);
        
        Chi2Bin curChi2Bin(method, firstBin, lastBin, ptLeadHist,
          (method == MultijetCrawlingBins::Method::MPF) ? balProfile : nullptr,
          sumProj, simBalSplineIt->second, std::pow(balProfileRebinned->GetBinError(binChi2), 2));


        // Add systematic variations for the newly constructed bin
        for (auto const &syst: dataVariations)
        {
            unsigned const systIndex = nuisanceDefs.Register(syst.first);
            curChi2Bin.AddDataSyst(systIndex,
              syst.second[0]->GetBinContent(binChi2), syst.second[1]->GetBinContent(binChi2));
        }

        for (auto const &syst: simVariations)
        {
            unsigned const systIndex = nuisanceDefs.Register(syst.first);
            auto const &splinePair = syst.second.at(splineIndex).second;
            curChi2Bin.AddSimSyst(systIndex, splinePair[0], splinePair[1]);
        }


        chi2Bins.emplace_back(curChi2Bin);
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


TGraphErrors MultijetCrawlingBins::ComputeResiduals(JetCorrBase const &corrector,
  Nuisances const &nuisances) const
{
    jetCache->Update(corrector);
    TGraphErrors graph(chi2Bins.size());

    for (unsigned i = 0; i < chi2Bins.size(); ++i)
    {
        auto const &chi2Bin = chi2Bins[i];
        double const simBalance = chi2Bin.MeanSimBalance(nuisances);
        graph.SetPoint(i, chi2Bin.MeanPt(), chi2Bin.MeanBalance(nuisances) / simBalance - 1.);
        graph.SetPointError(i, 0., chi2Bin.Uncertainty() / simBalance);
    }

    return graph;
}


unsigned MultijetCrawlingBins::GetDim() const
{
    return std::count(chi2BinMask.begin(), chi2BinMask.end(), true);
}


double MultijetCrawlingBins::Eval(JetCorrBase const &corrector, Nuisances const &nuisances) const
{
    jetCache->Update(corrector);
    double chi2 = 0.;
    
    for (unsigned i = 0; i < chi2Bins.size(); ++i)
    {
        if (chi2BinMask[i])
            chi2 += chi2Bins[i].Chi2(nuisances);
    }
    
    return chi2;
}


TH1D MultijetCrawlingBins::RecomputeBalanceData(JetCorrBase const &corrector,
  Nuisances const &nuisances) const
{
    // Read the binning with up to L2Res correction applied
    std::vector<double> binning;
    binning.reserve(chi2Bins.size() + 1);

    for (auto const &chi2Bin: chi2Bins)
        binning.emplace_back(chi2Bin.PtRange().first);

    binning.emplace_back(chi2Bins[chi2Bins.size() - 1].PtRange().second);


    // Apply the given L3Res correction to take into account the migration in pt of the leading jet
    for (auto &edge: binning)
        edge = corrector.Apply(edge);


    TH1D histBalance("MeanBalance", "", binning.size() - 1, binning.data());
    histBalance.SetDirectory(nullptr);

    jetCache->Update(corrector);

    for (unsigned i = 0; i < chi2Bins.size(); ++i)
    {
        histBalance.SetBinContent(i + 1, chi2Bins[i].MeanBalance(nuisances));
        histBalance.SetBinError(i + 1, chi2Bins[i].Uncertainty());
    }

    return histBalance;
}


TH1D MultijetCrawlingBins::RecomputeBalanceSim(JetCorrBase const &corrector,
  Nuisances const &nuisances) const
{
    std::vector<double> binning;
    binning.reserve(chi2Bins.size() + 1);

    for (auto const &chi2Bin: chi2Bins)
        binning.emplace_back(chi2Bin.PtRange().first);

    binning.emplace_back(chi2Bins[chi2Bins.size() - 1].PtRange().second);


    TH1D histBalance("MeanBalance", "", binning.size() - 1, binning.data());
    histBalance.SetDirectory(nullptr);

    // Update jet cache as this determines positions in pt at which the splines are evaluated
    jetCache->Update(corrector);

    for (unsigned i = 0; i < chi2Bins.size(); ++i)
        histBalance.SetBinContent(i + 1, chi2Bins[i].MeanSimBalance(nuisances));

    return histBalance;
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
