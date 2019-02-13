#pragma once

#include <FitBase.hpp>

#include <Morphing.hpp>
#include <Nuisances.hpp>

#include <TGraphErrors.h>
#include <TH1.h>
#include <TH2.h>
#include <TProfile.h>
#include <TSpline.h>

#include <array>
#include <map>
#include <memory>
#include <set>
#include <string>
#include <vector>
#include <utility>


/**
 * \class MultijetCrawlingBins
 * 
 * Measurement in multijet topology with the approach with "crawling bins"
 * 
 * The divergence between data and simulation is computed as the chi^2 difference between mean
 * values of a balance observable, in bins of pt of the leading jet. The computation follows the
 * approach described in [1]. The mean balance in data in a given bin in pt of the leading jet,
 * which is referred to as a chi^2 bin, is computed using the binned-sum method. Due to the effect
 * of jet corrections, the position of the chi^2 bin along the axis of the pt of the leading jet
 * changes. The mean balance in simulation is approximated with a continuous function of pt of the
 * leading jet. Its value in a chi^2 bin is then computed using the shifted position of the bin.
 * 
 * In order to obtain a more accurate approximation, each chi^2 bin typically consists of a number
 * of narrow bins in pt of the leading jet. They are used internally in the computation. Underlying
 * histograms are filled separately for different triggers, which are mapped to non-overlapping
 * ranges in pt of the leading jet. Corresponding groups of chi^2 bins are referred to as trigger
 * bins.
 *
 * All systematic variations found in the input file are applied (separately for data and
 * simulation). However, user can disable selected ones by providing their names to the constructor.
 *
 * The class can also construct histograms with mean values of the chosen balance observables for
 * the given jet correction and set of nuisance parameters. This is done with methods
 * RecomputeBalanceData and RecomputeBalanceSim. The residual deviations can be computed using
 * method ComputeResiduals. Using this is the preferred way to visualize the performance of the fit.
 * 
 * [1] https://indico.cern.ch/event/780845/#16-multijet-analysis-with-craw
 */
class MultijetCrawlingBins: public MeasurementBase
{
public:
    /// Supported methods of computation
    enum class Method
    {
        PtBal,
        MPF
    };

    using Spline = TSpline3;
    
private:
    /**
     * \class JetCache
     * 
     * Auxiliary class that implements caching for jet corrections and weights
     * 
     * Jet corrections are computed on a two-dimensional rectangular grid. It first axis represents
     * typical pt in bins in the pt of the leading jet. The second axis represents typical pt in
     * bins in pt of any other jet in the event. This class provides access to thus precomputed jet
     * corrections.
     * 
     * The pt balance observable is defined using a smooth threshold: jets are included in the
     * computation with a certain weight that changes from 0 to 1 between two reference points.
     * With the grid approximation, this translates into weights along the second axis, which are
     * also cached.
     * 
     * Bin indices exposed in the interface of this class always follow the ROOT convention, i.e.
     * start from 1.
     */
    class JetCache
    {
    public:
        /**
         * Constructor
         * 
         * \param meanPtLead  Ordered vector of typical pt in bins along the first axis.
         * \param meanPtJet  Ordered vector of typical pt in bins along the second axis.
         * \param thresholdStart, thresholdEnd  Reference values of pt that define the smooth
         *     pt threshold.
         */
        JetCache(std::vector<double> const &meanPtLead, std::vector<double> const &meanPtJet,
          double thresholdStart, double thresholdEnd);
        
    public:
        /// Returns mean pt of the leading jet in the given bin with applied correction
        double CorrectedMeanPtLead(unsigned bin) const;
        
        /// Returns correction for typical pt in the given bin along the first axis
        double CorrectionPtLead(unsigned bin) const;
        
        /// Returns correction for typical pt in the given bin along the second axis
        double CorrectionPtJet(unsigned bin) const;
        
        /**
         * Returns range of bins with non-trivial content along the second axis
         * 
         * The returned pair consist of the first bin along the second axis in which the weight is
         * non-zero, and the last bin along the axis. When iterating over the second axis, this
         * information allows to skip immediately bins with zero weights.
         */
        std::pair<unsigned, unsigned> PtJetBinRange() const;
        
        /// Updates all cached values for the given correction
        void Update(JetCorrBase const &corrector);
        
        /// Returns weight for the given bin along the second axis
        double Weight(unsigned bin) const;
        
    private:
        /**
         * Computes weight for the given corrected pt
         * 
         * The weight changes smoothly from 0 below thresholdStart to 1 above thresholdEnd.
         */
        double JetWeight(double pt) const;
        
    private:
        /// Typical values of pt along the two axes
        std::vector<double> meanPtLead, meanPtJet;
        
        /// Reference points defining the smooth threshold for pt balance
        double thresholdStart, thresholdEnd;
        
        /// Cached values of the correction for typical pt values along the two axes
        std::vector<double> ptLeadCorrections, ptJetCorrections;
        
        /// Cached values of weights for bins along the second axis
        std::vector<double> jetWeights;
        
        /**
         * Range of bins along the second axis that have non-zero contribution
         * 
         * The boundaries of the range are included.
         */
        unsigned firstPtJetBin, lastPtJetBin;
    };
    
    /**
     * \class Chi2Bin
     * 
     * Auxiliary class representing a single bin contributing to chi^2
     * 
     * A single chi^2 bin aggregates one or (typically) more bins in pt of the leading jet in the
     * underlying dense binning. This class computes mean values of the balance observable in data
     * and simulation in the current chi^2 bin, as well as the corresponding chi^2 value. The
     * algorithm is described in the documentation for class MultijetCrawlingBins. Values of jet
     * corrections are accessed exclusively through a JetCache object.
     *
     * Systematic variations in the mean value of the balance observable in data and simulation are
     * supported. They need to be registered with methods AddDataSyst and AddSimSyst.
     */
    class Chi2Bin
    {
    public:
        /**
         * Constructor
         * 
         * \param method  Computation method, i.e. the observable to be used for chi^2.
         * \param firstBin, lastBin  Range of bins in pt of the leading jet in the histograms given
         *     as other arguments, that contribute to the current chi^2 bin.
         * \param ptLeadHist  Histogram of event counts in bins of pt of the leading jet in data.
         * \param mpfProfile  Profile with mean values of the MPF observable.
         * \param sumProj  Histogram of jet projections in data. See description of data member
         *     with the same name.
         * \param simBalSpline  Spline that approximates mean value of the balance observable in
         *     simulation. See description of data member with the same name.
         * \param unc2  Squared uncertainty to be used in the computation of chi^2.
         */
        Chi2Bin(Method method, unsigned firstBin, unsigned lastBin,
          std::shared_ptr<TH1> ptLeadHist, std::shared_ptr<TProfile> mpfProfile,
          std::shared_ptr<TH2> sumProj, std::shared_ptr<Spline> simBalSpline, double unc2);
        
    public:
        /**
         * Add systematic variation in data
         *
         * \param nuisanceIndex  Index of nuisance parameter that controls this variation.
         * \param up  Reference up relative deviation.
         * \param down  Reference down relative deviation.
         *
         * The variation is applied to the mean value of the balance observable in data, and it is
         * evaluated for the chi^2 bin as a whole.
         */
        void AddDataSyst(unsigned nuisanceIndex, double up, double down);

        /**
         * Add systematic variation in simulation
         *
         * \param nuisanceIndex  Index of nuisance parameter that controls this variation.
         * \param up  Spline that defines the reference up relative deviation for the mean value of
         *     the balance observable in simulation. Parameterized as a function of the logarithm of
         *     pt of the leading jet.
         * \param down  Spline that defines the reference down variation, similarly to parameter up.
         */
        void AddSimSyst(unsigned nuisanceIndex, std::shared_ptr<Spline> up,
          std::shared_ptr<Spline> down);

        /// Computes value of chi^2 in this bin
        double Chi2(Nuisances const &nuisances) const;
        
        /// Computes mean value of the balance observable in data in this chi^2 bin
        double MeanBalance(Nuisances const &nuisances) const;

        /**
         * Computes mean value of pt of the leading jet in this chi^2 bin
         *
         * The cached jet correction is applied.
         */
        double MeanPt() const;
        
        /**
         * Computes mean value of the balance observable in simulation in this chi^2 bin
         *
         * The computation takes into account the shift in the position of this chi^2 bin along pt,
         * which caused by the cached jet correction.
         */
        double MeanSimBalance(Nuisances const &nuisances) const;
        
        /// Returns the range in pt of the leading jet for this chi^2 bin
        std::pair<double, double> PtRange() const;
        
        /// Computes mean value of the balance observable in simulation at given pt
        double SimBalance(double const ptLead, Nuisances const &nuisances) const;
        
        /// Updates JetCache object used in the computations
        void SetJetCache(JetCache const *jetCache);

        /// Returns statistical uncertainty in data
        double Uncertainty() const;
    
    private:
        /// Implements computation of mean value of the MPF observable in data
        double MeanMPF(Nuisances const &nuisances) const;
        
        /// Implements computation of mean value of the pt balance observable in data
        double MeanPtBal(Nuisances const &nuisances) const;
        
    private:
        /**
         * Range of bins along pt of the leading jet that are included in this chi^2 bin
         * 
         * Both boundaries are included.
         */
        unsigned firstBin, lastBin;
        
        /// Histogram of event counts in bins of pt of the leading jet
        std::shared_ptr<TH1> ptLeadHist;
        
        /**
         * Profile with mean values of the MPF observable
         * 
         * Not set when computing the pt balance observable.
         */
        std::shared_ptr<TProfile> mpfProfile;
        
        /**
         * Histogram of jet projections
         * 
         * The axes of the histogram are the pt of the leading jet and pt of any other jet in the
         * event. It is filled with projection of pt of a jet along the direction opposed to the
         * diretion of the pt of the leading jet, normalized by pt of the leading jet.
         */
        std::shared_ptr<TH2> sumProj;
        
        /**
         * Mean value of the balance observable in simulation
         * 
         * Parameterized as a function of the natural logarithm of pt.
         */
        std::shared_ptr<Spline> simBalSpline;
        
        /// Squared uncertainty to be used for chi^2
        double unc2;
        
        /// Pointer to method to compute mean balance in data using the chosen observable
        double (Chi2Bin::*meanBalanceCalc)(Nuisances const &) const;
        
        /// Non-owning pointer to a JetCache object
        JetCache const *jetCache;

        /**
         * Registered systematic variations in data
         *
         * The key of the map is the index of the nuisance parameter that corresponds to the
         * variation. The value is a PointMorph object that defines the relative deviation in the
         * full chi^2 bin.
         */
        std::map<unsigned, PointMorph> dataVariations;

        /**
         * Registered systematic variations in simulation
         *
         * The key of the map is the index of the nuisance parameter that corresponds to the
         * variation. The value is a pair of splines that give up and down relative deviations in
         * the mean value of the balance observable. As with simBalSpline, these splines are
         * parameterized with the natural logarithm of pt.
         */
        std::map<unsigned, std::array<std::shared_ptr<Spline>, 2>> simVariations;
    };
    
public:
    /**
     * Constructor
     *
     * \param fileName  Path to ROOT file with inputs.
     * \param method  Computation method.
     * \param nuisanceDefs  Object that will collect requested nuisance parameters.
     * \param systToExclude  Labels of systematic uncertainties that should not be included.
     */
    MultijetCrawlingBins(std::string const &fileName, Method method,
      NuisanceDefinitions &nuisanceDefs, std::set<std::string> systToExclude = {});
    
public:
    /**
     * Computes data-to-simulation residuals for given jet correction and nuisances
     *
     * The residuals are evaluated, in a given chi^2 bin, as the ratio between mean values of the
     * balance observable in data and in simulation, minus 1. The point corresponding to each chi^2
     * bin is assigned as the x coordinate the mean value of pt of the leading jet, which is
     * computed taking the jet correction into account. All chi^2 bins are included in the returned
     * graph, regardless of their mask statuses. The uncertainty for each point is set according to
     * the input uncertainties of the mean values of the balance observable in data.
     */
    TGraphErrors ComputeResiduals(JetCorrBase const &corrector, Nuisances const &nuisances) const;

    /**
     * Returns number of chi^2 bins
     * 
     * Implemented from MeasurementBase.
     */
    virtual unsigned GetDim() const override;
    
    /**
     * Computes chi^2 for the given jet corrector and set of nuisances
     * 
     * Implemented from MeasurementBase.
     */
    virtual double Eval(JetCorrBase const &corrector, Nuisances const &nuisances) const override;

    /**
     * Recompute mean balance observable in data for given jet correction and nuisances
     *
     * Return a TH1D histogram that represents mean values of the balance observable as a function
     * of pt of the leading jet. The binning is as for the chi^2 bins, but the given jet correction
     * is applied to it so that the migration in pt of the leading jet is taken into account. All
     * chi^2 bins are included, regardless of their mask statuses. The uncertainty in each bin is
     * set to the statistical uncertainty of the input mean values; as a result, it is not affected
     * by the given jet correction. The histogram is not associated with any ROOT directory.
     */
    TH1D RecomputeBalanceData(JetCorrBase const &corrector, Nuisances const &nuisances) const;
    
    /**
     * Recompute mean balance observable in simulation for given jet correction and nuisances
     *
     * Return a TH1D histogram that represents mean values of the balance observable as a function
     * of pt of the leading jet. The same binning as for chi^2 bins is used. All chi^2 bins are
     * included, regardless of their mask statuses. The uncertainties are set to zero. The histogram
     * is not associated with any ROOT directory.
     */
    TH1D RecomputeBalanceSim(JetCorrBase const &corrector, Nuisances const &nuisances) const;
    
    /**
     * Restricts computation to given range in pt of the leading jet
     * 
     * Given boundaries are rounded to the closest boundaries of chi^2 bins. Returns the actual
     * range that will be used in the computation.
     */
    std::pair<double, double> SetPtLeadRange(double minPt, double maxPt);
    
private:
    /// Method of computation
    Method method;
    
    /**
     * All chi^2 bins
     * 
     * The vector is sorted in the increasing order in pt.
     */
    std::vector<Chi2Bin> chi2Bins;
    
    /**
     * Masks for chi^2 bins
     * 
     * Only chi^2 bins at positions where the mask value is true should be used in the computation
     * of the overall chi^2.
     */
    std::vector<bool> chi2BinMask;
    
    /// An object to cache values of jet corrections
    mutable std::unique_ptr<JetCache> jetCache;
};
