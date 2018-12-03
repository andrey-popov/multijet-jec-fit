/**
 * A unit test to check evaluation of the multijet loss function.
 * 
 * This is currently a draft. Will put a small ROOT file into the repository and check computed
 * values of the loss function against a reference.
 */

#include <FitBase.hpp>
#include <MultijetBinnedSum.hpp>
#include <MultijetCrawlingBins.hpp>

#include <iostream>
#include <memory>
#include <string>


class JetCorr: public JetCorrBase
{
public:
    JetCorr();
    
public:
    virtual double Eval(double pt) const override;
};


JetCorr::JetCorr():
    JetCorrBase(1)
{}


double JetCorr::Eval(double pt) const
{
    double const b = 1.;
    double const ptmin = 15.;
    return 1. + parameters[0] * std::log(pt / ptmin) +
      parameters[0] / b * (std::pow(pt / ptmin, -b) - 1);
}


int main()
{
    JetCorr jetCorr;
    Nuisances dummyNuisances;
    std::unique_ptr<MeasurementBase> lossFunc;
    
    std::string inputDir("https://aapopov.web.cern.ch/aapopov/jec_inputs/prototype/");
    std::string inputFile;
    
    
    std::cout << "MultijetBinnedSum\n\n";
    inputFile = inputDir + "multijet_BinnedSum.root";
    
    std::cout << "Loss function for pt balancing with various jet corrections:\n";
    lossFunc.reset(new MultijetBinnedSum(inputFile, MultijetBinnedSum::Method::PtBal));
    
    for (auto const &p: {-2e-2, -1e-2, -5e-3, 0., 5e-3, 1e-2, 2e-2})
    {
        jetCorr.SetParams({p});
        std::cout << "  " << lossFunc->Eval(jetCorr, dummyNuisances) << std::endl;
    }
    
    std::cout << "\nLoss function for MPF with various jet corrections:\n";
    lossFunc.reset(new MultijetBinnedSum(inputFile, MultijetBinnedSum::Method::MPF));
    
    for (auto const &p: {-2e-2, -1e-2, -5e-3, 0., 5e-3, 1e-2, 2e-2})
    {
        jetCorr.SetParams({p});
        std::cout << "  " << lossFunc->Eval(jetCorr, dummyNuisances) << std::endl;
    }
    
    
    std::cout << "\n\nMultijetCrawlingBins\n\n";
    inputFile = "/gridgroup/cms/popov/Analyses/JetMET/2018.10.18_Crowling-bins/Analysis/multijet.root";
    
    std::cout << "Loss function for pt balancing with various jet corrections:\n";
    lossFunc.reset(new MultijetCrawlingBins(inputFile, MultijetCrawlingBins::Method::PtBal));
    
    for (auto const &p: {-2e-2, -1e-2, -5e-3, 0., 5e-3, 1e-2, 2e-2})
    {
        jetCorr.SetParams({p});
        std::cout << "  " << lossFunc->Eval(jetCorr, dummyNuisances) << std::endl;
    }
    
    return EXIT_SUCCESS;
}
