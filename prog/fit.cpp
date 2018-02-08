/**
 * Fits for jet correction combinning multiple analyses.
 */

#include <JetCorrDefinitions.hpp>
#include <FitBase.hpp>
#include <MultijetBinnedSum.hpp>
#include <PhotonJetBinnedSum.hpp>
#include <PhotonJetRun1.hpp>
#include <ZJetRun1.hpp>

#include <Minuit2/Minuit2Minimizer.h>
#include <Math/Functor.h>
#include <TMath.h>

#include <boost/algorithm/string.hpp>
#include <boost/program_options.hpp>

#include <cmath>
#include <fstream>
#include <iostream>
#include <list>
#include <memory>
#include <string>


int main(int argc, char **argv)
{
    using namespace std;
    namespace po = boost::program_options;
    
    
    // Parse arguments
    po::options_description options("Allowed options");
    options.add_options()
      ("help,h", "Prints help message")
      ("balance,b", po::value<string>()->default_value("PtBal"),
        "Type of balance variable, PtBal or MPF")
      ("photonjet-run1", po::value<string>(), "Input file for photon+jet analysis, Run 1 style")
      ("photonjet-binnedsum", po::value<string>(),
        "Input file for photon+jet analysis, binned sum")
      ("zjet-run1", po::value<string>(), "Input file for Z+jet analysis, Run 1 style")
      ("multijet-binnedsum", po::value<string>(), "Input file for multijet analysis, binned sum")
      ("output,o", po::value<string>()->default_value("fit.out"),
        "Name for output file with results of the fit");
    
    po::variables_map optionsMap;
    
    po::store(
      po::command_line_parser(argc, argv).options(options).run(),
      optionsMap);
    po::notify(optionsMap);
    
    if (optionsMap.count("help"))
    {
        cerr << "Fits for jet correction combinning multiple analyses.\n";
        cerr << "Usage: fit [options]\n";
        cerr << options << endl;
        return EXIT_FAILURE;
    }
    
    
    bool useMPF = false;
    string balanceVar(optionsMap["balance"].as<string>());
    boost::to_lower(balanceVar);
    
    if (balanceVar == "mpf")
        useMPF = true;
    else if (balanceVar != "ptbal")
    {
        cerr << "Do not recognize balance variable \"" <<
          optionsMap["balance"].as<string>() << "\".\n";
        return EXIT_FAILURE;
    }
    
    
    // Construct all requested measurements
    list<unique_ptr<MeasurementBase>> measurements;
    
    if (optionsMap.count("photonjet-run1"))
        measurements.emplace_back(new PhotonJetRun1(optionsMap["photonjet-run1"].as<string>(),
          (useMPF) ? PhotonJetRun1::Method::MPF : PhotonJetRun1::Method::PtBal));
    
    if (optionsMap.count("photonjet-binnedsum"))
        measurements.emplace_back(new PhotonJetBinnedSum(
          optionsMap["photonjet-binnedsum"].as<string>(),
          (useMPF) ? PhotonJetBinnedSum::Method::MPF : PhotonJetBinnedSum::Method::PtBal));
    
    if (optionsMap.count("zjet-run1"))
        measurements.emplace_back(new ZJetRun1(optionsMap["zjet-run1"].as<string>(),
          (useMPF) ? ZJetRun1::Method::MPF : ZJetRun1::Method::PtBal));
    
    if (optionsMap.count("multijet-binnedsum"))
        measurements.emplace_back(new MultijetBinnedSum(
          optionsMap["multijet-binnedsum"].as<string>(),
          (useMPF) ? MultijetBinnedSum::Method::MPF : MultijetBinnedSum::Method::PtBal));
    
    if (measurements.empty())
    {
        cerr << "No measurements requested.\n";
        return EXIT_FAILURE;
    }
    
    
    // Construct an object to evaluate the loss function
    auto jetCorr = make_unique<JetCorrStd2P>();
    CombLossFunction lossFunc(move(jetCorr));
    
    for (auto const &measurement: measurements)
        lossFunc.AddMeasurement(measurement.get());
    
    unsigned const nPars = lossFunc.GetNumParams();
    
    
    // Create minimizer
    ROOT::Minuit2::Minuit2Minimizer minimizer;
    ROOT::Math::Functor func(&lossFunc, &CombLossFunction::EvalRawInput, nPars);
    minimizer.SetFunction(func);
    minimizer.SetStrategy(2);   // high quality
    minimizer.SetErrorDef(1.);  // error level for a chi2 function
    minimizer.SetPrintLevel(3);
    
    
    // Initial point
    for (unsigned i = 0; i < nPars; ++i)
        minimizer.SetVariable(i, "p" + to_string(i), 0., 1e-2);
    
    
    // Run minimization
    minimizer.Minimize();
    
    
    // Print results
    cout << "\n\n\e[1mSummary\e[0m:\n";
    cout << "  Status: " << minimizer.Status() << '\n';
    cout << "  Covariance matrix status: " << minimizer.CovMatrixStatus() << '\n';
    cout << "  Minimal value: " << minimizer.MinValue() << '\n';
    cout << "  NDF: " << lossFunc.GetNDF() << '\n';
    
    double const pValue = TMath::Prob(minimizer.MinValue(), lossFunc.GetNDF());
    cout << "  p-value: " << pValue << '\n';
    
    double const *results = minimizer.X();
    double const *errors = minimizer.Errors();
    cout << "  Parameters:\n";
    
    for (unsigned i = 0; i < nPars; ++i)
        cout << "    " << minimizer.VariableName(i) << ":  " << results[i] << " +- " <<
          errors[i] << "\n";
    
    
    // Save fit results in a text file
    string const resFileName(optionsMap["output"].as<string>());
    ofstream resFile(resFileName);
    
    resFile << "# Fitted parameters\n";
    
    for (unsigned i = 0; i < nPars; ++i)
        resFile << results[i] << " ";
    
    resFile << "\n\n# Covariance matrix:\n";
    
    for (unsigned i = 0; i < nPars; ++i)
    {
        for (unsigned j = 0; j < nPars; ++j)
            resFile << minimizer.CovMatrix(i, j) << " ";
        
        resFile << '\n';
    }
    
    resFile << "\n# Minimal chi^2, NDF, p-value:\n";
    resFile << minimizer.MinValue() << " " << lossFunc.GetNDF() << " " << pValue << '\n';
    
    resFile.close();
    
    
    cout << "\nResults saved to file \"" << resFileName << "\".\n";
    
    
    return EXIT_SUCCESS;
}
