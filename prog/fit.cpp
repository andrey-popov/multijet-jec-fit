/**
 * Fits jet correction using photon+jet and multijet data.
 */

#include <JetCorrDefinitions.hpp>
#include <FitBase.hpp>
#include <MultijetBinnedSum.hpp>
#include <PhotonJetRun1.hpp>
#include <ZJetRun1.hpp>

#include <Minuit2/Minuit2Minimizer.h>
#include <Math/Functor.h>
#include <TMath.h>

#include <cmath>
#include <fstream>
#include <iostream>
#include <memory>
#include <string>


int main(int argc, char **argv)
{
    using namespace std;
    
    
    // Parse arguments
    if (argc < 4 or argc > 5)
    {
        cerr << "Usage: fit inputPhotonJet.root inputMultijet.root inputZJet.root [method]\n";
        return EXIT_FAILURE;
    }
    
    bool useMPF = false;
    
    if (argc >= 5)
    {
        string const methodLabel(argv[4]);
        
        if (methodLabel == "MPF")
            useMPF = true;
        else if (methodLabel != "PtBal")
        {
            cerr << "Do not recognize method \"" << methodLabel << "\".\n";
            return EXIT_FAILURE;
        }
    }
    
    
    // Construct an object to evaluate the loss function
    auto jetCorr = make_unique<JetCorrStd2P>();
    CombLossFunction lossFunc(move(jetCorr));
    
    PhotonJetRun1 measurementPhotonJet(argv[1],
      (useMPF) ? PhotonJetRun1::Method::MPF : PhotonJetRun1::Method::PtBal);
    lossFunc.AddMeasurement(&measurementPhotonJet);
    
    ZJetRun1 measurementZJet(argv[3],
      (useMPF) ? ZJetRun1::Method::MPF : ZJetRun1::Method::PtBal);
    lossFunc.AddMeasurement(&measurementZJet);
    
    MultijetBinnedSum measurementMultijet(argv[2],
      (useMPF) ? MultijetBinnedSum::Method::MPF : MultijetBinnedSum::Method::PtBal);
    lossFunc.AddMeasurement(&measurementMultijet);
    
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
        minimizer.SetVariable(i, "p"s + to_string(i), 0., 1e-2);
    
    
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
    string const resFileName("fit.out");
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
