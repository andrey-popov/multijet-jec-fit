#include <PythonWrapping.hpp>


ROOT::Math::Functor WrapLossFunction(CombLossFunction const *lossFunc)
{
    return ROOT::Math::Functor(lossFunc, &CombLossFunction::EvalRawInput,
      lossFunc->GetNumParams());
}
