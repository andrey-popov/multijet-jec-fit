#pragma once


#include <FitBase.hpp>

#include <Math/Functor.h>


ROOT::Math::Functor WrapLossFunction(CombLossFunction const *lossFunc);
