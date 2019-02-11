#include <Nuisances.hpp>

#include <algorithm>
#include <iterator>
#include <sstream>
#include <stdexcept>


unsigned NuisanceDefinitions::GetIndex(std::string const &name) const
{
    auto const res = indices.find(name);
    
    if (res == indices.end())
    {
        std::ostringstream message;
        message << "NuisanceDefinitions::GetIndex: Parameter with name \"" << name <<
          "\" has not been registered.";
        throw std::runtime_error(message.str());
    }

    return res->second;
}


std::string const &NuisanceDefinitions::GetName(unsigned index) const
{
    if (index >= registeredNames.size())
    {
        std::ostringstream message;
        message << "NuisanceDefinitions::GetName: Requesting parameter with index " << index <<
          " while only " << registeredNames.size() << " parameters have been registered.";
        throw std::runtime_error(message.str());
    }
    
    return registeredNames[index];
}


std::vector<std::string> const &NuisanceDefinitions::GetNames() const
{
    return registeredNames;
}


unsigned NuisanceDefinitions::GetNumParams() const
{
    return registeredNames.size();
}


bool NuisanceDefinitions::operator==(NuisanceDefinitions const &other) const
{
    return indices == other.indices;
}


bool NuisanceDefinitions::operator!=(NuisanceDefinitions const &other) const
{
    return not(operator==(other));
}


unsigned NuisanceDefinitions::Register(std::string const &name)
{
    auto const res = indices.find(name);
    
    if (res == indices.end())
    {
        registeredNames.emplace_back(name);
        indices[name] = registeredNames.size() - 1;
        return registeredNames.size() - 1;
    }
    else
    {
        return res->second;
    }
}



Nuisances::Nuisances(NuisanceDefinitions const &definitions_):
    definitions(definitions_)
{
    values.resize(definitions.GetNumParams(), 0.);
}


double Nuisances::Eval() const
{
    double sum = 0.;
    
    for (auto const &v: values)
        sum += v * v;
    
    return sum;
}


NuisanceDefinitions const &Nuisances::GetDefinitions() const
{
    return definitions;
}


unsigned Nuisances::GetNumParams() const
{
    return values.size();
}


void Nuisances::SetValues(Nuisances const &source)
{
    if (definitions != source.definitions)
    {
        std::ostringstream message;
        message << "Nuisances::SetValues: Mismatched set of nuisances in the source.";
        throw std::runtime_error(message.str());
    }
    
    values = source.values;
}


void Nuisances::SetValues(double const *values_)
{
    for (unsigned i = 0; i < values.size(); ++i)
        values[i] = values_[i];
}


double Nuisances::operator[](std::string const &name) const
{
    return values[definitions.GetIndex(name)];
}


double &Nuisances::operator[](std::string const &name)
{
    return values[definitions.GetIndex(name)];
}


double Nuisances::operator[](unsigned index) const
{
    if (index >= values.size())
    {
        std::ostringstream message;
        message << "Nuisances::operator[]: Requesting parameter with index " << index <<
          " while only " << values.size() << " parameters have been registered.";
        throw std::runtime_error(message.str());
    }
    
    return values[index];
}


double &Nuisances::operator[](unsigned index)
{
    if (index >= values.size())
    {
        std::ostringstream message;
        message << "Nuisances::operator[]: Requesting parameter with index " << index <<
          " while only " << values.size() << " parameters have been registered.";
        throw std::runtime_error(message.str());
    }
    
    return values[index];
}

