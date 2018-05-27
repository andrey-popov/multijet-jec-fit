#include <Nuisances.hpp>

#include <sstream>
#include <stdexcept>


Nuisances::Nuisances()
{}


double Nuisances::Eval() const
{
    double sum = 0.;
    
    for (auto const &v: values)
        sum += v * v;
    
    return sum;
}


unsigned Nuisances::GetNumParams() const
{
    return values.size();
}


unsigned Nuisances::Register(std::string const &name)
{
    auto const res = indices.find(name);
    
    if (res == indices.end())
    {
        values.push_back(0.);
        unsigned const index = values.size() - 1;
        indices[name] = index;
        return index;
    }
    else
    {
        return res->second;
    }
}


std::string const &Nuisances::GetName(unsigned index) const
{
    if (index >= values.size())
    {
        std::ostringstream message;
        message << "Nuisances::GetName: Requesting parameter with index " << index <<
          " while only " << values.size() << " parameters have been registered.";
        throw std::runtime_error(message.str());
    }
    
    for (auto const &it: indices)
    {
        if (it.second == index)
            return it.first;
    }
    
    // Not supposed to reach this point
    std::ostringstream message;
    message << "Nuisances::GetName: Failed to find name for parameter with index " << index << ".";
    throw std::logic_error(message.str());
}


void Nuisances::SetValues(Nuisances const &source)
{
    if (indices != source.indices)
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
    auto const res = indices.find(name);
    
    if (res == indices.end())
    {
        std::ostringstream message;
        message << "Nuisances::operator[]: Parameter with name \"" << name <<
          "\" has not been registered.";
        throw std::runtime_error(message.str());
    }
    
    return values[res->second];
}


double &Nuisances::operator[](std::string const &name)
{
    auto const res = indices.find(name);
    
    if (res == indices.end())
    {
        std::ostringstream message;
        message << "Nuisances::operator[]: Parameter with name \"" << name <<
          "\" has not been registered.";
        throw std::runtime_error(message.str());
    }
    
    return values[res->second];
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
