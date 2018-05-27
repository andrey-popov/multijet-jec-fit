#pragma once

#include <map>
#include <string>
#include <vector>


/**
 * \struct Nuisances
 * \brief A class to describe nuisance parameters
 * 
 * Nuisance parameters are identified by names and must be registered beforehand. They can be
 * accessed by names or by indices.
 */
class Nuisances
{
public:
    /// Trivial constructor
    Nuisances();
    
public:
    /**
     * \brief Evaluates chi^2 based on current values of the nuisance parameters
     * 
     * For all nuisance parameters, central values are 0 and input uncertainties are 1.
     */
    double Eval() const;
    
    /**
     * \brief Finds name of the nuisance parameter with given index
     * 
     * Loops over all parameters and therefore can be slow.
     */
    std::string const &GetName(unsigned index) const;
    
    /// Returns number of registered nuisance parameteres
    unsigned GetNumParams() const;
    
    /**
     * \brief Registers a new nuisance parameter
     * 
     * Registers a new parameter with the given name and returns its index. The value of the
     * parameter is set to zero. If a parameter with the same name has already been registered,
     * only returns its index.
     */
    unsigned Register(std::string const &name);
    
    /**
     * \brief Sets values of nuisance parameters to ones in the given copy
     * 
     * Only values but not the names are updated. If the names of registered nuisance parameters
     * are not identical, an exception is thrown.
     */
    void SetValues(Nuisances const &source);
    
    /**
     * \brief Sets values of nuisance parameters to ones from the given buffer
     * 
     * Reads GetNumParams() elements from the buffer.
     */
    void SetValues(double const *values);
    
    /// Returns value of the parameter with the given name
    double operator[](std::string const &name) const;
    
    /// Returns reference to the parameter with the given name
    double &operator[](std::string const &name);
    
    /// Returns value of the parameter with the given index
    double operator[](unsigned index) const;
    
    /// Returns reference to the parameter with the given index
    double &operator[](unsigned index);
    
private:
    /// Map from names of nuisance parameters to their indices in vector values
    std::map<std::string, unsigned> indices;
    
    /// Values of nuisance parameters
    std::vector<double> values;
};
