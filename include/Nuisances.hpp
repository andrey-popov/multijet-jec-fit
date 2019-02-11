#pragma once

#include <map>
#include <string>
#include <vector>



/**
 * \class NuisanceDefinitions
 * \brief A class that maintains a list of registered nuisance parameters
 *
 * This class allows to register, at runtime, a set of nuisance parameters with arbitrary (but
 * unique) names. Parameters can be identified by their registration indices.
 */
class NuisanceDefinitions
{
public:
    /// Trivial constructor
    NuisanceDefinitions() = default;

public:
    /**
     * \brief Returns index of nuisance parameter with given name
     *
     * Throws an exception if no parameter with this name has been registered.
     */
    unsigned GetIndex(std::string const &name) const;

    /**
     * \brief Returns name of the nuisance parameter with given index
     *
     * Equivalent to GetNames().at(index) but provides a better error message.
     */
    std::string const &GetName(unsigned index) const;
    
    /// Returns a vector of names of nuisance parameters registered so far
    std::vector<std::string> const &GetNames() const;

    /// Returns number of nuisance parameteres registered so far
    unsigned GetNumParams() const;
    
    /**
     * \brief Checks if the two sets of nuisance parameters are identical
     *
     * Returns true if the two sets contain the same parameters and in the same order.
     */
    bool operator==(NuisanceDefinitions const &other) const;

    /// Comparison operator complementary to operator==
    bool operator!=(NuisanceDefinitions const &other) const;

    /**
     * \brief Registers a new nuisance parameter
     * 
     * Registers a new parameter with the given name and returns its index. If a parameter with the
     * same name has already been registered, only returns its index.
     */
    unsigned Register(std::string const &name);

private:
    /// Names of registered nuisance parameters
    std::vector<std::string> registeredNames;
    
    /// Map from names of nuisance parameters to their indices
    std::map<std::string, unsigned> indices;
};



/**
 * \class Nuisances
 * \brief A vector of values of nuisance parameters
 *
 * Nuisance parameters are defined by an instance of class NuisanceDefinitions. In this class they
 * are assigned concrete values, which can be accessed by index (preferred) or by name.
 */
class Nuisances
{
public:
    /**
     * Constructor from a list of registered nuisance parameters
     *
     * All nuisance parameters are initialized at zero.
     */
    Nuisances(NuisanceDefinitions const &definitions);
    
public:
    /**
     * \brief Evaluates chi^2 based on current values of the nuisance parameters
     * 
     * For all nuisance parameters, central values are 0 and input uncertainties are 1.
     */
    double Eval() const;
    
    /// Returns the underlying object with definitions of nuisance parameters
    NuisanceDefinitions const &GetDefinitions() const;

    /// Returns number of nuisance parameteres
    unsigned GetNumParams() const;
    
    /**
     * \brief Sets values of nuisance parameters to ones in the given copy
     * 
     * If the names of registered nuisance parameters in this and source are not identical, an
     * exception is thrown.
     */
    void SetValues(Nuisances const &source);
    
    /**
     * \brief Sets values of nuisance parameters to ones from the given buffer
     * 
     * Reads GetNumParams() elements from the buffer.
     */
    void SetValues(double const *values);
    
    /**
     * Returns value of the parameter with the given name
     *
     * Access by index should be preferred to this method.
     */
    double operator[](std::string const &name) const;
    
    /**
     * Returns reference to the parameter with the given name
     *
     * Access by index should be preferred to this method.
     */
    double &operator[](std::string const &name);
    
    /// Returns value of the parameter with the given index
    double operator[](unsigned index) const;
    
    /// Returns reference to the parameter with the given index
    double &operator[](unsigned index);
    
private:
    /// Object that keeps track of names and indices of nuisance parameters
    NuisanceDefinitions const definitions;

    /// Values of nuisance parameters
    std::vector<double> values;
};

