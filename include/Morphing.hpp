#pragma once

#include <TH1.h>

#include <vector>


/**
 * \class PointMorph
 * \brief Performs three-point morphing
 *
 * For a triplet of reference points, this class implements a smooth interpolation between them and
 * a linear extrapolation.
 */
class PointMorph
{
public:
    /// Trivial default constructor
    PointMorph() = default;

    /// Constructor from central, up, and down reference points
    PointMorph(double central, double up, double down);

public:
    /**
     * \brief Computes interpolated/extrapolated value
     * 
     * Reference central, up, and down values are reproduced for x = 0, +1, -1.
     */
    double Eval(double x) const;

    /**
     * \brief Implementation of three-point interpolation and extrapolation
     *
     * Perform the morphing based on the three reference points. The reference values are
     * reproduced for x = 0, +1, -1.
     */
    static double Morph(double central, double up, double down, double x);

    /// Alias for method Eval
    double operator()(double x) const;
    
    /// Smooth step function
    static double SmoothStep(double x);

private:
    /// Reference points
    double central, up, down;
};


/**
 * \class HistMorph
 * \brief Performs three-point morphing of histograms
 * 
 * This class performs a smooth interpolation between three histograms and a linear extrapolation.
 * Each bin is treated independetly.
 */
class HistMorph
{
public:
    /// Trivial default constructor
    HistMorph() = default;
    
    /// Constructor from central, up, and down reference points
    HistMorph(std::vector<double> const &central, std::vector<double> const &up,
      std::vector<double> const &down);
    
    /**
     * \brief Constructor from up and down reference points
     * 
     * Central reference points are set to zero.
     */
    HistMorph(std::vector<double> const &up, std::vector<double> const &down);
    
    /// Constructor from central, up, and down reference points represented with histograms
    HistMorph(TH1 const &central, TH1 const &up, TH1 const &down);
    
    /**
     * \brief Constructor from up and down reference points represented with histograms
     * 
     * Central refernce points are set to zero.
     */
    HistMorph(TH1 const &up, TH1 const &down);
    
public:
    /**
     * \brief Computes interpolated/extrapolated value in the given bin
     * 
     * Reference central, up, and down values are reproduced for x = 0, +1, -1. Bin index is
     * zero-based.
     */
    double Eval(unsigned bin, double x) const;
    
private:
    /// Morphing objects for individual bins
    std::vector<PointMorph> bins;
};

