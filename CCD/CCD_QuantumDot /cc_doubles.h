#ifndef CC_DOUBLES_H
#define CC_DOUBLES_H
#include <vector>
#include <iostream>
#include "quantumstate.h"
#include "Coulomb_Functions.hpp"
#include <armadillo>

class CC_DOUBLES{
public:
    QuantumDot(int, double, int);
    int EnergyCutOff;
    int NumberOfParticles;
    double homega;
    void getQuantumDotStates();
    void diagonalizeHFMatrix();
    void getQuantumDotStatesNumber();
    //void applyHartreeFockMethod();

private:
    std::vector<QuantumState> m_shells;
    int m_sm = -1;
    const double m_s = 0.5;
    arma::mat m_C;              //Coefficient matrix
    int m_EnergyCutOff;
    //arma::mat m_DensityMatrix;
    arma::mat m_HF;             //Hartree-Fock matrix
    arma::vec eigval_previous;
    arma::vec eigval;
    arma::mat eigvec;
    arma::mat m_HOEnergies;
    char const* ResultsFile = "HF_energies";
    //void setUpStatesCartesian(int);
    void setUpStatesPolar(int, double h_omega, int);
    void setUpStatesPolarSorted(int, double h_omega, int);
    void setCoefficientMatrix(arma::mat);
    void computeHFmatrix(arma::mat);
    arma::mat computeDensityMatrix();
    //void CalculateNonIntEnergy();
    //double computeHartreeFockEnergyDifference();
    //void computeHartreeFockEnergy(arma::mat);
    //void writeToFile(double, int, int, double);
    void fillTwoBodyElements();
    double**** m_twoBodyElements;
    double**** create4dArray(int, int, int, int);
    //int m_numOfIterations;

    //Function for first iteration:

    //Function for TBME in FH basis




};



#endif // CC_DOUBLES_H
