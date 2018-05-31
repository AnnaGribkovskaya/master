#ifndef QUANTUMDOT_H
#define QUANTUMDOT_H
#include <vector>
#include <iostream>
#include "quantumstate.h"
#include "Coulomb_Functions.hpp"
#include <armadillo>

class QuantumDot{

public:
    QuantumDot(int, double, int);
    int EnergyCutOff;
    int NumberOfParticles;
    double homega;
    void getQuantumDotStates();
    void diagonalizeHFMatrix();
    void getQuantumDotStatesNumber();
    void applyHartreeFockMethod();

private:
    std::vector<QuantumState> m_shells;
    int m_sm = -1;
    const double m_s = 0.5;
    arma::mat m_C;              //Coefficient matrix
    int m_EnergyCutOff;
    //arma::mat m_DensityMatrix;
    arma::mat m_HF;             //Hartree-Fock matrix
    arma::mat m_CCD_t1_old;     //Matrix for T1 amplitudes
    arma::mat m_CCD_t1;
    //arma::mat m_CCD_t2_old;     //Matrix for T21 amplitudes
    //arma::mat m_CCD_t2;
    arma::mat m_fm;             //f-matrix for CC in polar



    arma::vec eigval_previous;
    arma::vec eigval;
    arma::mat eigvec;
    arma::mat m_HOEnergies;
    char const* ResultsFile = "HF_energies";
    void setUpStatesCartesian(int);
    void setUpStatesPolar(int, double h_omega, int);
    void setUpStatesPolarSorted(int, double h_omega, int);
    void setCoefficientMatrix(arma::mat);
    void computeHFmatrix(arma::mat);
    arma::mat computeDensityMatrix();
    void CalculateNonIntEnergy();
    double computeHartreeFockEnergyDifference();
    void computeHartreeFockEnergy(arma::mat);
    void writeToFile(double, int, int, double);
    void fillTwoBodyElements();
    double**** m_twoBodyElements;

    double**** create4dArray(int, int, int, int);
    int m_numOfIterations;



    double computeInitialCCDCorrEnergy();
    void computeInitialCCDAmplitudes();
    double computeCCDCorrEnergy();
    double**** m_CCD_t_old;
    double**** m_CCD_t;
    void computeCCDAmplitudes();
    void applyCCDMethod();
    void updateOldAmplitudes();
    void computeInitialT1Amplitudes();
    void setUpFmatrix();
    void computeT1Amplitudes();
    void computeT2Amplitudes();




};

#endif // QUANTUMDOT_H


