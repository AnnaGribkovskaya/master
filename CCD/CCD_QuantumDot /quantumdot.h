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
    void applyCoupledClusterDoubles();
    void applyCCD_Not_HF();

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
    double**** m_twoBodyElementsInHF;
    double**** create4dArray(int, int, int, int);
    double**** m_InitialAmplitudes;
    double**** m_Amplitudes;
    double**** m_Amplitudes_previous;
    void SetUpInitialAmplitudes();
    void UpdateTheInitialAmplitudes();
    void  UpdateTheAmplitudes();
    void SetUpTwoBobyMatrixForHartreeFock();
    int m_numOfIterations;
    void ComputeCorrelationEnergy();
    double  m_CorrelationEnergy;
    double m_InitialCorrelationEnergy;
    double m_ReferenceEnergyHF;
    void ComputeAmplitudes();
    void ComputeEpsilon();
    void ComputeInitialCorrelationEnergy();
    double**** m_Epsilon;
    void SetUpInitialAmplitudes_not_HF();
    double**** m_Amplitudes_notHF;
    double**** m_Amplitudes_previous_notHF;
    double**** m_InitialAmplitudes_notHF;
    void ComputeInitialCorrelationEnergy_not_HF();
    double m_InitialCorrelationEnergy_notHF;
    void UpdateTheInitialAmplitudes_not_HF();
    void ComputeAmplitudes_not_HF();


};

#endif // QUANTUMDOT_H


