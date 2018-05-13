#include "qdotchannelset.h"

qdotchannelset::qdotchannelset(){

}

void qdotchannelset::setUpChannels(generalSPclass * qsystem){
    qsys = qsystem;
    int Mlmax = 2*(qsys->getShellsStochastic() - 1);
    int Smax = 2;
    for(int Ml = -Mlmax; Ml <= Mlmax; Ml++){
          for(int S = -Smax; S <= Smax; S = S + 2){
            ChannelVariety.emplace_back(channel());
            for(int i = 0; i < qsys->getFermiLevel(); i++){
              for(int j = 0; j < qsys->getFermiLevel(); j++){
                qstate *QuantumState = qsys->sumState(i,j);
                if(Ml == QuantumState->m() && S == QuantumState->s() && i != j){
                  ChannelVariety.back().m_HoleHoleVec.emplace_back(channelindexpair());
                  ChannelVariety.back().m_HoleHoleVec.back().set(i, j);
                } delete QuantumState;
              }
            }
            for(int a = qsys->getFermiLevel(); a < qsys->getStatesStochastic(); a++){
              for(int b = qsys->getFermiLevel(); b < qsys->getStatesStochastic(); b++){
                qstate *QuantumState = qsys->sumState(a,b);
                if(Ml == QuantumState->m() && S == QuantumState->s() && a != b){
                  ChannelVariety.back().m_ParticleParticleVec.emplace_back(channelindexpair());
                  ChannelVariety.back().m_ParticleParticleVec.back().set(a, b);
                } delete QuantumState;
              }
            }
            for(int c = qsys->getFermiLevel(); c < qsys->getStatesStochastic(); c++){
              for(int k = 0; k < qsys->getFermiLevel(); k++){
                qstate *QuantumState = qsys->sumState(c,k);
                if(Ml == QuantumState->m() && S == QuantumState->s() && c != k){
                  ChannelVariety.back().m_ParticleHoleVec.emplace_back(channelindexpair());
                  ChannelVariety.back().m_ParticleHoleVec.back().set(c, k);
                } delete QuantumState;
              }
            }
            for(int c = qsys->getFermiLevel(); c < qsys->getStatesStochastic(); c++){
              for(int k = 0; k < qsys->getFermiLevel(); k++){
                qstate *QuantumState = qsys->substractState(c,k);
                if(Ml == QuantumState->m() && S == QuantumState->s() && c != k){
                  ChannelVariety.back().m_ParticleMinusHoleVec.emplace_back(channelindexpair());
                  ChannelVariety.back().m_ParticleMinusHoleVec.back().set(c, k);
                } delete QuantumState;
              }
            }
            for(int l = 0; l < qsys->getFermiLevel(); l++){
              for(int d = qsys->getFermiLevel(); d < qsys->getStatesStochastic(); d++){
                qstate *QuantumState = qsys->sumState(l,d);
                if(Ml == QuantumState->m() && S == QuantumState->s() && l != d){
                  ChannelVariety.back().m_HoleParticleVec.emplace_back(channelindexpair());
                  ChannelVariety.back().m_HoleParticleVec.back().set(l, d);
                } delete QuantumState;
              }
            }
            for(int l = 0; l < qsys->getFermiLevel(); l++){
              for(int d = qsys->getFermiLevel(); d < qsys->getStatesStochastic(); d++){
                qstate *QuantumState = qsys->substractState(l,d);
                if(Ml == QuantumState->m() && S == QuantumState->s() && l != d){
                  ChannelVariety.back().m_HoleMinusParticleVec.emplace_back(channelindexpair());
                  ChannelVariety.back().m_HoleMinusParticleVec.back().set(l, d);
                } delete QuantumState;
              }
            }
            /*
            for(int a = qsys->getFermiLevel(); a < qsys->getStatesStochastic(); a++){
              for(int b = qsys->getFermiLevel(); b < qsys->getStatesStochastic(); b++){
                qstate*  QuantumStateAB = qsys->sumState(a,b);
                if(Ml == QuantumStateAB->m() && S == QuantumStateAB->s() && a != b){
                  for(int i = 0; i < qsys->getFermiLevel(); i++){
                     qstate *QuantumState = qsys->sumSubstractState(a,b,i);
                     for(int j = 0; j < qsys->getFermiLevel(); j++){
                        qstate*  QuantumState1 = qsys->oneState(j);
                        if(   QuantumState->nx() == QuantumState1->nx()
                           && QuantumState->ny() == QuantumState1->ny()
                           && QuantumState->nz() == QuantumState1->nz()
                           && QuantumState->s()  == QuantumState1->s()) {
                            ChannelVariety.back().m_ParticlePlusParticleMinusHoleVec.emplace_back(channelindexpair());
                            ChannelVariety.back().m_ParticlePlusParticleMinusHoleVec.back().setThree(a, b, i);
                            ChannelVariety.back().m_HoleVec.emplace_back(channelindexpair());
                            ChannelVariety.back().m_HoleVec.back().setOne(j);
                         } delete QuantumState1;
                     } delete QuantumState;
                  }
                } delete QuantumStateAB;
              }
            }
          }
        }*/
      }
    }
}

void qdotchannelset::setUpChannelsQ3(generalSPclass * qsystem){
    qsys = qsystem;
    int Mlmax = 0; // AAAAAA dopisat' suda E na obolochke fermi
    int Smax = 1;
    for(int Ml = -Mlmax; Ml <= Mlmax; Ml++){
       for(int S = -Smax; S <= Smax; S++){
           ChannelVariety1.emplace_back(channel());
           for(int i = 0; i < qsys->getFermiLevel(); i++){
              qstate*  OneQS = qsys->oneState(i);
              if(Ml == OneQS->m() && S == OneQS->s()){
                for(int a = qsys->getFermiLevel(); a < qsys->getStatesStochastic(); a++){
                  for(int b = qsys->getFermiLevel(); b < qsys->getStatesStochastic(); b++){
                    for(int j = 0; j < qsys->getFermiLevel(); j++){
                      qstate *QuantumState = qsys->sumSubstractState(a,b,j);
                      if( QuantumState->m() == OneQS->m()
                       && QuantumState->s()  == OneQS->s()){
                          ChannelVariety1.back().m_ParticlePlusParticleMinusHoleVec.emplace_back(channelindexpair());
                          ChannelVariety1.back().m_ParticlePlusParticleMinusHoleVec.back().setThree(a, b, j);
                          ChannelVariety1.back().m_HoleVec.emplace_back(channelindexpair());
                          ChannelVariety1.back().m_HoleVec.back().setOne(i);
                      }
                    }
                  }
                }
              }
            }
        }
    }
}

void qdotchannelset::setUpChannelsQ4(generalSPclass * qsystem){
    qsys = qsystem;
    int Mlmax = qsys->getShellsStochastic() - 1;; // AAAAAA dopisat' suda E na obolochke fermi
    int Smax = 1;
    for(int Ml = -Mlmax; Ml <= Mlmax; Ml++){
       for(int S = -Smax; S <= Smax; S++){
           ChannelVariety2.emplace_back(channel());
           for(int a = 0; a < qsys->getStatesStochastic(); a++){
              qstate*  OneQS = qsys->oneState(a);
              if(Ml == OneQS->m() && S == OneQS->s()){
                for(int i = 0; i < qsys->getFermiLevel(); i++){
                  for(int j = 0; j < qsys->getFermiLevel(); j++){
                    for(int b = qsys->getFermiLevel(); b < qsys->getStatesStochastic(); b++){
                      qstate *QuantumState = qsys->sumSubstractState(i,j,b);
                      if( QuantumState->m() == OneQS->m()
                       && QuantumState->s()  == OneQS->s()){
                       ChannelVariety2.back().m_HolePlusHoleMinusParticleVec.emplace_back(channelindexpair());
                       ChannelVariety2.back().m_HolePlusHoleMinusParticleVec.back().setThree(i, j, b);
                       ChannelVariety2.back().m_ParticleVec.emplace_back(channelindexpair());
                       ChannelVariety2.back().m_ParticleVec.back().setOne(a);
                      }
                    }
                  }
               }
            }
        }
    }
}
