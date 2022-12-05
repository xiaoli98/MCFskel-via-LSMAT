//
// Created by boy on 04/12/22.
//

#ifndef MCFSKET_COLLAPSER_H
#define MCFSKET_COLLAPSER_H

#include "utils.hpp"

class Collapser{
private:
    MyMesh *m;
    double omega_L_0, omega_L;
    double omega_H_0, omega_H;
    double omega_M_0, omega_M;
    double zero_TH;
public:
    Collapser(){
        omega_L = omega_L_0 = 1;
        omega_H = omega_H_0 = 20;
        omega_M = omega_M_0 = 40;
        zero_TH = 1e-7;
    }
};


#endif //MCFSKET_COLLAPSER_H
