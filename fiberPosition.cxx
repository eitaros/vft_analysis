//#include <TH2.h>
//#include <TStyle.h>
//#include <TCanvas.h>

#include <vector>
#include <string>
#include <iostream>
#include <cmath>

using namespace std; 

int main(int argc, char* argv[]){
    //int mppc0=0, int ch0=0, int mppc1=0, int ch1=0
    const double inner_radius = 112.;
    const double inner_theta = 48;
    const double outer_radius = 126.;
    const double outer_theta = -53.5;
    const double dphi = 360./224.;
    const double zmax = 555.;

    double in_phi = atof(argv[1]);
    double z0, z1;
    double phi0, phi1;

    z0 = (M_PI*inner_radius*tan(inner_theta)/360.)*in_phi - (M_PI*inner_radius*tan(inner_theta)/360.)*dphi*0;
    z1 = (M_PI*outer_radius*tan(outer_theta)/360.)*in_phi - (M_PI*outer_radius*tan(outer_theta)/360.)*dphi*0;
    cout << z1 << endl;
    
    return 0;
}