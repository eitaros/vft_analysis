#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

#include <string>
#include <vector>
#include <iostream>

#include "vftana.C"

using namespace std;

void anatest(int runnum=784){
    vftana a(runnum,false);

    a.Show_MultiHitEvent(0,0,0);

    return;
}
