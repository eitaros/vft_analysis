#include <vector>
#include <string>
#include <iostream>
#include <algorithm>
#include <iterator>
#include <set>
#include <filesystem>

#include "vftana.C"

using namespace std; 

void relation(int runnum=784){
    vftana ana(runnum,true);

    for(int ich=0; ich<64; ich++){
        ana.Show_ChRelation(0,ich);
    }

    gSystem->Exit(1);
    return;
}