#include <TStyle.h>
#include <TROOT.h>
#include <TSystem.h>

#include <vector>
#include <string>
#include <iostream>
#include <algorithm>
#include <iterator>
#include <set>
#include <filesystem>
#include <cstdlib>
#include <fstream>
#include <map>
#include <unordered_map>

#include "vftana.C"

using namespace std; 

void doana(int runnum=784){
    vftana ana(runnum,true);

    for(int imppc=0; imppc<14; imppc++){
        for(int ich=0; ich<64; ich++){
            //ana.Show_rawTDC(imppc,ich);
            //ana.Show_rawTOT(imppc,ich);
            //if(imppc!=0) ana.Show_ChRelation(imppc,ich);
            //ana.Show_tdcmultiplicity(imppc,ich);
        }
        //ana.Show_2Dhitcorrelation(imppc);
    }

    /* 
    for(int ich=0; ich<64; ich++){
        ana.Show_ChRelation(0,ich);
    } 
    */

    //ana.Show_layerMultiplicity();
    
    //ana.Show_layerRelation(true);
    ana.Show_layerRelation(false);

    for(int i=0; i<6; i++){
        //ana.Show_mppcconnection(i,i+1,true);
        ana.Show_mppcconnection(i,i+1,false);
    }
    //ana.Show_mppcconnection(6,0,true);
    ana.Show_mppcconnection(6,0,false);
    //ana.Show_mppcconnection(1,2,true);

    if(false){
        vector<int> in = {13, 5, 3, 1, 0, 8, 10};
        vector<int> ou = {4, 2, 7, 9, 11, 12, 6};
        vector<int> m0 = {4, 6, 10, 12, 13, 11, 7, 5, 2, 0, 8, 14, 15, 9, 1, 3, 18, 16, 24, 30, 31, 25, 17, 19, 20, 22, 26, 28, 29, 27, 23, 21};
        vector<int> m1 = {53, 55, 59, 61, 60, 58, 54, 52, 51, 49, 57, 63, 62, 56, 48, 50, 35, 33, 41, 47, 46, 40, 32, 34, 37, 39, 43, 45, 44, 42, 38, 36};
        string chnum;
        char num[8];
        pair<int, int> result;
        ofstream ofs(("./text/convertcheck/" + to_string(runnum) +"_convert0728.txt").c_str());
        for(int i=0; i<7; i++){
            for(int j=0; j<32; j++){
                chnum.clear();
                result = ana.Test_Combert_mppctofiber(in.at(i), m0.at(j));
                snprintf(num, sizeof(num),"%02d",result.second);
                chnum += to_string(result.first) + "," + string(num) + "  ";

                result = ana.Test_Combert_mppctofiber(in.at(i), m1.at(j));
                snprintf(num, sizeof(num),"%02d",result.second);
                chnum += to_string(result.first) + "," + string(num) + "  ";

                result = ana.Test_Combert_mppctofiber(ou.at(i), m0.at(j));
                snprintf(num, sizeof(num),"%02d",result.second);
                chnum += to_string(result.first) + "," + string(num) + "  ";

                result = ana.Test_Combert_mppctofiber(ou.at(i), m1.at(j));
                snprintf(num, sizeof(num),"%02d",result.second);
                chnum += to_string(result.first) + "," + string(num) + "  ";

                cout << chnum << endl;
                ofs << chnum << endl;
            }
        }
        ofs.close();
    }
    
    gSystem->Exit(1);
    return;
}