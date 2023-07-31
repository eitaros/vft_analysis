#include <vector>
#include <string>
#include <iostream>
#include <algorithm>
#include <iterator>
#include <set>
#include <filesystem>
#include <fstream>

#include <TCanvas.h>
#include <TGraph.h>
#include <TStyle.h>
#include <TROOT.h>

using namespace std; 

int relation_entries(int runnum=785, int mpppcnum=0){
    gStyle->SetLabelSize(0.045, "XY");
    gStyle->SetTitleSize(0.05, "XY");
    gStyle->SetTitleOffset(0.9, "X");
    gStyle->SetTitleOffset(1.05, "Y");

    const int N = 64;
    double X[N];
    double Y[N];
    
    string ifname = "./text/chrelation/" + to_string(runnum) + "_" + to_string(mpppcnum) + ".txt";
    ifstream ifs(ifname.c_str());
    if(!ifs) return -1;
    double val, ich;
    int howmany = 0;
    while(ifs>>ich>>val){
        if(howmany >= N) return -2;
        X[howmany] = ich;
        Y[howmany] = val;
        //cout << howmany << " " << num << endl;
        howmany++;
    }  

    TGraph *gr = new TGraph(N,X,Y);
    string grtitle = to_string(runnum) + " mppc "+ to_string(mpppcnum) + " relation entries;mppc channel;entries/hits";
    gr->SetTitle(grtitle.c_str());
    TCanvas *c1 = new TCanvas();
    gPad->SetGrid();

    gr->SetMarkerColor(4);
    gr->SetMarkerStyle(34);
    gr->SetMarkerSize(2);

    gr->Draw("AP");

    string picname = "./pic/chrelation/" + to_string(runnum) + "_" + to_string(mpppcnum) + ".png";
    c1->Print(picname.c_str());
    
    return 0;
}

int main(int argc, char* argv[]){
    int runnum, mppcnum;
    if(argc!=3){
        runnum = 785;
        mppcnum = 0;
    }else{
        runnum = atoi(argv[1]);
        mppcnum = atoi(argv[2]);
    }

    if((runnum!=784 && runnum!=785) || (mppcnum<0||mppcnum>13)){
        return -1;
    }

    relation_entries(runnum,mppcnum);

    return 0;

}