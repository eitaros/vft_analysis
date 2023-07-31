#include <vector>
#include <string>
#include <iostream>
#include <algorithm>
#include <iterator>
#include <set>
#include <filesystem>

using namespace std; 

int relation_entries(int mpppcnum=0){
    gStyle->SetLabelSize(0.045, "XY");
    gStyle->SetTitleSize(0.05, "XY");
    gStyle->SetTitleOffset(0.9, "X");
    gStyle->SetTitleOffset(1.05, "Y");

    const int N = 64;
    double X[N];
    double Y[N];
    
    string ifname = "../text/chrelation/785_" + to_string(mpppcnum) + ".txt";
    ifstream ifs(ifname.c_str());
    if(!ifs) return -1;
    string num;
    int howmany = 0;
    while(getline(ifs,num)){
        if(howmany >= N) return -2;
        X[howmany] = (double)howmany;
        Y[howmany] = stod(num);
        //cout << howmany << " " << num << endl;
        howmany++;
    }  

    TGraph *gr = new TGraph(N,X,Y);
    string grtitle = "785 mppc "+ to_string(mpppcnum) + " relation entries;mppc channel;entries/hits";
    gr->SetTitle(grtitle.c_str());
    TCanvas *c1 = new TCanvas();

    gr->SetMarkerColor(4);
    gr->SetMarkerStyle(34);
    gr->SetMarkerSize(2);

    gr->Draw("AP");

    string picname = "../pic/chrelation/785_" + to_string(mpppcnum) + ".png";
    c1->Print(picname.c_str());
    
    return 0;
}