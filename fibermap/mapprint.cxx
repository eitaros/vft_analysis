#include <vector>
#include <string>
#include <iostream>
#include <algorithm>
#include <iterator>
#include <fstream>
#include <sstream>

using namespace std; 

#define udswitch true //falseで0~31と32~63反転

int changechnum(string chname){
    int num = -1;
    if(udswitch){
        num = ((int)chname.at(0) - 65)*8 + (int)(chname.at(1) - '0');
    } else {
        num = ((int)chname.at(0) - 65)*8 + (int)(chname.at(1) - '0');
        if(num>32){
            num -= 32;
        } else {
            num += 32;
        }
    }
    
    return num;
}

int changechnum2(string chname){
    int num = -1;
    if(udswitch){
        num = ((int)chname.at(0) - 65)*8 + (int)(chname.at(1) - '1');
    } else {
        num = ((int)chname.at(0) - 65)*8 + (int)(chname.at(1) - '1');
        if(num>31){
            num -= 32;
        } else {
            num += 32;
        }
    }
    
    return num;
}

int mapprint(){
    string ifname = "../honban0606.csv";
    ifstream ifs(ifname.c_str());
    if(!ifs){
        cout << "cant open csv" << endl;
        return 0;
    }

    vector<int> vmech(64, -1);
    vector<int> vmechr(64, -1);
    string input, tmp;
    int ho = 0;
    int ta = 0;
    int ch, chr, chnum;
    while(getline(ifs, input)){
        tmp = "";
        ta = 0;
        ch = -1;
        chr = -1;
        istringstream stream(input);
        while(getline(stream, tmp, ',')){
            if(ta==0 && ho!=0){
                ch = stoi(tmp);
            }
            if(ta==1 && ho!=0){
                chr = stoi(tmp);
            }
            if(ta==3 && ho!=0){
                chnum = changechnum2(tmp);

                vmech.at(chnum) = ch;
                vmechr.at(chnum) = chr;
            }
            ta++;
        }
        ho++;
    }
    ifs.close();


    string ofname;
    ofstream of;
    char fnum[8], numnum[8];
    for(int imppc=0; imppc<14; imppc++){
        snprintf(fnum,sizeof(fnum),"%02d",imppc);
        ofname = "map" + string(fnum) + ".txt";
        of.open(ofname.c_str());
        if(!of) cout << "hoge ";
        for(int i=0; i<8; i++){
            for(int j=0; j<8; j++){
                if(imppc==13 && (i*8+j)>31){
                    of << vmechr.at(i*8 + j) << " ";
                }else{
                    of << vmech.at(i*8 + j) << " ";
                }
            }
            of << endl;
        }
        of.close();
    }


    return 0;
}

int main(){
    mapprint();
    return 0;
}
