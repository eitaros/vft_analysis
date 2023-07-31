//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Tue Jul 11 18:12:12 2023 by ROOT version 6.26/06
// from TTree tree/vft raw data
// found on file: datas/ExVFT_00784.root
//////////////////////////////////////////////////////////

#ifndef vftana_h
#define vftana_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TStyle.h>


// Header file for the classes stored in the TTree if any.
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

using namespace std; 

/* 
53 55 59 61 60 58 54 52 
51 49 57 63 62 56 48 50 
35 33 41 47 46 40 32 34 
37 39 43 45 44 42 38 36 
04 06 10 12 13 11 07 05 
02 00 08 14 15 09 01 03 
18 16 24 30 31 25 17 19 
20 22 26 28 29 27 23 21
*/

class vftana {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
   Int_t           runnum;
   Int_t           evnum;
   vector<short>   *vft_layer;
   vector<short>   *vft_channel;
   vector<vector<double> > *vft_leading;
   vector<vector<double> > *vft_trailing;

   // List of branches
   TBranch        *b_runnum;   //!
   TBranch        *b_evnum;   //!
   TBranch        *b_vft_layer;   //!
   TBranch        *b_vft_channel;   //!
   TBranch        *b_vft_leading;   //!
   TBranch        *b_vft_trailing;   //!

   //
   int runnumber;
   bool saveflag;
   vector<double> totpair_leading;
   vector<double> totpair_trailing;
   vector<unordered_map<int, int> > mppcmap0;
   vector<unordered_map<int, int> > mppcmap1;
   unordered_map<int, int> innermppc;
   unordered_map<int, int> outermppc;
   vector<bool>          tdcisgoodhit;
   vector<vector<bool> > tdcisgoodvalue;
   vector<short> fiber_layer;
   vector<short> fiber_fiber;

   vftana(int id=784, bool flag=false);
   virtual ~vftana();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
   
   //initialise
   virtual void     Load_multihitevent();
   virtual void     Load_mppcmap(int mppcnum=0);
   virtual void     Initialize_values();
   virtual void     Load_layermap();
   //test
   virtual void     Test(int num=1);
   //correction, cut
   virtual void     Select_TOTpair(vector<double> l, vector<double> t);
   virtual int      Cut_byTDCValue(int cutmin=760, int cutmax=820);
   virtual void     Combert_mppctofiber();
   virtual pair<int, int> Test_Combert_mppctofiber(int mppcnum=0, int chnum=0);

   //raw 
   virtual void     Show_rawdatas(int evnum=0);
   virtual void     Show_rawTDC(int mppcnum=0, int chnum=0);
   virtual void     Show_rawTOT(int mppcnum=0, int chnum=0);
   //multiplicity
   virtual void     Show_tdcmultiplicity(int mppcnum=0, int chnum=0);
   virtual void     Show_layerMultiplicity();
   //for vft position
   virtual void     Show_ChRelation(int mppcnum=0, int chnum=0);
   virtual void     Show_2Dhitcorrelation(int mppcnum=0, bool reverseflag=false);
   virtual void     Show_layerRelation(bool isinner=true);
   virtual void     Show_mppcconnection(int mppc0=0,int mppc1=1, bool isinner=true);
};

#endif

#ifdef vftana_cxx
vftana::vftana(int id, bool flag) : fChain(0) 
{
   runnumber = id;
   saveflag = flag;
   TTree *tree = 0;
   char datafilenamechar[128];
   snprintf(datafilenamechar,sizeof(datafilenamechar),"./datas/ExVFT_%05d.root",runnumber);
   string datafilename = datafilenamechar;
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject(datafilename.c_str());
      if (!f || !f->IsOpen()) {
         f = new TFile(datafilename.c_str());
         //cout << "cant open data file" << endl;
      }
      f->GetObject("tree",tree);

   }
   Init(tree);

   TStyle *picst = new TStyle("vft_pic_style","vft pic style");

   picst->SetLabelSize(0.045, "XY");
   picst->SetTitleSize(0.05, "XY");
   picst->SetTitleOffset(0.9, "X");
   picst->SetTitleOffset(1.05, "Y");
   picst->cd(); 

   mppcmap0 = vector<unordered_map<int, int> >(14);
   mppcmap1 = vector<unordered_map<int, int> >(14);
   fiber_fiber.clear();
   fiber_layer.clear();

   Load_multihitevent();
   for(int imppc=0; imppc<14; imppc++){
      Load_mppcmap(imppc);
   }
   Load_layermap();
}

vftana::~vftana()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t vftana::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t vftana::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void vftana::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set object pointer
   vft_layer = 0;
   vft_channel = 0;
   vft_leading = 0;
   vft_trailing = 0;
   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("runnum", &runnum, &b_runnum);
   fChain->SetBranchAddress("evnum", &evnum, &b_evnum);
   fChain->SetBranchAddress("vft_layer", &vft_layer, &b_vft_layer);
   fChain->SetBranchAddress("vft_channel", &vft_channel, &b_vft_channel);
   fChain->SetBranchAddress("vft_leading", &vft_leading, &b_vft_leading);
   fChain->SetBranchAddress("vft_trailing", &vft_trailing, &b_vft_trailing);
   Notify();
}

Bool_t vftana::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void vftana::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t vftana::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}

void vftana::Load_multihitevent(){
   /* multihitevent.clear();

   string listfilename = "./config/" + to_string(runnumber) + "_multihitevent.txt";
   ifstream ifs(listfilename.c_str());
   if(!ifs){
      cout << "cant load multihiteventfile" << endl;
      multihitevent.clear();
      return;
   }
   int num;
   while(ifs >> num){
      multihitevent.insert(num);
   }
   ifs.close();

   cout << "load multihiteventfile" 
        << endl;
   
   cout << "multihit "
        << multihitevent.size() 
        << " : " 
        << fChain->GetEntriesFast() 
        << " " 
        << 100*multihitevent.size()/fChain->GetEntriesFast() 
        << "%" 
        << endl;
    */
   return;
}

void vftana::Load_mppcmap(int mppcnum){
   char ifnum[8];
   snprintf(ifnum,sizeof(ifnum),"%02d",mppcnum);
   string ifname = "./fibermap/map" + string(ifnum) + ".txt";
   ifstream ifs(ifname);
   if(!ifs) {
      cout << "cant load map" << endl;
      return;
   }
   
   vector<int> input;
   string line, word;
   int num, iline=0;
   while(getline(ifs,line)){
      istringstream sstream(line);
      while(getline(sstream,word,' ')){
         if(word == "\n") break;
         input.push_back(stoi(word));
      }
   }
   if(input.size()!=64){
      cout << "loading error" << endl;
      return;
   }
   ifs.close();

   for(int i=0; i<64; i++){
      if(i<32) mppcmap1.at(mppcnum).emplace(input.at(i), i);
      else     mppcmap0.at(mppcnum).emplace(input.at(i), i-32);
   }

   return;
}

void vftana::Load_layermap(){
   //const vector<int>    innermppcvec = {13, 3, 5, 1, 0, 8, 10};
   const vector<int>    innermppcvec = {10,8,0,1,3,5,13};
   //const vector<int>    outermppcvec = {4, 2, 7, 9, 11, 12, 6};
   const vector<int>    outermppcvec = {6, 12, 11, 9, 7, 2, 4};


   for(int i=0; i<7; i++){
      innermppc.insert(make_pair(innermppcvec.at(i),i));
      outermppc.insert(make_pair(outermppcvec.at(i),i));
   }
   
   return;
}

#endif // #ifdef vftana_cxx
