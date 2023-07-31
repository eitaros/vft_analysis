#define vftana_cxx
#include "vftana.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TLegend.h>


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

void vftana::Initialize_values(){
   totpair_leading.clear();
   totpair_trailing.clear();

   return;
}

void vftana::Select_TOTpair(vector<double> l, vector<double> t){
   totpair_leading.clear();
   totpair_trailing.clear();

   if(l.size()==0 || t.size()==0) return;
   
   vector<double> lcopy, tcopy;
   copy(l.begin(),l.end(),back_inserter(lcopy));
   copy(t.begin(),t.end(),back_inserter(tcopy));

   sort(lcopy.rbegin(),lcopy.rend());
   sort(tcopy.rbegin(),tcopy.rend());

   //pairing
   vector<vector<int> > indexpair;
   vector<int> p;
   int lsize = lcopy.size();
   int tsize = tcopy.size();
   double nowleading, nowtrailing;
   int iloop = 0;
   for(int i=0; i<lsize; i++){
      p.clear();
      nowleading = lcopy.at(i);
      
      while(iloop < tsize){
         nowtrailing = tcopy.at(iloop);

         if(nowleading > nowtrailing){
            p.push_back(i);
            p.push_back(iloop);
            indexpair.push_back(p);
            break;
         }else{
            iloop++;
         }
      }

      if(iloop == tsize) break;
   }

   //solve same index pair , collect pair
   for(int i=0; i<indexpair.size(); i++){
      if(i!=(indexpair.size()-1) && (indexpair.at(i).at(1) == indexpair.at(i+1).at(1))) continue;

      totpair_leading.push_back(lcopy.at(indexpair.at(i).at(0)));
      totpair_trailing.push_back(tcopy.at(indexpair.at(i).at(1)));
   
   }

   return;
}

void vftana::Combert_mppctofiber(){
   fiber_fiber.clear();
   fiber_layer.clear();
   fiber_layer.resize(vft_layer->size(),-1);
   fiber_fiber.resize(vft_layer->size(),-1);

   const short fiberoffset = -5;
   for(int ihit=0; ihit<vft_layer->size(); ihit++){
      if(innermppc.count(vft_layer->at(ihit))){
         if(vft_channel->at(ihit)<32){
            fiber_layer.at(ihit) = 0;
            fiber_fiber.at(ihit) = (short)mppcmap0.at(vft_layer->at(ihit)).at(vft_channel->at(ihit)) + (short)32*innermppc.at(vft_layer->at(ihit));
         }else{
            fiber_layer.at(ihit) = 1;
            fiber_fiber.at(ihit) = fiberoffset + (short)mppcmap1.at(vft_layer->at(ihit)).at(vft_channel->at(ihit)) + (short)32*innermppc.at(vft_layer->at(ihit));
            if(fiber_fiber.at(ihit) < 0) fiber_fiber.at(ihit) += 224;
         }
      }else if(outermppc.count(vft_layer->at(ihit))){
         if(vft_channel->at(ihit)<32){
            fiber_layer.at(ihit) = 2;
            fiber_fiber.at(ihit) = (short)mppcmap0.at(vft_layer->at(ihit)).at(vft_channel->at(ihit)) + (short)32*outermppc.at(vft_layer->at(ihit));
         }else{
            fiber_layer.at(ihit) = 3;
            fiber_fiber.at(ihit) = fiberoffset + (short)mppcmap1.at(vft_layer->at(ihit)).at(vft_channel->at(ihit)) + (short)32*outermppc.at(vft_layer->at(ihit));
            if(fiber_fiber.at(ihit) < 0) fiber_fiber.at(ihit) += 224;
         }
      }
   }
   
   return;
}

pair<int, int> vftana::Test_Combert_mppctofiber(int mppcnum, int chnum){
   if(mppcmap0.at(0).size()==0){
      for(int imppc=0; imppc<14; imppc++){
         Load_mppcmap(imppc);
      }
   }
   if(innermppc.size()==0 || outermppc.size()==0){
      Load_layermap();
   }

   int lnum=-1, fnum=-1;
   const int fiberoffset = -5;
   if(innermppc.count(mppcnum)){
      if(chnum<32){
         lnum = 0;
         fnum = mppcmap0.at(mppcnum).at(chnum) + 32*innermppc.at(mppcnum);
      }else{
         lnum = 1;
         fnum = fiberoffset + mppcmap1.at(mppcnum).at(chnum) + 32*innermppc.at(mppcnum);
         if(fnum < 0) fnum += 224;
      }
   }else if(outermppc.count(mppcnum)){
      if(chnum<32){
         lnum = 2;
         fnum = mppcmap0.at(mppcnum).at(chnum) + 32*outermppc.at(mppcnum);
      }else{
         lnum = 3;
         fnum = fiberoffset + mppcmap1.at(mppcnum).at(chnum) + 32*outermppc.at(mppcnum);
         if(fnum < 0) fnum += 224;
      }
   }else{
      fnum = -1;
      lnum = -1;
   }
   
   pair<int, int> result;
   //result = make_pair(fnum,lnum);
   result.first = lnum;
   result.second = fnum;
   return result;
}

int vftana::Cut_byTDCValue(int cutmin, int cutmax){
   tdcisgoodhit.clear();
   tdcisgoodvalue.clear();

   vector<bool> input;
   bool hitisgood;
   for(int ihit=0; ihit<vft_layer->size(); ihit++){
      hitisgood = false;
      input.clear();

      for(int itdc=0; itdc<vft_leading->at(ihit).size(); itdc++){
         if(vft_leading->at(ihit).at(itdc) > cutmin && vft_leading->at(ihit).at(itdc) < cutmax){
            input.push_back(true);
            hitisgood = true;
         }else{
            input.push_back(false);
         }
      }

      if(input.size()!=vft_leading->at(ihit).size()) return -1;
      tdcisgoodvalue.push_back(input);
      tdcisgoodhit.push_back(hitisgood);
   }

   if(tdcisgoodhit.size()!=vft_leading->size()) return -2;
   if(tdcisgoodvalue.size()!=vft_leading->size()) return -3;

   return 1;
}

void vftana::Loop()
{
//   In a ROOT session, you can do:
//      root> .L vftana.C
//      root> vftana t
//      root> t.GetEntry(12); // Fill t data members with entry number 12
//      root> t.Show();       // Show values of entry 12
//      root> t.Show(16);     // Read and show values of entry 16
//      root> t.Loop();       // Loop on all entries
//

//     This is the loop skeleton where:
//    jentry is the global entry number in the chain
//    ientry is the entry number in the current Tree
//  Note that the argument to GetEntry must be:
//    jentry for TChain::GetEntry
//    ientry for TTree::GetEntry and TBranch::GetEntry
//
//       To read only selected branches, Insert statements like:
// METHOD1:
//    fChain->SetBranchStatus("*",0);  // disable all branches
//    fChain->SetBranchStatus("branchname",1);  // activate branchname
// METHOD2: replace line
//    fChain->GetEntry(jentry);       //read all branches
//by  b_branchname->GetEntry(ientry); //read only this branch
   //
   if (fChain == 0) return;

   Long64_t nentries = fChain->GetEntriesFast();

   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      // if (Cut(ientry) < 0) continue;
   }
   //
}

void vftana::Test(int num){
   int mulflag = 0;
   int ev = -1;
   bool flag;
   int ltdiffer = 0;
   bool flag2, flag3;
   int mppc0hit = 0;
   
   //
   if (fChain == 0) return;

   Long64_t nentries = fChain->GetEntriesFast();

   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      // if (Cut(ientry) < 0) continue;
      Initialize_values();

      Combert_mppctofiber();
   }
   //

   //cout << mppc0hit << endl;

   return;
}

void vftana::Show_rawdatas(int evnum){
   if(evnum > (int)fChain->GetEntriesFast()){
      cout << "bad event number" << endl;
      return;
   }

   if (fChain == 0) return;
   Long64_t ientry = LoadTree((Long64_t)evnum);
   if (ientry < 0) return;
   Long64_t nb = fChain->GetEntry((Long64_t)evnum);
   
   if(vft_layer->size()==0) return;

   if(saveflag){
      string ofname = "./text/" + to_string(runnumber) + "_rawdatas.txt";
      ofstream ofs(ofname.c_str());
      ofs << evnum << " ======================" << endl;
      int noftdc;
      for(int ihit = 0; ihit<vft_layer->size(); ihit++){
         ofs << vft_layer << "  "
             << vft_channel << "  ";
         if(vft_leading->at(ihit).size()>vft_trailing->at(ihit).size()) noftdc = vft_leading->at(ihit).size();
         else vft_trailing->at(ihit).size();
         for(int itdc=0; itdc<noftdc; itdc++){
            if(itdc==0) ofs << "        ";
            if(vft_leading->at(ihit).size()<itdc) ofs << vft_leading->at(ihit).at(itdc) << "  ";
            else ofs << "     ";
            if(vft_trailing->at(ihit).size()<itdc) ofs << vft_trailing->at(ihit).at(itdc) << "  ";
            else ofs << "     ";
            ofs << endl;
         }
         ofs << endl;
      }
      ofs.close();
   }else{
      cout << evnum << " ======================" << endl;
      int noftdc;
      for(int ihit = 0; ihit<vft_layer->size(); ihit++){
         cout << vft_layer->at(ihit) << "  "
             << vft_channel->at(ihit) << "  ";
         if(vft_leading->at(ihit).size()>vft_trailing->at(ihit).size()) noftdc = vft_leading->at(ihit).size();
         else noftdc = vft_trailing->at(ihit).size();
         for(int itdc=0; itdc<noftdc; itdc++){
            if(itdc!=0) cout << "          ";
            else cout << "  ";
            if(vft_leading->at(ihit).size()>itdc) cout << vft_leading->at(ihit).at(itdc) << "  ";
            else cout << "     ";
            if(vft_trailing->at(ihit).size()>itdc) cout << vft_trailing->at(ihit).at(itdc) << "  ";
            else cout << "     ";
            cout << endl;
         }
         //cout << endl;
      }

      /* for(int i=0; i<vft_layer->size(); i++){
         cout << vft_layer->at(i) << endl;
      } */
   }
   
   return;
}

void vftana::Show_rawTDC(int mppcnum, int chnum){
   if(mppcnum<0 || mppcnum>13){
      cout << "bad mppcnum" <<endl;
      return;
   }
   if(chnum<0 || chnum>63){
      cout << "bad chnum" <<endl;
      return;
   }

   char hnamechar[64];
   snprintf(hnamechar,sizeof(hnamechar),"%02d_%02d",mppcnum,chnum);
   string hname, htitle;

   double tdcmax = 1000;
   double tdcmin = 0;
   double tdcbin = tdcmax - tdcmin;
   hname = "hleading_" + string(hnamechar);
   htitle = "hleading_" + string(hnamechar) + ";TDCchannel;count";
   TH1D *htdcl = new TH1D(hname.c_str(),htitle.c_str(),tdcbin,tdcmin,tdcmax);
   hname = "htrailing_" + string(hnamechar);
   htitle = "htrailing_" + string(hnamechar) + ";TDCchannel;count";
   TH1D *htdct = new TH1D(hname.c_str(),htitle.c_str(),tdcbin,tdcmin,tdcmax);
   
   //
   if (fChain == 0) return;

   Long64_t nentries = fChain->GetEntriesFast();

   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      // if (Cut(ientry) < 0) continue;
      Initialize_values();

      for(int ihit=0; ihit<vft_layer->size(); ihit++){
         if(vft_layer->at(ihit) == mppcnum && vft_channel->at(ihit) == chnum){
            for(int il=0; il<vft_leading->at(ihit).size(); il++){
               htdcl->Fill(vft_leading->at(ihit).at(il));
            }
            for(int it=0; it<vft_trailing->at(ihit).size(); it++){
               htdct->Fill(vft_trailing->at(ihit).at(it));
            }
            
            break;
         }
      }
   }
   //

   TCanvas *ctdc = new TCanvas("ctdc","ctdc");
   ctdc->Divide(2);
   ctdc->cd(1);
   htdcl->Draw();
   ctdc->cd(2);
   htdct->Draw();

   ctdc->UseCurrentStyle();
   
   if(saveflag){
      string picname = "./pic/rawhist/tdc/" + to_string(runnumber) + "_" + string(hnamechar) + ".png";
      ctdc->Print(picname.c_str());

      delete htdcl;
      delete htdct;
      delete ctdc;
   }

   return;
}

void vftana::Show_rawTOT(int mppcnum, int chnum){
   if(mppcnum<0 || mppcnum>13){
      cout << "bad mppcnum" <<endl;
      return;
   }
   if(chnum<0 || chnum>63){
      cout << "bad chnum" <<endl;
      return;
   }

   char hnamechar[64];
   snprintf(hnamechar,sizeof(hnamechar),"%02d_%02d",mppcnum,chnum);
   string hname, htitle;

   double totmax = 300;
   double totmin = 0;
   double totbin = totmax - totmin;
   hname = "htot_" + string(hnamechar);
   htitle = "tot_" + string(hnamechar) + ";TDCchannel;count";
   TH1D *htot = new TH1D(hname.c_str(),htitle.c_str(),totbin,totmin,totmax);

   int debug=0;
   //
   if (fChain == 0) return;

   Long64_t nentries = fChain->GetEntriesFast();

   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      // if (Cut(ientry) < 0) continue;
      Initialize_values();

      for(int ihit=0; ihit<vft_layer->size(); ihit++){
         if(vft_layer->at(ihit) == mppcnum && vft_channel->at(ihit) == chnum){
            Select_TOTpair(vft_leading->at(ihit), vft_trailing->at(ihit));
            
            if(totpair_leading.size()!=totpair_trailing.size()){
               debug++;
               break;
            }
            
            for(int itdc=0; itdc<totpair_leading.size(); itdc++){
               htot->Fill(totpair_leading.at(itdc) - totpair_trailing.at(itdc));
            }

            break;
         }
      }
   }
   //

   TCanvas *ctot = new TCanvas("ctot","ctot");
   htot->Draw();
   ctot->UseCurrentStyle();

   if(saveflag){
      string picname = "./pic/rawhist/tot/" + to_string(runnumber) + "_" + string(hnamechar) + ".png";
      ctot->Print(picname.c_str());

      delete htot;
      delete ctot;
   }

   cerr << debug << endl;
   
   return;
}

void vftana::Show_tdcmultiplicity(int mppcnum, int chnum){
   if(mppcnum<0 || mppcnum>13){
      cout << "bad mppcnum" <<endl;
      return;
   }
   if(chnum<0 || chnum>63){
      cout << "bad chnum" <<endl;
      return;
   }

   char hnamechar[64];
   snprintf(hnamechar,sizeof(hnamechar),"%02d_%02d",mppcnum,chnum);
   string hname, htitle;
   double hmax = 10;
   double hmin = 0;
   double hbin = hmax - hmin;
   hname = "h_multi_l_" + string(hnamechar);
   htitle = "multi_l_" + string(hnamechar) + ";nofhit;count";
   TH1D *hmultil = new TH1D(hname.c_str(),htitle.c_str(),hbin,hmin,hmax);
   hname = "h_multi_t_" + string(hnamechar);
   htitle = "multi_t_" + string(hnamechar) + ";nofhit;count";
   TH1D *hmultit = new TH1D(hname.c_str(),htitle.c_str(),hbin,hmin,hmax);

   //
   if (fChain == 0) return;

   Long64_t nentries = fChain->GetEntriesFast();

   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      // if (Cut(ientry) < 0) continue;
      Initialize_values();

      for(int ihit=0; ihit<vft_layer->size(); ihit++){
         if(vft_layer->at(ihit)>mppcnum) break;
         if(vft_layer->at(ihit) == mppcnum && vft_channel->at(ihit) == chnum){
            hmultil->Fill(vft_leading->at(ihit).size());
            hmultit->Fill(vft_trailing->at(ihit).size());
         }
      }
   }
   //

   TCanvas *cmulti = new TCanvas("cmulti,cmulti");
   gPad->SetLogy();
   hmultil->SetLineColor(4);
   hmultit->SetLineColor(2);
   hmultil->Draw();
   hmultit->Draw("same");
   TLegend *leg = new TLegend(0.8, 0.68, 0.99, 0.78);
   leg->AddEntry(hmultil,"leading","l");
   leg->AddEntry(hmultit,"trailing","l");
   leg->Draw("same");

   if(saveflag){
      string picname = "./pic/multihit/hists/" + to_string(runnumber) + "_" + string(hnamechar) + ".png";
      cmulti->Print(picname.c_str());

      delete cmulti;
      delete hmultil;
      delete hmultit;
      delete leg;
   }

   return;
}

void vftana::Show_layerMultiplicity(){
   int hmax = 100;
   int hmin = 0;
   int hbin = hmax - hmin;
   string hname, htitle;
   TH1I *hlmulti[4], *hlmulti_cut0[4];
   for(int i=0; i<4; i++){
      hname = "hlayer" + to_string(i) + "multi";
      htitle = "layer" + to_string(i) + "multiplicity;multiplicity;count";
      hlmulti[i] = new TH1I(hname.c_str(),htitle.c_str(),hbin,hmin,hmax);
      hlmulti[i]->SetLineColor(4);
      hname = "hlayer" + to_string(i) + "multi_cut";
      htitle = "layer" + to_string(i) + "multiplicity_cut;multiplicity;count";
      hlmulti_cut0[i] = new TH1I(hname.c_str(),htitle.c_str(),hbin,hmin,hmax);
      hlmulti_cut0[i]->SetLineColor(2);
   }
   int nofhitlayer[4] = {0,0,0,0};
   int nofhitlayer_cut0[4] = {0,0,0,0};
   bool fillflag_cut0;
   int cutresult;

   //
   if (fChain == 0) return;

   Long64_t nentries = fChain->GetEntriesFast();

   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      // if (Cut(ientry) < 0) continue;
      Initialize_values();
      for(int i=0; i<4; i++){
         nofhitlayer[i] = 0;
         nofhitlayer_cut0[i] = 0;
      }

      cutresult = Cut_byTDCValue(760, 820);
      if(cutresult<0){
         cout << "cut failed" << endl;
         break;
      }

      for(int ihit=0; ihit<vft_layer->size(); ihit++){
         /* if(innermppc.count(vft_layer->at(ihit))){
            nofhitlayer[vft_channel->at(ihit)/32] += 1;
         }else{
            nofhitlayer[vft_channel->at(ihit)/32 + 2] += 1;
         } */
         if(innermppc.count(vft_layer->at(ihit)) && mppcmap0.at(vft_layer->at(ihit)).count(vft_channel->at(ihit))) nofhitlayer[0] += 1;
         else if(innermppc.count(vft_layer->at(ihit)) && mppcmap1.at(vft_layer->at(ihit)).count(vft_channel->at(ihit))) nofhitlayer[1] += 1;
         else if(outermppc.count(vft_layer->at(ihit)) && mppcmap0.at(vft_layer->at(ihit)).count(vft_channel->at(ihit))) nofhitlayer[2] += 1;
         else if(outermppc.count(vft_layer->at(ihit)) && mppcmap1.at(vft_layer->at(ihit)).count(vft_channel->at(ihit))) nofhitlayer[3] += 1;
         
         if(tdcisgoodhit.at(ihit)){
            if(innermppc.count(vft_layer->at(ihit)) && mppcmap0.at(vft_layer->at(ihit)).count(vft_channel->at(ihit))) nofhitlayer_cut0[0] += 1;
            else if(innermppc.count(vft_layer->at(ihit)) && mppcmap1.at(vft_layer->at(ihit)).count(vft_channel->at(ihit))) nofhitlayer_cut0[1] += 1;
            else if(outermppc.count(vft_layer->at(ihit)) && mppcmap0.at(vft_layer->at(ihit)).count(vft_channel->at(ihit))) nofhitlayer_cut0[2] += 1;
            else if(outermppc.count(vft_layer->at(ihit)) && mppcmap1.at(vft_layer->at(ihit)).count(vft_channel->at(ihit))) nofhitlayer_cut0[3] += 1;
         }
      }

      for(int i=0; i<4; i++){
         hlmulti[i]->Fill(nofhitlayer[i]);
         hlmulti_cut0[i]->Fill(nofhitlayer_cut0[i]);
      }
   }
   //

   TCanvas *clmulti = new TCanvas("clmulti","clmulti");
   clmulti->UseCurrentStyle();
   clmulti->Divide(2,2);
   for(int i=0; i<4; i++){
      clmulti->cd(i+1);
      gPad->SetLogy(1);
      hlmulti[i]->Draw();
      hlmulti_cut0[i]->Draw("same");
   }
   //clmulti->UseCurrentStyle();

   if(saveflag){
      string picname = "./pic/multihit/" + to_string(runnumber) + "layermulti.png";
      clmulti->Print(picname.c_str());

      delete clmulti;
      for(int i=0; i<4; i++){
         delete hlmulti[i];
      }
   }
   
   return;
}

void vftana::Show_ChRelation(int mppcnum, int chnum){
   if(mppcnum<0 || mppcnum>13){
      cout << "bad mppcnum" <<endl;
      return;
   }
   if(chnum<0 || chnum>63){
      cout << "bad chnum" <<endl;
      return;
   }

   char hnamechar[64];
   snprintf(hnamechar,sizeof(hnamechar),"%02d_%02d",mppcnum,chnum);
   string hname, htitle;

   double histmax = 64;
   double histmin = 0;
   double histbin = histmax - histmin;
   hname = "hrelation_" + string(hnamechar);
   htitle = "hitwith_" + string(hnamechar) + ";MPPCchannel;count";
   TH1I *hrelation = new TH1I(hname.c_str(),htitle.c_str(),histbin,histmin,histmax);

   bool fillflag = false;
   int mppcnofhit = 0;
   int targetchhit = 0;

   //
   if (fChain == 0) return;

   Long64_t nentries = fChain->GetEntriesFast();

   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      // if (Cut(ientry) < 0) continue;
      Initialize_values();
      fillflag = false;
      mppcnofhit = 0;

      if(vft_layer->size()!=4) continue;

      for(int ihit=0; ihit<vft_layer->size(); ihit++){
         if(vft_layer->at(ihit)>mppcnum) break;
         if(vft_layer->at(ihit) == mppcnum){
            mppcnofhit++;
            if(vft_channel->at(ihit) == chnum){
               targetchhit++;
               fillflag = true;
            }
         }
      }

      if(fillflag){
         for(int ihit=0; ihit<vft_layer->size(); ihit++){
            if(vft_layer->at(ihit)>mppcnum) break;
            if(vft_layer->at(ihit) == mppcnum && vft_channel->at(ihit) != chnum){
               hrelation->Fill(vft_channel->at(ihit));
            }
         }
      }
   }
   //

   TCanvas *crelation = new TCanvas("crelation","crelation");
   hrelation->Draw();
   crelation->UseCurrentStyle();

   if(saveflag){
      string picnme = "./pic/chrelation/eachhist/" + to_string(runnumber) + "_" + string(hnamechar) + ".png";
      crelation->Print(picnme.c_str());

      string textname = "./text/chrelation/" + to_string(runnumber) + "_" + to_string(mppcnum) + ".txt";
      ofstream ofs(textname.c_str(),std::ios::app);
      ofs << chnum << " " << (double)hrelation->GetEntries()/(double)targetchhit << endl;
      ofs.close();


      delete hrelation;
      delete crelation;
   }

   
   return;
}

void vftana::Show_2Dhitcorrelation(int mppcnum, bool reverseflag){
   if(mppcnum<0 || mppcnum>13){
      cout << "bad mppcnum" <<endl;
      return;
   }
   if(mppcmap0.at(mppcnum).size()!=32 || mppcmap1.at(mppcnum).size()!=32){
      cout << "map load fialed" <<endl;
      return;
   }

   char hnamechar[64];
   snprintf(hnamechar,sizeof(hnamechar),"%02d",mppcnum);
   string hname, htitle;
   int hmax = 64;
   int hmin = 0;
   int hbin = hmax - hmin;
   hname = "hhitcorrelation" + string(hnamechar);
   htitle = "hitcorrelation" + string(hnamechar);
   TH2I *hhitcor = new TH2I(hname.c_str(),htitle.c_str(),32,0,32,32,0,32);

   int x=-1, y=-1, mppcnofhit=0;

   //
   if (fChain == 0) return;

   Long64_t nentries = fChain->GetEntriesFast();

   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      // if (Cut(ientry) < 0) continue;
      Initialize_values();
      mppcnofhit = 0;
      x = -1;
      y = -1;

      if(vft_layer->size()!=4) continue;

      for(int ihit=0; ihit<vft_layer->size(); ihit++){
         if(vft_layer->at(ihit)>mppcnum) break;
         if(vft_layer->at(ihit) == mppcnum){
            mppcnofhit++; 
         }
      }

      if(mppcnofhit==2){
         for(int ihit=0; ihit<vft_layer->size(); ihit++){
            if(vft_layer->at(ihit)>mppcnum) break;
            if(vft_layer->at(ihit) == mppcnum){
               if(vft_channel->at(ihit)<32) x = vft_channel->at(ihit);
               else y = vft_channel->at(ihit);
            }
         }
         if(x!=-1 && y!=-1){
            hhitcor->Fill(mppcmap0.at(mppcnum).at(x), mppcmap1.at(mppcnum).at(y));
         }
      }
   }
   //

   TCanvas *chitcor = new TCanvas("chitcor","chitcor");
   hhitcor->SetStats(0);
   hhitcor->Draw("colz");

   if(saveflag){
      string picname = "./pic/chrelation/" + to_string(runnumber) + "_" + string(hnamechar) + ".png";
      chitcor->Print(picname.c_str());

      delete chitcor;
      delete hhitcor;
   }
   
   return;
}

void vftana::Show_layerRelation(bool isinner){
   string hname, htitle;
   int hmax = 224;
   int hmin = 0;
   int hbin = hmax - hmin;
   if(isinner){
      hname = "hinnerlayer";
      htitle = "inner layer fiber correlation;U fiber;U' fiber";
   }else{
      hname = "houterlayer";
      htitle = "outer layer fiber correlation;V fiber;V' fiber";
   }
   TH2I *hfibercor = new TH2I(hname.c_str(),htitle.c_str(),hbin,hmin,hmax,hbin,hmin,hmax);

   int noflayerhit, x, y;

   //
   if (fChain == 0) return;

   Long64_t nentries = fChain->GetEntriesFast();

   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      // if (Cut(ientry) < 0) continue;
      Initialize_values();
      Combert_mppctofiber();
      noflayerhit = 0;
      x = -1;
      y = -1;
      
      for(int ihit=0; ihit<vft_layer->size(); ihit++){
         if(isinner && innermppc.count(vft_layer->at(ihit))) noflayerhit++;
         if(!isinner && outermppc.count(vft_layer->at(ihit))) noflayerhit++;
      }

      if(noflayerhit == 2){
         if(true){
            for(int ihit=0; ihit<vft_layer->size(); ihit++){
               if(isinner && innermppc.count(vft_layer->at(ihit))){
                  if(vft_channel->at(ihit)<32) x = mppcmap0.at(vft_layer->at(ihit)).at(vft_channel->at(ihit)) + 32*innermppc[vft_layer->at(ihit)];
                  else                         y = mppcmap1.at(vft_layer->at(ihit)).at(vft_channel->at(ihit)) + 32*innermppc[vft_layer->at(ihit)];
               }
               if(!isinner && outermppc.count(vft_layer->at(ihit))){
                  if(vft_channel->at(ihit)<32) x = mppcmap0.at(vft_layer->at(ihit)).at(vft_channel->at(ihit)) + 32*outermppc[vft_layer->at(ihit)];
                  else                         y = mppcmap1.at(vft_layer->at(ihit)).at(vft_channel->at(ihit)) + 32*outermppc[vft_layer->at(ihit)];
               }
            }
         }else{
            for(int ihit=0; ihit<vft_layer->size(); ihit++){
               if(fiber_layer.at(ihit)%2==0) x = fiber_fiber.at(ihit);
               if(fiber_layer.at(ihit)%2==1) y = fiber_fiber.at(ihit);
            }
         }
         if(x!=-1 && y!=-1) hfibercor->Fill(x,y);
      }
   }
   //

   TCanvas *cfibercor = new TCanvas("cfibercor","cfibercor");
   cfibercor->UseCurrentStyle();
   hfibercor->SetStats(0);
   hfibercor->Draw("colz");

   if(saveflag){
      string picname;
      if(isinner) picname = "./pic/chrelation/" + to_string(runnumber) + "_inner.png";
      else        picname = "./pic/chrelation/" + to_string(runnumber) + "_outer.png";
      cfibercor->Print(picname.c_str());

      delete cfibercor;
      delete hfibercor;
   }

   
   return;
}

void vftana::Show_mppcconnection(int mppc0, int mppc1, bool isinner){
   if(!((mppc0-mppc1)==-1 || (mppc0==6&&mppc1==0))){
      cout << "bad argument" << endl;
      return;
   }
   int hmax = 224;
   int hmin = 0;
   int hbin = hmax - hmin;
   string hname = (isinner? "inner_":"outer_") + to_string(mppc0) + "_" + to_string(mppc1);
   string htitle = to_string(runnumber) + "_" + hname + (isinner? ";U fiber;U' fiber":";V fiber;V' fiber");
   TH2I *h2mppccor = new TH2I(hname.c_str(),htitle.c_str(),hbin,hmin,hmax,hbin,hmin,hmax);
   
   int nofhit=0, x, y;
   int nlx, nfx, nly, nfy;
   string connectchinfo;
   ofstream osf(("./text/convertcheck/"+to_string(runnumber)+"_1_2ch.txt").c_str());

   //
   if (fChain == 0) return;

   Long64_t nentries = fChain->GetEntriesFast();

   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      // if (Cut(ientry) < 0) continue;
      Initialize_values();
      //if(vft_layer->size()!=4) continue;

      Combert_mppctofiber();
      nofhit=0;
      x = -1; y = -1;
      nlx = -1; nly = -1; nfx = -1; nfy = -1;
      connectchinfo.clear();
      
      if(vft_layer->size()!=fiber_layer.size()) continue;
      for(int ihit=0; ihit<fiber_layer.size(); ihit++){
         if(fiber_layer.at(ihit)==-1 || fiber_fiber.at(ihit)==-1) continue;

         if(isinner){
            if(!innermppc.count(vft_layer->at(ihit))) continue;
            if((innermppc.at(vft_layer->at(ihit))==mppc0 || innermppc.at(vft_layer->at(ihit))==mppc1)
               && fiber_layer.at(ihit)==0){
               nofhit++;
               x = fiber_fiber.at(ihit);
               nlx = vft_layer->at(ihit);
               nfx = vft_channel->at(ihit);
            }else if((innermppc.at(vft_layer->at(ihit))==mppc0 || innermppc.at(vft_layer->at(ihit))==mppc1)
                     && fiber_layer.at(ihit)==1){
               nofhit++;
               y = fiber_fiber.at(ihit);
               nly = vft_layer->at(ihit);
               nfy = vft_channel->at(ihit);
            }
         }else{
            if(!outermppc.count(vft_layer->at(ihit))) continue;
            if((outermppc.at(vft_layer->at(ihit))==mppc0 || outermppc.at(vft_layer->at(ihit))==mppc1)
               && fiber_layer.at(ihit)==2){
               nofhit++;
               x = fiber_fiber.at(ihit);
            }else if((outermppc.at(vft_layer->at(ihit))==mppc0 || outermppc.at(vft_layer->at(ihit))==mppc1)
                     && fiber_layer.at(ihit)==3){
               nofhit++;
               y = fiber_fiber.at(ihit);
            }
         }
      }
      if(nofhit==2 && x!=-1 && y!=-1){
         h2mppccor->Fill(x,y);
         if(x>63 && y<(64-5)) {
            connectchinfo += "X: " + to_string(x) + " " + to_string(nlx) + " " + to_string(nfx) + " "+ ",  Y: " + to_string(y) + " " + to_string(nly) + " " + to_string(nfy);
            //osf << connectchinfo << endl;
            //cout << connectchinfo << endl;
         }
      }
   }
   //
   
   osf.close();

   TCanvas *c2mppccor = new TCanvas("c2mppccor","c2mppccor");
   h2mppccor->SetStats(0);
   h2mppccor->Draw("colz");

   if(saveflag){
      string picname;
      if(isinner) picname = "./pic/chrelation/2mppc/" + to_string(runnumber) + "_inner" + to_string(mppc0) + "_" + to_string(mppc1) + ".png";
      else picname = "./pic/chrelation/2mppc/" + to_string(runnumber) + "_outer" + to_string(mppc0) + "_" + to_string(mppc1) + ".png";
      c2mppccor->Print(picname.c_str());

      delete c2mppccor;
      delete h2mppccor;
   }


   return;
}