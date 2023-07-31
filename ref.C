void ref(){
  TString filename="./datas/ExVFT_00784.root";
  TFile *f=new TFile(filename);
  TTree *tree=(TTree*)f->Get("tree");
  TTreeReader fReader;  
  fReader.SetTree(tree);
  TTreeReaderArray<short> vft_layer   = {fReader, "vft_layer"};
  TTreeReaderArray<short> vft_channel = {fReader, "vft_channel"};
  TTreeReaderArray<short> vft_leading = {fReader, "vft_leading"};
  TTreeReaderArray<short> vft_trailing= {fReader, "vft_trailing"};
  int nev=fReader.GetEntries();
  //  nev=1000;
  TH1* hprof=new TH1I("hprof","hit profile",64,-0.5,63);
  TH1* htdc=new TH1I("htdc","tdc",1e3,0,1e3);
  for(int i=0;i<nev;i++){
    fReader.SetLocalEntry(i);
    for(int ihit=0;ihit<vft_layer.GetSize();ihit++){
      if(vft_layer.At(ihit)==0){
	      hprof->Fill(vft_channel.At(ihit));
	      if(vft_channel.At(ihit)==0){
	        htdc->Fill(vft_leading.At(ihit));
	      }
      }
    }
  }
  TCanvas *c1=new TCanvas();
  c1->Divide(2,2);
  c1->cd(1); hprof->Draw();
  c1->cd(2); htdc->Draw();
  c1->cd(3);
  tree->Draw("vft_leading>>h(1000,0,1000)","vft_layer==0&&vft_channel==0");
}

