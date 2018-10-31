#include<vector>
#include<math.h>
#include<TRandom3.h>
#include<TMath.h>
#include<TROOT.h>
#include<TCanvas.h>
#include<TTree.h>
#include<TBranch.h>
#include<TF1.h>
#include<TF2.h>
#include<TGraph.h>
#include<TFile.h>
#include<iostream>
#include<fstream>
#include<TH1F.h>
#include<TH1D.h>
#include<TH2D.h>
#include<TH3D.h>
#include<TStyle.h>
#include<TMinuit.h>
#include"pmts_profileid18.hh"
#include"pmts_profileid19.hh"
#include"pmts_profileid20.hh"
#include"pmts_profileid21.hh"
#include"pmts_profileid22.hh"
#include"pmts_profileid23.hh"
#include"pmts_profileid24.hh"
#include"pmts_profileid25.hh"
#include"run_number_profile_id_array.hh"
#include<TBits.h>
#include <math.h>
#include<string>
#include <algorithm>
#include<iomanip>

//cambiamento per git!! 

double pmtPos[2241][6];

void init_pmtPos(int profileN){

    for(int i=0; i<2241; i++){
        for(int j=0; j<6; j++){
            pmtPos[i][j] = 0.0;
        }
    }
    for(int i=0; i<2052; i++){
        for(int j=0; j<6; j++){
     if(profileN==18) pmtPos[ int(pmt_profile_18[i][0]) ][j] = pmt_profile_18[i][j+1];
     if(profileN==19) pmtPos[ int(pmt_profile_19[i][0]) ][j] = pmt_profile_19[i][j+1];
     if(profileN==20) pmtPos[ int(pmt_profile_20[i][0]) ][j] = pmt_profile_20[i][j+1];
     if(profileN==21) pmtPos[ int(pmt_profile_21[i][0]) ][j] = pmt_profile_21[i][j+1];
     if(profileN==22) pmtPos[ int(pmt_profile_22[i][0]) ][j] = pmt_profile_22[i][j+1];
     if(profileN==23) pmtPos[ int(pmt_profile_23[i][0]) ][j] = pmt_profile_23[i][j+1];
     if(profileN==24) pmtPos[ int(pmt_profile_24[i][0]) ][j] = pmt_profile_24[i][j+1];
     if(profileN==25) pmtPos[ int(pmt_profile_25[i][0]) ][j] = pmt_profile_25[i][j+1];
      }
    }
}
double get_distance(int lg, double x, double y, double z){
    double result = sqrt( (pmtPos[lg][0]-x)*(pmtPos[lg][0]-x)+(pmtPos[lg][1]-y)*(pmtPos[lg][1]-y)+(pmtPos[lg][2]-z)*(pmtPos[lg][2]-z) );
    return result;
}
double get_angle_direction(int lg, double x, double y, double z, double xdir, double ydir, double zdir){
    double pmt_dist = get_distance(lg, x, y, z);
    double pmt_event_x = (pmtPos[lg][0]-x)/pmt_dist;
    double pmt_event_y = (pmtPos[lg][1]-y)/pmt_dist;
    double pmt_event_z = (pmtPos[lg][2]-z)/pmt_dist;
    double angle = (xdir*pmt_event_x + ydir*pmt_event_y + zdir*pmt_event_z);
    return angle;
}


void make_time_histo(float neff = 1.52, int n_aBins = 10, int n_hits_per_time = 10000){

cout<<"LEGGO I DATI VERI-----------------------------------"<<endl;

        TH1D *time_tof = new TH1D("time_ToF","time_ToF", 500,0.0, 500.0);
        TH1D *time_tof_noTread = new TH1D("time_ToF_noTread","time_ToF_noTread", 500,0.0, 500.0);
        TH1D *time_Tread = new TH1D("time_Tread","time_Tread", 500,0.0, 500.0);
        TH1D *time_riord = new TH1D("time_riord","time_riord", 500,0.0, 500.0);


	vector<double> time_vec;
	vector<int> run_vec;
	vector<float> phi_vec;
	vector<float> theta_vec;
       TFile* f = new TFile("mioTree_nusol_data_14Sep_reducedToMc.root", "READ"); //load Data file to read run and ev_num because the MC-files are named after them
 
   TTree* tree = (TTree*)f->Get("mioTree");
    int n_events = tree->GetEntries();
    cout << "n_events " << n_events << endl;
    float X_read;
    float Y_read;
    float Z_read;
    float T_read;
    float Xdir;
    float Ydir;
    float Zdir;
    float azimuth_read;
    float altitude_read;
    int nhits_decoded_read;
    int nhits_clustered_read;
    int lg_read[1000];
    float TimeRaw_read[1000];
    int hit_clustered_read[1000];
    int run_read;
    int evnum_read;
    float n_norm_pmts_read;
    float n_live_pmts_geo_read;

    tree->SetBranchAddress("Pos_Reco_x",&X_read);
    tree->SetBranchAddress("Pos_Reco_y",&Y_read);
    tree->SetBranchAddress("Pos_Reco_z",&Z_read);
    tree->SetBranchAddress("Pos_Reco_t",&T_read);
    tree->SetBranchAddress("nhits_decoded",&nhits_decoded_read);
    tree->SetBranchAddress("nhits_clustered",&nhits_clustered_read);
    tree->SetBranchAddress("TimeRaw",TimeRaw_read);
    tree->SetBranchAddress("hit_clustered",hit_clustered_read);
    tree->SetBranchAddress("lg",lg_read);
    tree->SetBranchAddress("azimuth",&azimuth_read);
    tree->SetBranchAddress("altitude",&altitude_read);
    tree->SetBranchAddress("n_norm_pmts",&n_norm_pmts_read);
    tree->SetBranchAddress("n_live_pmts_geo",&n_live_pmts_geo_read);
    tree->SetBranchAddress("run",&run_read);
    tree->SetBranchAddress("evnum",&evnum_read);
	int pf_last=-1;
   
 
   
    for(int iev=0; iev<n_events; iev++) {
        if(iev%(n_events/100)==0) cout<<iev/(n_events/100)<<endl;
        tree->GetEntry( iev );
		if(run_read < 17329 || run_read > 26568) continue;
        if(n_norm_pmts_read < 260. || n_norm_pmts_read > 320.) continue;
        int pf=0;
        for(int j=2; j<11; j++){
            if(run_read <= RunNumber_ProfileID[j][0]){
                pf=RunNumber_ProfileID[j][1];
                break;
            }
        }
        if(pf!=pf_last) {
            init_pmtPos(pf);
            cout<<"Cambiato il profilo! " << pf <<endl;
        }
        pf_last=pf;

	double phi = TMath::Pi()*(308.0 - azimuth_read)/180.0;
        double theta = acos(cos(TMath::Pi()*( 90.0 - altitude_read  )/180.0));
        Xdir = -sin(theta)*cos(phi);
        Ydir = -sin(theta)*sin(phi);
        Zdir = -cos(theta);

        run_vec.push_back(run_read);
	phi_vec.push_back(TMath::ATan2(Ydir,Xdir));
	theta_vec.push_back(acos(Zdir));

        vector<double> time_vec;
        int ccc=0;
	double time0=TimeRaw_read[0];
        for(int d=0; d<nhits_decoded_read; d++){
            if(lg_read[d] <= 0 || hit_clustered_read[d] != 1 ) continue;
            double pmt_dist = get_distance(lg_read[d],X_read,Y_read,Z_read);
            double angle = get_angle_direction(lg_read[d],X_read,Y_read,Z_read , Xdir, Ydir, Zdir);
            double time_corrected = TimeRaw_read[d] - time0 + T_read - neff*(pmt_dist)/0.3;
            double time_noTread = TimeRaw_read[d] - time0  - neff*(pmt_dist)/0.3+62;
         time_vec.push_back(TimeRaw_read[d]-neff*(pmt_dist)/0.3);
         time_tof->Fill(time_corrected);
         time_tof_noTread->Fill(time_noTread);
         time_Tread->Fill(T_read);
         ccc++;
         }
         
         std::sort(time_vec.begin(),time_vec.end());
         
        if(time_vec.size()!=ccc) cout<<"ERROR "<<endl; 
        
         for(int i=0; i<time_vec.size(); i++){
            time_riord->Fill(time_vec[i]-time_vec[0]+32);
         }
    }
  
   f->Close();


//LEGGO I DATI MC

cout<<"LEGGO I DATI MC be7---------------------------------"<<endl;

        TH1D *MCtime_tof = new TH1D("MCtime_ToF","MCtime_ToF", 400,0.0, 400.0);
        TH1D *MCtime_tof_noTread = new TH1D("MCtime_ToF_noTread","MCtime_ToF_noTread", 400,0.0, 400.0);
        TH1D *MCtime_Tread = new TH1D("MCtime_Tread","MCtime_Tread", 400,0.0, 400.0);
        TH1D *MCt_time_tof = new TH1D("MCt_time_ToF","MCt_time_ToF", 400,0.0, 400.0);
        TH1D *MCt_time_tof_noTread = new TH1D("MCt_time_ToF_noTread","MCt_time_ToF_noTread", 400,0.0, 400.0);
        TH1D *MC_time_riord = new TH1D("MC_time_riord","MC_time_riord", 400,0.0, 400.0);

    TFile* fmc = new TFile("miotree_nusol_be7MC_14Sep.root", "READ"); 
    
    TTree* treemc = (TTree*)fmc->Get("mioTree");
    n_events = treemc->GetEntries();
    cout << "n_events " << n_events << endl;
    /*
    float X_read;
    float Y_read;
    float Z_read;
    float T_read;
    float Xdir;
    float Ydir;
    float Zdir;
    float azimuth_read;
    float altitude_read;
    int nhits_decoded_read;
    int nhits_clustered_read;
    int lg_read[1000];
    float TimeRaw_read[1000];
    int hit_clustered_read[1000];
    int run_read;
    int evnum_read;
    float n_norm_pmts_read;
    float n_live_pmts_geo_read;
*/

    float mc_x;
    float mc_y;
    float mc_z;
    float mc_e;

    treemc->SetBranchAddress("Pos_Reco_x",&X_read);
    treemc->SetBranchAddress("Pos_Reco_y",&Y_read);
    treemc->SetBranchAddress("Pos_Reco_z",&Z_read);
    treemc->SetBranchAddress("Pos_Reco_t",&T_read);
    treemc->SetBranchAddress("nhits_decoded",&nhits_decoded_read);
    treemc->SetBranchAddress("nhits_clustered",&nhits_clustered_read);
    treemc->SetBranchAddress("TimeRaw",TimeRaw_read);
    treemc->SetBranchAddress("hit_clustered",hit_clustered_read);
    treemc->SetBranchAddress("lg",lg_read);
    treemc->SetBranchAddress("n_norm_pmts",&n_norm_pmts_read);
    treemc->SetBranchAddress("n_live_pmts_geo",&n_live_pmts_geo_read);
    treemc->SetBranchAddress("run",&run_read);
    treemc->SetBranchAddress("evnum",&evnum_read);
    treemc->SetBranchAddress("mc_dir_x",&Xdir);
    treemc->SetBranchAddress("mc_dir_y",&Ydir);
    treemc->SetBranchAddress("mc_dir_z",&Zdir);
    treemc->SetBranchAddress("mc_pos_x",&mc_x);
    treemc->SetBranchAddress("mc_pos_y",&mc_y);
    treemc->SetBranchAddress("mc_pos_z",&mc_z);
    treemc->SetBranchAddress("mc_energy",&mc_e);
	
	pf_last=-1;
   
    
    for(int iev=0; iev<n_events; iev++) {
        if(iev%(n_events/100)==0) cout<<iev/(n_events/100)<<endl;
        treemc->GetEntry( iev );
		if(run_read < 17329 || run_read > 26568) continue;
        if(n_norm_pmts_read < 260. || n_norm_pmts_read > 320.) continue;
        int pf=0;
        for(int j=2; j<11; j++){
            if(run_read <= RunNumber_ProfileID[j][0]){
                pf=RunNumber_ProfileID[j][1];
                break;
            }
        }
        if(pf!=pf_last) {
            init_pmtPos(pf);
            cout<<"Cambiato il profilo! " << pf <<endl;
        }
        pf_last=pf;

        run_vec.push_back(run_read);
	phi_vec.push_back(TMath::ATan2(Ydir,Xdir));
	theta_vec.push_back(acos(Zdir));
        vector<double> time_v;

	double time0=TimeRaw_read[0];
        for(int d=0; d<nhits_decoded_read; d++){
            if(lg_read[d] <= 0 || hit_clustered_read[d] != 1 ) continue;
            double pmt_dist = get_distance(lg_read[d],X_read,Y_read,Z_read);
            double angle = get_angle_direction(lg_read[d],X_read,Y_read,Z_read , Xdir, Ydir, Zdir);
            double time_corrected = TimeRaw_read[d] - time0 + T_read - (neff-0.02)*(pmt_dist)/0.3;
            double time_noTread = TimeRaw_read[d] - time0  - (neff-0.02)*(pmt_dist)/0.3+62;

         MCtime_tof->Fill(time_corrected);
         MCtime_tof_noTread->Fill(time_noTread);
         MCtime_Tread->Fill(T_read);
         
             //tempi calcolati con MCtruth
         pmt_dist = get_distance(lg_read[d],mc_x,mc_y,mc_z);
         angle = get_angle_direction(lg_read[d],mc_x,mc_y,mc_z, Xdir, Ydir, Zdir);
         time_corrected = TimeRaw_read[d] - time0 + T_read - (neff-0.02)*(pmt_dist)/0.3;
         time_noTread = TimeRaw_read[d] - time0  - (neff-0.02)*(pmt_dist)/0.3+62;

         MCt_time_tof->Fill(time_corrected);
         MCt_time_tof_noTread->Fill(time_noTread);
         time_v.push_back(TimeRaw_read[d]-neff*(pmt_dist)/0.3);
  
        }
     
        std::sort(time_v.begin(),time_v.end());
   //     if(time_vec.size()!=ccc) cout<<"ERROR "<<endl;

         for(int i=0; i<time_v.size(); i++){
           MC_time_riord->Fill(time_v[i]-time_v[0]+32);
         }

   }
 
   fmc->Close();


//PLOT

//normalizzo
time_tof_noTread->Scale(1/time_tof_noTread->Integral());
time_tof->Scale(1/time_tof->Integral());
time_Tread->Scale(1/time_Tread->Integral());
time_riord->Scale(1/time_riord->Integral());

MCtime_tof_noTread->Scale(1/MCtime_tof_noTread->Integral());
MCtime_tof->Scale(1/MCtime_tof->Integral());
MCtime_Tread->Scale(1/MCtime_Tread->Integral());

MCt_time_tof_noTread->Scale(1/MCt_time_tof_noTread->Integral());
MCt_time_tof->Scale(1/MCt_time_tof->Integral());
MC_time_riord->Scale(1/MC_time_riord->Integral());

/*
   TH1D *sub_time_tof = new TH1D("sub_time_tof","sub_time_tof", 150,-30.0, 120.0);
   TH1D *sub_time_tof_noTread = new TH1D("sub_time_tof_noTread","sub_time_tof_noTread", 150,-30.0, 120.0);
   TH1D *sub_time_Tread = new TH1D("sub_time_Tread","sub_time_Tread", 150,-30.0, 120.0);

    TH1D *subMCt_time_tof=new TH1D("subMCt_time_tof","subMCt_time_tof",150, -30.0,120.0);

for(int i=1; i<=150; i++){
  float val;
  if(time_tof->GetBinContent(i)!=0) val=(time_tof->GetBinContent(i)-MCtime_tof->GetBinContent(i))/time_tof->GetBinContent(i);
  else val=0;
sub_time_tof->SetBinContent(i,val);
//cout<<i<<"contenuto time_tof "<<time_tof->GetBinContent(i)<<" "<<MCtime_tof->GetBinContent(i)<<" "<<val<<endl;
if(time_tof_noTread->GetBinContent(i)!=0) val=(time_tof_noTread->GetBinContent(i)-MCtime_tof_noTread->GetBinContent(i))/time_tof_noTread->GetBinContent(i);
else val=0;
sub_time_tof_noTread->SetBinContent(i,val);

if(time_tof->GetBinContent(i)!=0) val=(time_tof->GetBinContent(i)-MCt_time_tof->GetBinContent(i))/time_tof->GetBinContent(i);
else val=0;
subMCt_time_tof->SetBinContent(i,val);

if(time_Tread->GetBinContent(i)!=0) val=(time_Tread->GetBinContent(i)-MCtime_Tread->GetBinContent(i))/time_Tread->GetBinContent(i);
else val=0;
sub_time_Tread->SetBinContent(i,val);


}
*/
 
   TCanvas *c0 =  new TCanvas();
   time_Tread->SetLineColor(2);
   MCtime_Tread->SetLineColor(1);
   time_Tread->Draw();
   MCtime_Tread->Draw("SAME");

   time_tof_noTread->SetLineColor(3);
   time_tof_noTread->Draw("SAME");
   MCtime_tof_noTread->SetLineColor(1);
   MCtime_tof_noTread->Draw("SAME");
  // MCt_time_tof_noTread->SetLineColor(3);
  // MCt_time_tof_noTread->Draw("SAME");
   
   time_tof->SetLineColor(4);
   time_tof->Draw("SAME");
   MCtime_tof->SetLineColor(1);
   MCt_time_tof->SetLineColor(1);
 
   time_riord->SetLineColor(6);
   time_riord->Draw("SAME");

   MC_time_riord->SetLineColor(1);
   MC_time_riord->Draw("SAME");
/*
   TCanvas *c1 =  new TCanvas();

   sub_time_tof->SetLineColor(2);
   sub_time_tof_noTread->SetLineColor(3);
   sub_time_Tread->SetLineColor(1);

   sub_time_tof->Draw();
   sub_time_tof_noTread->Draw("SAME");
   sub_time_Tread->Draw("SAME");
*/

   TFile out("isto_time_nusol_simili.root","RECREATE");
  
   MCtime_Tread->Write();
   time_tof_noTread->Write();
   MCtime_tof_noTread->Write();
   MCt_time_tof_noTread->Write();

   MCt_time_tof->Write();
  MCtime_tof->Write();
   time_tof->Write();
  MC_time_riord->Write();
   time_riord->Write();
/*
   sub_time_tof->Write();
   sub_time_tof_noTread->Write();
   sub_time_Tread->Write();
   subMCt_time_tof->Write();
*/
  out.Close();

}



