
// test_analyzer.cpp
// Charles E. Hyde, 05 November 2021
//
// Analyze the output of a Delphes simulation of DVCS event file in CORE
//  Tuned for hepmc/delphes files of Andrey Kim
//  beam particles in hempc input file and delphes output:
//      Particles bank with Particles->Status==4

#include <stdio.h>

#include <TTree.h>
#include <TFile.h>
#include <TSystem.h>
#include <TCanvas.h>
#include <TH1.h>
#include <TH2D.h>
#include <TChain.h>
#include <TF1.h>
#include <TLorentzVector.h>
#include <TRandom3.h>
#include <TStyle.h>
#include <TLatex.h>
#include <TLegend.h>
#include <TDatabasePDG.h>

#ifdef __CLING__
R__LOAD_LIBRARY(libDelphes)
#include "classes/DelphesClasses.h"
#include "external/ExRootAnalysis/ExRootProgressBar.h"
#include "external/ExRootAnalysis/ExRootTreeBranch.h"
#include "external/ExRootAnalysis/ExRootTreeReader.h"
#endif

// Access physical constants, masses in pdg.lbl.gov database implemented in root
auto dbPDG = TDatabasePDG::Instance();
const double umass = 0.931494; // GeV, atomic mass unit u
const double alpha_mass = 4.0*umass+2.4249e-3;

const double Crossing_angle = 0.035;  // switch future input files to positive


int test_analyzer(const char *fname="output_CORE"){
    
    dbPDG->TDatabasePDG::AddParticle("4He","alpha", alpha_mass, kTRUE,0.0,6.0,"nuclei", 1000020040,-1,0);
    TVector3 k3Beam;
    TVector3 P3Beam;
    
    char inputfile[128];
    sprintf(inputfile,"%s.root",fname);
    printf("%s\n",inputfile);
    gSystem->Load("libDelphes");
    TRandom3 ran3;
    TLatex tt;

    // Create chain of root trees
    TChain chain("Delphes");
    chain.Add(inputfile);

    
    // Create object of class ExRootTreeReader
    ExRootTreeReader *treeReader = new ExRootTreeReader(&chain);
    Long64_t numberOfEntries = treeReader->GetEntries();
    
    // Get pointers to branches used in this analysis
    TClonesArray *branchJet = treeReader->UseBranch("Jet");
    TClonesArray *branchElectron = treeReader->UseBranch("Electron");
    TClonesArray *branchParticle = treeReader->UseBranch("Particle");
    TClonesArray *branchPhoton = treeReader->UseBranch("Photon");
    TClonesArray *branchEFlow = treeReader->UseBranch("EFlowTrack");
    TClonesArray *branchTrack = treeReader->UseBranch("Track");

    GenParticle * thisParticle;
    Photon * thisPhoton;
    Electron * thisElectron;
    Track * thisEFlow;
    Track * thisTrack;
    // Book histograms
//    TH1 *histJetPT = new TH1F("jet_pt", "jet P_{T}", 100, 0.0, 100.0);
    TH1 *histMass = new TH1F("mass", "M_{inv}(e_{1}, e_{2})", 100, 0.0, 200.0);
    TH1 *h_MissingMassSq = new TH1F("h_MissingMassSq","M_{X}^{2} (GeV2)",100.,-5.0,5.0);
    TH1 *h_photons = new TH1F("h_photons","Photon Multiplicity",10,-0.5,9.5);
    TH1 *h_electrons = new TH1F("h_electrons","Electron Multiplicity",10,-0.5,9.5);
    TH1 *h_tracks = new TH1F("h_tracks","Reconstructed Track Multiplicity",10,-0.5,9.5);
    TH1 *h_qprime = new TH1F("h_qprime","(Q2-Q2_{Gen})/Q2_{Gen}",250,-1.0,1.0);
    TH1 *h_PT = new TH1F("h_PT","Photon PT ; GeV ; Events", 200,0.0,20.);
    TH1 *h_E = new TH1F("h_E","Photon Energy ; GeV ; Events", 200,0.0,20.);
    TH1 *h_error_g = new TH1F("h_error_g","Photon Reconstruction ; #Delta E (GeV) ; Events", 200,-0.5,0.5);
    TH1 *h_eta_g = new TH1F("h_eta_g","Photon Pseudo Rapidity ; #eta ; Events", 200,-4.0,4.0);
    TH2 *h2_E_vs_eta_e = new TH2F("h2_E_vs_eta_e","Detected Electrons; #eta ; E_{e} (GeV)",
                                  80.,-4.0,4.0,40.,0.0,20.0);
    TH1 *h_error_e = new TH1F("h_error_e","Electron Reconstruction; |#Delta k3|/k3  ; Events",
                              200.,-0.5,0.5);
    TH1 *h_DeltaPerp = new TH1F("h_DeltaPerp",
                                "#alpha(e,e'#gamma)#alpha; -#Delta_{#perp}^{2} GeV^{2} ; Events",
                                50,0.0,0.25);
    TH1 *h_DeltaPerpGen = new TH1F("h_DeltaPerpGen",
                                "#alpha(e,e'#gamma)#alpha; -#Delta_{#perp}^{2} GeV^{2} ; Events",
                                50,0.0,0.25);
    TH1 *h_norm = new TH1F("h_norm","LightCone Mass",200,-4.0,4.0);
    TH1 *h_DeltaPerp_err = new TH1F("h_DeltaPerp_err","-(#Delta_{#perp} -#Delta_{#perp Gen}^{2}); GeV^{2}    ",
                                    500,-0.1,0.4);
    TH1 *h_DeltaPerp2_err = new TH1F("h_DeltaPerp2_err",
                                     " #Delta_{#perp}^{2} -#Delta_{#perp Gen}^{2}; GeV^{2}    ",
                                    600,-0.15,0.15);
    TH1 *h_Delta2_err = new TH1F("h_Delta2_err",
                                     " #Delta^{2} -#Delta_{Gen}^{2}; GeV^{2}    ",
                                    600,-0.15,0.15);

    const int nbin_xA = 20;
    double xA_fact = pow(10.,-1./5.);
    double xA_BinLowEdge[nbin_xA+1];
    double xxA=1.0;
    printf("xBj low edge" );
    for (int ibin=0; ibin<=nbin_xA; ibin++) {
        printf("%8.3g, ", xxA );
        xA_BinLowEdge[nbin_xA-ibin] = xxA;
        xxA *= xA_fact;
    }
    printf(" \n");
    const double Q2Min = 1.0, Q2Max=200.0;
    const int nbin_Q2 = 10;
    double Q2_BinLowEdge[nbin_Q2+1];
    double Q2_fact = pow(Q2Max/Q2Min,1.0/(double)nbin_Q2);
    double Q2val = Q2Min;
    for (int ibin=0; ibin<=nbin_Q2+1; ibin++) {
        Q2_BinLowEdge[ibin] = Q2val;
        Q2val*=Q2_fact;
    }

    TH2 *h2_Q2_xBj = new TH2F("h2_Q2_xBj","Q^{2} vs x_{Bj};x_{Bj}    ;Q^{2}      ",
                              nbin_xA,xA_BinLowEdge,nbin_Q2,Q2_BinLowEdge);
    //h2_Q2_xBj->SetBins(nbin_xA,xA_BinLowEdge,nbin_Q2,Q2_BinLowEdge);

    TLorentzVector k4Beam, k4Scat, P4Beam, P4Scat, q4prime, Delta, P4Miss, k4ScatCheck;
    TLorentzVector k4BeamGen, k4ScatGen, k4EFlow, P4BeamGen, P4ScatGen, q4primeGen, P4MissGen;
    TLorentzVector Delta4Gen, Zero, q4Virtual, n4_q, n4Tilde_q, Delta_perp, Delta_perpGen;
    TVector3 q3prime,q3primeGen, q3Diff, k3Scat, k3ScatGen, Delta3Scat;
    double DeltaSqGen, DeltaSq, Q2, Q2Gen, kScat_E, yGen, relErr, radiation;
    double aEM=0.01, bEM=0.02, esmear,gsmear;
    double q_dot_P, deltaQ, sqrt_one_d, norm;
    
    //k4BeamGen.SetVectM(k3Beam,dbPDG->GetParticle(ElPID)->Mass());
    //P4BeamGen.SetVectM(P3Beam,dbPDG->GetParticle(IonPID)->Mass());
    // Loop over all events
    Int_t jentry=0;
    Int_t nElectron=0;
    Int_t nJet=0;
    Int_t nPhoton=0;
    // Delphes status codes, incident beam or stable produced particle
    Int_t iBeamStatus = 4;
    Int_t iStableStatus = 1;
    // PDG codes
    const int PDGelectron = 11;
    const int PDGphoton = 22;
    const int PDGproton = 2212;
    const int PDGbaryons= 2000;
    for(jentry = 0; jentry < numberOfEntries; ++jentry)
    {
        // Load selected branches with data from specified event
        treeReader->ReadEntry(jentry);
        double Egamma = 0.0;
        int igamma = -1;
        for(int ipart=0; ipart<branchParticle->GetEntries(); ipart++){
            thisParticle = (GenParticle *) branchParticle->At(ipart);
            if (thisParticle->PID==PDGelectron) {
                if (thisParticle->Status==iBeamStatus){
                    k4BeamGen = thisParticle->P4();
                } else {
                k4ScatGen = thisParticle->P4();
                k3ScatGen = k4ScatGen.Vect();
                }
            } else if( (thisParticle->PID>PDGbaryons)&&(thisParticle->Status==iBeamStatus) ){
                P4BeamGen =thisParticle->P4();
            } else if (thisParticle->PID==PDGphoton) {
                q4primeGen = thisParticle->P4();
                if (q4primeGen.E()>Egamma){
                    igamma = ipart;
                    q3primeGen = q4primeGen.Vect();
                }
            }
        }
        nElectron += branchElectron->GetEntries();
        h_electrons->Fill((float)branchElectron->GetEntries());
        nJet += branchJet->GetEntries();
        nPhoton += branchPhoton->GetEntries();
        h_photons->Fill((float)branchPhoton->GetEntries());
        h_tracks->Fill((float)branchTrack->GetEntries());
        for (Int_t iphot=0; iphot<branchPhoton->GetEntries(); iphot++) {
            thisPhoton = (Photon *) branchPhoton->At(iphot);
            q4prime = thisPhoton->P4();
            h_PT->Fill(thisPhoton->PT);
            h_E->Fill(q4prime.E());
            h_eta_g->Fill(q4prime.Eta());
            if (q4prime.E()>Egamma){
                Egamma = q4prime.E();
                igamma = iphot;
            }
            /*
            thisParticle = (GenParticle *)thisPhoton->Particle.GetObject();
            if (thisParticle->PID==PDGphoton) {
                q4primeGen = thisParticle->P4();
                h_error_g->Fill(q4prime.E()- q4primeGen.E());
            }
             */
            for(Int_t jpart=0; jpart<branchParticle->GetEntries(); jpart++){
                thisParticle = (GenParticle *) branchParticle->At(jpart);
                if (thisParticle->PID==PDGphoton){
                    q4primeGen = thisParticle->P4();
                    h_error_g->Fill((q4prime.E()- q4primeGen.E())/q4primeGen.E());
                }
            }
        }
        for (Int_t iElec=0; iElec<branchElectron->GetEntries(); iElec++) {
            thisElectron = (Electron *) branchElectron->At(iElec);
            k4Scat = thisElectron->P4();
            h2_E_vs_eta_e->Fill(k4Scat.Eta(),k4Scat.E());
            k3Scat = k4Scat.Vect();
            Delta3Scat = k3Scat-k3ScatGen;
 //           kScat_err = sqrt(Delta3Scat.Dot(Delta3Scat))/k3ScatGen.P();
            h_error_e->Fill(sqrt(Delta3Scat.Dot(Delta3Scat))/k3ScatGen.Mag());
        }
        if (igamma>-1 && branchElectron->GetEntries()==1) {
            thisPhoton = (Photon *) branchPhoton->At(igamma);
            q4prime = thisPhoton->P4();
            // Find Lorentz Invariant Delta_perp
            //
            //  Make the Event Light cone vectors n4_q and n4Tilde_q
            //
            q4Virtual = k4BeamGen-k4ScatGen;
            q_dot_P   = q4Virtual.Dot(P4BeamGen);
            deltaQ    = -q4Virtual.M2()*P4BeamGen.M2()/(q_dot_P*q_dot_P);
            sqrt_one_d= sqrt(1.0+deltaQ);
            norm      = 1.0/(sqrt(q_dot_P)*sqrt_one_d);
            n4_q  = (1.0+sqrt_one_d)* q4Virtual;
            n4_q += ( -q4Virtual.M2()/q_dot_P ) * P4BeamGen;
            n4_q *= norm/2.0;
            n4Tilde_q  = (-P4BeamGen.M2()/((1.0+sqrt_one_d)*q_dot_P)) *q4Virtual ;
            n4Tilde_q += P4BeamGen;
            n4Tilde_q *= norm;
            if (abs(n4Tilde_q.Dot(n4_q)-1.0)>0.001 || abs(n4Tilde_q.M2())>0.001 || abs(n4_q.M2())>0.001) {
                n4_q.Print();
                n4Tilde_q.Print();
                return jentry;
            }
            if (jentry%10000==0){
                printf("%d nqSq, nqTildeSq, nq(nqTilde) = %10.3g, %10.3g, %10.3g \n",
                       jentry, n4_q.M2(), n4Tilde_q.M2(), n4Tilde_q.Dot(n4_q));
            }
            Delta4Gen = q4Virtual - q4primeGen;
            Delta_perpGen = Delta4Gen-(Delta4Gen.Dot(n4_q))*n4Tilde_q - (Delta4Gen.Dot(n4Tilde_q))*n4_q;
            h_DeltaPerpGen->Fill(-Delta_perpGen.M2());
            //h_norm->Fill(n4Tilde_q.Dot(n4_q));
            // Now computer reconstructed Delta
            thisPhoton = (Photon *) branchPhoton->At(igamma);
            q4prime = thisPhoton->P4();
            // no beam smearing yet
            k4Beam = k4BeamGen;
            Delta = k4Beam-k4Scat-q4prime;
            Delta_perp = Delta-(Delta.Dot(n4_q))*n4Tilde_q - (Delta.Dot(n4Tilde_q))*n4_q;
            h_DeltaPerp->Fill(-Delta_perp.M2());
            // Square of the difference or difference of the squares?
            // Either is valid
            //h_DeltaPerp2_err->Fill(2.0*Delta_perp.Dot(Delta_perpGen)  -Delta_perp.M2()-Delta_perpGen.M2());
            h_DeltaPerp2_err->Fill(Delta_perp.M2()-Delta_perpGen.M2());
            h_Delta2_err->Fill(Delta.M2()-Delta4Gen.M2());
        }
    }
    printf("Electrons, Photons Jets found in branches = (%d, %d, %d) \n",
           nElectron,nPhoton,nJet);
   
    TFile *froot= new TFile("histOut.root","recreate");

    TCanvas *c1 = new TCanvas("c1","Delphes",75,50,900,600);
    c1->SetLogy();
    c1->Divide(3,1);
    c1->cd(1);
    h_electrons->Draw();
    gPad->SetLogy();
//    gPad->SetLogy();  gPad->SetLogx();
    c1->cd(2);
    h_photons->Draw();
    gPad->SetLogy();
    c1->cd(3);
    h_tracks->Draw();
    TCanvas *c2 = new TCanvas("c2","Photons",100,75,900,600);
    c2->Divide(3,1);
    c2->cd(1);
    h_E->Draw();
    gPad->SetLogy();
    c2->cd(2);
    h_error_g->Draw();
    gPad->SetLogy();
    gStyle->SetOptStat(1111);
    TF1 *fGauss = new TF1("fGauss","gaus(x)",-0.07,0.07);
    fGauss->SetParameters(jentry/100.,0.0,0.012);
    h_error_g->Fit("fGauss","I","",-0.07,0.07);
    gStyle->SetOptFit(1111);
    c2->cd(3);
    h_eta_g->Draw();
    gPad->SetLogy();
    TCanvas *c3 = new TCanvas("c3","Electons",125,100,900,600);
    c3->Divide(2,1);
    c3->cd(1);
    h2_E_vs_eta_e->Draw("cont1z");
    gPad->SetLogz();
    gPad->SetRightMargin(0.15);
    c3->cd(2);
    h_error_e->Draw();
    gPad->SetLogy();
    
    TCanvas *c4 = new TCanvas("c4","#Delta_{#perp}",105,125,900,600);
    c4->Divide(3,1);
    c4->cd(1);
    gPad->SetLeftMargin(0.15);
    gPad->SetRightMargin(0.05);
    TLegend *leg4 = new TLegend(0.125,0.15,0.4,0.3,"CORE");
    gStyle->SetTitleSize(0.05);
    gStyle->SetLabelSize(0.05);
    h_DeltaPerpGen->SetLineWidth(2);
    h_DeltaPerpGen->SetLineColor(kRed);
    h_DeltaPerpGen->SetTitleSize(0.10);
    h_DeltaPerpGen->Draw();
    h_DeltaPerpGen->GetXaxis()->SetTitleSize(0.05);
    h_DeltaPerpGen->GetXaxis()->SetLabelSize(0.05);
    h_DeltaPerpGen->GetYaxis()->SetTitleSize(0.05);
    h_DeltaPerpGen->GetYaxis()->SetLabelSize(0.05);
    leg4->AddEntry(h_DeltaPerpGen,"TOPEG","L");
    gStyle->SetOptStat(11);
    gPad->SetLogy();
    h_DeltaPerp->SetLineWidth(2);
    h_DeltaPerp->SetLineColor(kBlack);
    h_DeltaPerp->SetMarkerStyle(21);
    h_DeltaPerp->SetMarkerColor(kBlack);
    h_DeltaPerp->SetMarkerSize(0.5);
    h_DeltaPerp->Draw("sameEP");
    leg4->AddEntry(h_DeltaPerp,"DELPHES","EP");
    /*
    tt.SetTextAlign(12);
    tt.DrawLatex(0.025,10.,"Resolution #oplus Bin-Migration");
     */
    leg4->Draw();
    c4->cd(2);
    h_DeltaPerp2_err->Draw();
    c4->cd(3);
    h_Delta2_err->Draw();
    
    h_E->Write();
    h_error_g->Write();

    h2_E_vs_eta_e->Write();
    h_error_e->Write();
    
    h_DeltaPerpGen->Write();
    h_DeltaPerp->Write();
    h_DeltaPerp2_err->Write();

    froot->Close();
    return jentry;
}
