
#include "MbdQA.h"
//MBD
#include <mbd/MbdPmtContainer.h>
#include <mbd/MbdGeom.h>
#include <mbd/MbdPmtHit.h>
#include <phool/getClass.h>
#include <mbd/MbdPmtContainerV1.h>
#include <mbd/MbdOut.h>
#include <mbd/MbdDefs.h>

//Vertx
#include <calobase/TowerInfoContainer.h>
#include <calobase/TowerInfo.h>
#include <globalvertex/GlobalVertexMap.h>
#include <globalvertex/GlobalVertex.h>
#include <globalvertex/MbdVertexMap.h>
#include <globalvertex/MbdVertex.h>


#include <phool/PHCompositeNode.h>

//Fun4All
#include <fun4all/Fun4AllReturnCodes.h>
#include <fun4all/Fun4AllServer.h>
#include <fun4all/Fun4AllHistoManager.h>
#include <phool/PHCompositeNode.h>
#include <phool/getClass.h>
#include <phool/phool.h>
#include <phool/recoConsts.h>
#include <ffaobjects/EventHeader.h>
#include <ffarawobjects/Gl1Packet.h>
#include <ffaobjects/EventHeaderv1.h>

//ROOT stuff
#include <TH1F.h>
#include <TH2F.h>
#include <TH3F.h>
#include <TFile.h>
#include <TLorentzVector.h>
#include <TTree.h>
#include <TGraphErrors.h>
#include <TCanvas.h>
#include <TString.h>
#include <TSystem.h>

#include <iostream>
#include <cmath>
#include <cstdio>

// Centrality
#include <centrality/CentralityInfo.h>

// Minimum Bias
#include <calotrigger/MinimumBiasInfo.h>

using namespace std;
using namespace MbdDefs;


//---------------------------------------------->
MbdQA::MbdQA(const std::string &name):
 SubsysReco(name)
    ,T(nullptr)
    ,T2(nullptr)
    ,_tres( 0.05 ) //time in each mbd, per event
    ,_savefile(name)
{
  cout << "MbdQA::MbdQA(const std::string &name) :)  Calling constuctor" << std::endl;
}

//---------------------------------------------------------->
MbdQA::~MbdQA()
{
  std::cout << "MbdQA::~MbdQA() :) Calling deconstructors" << std::endl;
}

//--------------------------------------------------------------->
int MbdQA::Init(PHCompositeNode *topNode)
{
  cout << PHWHERE << " Saving to file " << Outfile << endl;
  _savefile = new TFile(_savefname.c_str(),"RECREATE");
  
  T = new TTree("T","MbdQA");
  T->Branch("run",&f_run);
  T->Branch("ch",&f_ch);
  T->Branch("f_bimp",&f_bimp);
  T->Branch("f_ncoll",&f_ncoll);
  T->Branch("f_npart",&f_npart);
  T->Branch("qmean",&f_qmean);
  T->Branch("qmerr",&f_qmerr);
  
  T->Branch("f_vx",&f_vx);
  T->Branch("f_vy",&f_vy);
  T->Branch("f_vz",&f_vz);
  T ->Branch("f_bz",&f_bz);
  T-> Branch("f_vt",&f_vt);
  T->Branch("f_vr",&f_vr);
  
  T->Branch("bq",&bq);
  T->Branch("bqs",&f_mbdq[0]);
  T->Branch("bqn",&f_mbdq[1]);
  T->Branch("bns",&f_mbdn[0]);
  T->Branch("bnn",&f_mbdn[1]);
  
  
  T->Branch("mbdt0",&f_mbdt0);
  T->Branch("bts",&f_mbdt[0]);   
  T->Branch("btn",&f_mbdt[1]);
  T->Branch("bts_EARhit",&f_mbdt_EARhit[0]);   
  T->Branch("btn_EARhit",&f_mbdt_EARhit[1]);
  
  
 //the output from PMTs
  T->Branch("bq",&bq); //get_q
  T->Branch("bt",&bt); //get_time
  T->Branch("btt",&btt);// get_tt 
  T->Branch("btq",&btq);//get_tq
  
  
  T->Branch("phi",&phi);
  //T->Branch("cent",&cent);
 // T-> Branch("evtPlaneAngle",&evtPlaneAngle);
 
  
  T2 = new TTree("T2","MbdQA Event Tree");
  T2->Branch("run",&f_run,"run/I");
  T2->Branch("evt",&f_evt,"evt/I");
  T2->Branch("cross",&f_cross,"cross/S");
  T2->Branch("rtrig",&f_rtrig,"rtrig/l");
  T2->Branch("ltrig",&f_ltrig,"ltrig/l");
  T2->Branch("strig",&f_strig,"strig/l");
  T2->Branch("bz",&f_bz,"bz/F");


  TString name, title;
  //Loop over all PMTs
  for (int ipmt=0; ipmt<128; ipmt++)
  {
    name = "h_q"; name += ipmt;
    title = "mbd charge, ch "; title += ipmt;
    h_mbdq[ipmt] = new TH1F(name,title,1200,0,60);

    // TGraph to track mean of mbdq distribution
    name = "g_mbdq"; name += ipmt;
    title = "mbdq, ch "; title += ipmt;
    g_mbdq[ipmt] = new TGraphErrors();
    g_mbdq[ipmt]->SetName(name);
    g_mbdq[ipmt]->SetTitle(name);

  }
  //Loop over each south and north arm
  for (int iarm=0; iarm<2; iarm++)
  {
    name = "h_mbdqtot"; name += iarm;
    title = "mbd charge, arm "; title += iarm;
    h_mbdqtot[iarm] = new TH1F(name,title,1400,0,1400); 

    name = "h2_mbdqtot_bz"; name += iarm;
    title = "mbd charge, arm "; title += iarm; title += " vs z";
    h2_mbdqtot_bz[iarm] = new TH2F(name,title,20,-50.,50.,1000,0.,3000.);

    
    name = "hevt_mbdt"; name += iarm;
    title = "mbd times, arm "; title += iarm;
    hevt_mbdt[iarm] = new TH1F(name,title,200,7.5,11.5);
    hevt_mbdt[iarm]->SetLineColor(4);
  }
  hstack_bqs&bqn = new THStack("hstack_bqs&bqn", "The total charge on the south and North MBD; bqn & bqs [A.U.]");
  
 
  h_mbdqsum = new TH1F("h_mbdqsum","BBC/MBD north + south charge sum",3000,0.,3000.);
  h2_mbdqsum = new TH2F("h2_mbdqsum","north MBDQ vs South MBDQ",1400,0.,1400.,1400,0.,1400.);

  h2_mbdqsum_bz = new TH2F("h2_mbdqsum_bz","MBDQsum vs z",20,-50.,50.,1000,0.,3000.);
  
  h2_mbdn = new TH2F("h2_mbdn"," number of event that have hits on MBDSouth & MBDNorth; South MBD Nhits; North MBD Nhits ", 70,0,70,70,0,70);

  h_zdce = new TH1F("h_zdce","ZDC Energy",820,-50,4050);
  h_zdcse = new TH1F("h_zdcse","ZDC.S Energy",500,0,250);
  h_zdcne = new TH1F("h_zdcne","ZDC.N Energy",500,0,250);
  h_zdctimecut = new TH1F("h_zdctimecut", "zdctimecut", 50, -17.5 , 32.5);

  h_emcale = new TH1F("h_emcale","Emcal Energy",3000,-100,2900);
  h_emcaltimecut = new TH1F("h_emcaltimecut", "emcaltimecut", 50, -17.5 , 32.5);

  h_ohcale = new TH1F("h_ohcale","OHCAL Energy",1000,-100,900);
  h_ohcaltimecut = new TH1F("h_ohcaltimecut", "ohcaltimecut", 50, -17.5 , 32.5);

  h_ihcale = new TH1F("h_ihcale","IHCAL Energy",1000,-100,900);
  h_ihcaltimecut = new TH1F("h_ihcaltimecut", "ihcaltimecut", 50, -17.5 , 32.5);

  h_bz = new TH1F("h_bz","MBD z-vertex",1200,-300,300);
  h_bz->SetXTitle("z_{VTX} [cm]");

  for (int itrig=0; itrig<5; itrig++)
  {
    name = "h_bz"; name += itrig;
    title = "MBD z-vertex, trig "; title += itrig;
    h_bztrig[itrig] = new TH1F(name,title,1200,-300,300);
    h_bztrig[itrig]->SetXTitle("z_{VTX} [cm]");
  }

  h_bpmt_bad = new TH1F("h_bpmt_bad","PMT for BAD MBD z-vertex",128,0,128);

  
  for (int ipmt=0; ipmt<128; ipmt++)
  {
    name = "h2_slew"; name += ipmt;
    title = "slew curve, ch "; title += ipmt;
    h2_slew[ipmt] = new TH2F(name,title,4000,0.,100,1100,-5,6);
  }
  h2_tq = new TH2F("h2_tq","ch vs tq",900,-150,150,128,-0.5,128-0.5);
  h2_tt = new TH2F("h2_tt","ch vs tt",900,-150,150,128,-0.5,128-0.5);

  gaussian = new TF1("gaussian","gaus",0,20);
  gaussian->FixParameter(2,0.05);   // set sigma to 50 ps

  c_mbdt = new TCanvas("c_mbdt","MBD Times",1200,800);
  c_mbdt->Divide(1,2);

  // MBD triggers
  mbdtrigbits.push_back(0x0400);    // MBDNS>=1
  mbdtrigbits.push_back(0x0800);    // MBDNS>=2
  mbdtrigbits.push_back(0x1000);    // MBD60
  mbdtrigbits.push_back(0x2000);    // MBD30
  mbdtrigbits.push_back(0x4000);    // MBD10


  
  return Fun4AllReturnCodes::EVENT_OK;
}

//-------------------------------------------------------------------------->
int MbdQA::InitRun(PHCompositeNode *topNode)
{  
  //cout << "MbdQA::InitRun(PHCompositeNode *topNode) Initializing for Run XXX" << endl;
  recoConsts *rc = recoConsts::instance();
  f_run = rc->get_IntFlag("RUNNUMBER");

  GetNodes(topNode);

  return Fun4AllReturnCodes::EVENT_OK;
}

//------------------------------------------------------------------------------>

//Call user instructions for every event
int MbdQA::process_event(PHCompositeNode *topNode)
{
  // Get the raw gl1 data from event combined DST
 Gl1Packet* _gl1raw = findNode::getClass<Gl1Packet>(topNode, "GL1Packet");

  if(!_gl1raw && f_evt<4) cout << PHWHERE << " Gl1Packet node not found on node tree" << endl;

  nprocessed++;

  f_evt = _evtheader->get_EvtSequence();
  if (f_evt%1000==0) cout << PHWHERE << "Event " << f_evt << "\t" << ", nprocessed = " << nprocessed << endl;

  // Only use MBDNS triggered events
  if ( _gl1raw != nullptr )
  {
    const uint64_t MBDTRIGS = 0x7c00;  // MBDNS trigger bits
    const uint64_t ZDCNS = 0x8;        // ZDCNS trigger bits

    f_cross = _gl1raw->getBunchNumber();
    f_rtrig = _gl1raw->getTriggerVector();
    f_ltrig = _gl1raw->getLiveVector();
    f_strig = _gl1raw->getScaledVector();

    if ( (f_strig&MBDTRIGS) == 0 )
    {
      cout << PHWHERE << " scaled vector and MBDTRIGS = 0 " << endl;
      return Fun4AllReturnCodes::ABORTEVENT;
    }

    if ( nprocessed<10 )
    {
      cout << "Trig " << nprocessed << endl;
      cout << hex;
      cout << "Raw\t" << f_rtrig << endl;
      cout << "Live\t" << f_ltrig << endl;
      cout << "Scaled\t" << f_strig <<endl;
      cout << "AND\t" << (f_strig&MBDTRIGS) << endl;
      cout << dec;

    }

   }

  // Initialize Variables
  f_mbdn[0] = 0;
  f_mbdn[1] = 0;
  f_mbdq[0] = 0.;
  f_mbdq[1] = 0.;
  f_mbdt[0] = -9999.;
  f_mbdt[1] = -9999.;
  f_mbdte[0] = -9999.;
  f_mbdte[1] = -9999.;
  f_bz = NAN;
  f_mbdt0 = NAN;
  hevt_mbdt[0]->Reset();
  hevt_mbdt[1]->Reset();

  
  

  CheckDST(topNode);

  T2->Fill();
  return 0;
}

//-------------------------------------------------------------->

void MbdQA::GetNodes(PHCompositeNode *topNode)
{


  //MbdPmt information
  _mbdpmts = findNode::getClass<MbdPmtContainer>(topNode,"MbdPmtContainer"); 
    if(!_mbdpmts && f_evt<4) 
      {
        cout << PHWHERE << "Aborting : MbdPmtContainer node not found on node tree" << endl;
 
      return Fun4AllReturnCodes::ABORTEVENT;
     }

  //MbdOut 
  _mbdout = findNode::getClass<MbdOut>(topNode,"MbdOut"); 
  if(!_mbdout && f_evt<4) 
    { 
      cout << PHWHERE << " Aborting : MbdOut node not found on node tree" << endl;
 
      return Fun4AllReturnCodes::ABORTEVENT;
    }

  //Get DST object
  _evtheader = findNode::getClass<EventHeader>(topNode,"EventHeader"); 
  if (!_evtheader && f_evt<10) 
    {
      cout << PHWHERE << " Aborting : EventHeader node not found on node tree" << endl;
 
      return Fun4AllReturnCodes::ABORTEVENT;
    }
     // Get Truth Centrality Info
     f_bimp = _evtheader->get_ImpactParameter();
     f_ncoll = _evtheader->get_ncoll();
     f_npart = _evtheader->get_npart();

  _centInfo = findNode::getClass<CentralityInfo>(topNode,"CentralityInfo");
  if(!centInfo)
    {
      cout << PHWHERE << "Aborting: CentralityInfo node not found on node tree" << std::endl;

      return Fun4AllReturnCodes::ABORTEVENT;
    }
  cent = centInfo->get_centile(CentralityInfo::PROP::mbd_NS);

/*
  MinimumBiasInfo *minBiasInfo = findNode::getClass<MinimumBiasInfo>(topNode,"MinimumBiasInfo");
  Bool_t isMinBias = (minBiasInfo) ? minBiasInfo->isAuAuMinimumBias() : false;
  if(!minBiasInfo || !isMinBias)
  {
    std::cout << PHWHERE << "caloTreeGen::process_event: is not MinimumBias" << std::endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }
  */
  
  //MBD pmg geometry 
  _mbdgeom = findNode::getClass<MbdGeom>(topNode, "MbdGeom");
  if(!mbdgeom)
    {
      cout << "Aborting: MbdGeom node not found on node tree" << endl;
      return Fun4AllReturnCodes::ABORTEVENT;
    }

   f_npmt = _mbdpmts -> get_npmt();        
  cout << "nPMTS =  "<< f_npmt << endl;
  
  f_bz = _mbdout->get_zvtx();
  
  evtPlaneAngle = evtHeader -> get_EventPlaneAngle();
  cout <<" evtPlaneAngle = " <<endl;
  
  f_bimp = evtHeader->get_ImpactParameter();
  if(f_bimp >= 0.00  && f_bimp < 4.88) cent = 0.05;
  if(f_bimp >= 4.88  && f_bimp < 6.81) cent = 0.15;
  if(f_bimp >= 6.81  && f_bimp < 8.40) cent = 0.25;
  if(f_bimp >= 8.40  && f_bimp < 9.71) cent = 0.35;
  if(f_bimp >= 9.71  && f_bimp < 10.81) cent = 0.45;
  if(f_bimp >= 10.81 && f_bimp < 11.84) cent = 0.55;
  if(f_bimp >= 11.84 && f_bimp < 12.76) cent = 0.65;
  if(f_bimp >= 12.76 && f_bimp < 13.68) cent = 0.75;
  if(f_bimp >= 13.68 && f_bimp < 14.65) cent = 0.85;

  for(int i = 0; i < nPMTs; i++)
    {
      MbdPmtHit* mbdpmt = _mbdpmts -> get_pmt(i);        // grab ith pmt (changed from pmtID)
      
      //store quantities of interest
      bq.push_back(mbdpmt -> get_q());
      bt.push_back(mbdpmt -> get_time());
      f_vx.push_back(mbdgeom->get_x(i));
      f_vy.push_back(mbdgeom->get_y(i));
      f_vz.push_back(mbdgeom->get_z(i));
      f_vr.push_back(mbdgeom->get_r(i));
      phi.push_back(mbdgeom->get_phi(i));

    }
  
  T -> Fill();
  
   
 // return Fun4AllReturnCodes::EVENT_OK;
}

//------------------------------------------------------------>
int MbdQA::ResetEvent(PHCompositeNode *topNode)
{
  bq.clear();
  bt.clear();
  f_vx.clear();
  f_vy.clear();
  f_vz.clear();
  f_vr.clear();
  phi.clear();
  return Fun4AllReturnCodes::EVENT_OK;
}

//------------------------------------unused------------------------------------>
/*int MbdQA::EndRun(const int runnumber)
{
  std::cout << "MbdQA::EndRun(const int runnumber) Ending Run for Run " << runnumber << std::endl;
  return Fun4AllReturnCodes::EVENT_OK;
}
*/
//____________________________________________________________________________..
int MbdQA::End(PHCompositeNode *topNode)
{
  cout << "MbdQA::End(PHCompositeNode *topNode) Filing the histograms" << endl;
  _savefile->cd();

   Double_t nevents = h_mbdqsum->Integral();
  h_mbdqsum->Fill(-1000,nevents); // underflow bin keeps track of nevents
  h_mbdqtot[0]->Fill(-1000,nevents); // underflow bin keeps track of nevents
  h_mbdqtot[1]->Fill(-1000,nevents); // underflow bin keeps track of nevents
 
  
  
  Double_t norm = 1.0/nevents;
  h_mbdqsum->Scale( norm );
  h2_mbdqsum->Scale( norm );

  h_mbdqtot[0]->Scale( norm );
  h_mbdqtot[1]->Scale( norm );
  

  for (int ipmt=0; ipmt<MBD_N_PMT; ipmt++)
  {
    h_mbdq[ipmt]->Scale( norm );

    // Fill info on q distribution
    f_ch = ipmt;
    f_qmean = h_mbdq[ipmt]->GetMean();
    f_qmerr = h_mbdq[ipmt]->GetMeanError();
    T->Fill();
    //cout << f_run << "\t" << f_ch << "\t" << f_qmean << "\t" << f_qmerr << endl;

    g_mbdq[ipmt]->SetPoint(0,f_run,f_qmean);
    g_mbdq[ipmt]->SetPointError(0,0,f_qmerr);
    g_mbdq[ipmt]->Write();
  }

  cout << "Nevents processed integral " << nprocessed << "\t" << nevents << "\t" << nevents/nprocessed << endl;
  _savefile->Write();
  _savefile->Close();
  return Fun4AllReturnCodes::EVENT_OK;
}


void MbdQA::CheckDST(PHCompositeNode *topNode)
{
   cout << "MbdQA::CheckDST(PHCompositeNode *topNode) do not panic you are in CheckDST" << endl;
  
  // bqs = _mbdout->get_q(0)/120.;
  // bqn = _mbdout->get_q(1)/120.;
  bqs = _mbdout->get_q(0);
  bqn = _mbdout->get_q(1);
  
  
  bts = _mbdout->get_time(0);
  btn = _mbdout->get_time(1);
  
  f_bz = _mbdout->get_zvtx();
  h_bz->Fill( f_bz );
  
  bns = _mbdout->get_npmt(0);
  bnn = _mbdout->get_npmt(1);
  cout<< " bns = "<< bns <<"\t" << "bnn = " << bnn << endl;
   
  
  T->Fill();
  // Check the MbdVertexMap
  MbdVertexMap *mbdvtxmap = findNode::getClass<MbdVertexMap>(topNode, "MbdVertexMap");
  if ( mbdvtxmap && !mbdvtxmap->empty() )
  {
    MbdVertex *vtx = mbdvtxmap->begin()->second;
    if ( vtx )
    {
      float vtx_z = vtx->get_z();
      if ( !isnan(f_bz) && vtx_z != f_bz )
      {
        cout << "ERROR, vertices do not match " << vtx_z << "\t" << f_bz << endl;
      }
      else
      {
        static int counter = 0;
        if ( counter<3 )
        {
          cout << "GOOD, mbd vertices match " << vtx_z << "\t" << f_bz << endl;
          counter++;
        }
      }
    }
  }
  else
  {
    static int counter = 0;
    if ( counter < 4 )
    {
      if ( mbdvtxmap && mbdvtxmap->empty() )
      {
        cout << "MbdVertexMap is empty" << endl;
      }
      else
      {
        cout << "MbdVertexMap not found" << endl;
      }
      counter++;
    }
  }

  // Check the GlobalVertex
  GlobalVertexMap *vertexmap = findNode::getClass<GlobalVertexMap>(topNode, "GlobalVertexMap");
  if ( vertexmap && !vertexmap->empty() )
  {
    GlobalVertex *vtx = vertexmap->begin()->second;
    if ( vtx )
    {
      float vtx_z = vtx->get_z();
      if ( vtx_z != f_bz )
      {
        cout << "ERROR, vertices do not match " << vtx_z << "\t" << f_bz << endl;
      }
    }
  }
  else
  {
    static int counter = 0;
    if ( counter < 4 )
    {
      if ( vertexmap && vertexmap->empty() )
      {
        cout << "GlobalVertexMap is empty" << endl;
      }
      else
      {
        cout << "GlobalVertexMap not found" << endl;
      }
      counter++;
    }
  }

  //cout << "bz " << f_bz << endl;
  if ( fabs(f_bz)>3000. ) return;
  //if ( bqn<10 && bqs>150 ) return;

  //Float_t r = (4.4049/4.05882);
  Float_t r = 1.0;

  h_mbdqsum->Fill( (bqn+bqs)*r );
  h_mbdqtot[0]->Fill( bqs*r );
  h_mbdqtot[1]->Fill( bqn*r );
  hstack_bqs&bqn->Add( h_mbdqtot[0]);
  hstack_bqs&bqn->Add( h_mbdqtot[1]);
  
  h2_mbdqsum->Fill( bqn*r, bqs*r );

  h2_mbdqtot_bz[0]->Fill( f_bz, bqs );
  h2_mbdqtot_bz[1]->Fill( f_bz, bqn );
  h2_mbdqsum_bz->Fill( f_bz, bqn+bqs );
  
  h2_mbdn->Fill( bns, bnn );

  cout << "nPMTs " << _mbdpmts->get_npmt() << endl;
  
  // Fill info from each PMT
  //for (int ipmt=0; ipmt<_mbdpmts->get_npmt(); ipmt++)
  for (int ipmt=0; ipmt<128; ipmt++)
  {
    int arm = ipmt/64;
    //Float_t q = _mbdpmts->get_pmt(ipmt)->get_q()/120.;
     bq = _mbdpmts->get_pmt(ipmt)->get_q();
     bt = _mbdpmts->get_pmt(ipmt)->get_time();
     btt = _mbdpmts->get_pmt(ipmt)->get_tt();
     btq = _mbdpmts->get_pmt(ipmt)->get_tq();
     phi = mbdgeom->get_phi(ipmt);   // get phi angle

    if ( fabs(t) < 25. )
    {
      h_mbdq[ipmt]->Fill( q*r );
    }

    h2_tt->Fill( t, ipmt );
    h2_tq->Fill( tq, ipmt );

    h2_slew[ipmt]->Fill( q, t - mbdt[arm] );
    //cout << ipmt << ":\t" << q << "\t" << t << endl;
 
  }

  h_bz->Fill( f_bz );
  for (size_t itrig=0; itrig<mbdtrigbits.size(); itrig++)
  {
    if ( f_ltrig&mbdtrigbits[itrig] )
    {
      h_bztrig[itrig]->Fill( f_bz );
    }
  }

  // Analyze other subsystems
  //process_gl1( topNode );
  //process_zdc( topNode );
  /*
  process_emcal( topNode );
  process_ohcal( topNode );
  process_ihcal( topNode );
  */
}

int MbdQA::Getpeaktime(TH1 * h)
{
  int getmaxtime, tcut = -1;

  for(int bin = 1; bin < h->GetNbinsX()+1; bin++)
  {
    double c = h->GetBinContent(bin);
    double max = h->GetMaximum();
    int bincenter = h->GetBinCenter(bin);
    if(max == c)
    {
      getmaxtime = bincenter;
      if(getmaxtime != -1) tcut = getmaxtime;
    }
  }

  return tcut;
}

void MbdQA::process_gl1( PHCompositeNode *topNode )
{
  TowerInfoContainer* zdctowers = findNode::getClass<TowerInfoContainer>(topNode, "TOWERINFO_CALIB_ZDC");
  if (zdctowers)
  {
    int max_zdc_t = Getpeaktime( h_zdctimecut );
    int range = 1;
    double zdc_etot = 0.;
    double zdc_e[2] {0.,0.};

    int size = zdctowers->size(); //online towers should be the same!
    for (int ich = 0; ich < size; ich++)
    {
      TowerInfo* zdctower = zdctowers->get_tower_at_channel(ich);
      float zdce = zdctower->get_energy();
      int zdct = zdctowers->get_tower_at_channel(ich)->get_time();
      h_zdctimecut->Fill( zdct );

      if ( (zdct  < (max_zdc_t - range)) ||  (zdct  > (max_zdc_t + range)) )
      {
        continue;
      }

      //
      if (ich==0||ich==2||ich==4)
      {
        zdc_e[0] += zdce;
      }
      else if (ich == 8 || ich == 10 || ich == 12)
      {
        zdc_e[1] += zdce;
      }

    }

    h_zdcse->Fill( zdc_e[0] );
    h_zdcne->Fill( zdc_e[1] );
    zdc_etot = zdc_e[0] + zdc_e[1];
    h_zdce->Fill( zdc_etot );
  }
}


void MbdQA::process_zdc( PHCompositeNode *topNode )
{
  TowerInfoContainer* zdctowers = findNode::getClass<TowerInfoContainer>(topNode, "TOWERINFO_CALIB_ZDC");
  if (zdctowers)
  {
    int max_zdc_t = Getpeaktime( h_zdctimecut );
    int range = 1;
    double zdc_etot = 0.;

    int size = zdctowers->size(); //online towers should be the same!
    for (int ich = 0; ich < size; ich++)
    {
      TowerInfo* zdctower = zdctowers->get_tower_at_channel(ich);
      float zdce = zdctower->get_energy();
      int zdct = zdctowers->get_tower_at_channel(ich)->get_time();
      h_zdctimecut->Fill( zdct );

      if (ich == 0 || ich == 2 || ich == 4 || ich == 8 || ich == 10 || ich == 12)
      {
        if( zdct  > (max_zdc_t - range) &&  zdct  < (max_zdc_t + range))
        {
          zdc_etot += zdce;
        }
      }
    }

    h_zdce->Fill( zdc_etot);
  }
}


void MbdQA::process_emcal( PHCompositeNode *topNode )
{
  TowerInfoContainer* towers = findNode::getClass<TowerInfoContainer>(topNode, "TOWERINFO_CALIB_CEMC");
  if (towers)
  {
    int max_emcal_t = Getpeaktime( h_emcaltimecut );
    int range = 1;
    double etot = 0.;

    int size = towers->size(); //online towers should be the same!
    for (int channel = 0; channel < size;channel++)
    {
      TowerInfo* tower = towers->get_tower_at_channel(channel);
      float energy = tower->get_energy();
      int t = towers->get_tower_at_channel(channel)->get_time();
      h_emcaltimecut->Fill(t);

      if( abs(t - max_emcal_t) <  range )
      {
        etot += energy;
      }
    }

    h_emcale->Fill(etot);
  }
}

void MbdQA::process_ohcal( PHCompositeNode *topNode )
{
  TowerInfoContainer* towers = findNode::getClass<TowerInfoContainer>(topNode, "TOWERINFO_CALIB_HCALOUT");
  if (towers)
  {
    int max_hcal_t = Getpeaktime( h_ohcaltimecut );
    int range = 1;
    double etot = 0.;

    int size = towers->size(); //online towers should be the same!
    for (int channel = 0; channel < size;channel++)
    {
      TowerInfo* tower = towers->get_tower_at_channel(channel);
      float energy = tower->get_energy();
      int t = towers->get_tower_at_channel(channel)->get_time();
      h_ohcaltimecut->Fill(t);

      if( abs(t - max_hcal_t) <  range )
      {
        etot += energy;
      }
    }

    h_ohcale->Fill(etot);
  }
}

void MbdQA::process_ihcal( PHCompositeNode *topNode )
{
  TowerInfoContainer* towers = findNode::getClass<TowerInfoContainer>(topNode, "TOWERINFO_CALIB_HCALIN");
  if (towers)
  {
    int max_hcal_t = Getpeaktime( h_ihcaltimecut );
    int range = 1;
    double etot = 0.;

    int size = towers->size(); //online towers should be the same!
    for (int channel = 0; channel < size;channel++)
    {
      TowerInfo* tower = towers->get_tower_at_channel(channel);
      float energy = tower->get_energy();
      int t = towers->get_tower_at_channel(channel)->get_time();
      h_ihcaltimecut->Fill(t);

      if( abs(t - max_hcal_t) <  range )
      {
        etot += energy;
      }
    }

    h_ihcale->Fill(etot);
  }
}

//------------------------------------------------------->
int MbdQA::Reset(PHCompositeNode *topNode)
{
  
  pmtcharge.clear();
  pmttime.clear();
  f_vx.clear();
  f_vy.clear();
  f_vz.clear();
  pmtr.clear();
  pmtphi.clear();

  
  return Fun4AllReturnCodes::EVENT_OK;
}

//----------------------------------------------------------------------->
void MbdQA::Print(const std::string &what) const
{
  std::cout << "MbdQA::Print(const std::string &what) const Printing info for " << what << std::endl;
}

