#ifndef __MBDQA_H__
#define __MBDQA_H__

#include <fun4all/SubsysReco.h>
#include <string>
#include <vector>
#include <map>
#include <TFile.h>
//#include <THStack.h>

//Forward declerations
class PHCompositeNode;
class EventHeader;
class MbdOut;
class MbdPmtContainer;
class MbdGeom;
class CentralityInfo;
class MbdGeom;

class Gl1Packet;
class TFile;
class TTree;
class TDatabasePDG;
class TRandom3;
class TH1;
class TH2;
class TGraphErrors;
class TF1;
class TCanvas;
class THStack;

//Brief: basic TTree and histogram creation for sim evaluation
class MbdQA: public SubsysReco
{
  public: 

    //Default constructor
    MbdQA(const std::string &name="MbdQA");
     
    // destructor
     ~MbdQA() override;

    //Initialization, called for initialization
    int Init(PHCompositeNode *) override;

    //Initialization at start of every run
    int InitRun(PHCompositeNode *)  override;

    //Process Event, called for each event
    int process_event(PHCompositeNode *) override;

    ///Clean up internals after each event.
    int ResetEvent(PHCompositeNode *topNode) override;

    /// Called at the end of each run.
    int EndRun(const int runnumber)  override;

    //End, write and close files
    int End(PHCompositeNode *) override;

    // Reset
    int Reset(PHCompositeNode * /*topNode*/) override;

    void Print(const std::string &what = "ALL") const override;

    //Change output filename
    void set_savefile(const char *f) { _savefname = f; }

    void set_run(const int r) { f_run = r; }

  private:

    //
    void CheckDST(PHCompositeNode *topNode);

    int  Getpeaktime(TH1 *h);
    void process_gl1( PHCompositeNode *topNode );
    void process_zdc( PHCompositeNode *topNode );
    void process_emcal( PHCompositeNode *topNode );
    void process_ohcal( PHCompositeNode *topNode );
    void process_ihcal( PHCompositeNode *topNode );

    //output filename
    std::string _savefname; // "MbdQA.root"
    TFile* _savefile;

    int nprocessed{0};     // num events processed
    Int_t f_npmt = 128;
    int nPMTs=128;

    Float_t cent;
    Double_t m_cent;

    //Output
    TTree* T{ nullptr };
    TTree* T2{ nullptr };  // for event by event info
    Int_t    f_run;
    Int_t    f_evt = 0 ;
    Short_t  f_cross;
    Int_t    f_ch;
    Float_t  f_qmean;
    Float_t  f_qmerr;

    ULong64_t f_rtrig{0};
    ULong64_t f_ltrig{0};
    ULong64_t f_strig{0};

    Float_t  f_bimp;    // impact parameter
    Int_t    f_ncoll;   // number n-n collisions
    Int_t    f_npart;   // number participants

    // MbdGeom sim
    Float_t  f_vx;      // true x vertex of collision point
    Float_t  f_vy;      // true y vertex
    Float_t  f_vz;      // true z vertex
    Float_t  f_vt;       // true start time
    Float_t  f_vr;

    // MbdGeom real
     std::vector<uint64_t> pmt_x[128]{};
     std::vector<uint64_t> pmt_y[128]{};
     std::vector<uint64_t> pmt_z[128]{};
     std::vector<uint64_t> pmt_r[128]{};
     std::vector<uint64_t> pmt_phi[128]{};


    Short_t  f_mbdn[2];  // num hits for each arm (north and south)
    float  f_bq;       // total charge (sum)
    Float_t  f_mbdq[2];  // total charge (currently npe) in each arm
    Float_t  f_mbdt[2];  // time in arm
    // Float_t  f_mbdte[2]; // earliest hit time in arm
    Float_t  f_bz;       // z-vertex
    Float_t  f_mbdt0;    // start time
    Float_t bqs;         //total charge on south arm
    Float_t bqn ;       //total charge on north arm
    Float_t bts;         //time on south arm
    Float_t btn;        //time on north arm
    Float_t bq;         //total charge
    Float_t bt;          //time
    Float_t btt;       //time from the time channel
    Float_t btq ;      //time from the charge channel
    Float_t phi;       //  phi angle
    Float_t f_mbdt_EARhit[2]; //earliest hit time in arm
    Float_t bns;
    Float_t bnn;

   // Float_t  evtPlaneAngle;


    TH1 *h_mbdq[128];   // q in each tube
    TGraphErrors *g_mbdq[128];   // q in each tube
    TH1 *h_mbdqtot[2];  // total q in mbd arms
    TH1 *h_mbdqsum;     // total q in mbd arms
    TH2 *h2_mbdqsum;    // north q vs south q
    TH1 *h2_mbdqtot_bz[2];  // total q in mbd arms, vs bz
    TH2 *h2_mbdqsum_bz;    // north q vs south q, vs bz
    TH2 *h2_mbdn; //number of event that have hits on South & North 

    TH1 *h_emcale{nullptr};
    TH1 *h_emcaltimecut{nullptr};
    TH1 *h_ohcale{nullptr};
    TH1 *h_ohcaltimecut{nullptr};
    TH1 *h_ihcale{nullptr};
    TH1 *h_ihcaltimecut{nullptr};
    TH1 *h_zdce{nullptr};
    TH1 *h_zdcse{nullptr};
    TH1 *h_zdcne{nullptr};
    TH1 *h_zdctimecut{nullptr};

    TH1 *h_bz{nullptr};
    TH1 *h_bztrig[5]{nullptr}; // 5 MBDNS triggers
    TH2 *h2_tt{nullptr};
    TH2 *h2_tq{nullptr};
    TH2 *h2_slew[128]{nullptr};
    TH1 *h_cross{nullptr};

    std::vector<uint64_t> mbdtrigbits;

    TCanvas *c_mbdt;    // Canvas to 
    TH1 *hevt_mbdt[2]{nullptr};  // time in each mbd, per event
    //TF1 *gaussian{nullptr};
    THStack *hstack_bqs_bqn{nullptr}; //The total charge on the south and North MBD; bqn & bqs [A.U.]

    TH1 *h_bpmt_bad{nullptr};

    std::map<int,int> _pids;  // PIDs of tracks in the MBD

    //
    TRandom3*     _rndm;
    Float_t       _tres;    // time resolution of one channel

    //Get all the nodes
    void GetNodes(PHCompositeNode *);

    //Node pointers
    EventHeader* _evtheader;
    MbdOut* _mbdout;
    MbdPmtContainer* _mbdpmts;
    Gl1Packet* _gl1raw;
   // CentralityInfo* _centInfo;
    MbdGeom* _mbdgeom;

};

#endif //* __MBDQA_H__ *//

