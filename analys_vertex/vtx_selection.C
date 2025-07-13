#include "TMathBase.h"
#include "TMath.h"
#include "roottuple_var.h"


#include "TObject.h"
#include "TVector3.h"
#include <vector>
#include "TLorentzVector.h"

Double_t const M_Pion = 139.57039; // MeV

Int_t vtx_selection(){ // returns selected vertex id (if only one is found)

    Double_t const zvmax  = 4.0; // cm
    Double_t const rtvmax = 2.0; // cm

    Int_t nvip = 0;
    Int_t kvip[3];
    for (Int_t kv = 0; kv < nv; kv++){ //loop on vertices
        
        Double_t rtv = TMath::Sqrt((xv[kv]-bx)*(xv[kv]-bx)+(yv[kv]-by)*(yv[kv]-by));

        if (rtv<rtvmax && TMath::Abs(zv[kv]-bz) < zvmax){
            if (nvip<3) kvip[nvip] = kv; // save vertex index 
            nvip++;
        }
    }
    if (nvip != 1) return -666;
    return kvip[nvip];
}


TVector2 trkv_selection(Int_t vtx_id){ // use with the vertex index

    Int_t ntv_vtxid = 0;
    
    std::vector<TVector3> trkv_momenta;
    std::vector<Int_t> trkv_charge;
    std::vector<Int_t> trkv_index;


    for (Int_t ktv = 0; ktv < ntv ; ktv++){ //loop over tracks connected to vertices

        if (iv[ktv] == vtx_id){ // these tracks are the one connected to my vertex
            ntv_vtxid++;        // just to check that two tracks only have been found
            if (curv[ktv]>0) {
                trkv_charge.push_back(1);
            }
	    else {
                trkv_charge.push_back(-1);
            }
            TVector3 momentum(pxtv[ktv],pytv[ktv],pztv[ktv]);
            trkv_momenta.push_back(momentum);
            trkv_index.push_back(ktv);
        }

    }

    // selection of vertex with only two track connected and opposite charge
    if (ntv_vtxid == 2 && trkv_charge[0]*trkv_charge[1] < 0 ) { 
        // some further selection on momenta? 
    
        TVector2 trkv_sel;
        if (trkv_charge[0]>0) {
            trkv_sel.SetX(trkv_index[0]); // pi+ track index in NTV block 
            trkv_sel.SetY(trkv_index[1]); // pi- track index in NTV block
        } else {
            trkv_sel.SetX(trkv_index[1]); //pi- in the NTV track block 
            trkv_sel.SetY(trkv_index[0]); //pi+ in the NTV track block  
        }
        return trkv_sel;
        //TVector3 pionp = (trkv_charge[0]>0 ? trkv_momenta[0] : trkv_momenta[1]); // pi+ 
        //TVector3 pionm = (trkv_charge[0]<0 ? trkv_momenta[0] : trkv_momenta[1]); // pi-
        //TLorentzVector pip;
        //pip.SetVect(pionp);
        //pip.SetVectM(pionp,M_Pion)
    }   

}





