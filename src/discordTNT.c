#include "ergm_MHproposal.h"
#include "ergm_edgelist.h"
#include "ergm_changestat.h"
#include "ergm_rlebdm.h"
#include "ergm_MHstorage.h"
#include "ergm_unsorted_edgelist.h"
#include "changestats_lasttoggle.h"
#include "tergm_model.h"

typedef struct {
  UnsrtEL *discordantEdges;
  UnsrtEL *discordantNonEdges;
  
  int in_discord;
  
} discordTNTStorage; 

MH_I_FN(Mi_discordTNT) {
  MHp->ntoggles = 1;
  
  ALLOC_STORAGE(1, discordTNTStorage, sto);
  sto->discordantEdges = UnsrtELInitialize(0, NULL, NULL, FALSE);
  sto->discordantNonEdges = UnsrtELInitialize(0, NULL, NULL, FALSE);
  // we ignore discord for this initialization (assuming a TICK will precede any proposals)
}

MH_X_FN(Mx_discordTNT) {
  GET_STORAGE(discordTNTStorage, sto);
  
  if(type == TICK) {
    // "clear" the discordant dyads
    sto->discordantEdges->nedges = 0;
    sto->discordantNonEdges->nedges = 0;
  }
}

MH_P_FN(MH_discordTNT) {
  GET_AUX_STORAGE(StoreTimeAndLasttoggle, dur_inf);
  
  GET_STORAGE(discordTNTStorage, sto);
  
  int in_discord;
  int in_network;
  
  int nedges = EDGECOUNT(nwp);
  int nddyads = kh_size(dur_inf->discord);
  
  if(nddyads == 0 || unif_rand() < 0.5) {
    // propose from network
    if(nedges == 0 || unif_rand() < 0.5) {
      // propose toggling a random dyad
      GetRandDyad(Mtail, Mhead, nwp);
      in_network = IS_OUTEDGE(Mtail[0], Mhead[0]);
      in_discord = kh_get(DyadMapInt, dur_inf->discord, THKey(dur_inf->discord, Mtail[0], Mhead[0])) != kh_none;
    } else {
      // propose toggling off an edge in network
      GetRandEdge(Mtail, Mhead, nwp);
      in_network = TRUE;
      in_discord = kh_get(DyadMapInt, dur_inf->discord, THKey(dur_inf->discord, Mtail[0], Mhead[0])) != kh_none;      
    }

    if(in_discord) {
      // need to resample to know index
      if(in_network) {
        UnsrtELGetRand(Mtail, Mhead, sto->discordantEdges);
      } else {
        UnsrtELGetRand(Mtail, Mhead, sto->discordantNonEdges);
      }
    }
  } else {
    // propose from discord
    if(unif_rand() < sto->discordantEdges->nedges/((double) nddyads)) {
      UnsrtELGetRand(Mtail, Mhead, sto->discordantEdges);
      in_network = TRUE;
    } else {
      UnsrtELGetRand(Mtail, Mhead, sto->discordantNonEdges);
      in_network = FALSE;
    }
    in_discord = TRUE;
  }
  
  // compute logratio
  
  int ndyads = DYADCOUNT(nwp);
  
  // these ignore overall factor of 1/2 which cancels out in ratio
  double forward_discord = in_discord ? 1.0/nddyads : 0;
  double backward_discord = in_discord ? 0 : 1.0/(1 + nddyads);
  
  double forward_network = in_network ? (0.5/nedges + 0.5/ndyads) : (nedges == 0 ? 1.0/ndyads : 0.5/ndyads);
  double backward_network = in_network ? (nedges == 1 ? 1.0/ndyads : 0.5/ndyads) : (0.5/(nedges + 1) + 0.5/ndyads);
  
  if(nddyads == 0) forward_network *= 2;
  if(nddyads == 1 && in_discord) backward_network *= 2;
  
  double forward = forward_discord + forward_network;
  double backward = backward_discord + backward_network;

  MHp->logratio = log(backward/forward);

  sto->in_discord = in_discord;
}

MH_U_FN(Mu_discordTNT) {
  // add or remove the dyad from the appropriate discordance edgelist
  
  GET_STORAGE(discordTNTStorage, sto);
  
  if(sto->in_discord) {
    // currently in discord; take it out
    if(edgeflag) {        
      UnsrtELDelete(tail, head, sto->discordantEdges);
    } else {
      UnsrtELDelete(tail, head, sto->discordantNonEdges);  
    }
  } else {
    // not currently in discord; add it in
    if(edgeflag) {
      UnsrtELInsert(tail, head, sto->discordantNonEdges);
    } else {
      UnsrtELInsert(tail, head, sto->discordantEdges);        
    }
  }
}

MH_F_FN(Mf_discordTNT) {
  GET_STORAGE(discordTNTStorage, sto);
  
  UnsrtELDestroy(sto->discordantNonEdges);
  UnsrtELDestroy(sto->discordantEdges);
}



/********************
   discordStratTNT
********************/

typedef struct {
  UnsrtEL **nonDiscordantELs;
  UnsrtEL **discordantELs;
  UnsrtEL **discordantNonELs;
  
  double **nodesbycode;
  
  double *pmat;
  double *nodecountsbycode;
  double *dyadcounts;
  
  double *tailtypes;
  double *headtypes;
  
  int mixingtype;
  
  int nmixtypes;
  
  int in_discord;
} discordStratTNTStorage; 


MH_I_FN(Mi_discordStratTNT) {
  // process the inputs and initialize all the edgelists in storage; set MHp->ntoggles to 1
  MHp->ntoggles = 1;
  
  int nmixtypes = MHp->inputs[0];
    
  int nattrcodes = MHp->inputs[1 + 3*nmixtypes];
  
  double *vattr = MHp->inputs + 1 + 3*nmixtypes + 1 + nattrcodes + N_NODES;
  
  ALLOC_STORAGE(1, discordStratTNTStorage, sto);
  
  sto->nonDiscordantELs = (UnsrtEL **)Calloc(nmixtypes, UnsrtEL *);
  sto->discordantELs = (UnsrtEL **)Calloc(nmixtypes, UnsrtEL *);
  sto->discordantNonELs = (UnsrtEL **)Calloc(nmixtypes, UnsrtEL *);
  
  for(int i = 0; i < nmixtypes; i++) {
    sto->nonDiscordantELs[i] = UnsrtELInitialize(0, NULL, NULL, FALSE);
    sto->discordantELs[i] = UnsrtELInitialize(0, NULL, NULL, FALSE);
    sto->discordantNonELs[i] = UnsrtELInitialize(0, NULL, NULL, FALSE);
  }
    
  double *inputindmat = MHp->inputs + 1 + 3*nmixtypes + 1 + nattrcodes + N_NODES + N_NODES + nmixtypes;  
  
  double **indmat = (double **)Calloc(nattrcodes, double *);
  indmat[0] = inputindmat;
  for(int i = 1; i < nattrcodes; i++) {
    indmat[i] = indmat[i - 1] + nattrcodes;
  }
  
  // we are treating all edges as nondiscordant, 
  // assuming a TICK will precede any proposals
  Vertex head;
  Edge e;
  for(Vertex tail = 1; tail <= N_NODES; tail++) {
    STEP_THROUGH_OUTEDGES(tail, e, head) {
      int index = indmat[(int)vattr[tail - 1]][(int)vattr[head - 1]];
      if(index >= 0) {
        UnsrtELInsert(tail, head, sto->nonDiscordantELs[index]);
      }
    }
  }
  Free(indmat);
  
  sto->nodecountsbycode = MHp->inputs + 1 + 3*nmixtypes + 1;
  
  sto->nodesbycode = (double **)Calloc(nattrcodes, double *);
  sto->nodesbycode[0] = MHp->inputs + 1 + 3*nmixtypes + 1 + nattrcodes;
  for(int i = 1; i < nattrcodes; i++) {
    sto->nodesbycode[i] = sto->nodesbycode[i - 1] + (int)sto->nodecountsbycode[i - 1];
  }
  
  sto->nmixtypes = nmixtypes;
  sto->pmat = MHp->inputs + 1 + 2*nmixtypes;
  
  sto->dyadcounts = MHp->inputs + 1 + 3*nmixtypes + 1 + nattrcodes + N_NODES + N_NODES;
  
  sto->tailtypes = MHp->inputs + 1;
  sto->headtypes = MHp->inputs + 1 + nmixtypes;
} 

MH_X_FN(Mx_discordStratTNT) {
  GET_STORAGE(discordStratTNTStorage, sto);
    
  if(type == TICK) {
    // transfer discordant edges to nondiscordant edges
    // clear discordant edges and discordant nonedges
    for(int i = 0; i < sto->nmixtypes; i++) {
      Vertex *tails = sto->discordantELs[i]->tails;
      Vertex *heads = sto->discordantELs[i]->heads;
      int nedges = sto->discordantELs[i]->nedges;
      
      for(int j = 0; j < nedges; j++) {
        UnsrtELInsert(tails[j], heads[j], sto->nonDiscordantELs[i]);
      }
      
      sto->discordantELs[i]->nedges = 0;
      sto->discordantNonELs[i]->nedges = 0;      
    }
  }
}

MH_P_FN(MH_discordStratTNT) {
  GET_AUX_STORAGE(StoreTimeAndLasttoggle, dur_inf);
  GET_STORAGE(discordStratTNTStorage, sto);
      
  double ur = unif_rand();
  
  // find the first mixing type i with (cumulative) probability larger than ur
  int i = 0;
  while(ur > sto->pmat[i]) {
    i++;
  }
  
  // record the mixing type of the toggle, in case it's needed in the U function later
  sto->mixingtype = i;    
  
  int tailtype = sto->tailtypes[i];
  int headtype = sto->headtypes[i];
    
  // number of edges of this mixing type
  int nedgestype = sto->nonDiscordantELs[i]->nedges + sto->discordantELs[i]->nedges;

  // number of dyads of this mixing type
  int ndyadstype = sto->dyadcounts[i];
  
  // number of discordant dyads of this mixing type
  int nddyadstype = sto->discordantNonELs[i]->nedges + sto->discordantELs[i]->nedges;
  
  // flags
  int in_discord;
  int in_network;
  
  if(nddyadstype == 0 || unif_rand() < 0.5) {
    // propose from network
    if(nedgestype == 0 || unif_rand() < 0.5) {
      // propose toggling a random dyad of the specified mixing type
      
      int tailindex = sto->nodecountsbycode[tailtype]*unif_rand();
      int headindex;
      if(tailtype == headtype) {
        // need to avoid sampling a loop
        headindex = (sto->nodecountsbycode[headtype] - 1)*unif_rand();
        if(headindex == tailindex) {
          headindex = sto->nodecountsbycode[headtype] - 1;
        }
      } else {
        // any old head will do
        headindex = sto->nodecountsbycode[headtype]*unif_rand();
      }
            
      Vertex tail = sto->nodesbycode[tailtype][tailindex];
      Vertex head = sto->nodesbycode[headtype][headindex];
      
      if(tail > head && !DIRECTED) {
        Vertex tmp = tail;
        tail = head;
        head = tmp;
      }
      
      Mtail[0] = tail;
      Mhead[0] = head;

      in_network = IS_OUTEDGE(Mtail[0], Mhead[0]);
      in_discord = kh_get(DyadMapInt, dur_inf->discord, THKey(dur_inf->discord, Mtail[0], Mhead[0])) != kh_none;
      
      // if it resides in any of the edgelists we store, we need to resample
      if(in_network) {
        if(in_discord) {
          UnsrtELGetRand(Mtail, Mhead, sto->discordantELs[i]);
        } else {
          UnsrtELGetRand(Mtail, Mhead, sto->nonDiscordantELs[i]);
        }
      } else if(in_discord) {
        UnsrtELGetRand(Mtail, Mhead, sto->discordantNonELs[i]);
      }
    } else {
      // propose toggling off an edge of the specified mixing type
      if(unif_rand() < sto->nonDiscordantELs[i]->nedges/((double) nedgestype)) {
        UnsrtELGetRand(Mtail, Mhead, sto->nonDiscordantELs[i]);
        in_discord = FALSE;
      } else {
        UnsrtELGetRand(Mtail, Mhead, sto->discordantELs[i]);
        in_discord = TRUE;
      }
      in_network = TRUE;
    }
  } else {
    // propose from discord
    if(unif_rand() < sto->discordantELs[i]->nedges/((double) nddyadstype)) {
      UnsrtELGetRand(Mtail, Mhead, sto->discordantELs[i]);
      in_network = FALSE;
    } else {
      UnsrtELGetRand(Mtail, Mhead, sto->discordantNonELs[i]);
      in_network = TRUE;
    }
    in_discord = TRUE;
  }
  
  // compute logratio
  
  // these ignore overall factor of 1/2 which cancels out in ratio
  // and overall factor of prob(mixing type), which also cancels out in ratio
  double forward_discord = in_discord ? 1.0/nddyadstype : 0;
  double backward_discord = in_discord ? 0 : 1.0/(1 + nddyadstype);
  
  double forward_network = in_network ? (0.5/nedgestype + 0.5/ndyadstype) : (nedgestype == 0 ? 1.0/ndyadstype : 0.5/ndyadstype);
  double backward_network = in_network ? (nedgestype == 1 ? 1.0/ndyadstype : 0.5/ndyadstype) : (0.5/(nedgestype + 1) + 0.5/ndyadstype);
  
  if(nddyadstype == 0) forward_network *= 2;
  if(nddyadstype == 1 && in_discord) backward_network *= 2;
  
  double forward = forward_discord + forward_network;
  double backward = backward_discord + backward_network;

  MHp->logratio = log(backward/forward);

  sto->in_discord = in_discord;
}

MH_U_FN(Mu_discordStratTNT) {
  // add or remove edge from appropriate edgelist
  GET_STORAGE(discordStratTNTStorage, sto);
  
  if(edgeflag) {
    // we are removing an existing edge
    if(sto->in_discord) {
      UnsrtELDelete(tail, head, sto->discordantELs[sto->mixingtype]);
    } else {
      UnsrtELDelete(tail, head, sto->nonDiscordantELs[sto->mixingtype]);
      UnsrtELInsert(tail, head, sto->discordantNonELs[sto->mixingtype]);
    }
  } else {
    // we are adding a new edge
    if(sto->in_discord) {
      UnsrtELInsert(tail, head, sto->nonDiscordantELs[sto->mixingtype]);
      UnsrtELDelete(tail, head, sto->discordantNonELs[sto->mixingtype]);
    } else {
      UnsrtELInsert(tail, head, sto->discordantELs[sto->mixingtype]);
    }
  }
}

MH_F_FN(Mf_discordStratTNT) {
  // Free all the things
  GET_STORAGE(discordStratTNTStorage, sto);
  
  for(int i = 0; i < sto->nmixtypes; i++) {
    UnsrtELDestroy(sto->nonDiscordantELs[i]);
    UnsrtELDestroy(sto->discordantELs[i]);
    UnsrtELDestroy(sto->discordantNonELs[i]);    
  }

  Free(sto->nonDiscordantELs);
  Free(sto->discordantELs);
  Free(sto->discordantNonELs);

  Free(sto->nodesbycode);

  // MHp->storage itself should be Freed by MHProposalDestroy
}


