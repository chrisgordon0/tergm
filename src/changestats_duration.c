/*  File src/changestats_duration.c in package tergm, part of the Statnet suite
 *  of packages for network analysis, https://statnet.org .
 *
 *  This software is distributed under the GPL-3 license.  It is free,
 *  open source, and has the attribution requirements (GPL Section 7) at
 *  https://statnet.org/attribution
 *
 *  Copyright 2008-2019 Statnet Commons
 */
#include "changestats_duration.h"

#define CSD_TRANSFORM_ET(et)                    \
  double ett=0, ett1=1;                        \
  switch(transform) {                        \
    case 0: ett = et; ett1 = et+1; break;                \
    case 1: ett = log(et); ett1 = log(et+1); break;        \
    default: error("Unrecognized dyad age transformation code."); \
  }                                \
  (void) ett; (void) ett1; // Get rid of unused variable warnings, since either ett or ett1 may be unused.

/*****************
 edges_ageinterval

 This is essentially the edges statistic, which only counts dyads with "age"
 (time steps spent in the current state) in the interval [inputparams0,inputparams1).
*****************/

X_CHANGESTAT_FN(x_edges_ageinterval_mon){
  ZERO_ALL_CHANGESTATS();

  if(type == TICK) {
    GET_AUX_STORAGE(StoreTimeAndLasttoggle, dur_inf);
  
    for(Edge k=1; k <= N_EDGES; k++){
      Vertex tail, head;
      FindithEdge(&tail, &head, k, nwp);
      int age = ElapsedTime(tail,head,dur_inf) + 1; // lasttoggle auxiliary has *not* yet been updated via TICK
      for(unsigned int j=0; j<N_CHANGE_STATS; j++){
        unsigned int from = INPUT_PARAM[j*2], to = INPUT_PARAM[j*2+1];
        if(age+1 == from) CHANGE_STAT[j]++; // The tie "ages" into the interval.
        if(to!=0 && age+1 == to) CHANGE_STAT[j]--; // The tie "ages" out of the interval.
      }
    }
  }
}

C_CHANGESTAT_FN(c_edges_ageinterval_mon){
  GET_AUX_STORAGE(StoreTimeAndLasttoggle, dur_inf);

  int age = ElapsedTimeToggle(tail,head,dur_inf,edgeflag);
  // Only count if the age is in [from,to). ( to=0 ==> to=Inf )

  for(unsigned int j=0; j<N_CHANGE_STATS; j++){
    unsigned int from = INPUT_PARAM[j*2], to = INPUT_PARAM[j*2+1];
    if(edgeflag){ // If already an edge, we are dissolving.
      if(from<=age+1 && (to==0 || age+1<to)) CHANGE_STAT[j]--; // Statistic only changes if it's in the interval.
    }else{ // If not already an edge, we are forming.
      if(from<=age+1 && (to==0 || age+1<to)) CHANGE_STAT[j]++; // Statistic only changes if it's in the interval.
    }
  }
}

S_CHANGESTAT_FN(s_edges_ageinterval_mon){
  GET_AUX_STORAGE(StoreTimeAndLasttoggle, dur_inf);
  
  ZERO_ALL_CHANGESTATS(i);
  for (Edge k=1; k <= N_EDGES; k++){
    Vertex tail, head;
    FindithEdge(&tail, &head, k, nwp);
    int age = ElapsedTime(tail,head,dur_inf) + 1;
    for(unsigned int j=0; j<N_CHANGE_STATS; j++){
      unsigned int from = INPUT_PARAM[j*2], to = INPUT_PARAM[j*2+1];
      if(from<=age && (to==0 || age<to)) CHANGE_STAT[0]++;
    }
  }
}

/*****************
 edge_ages

 Sum of ages of all extant ties.

*****************/

X_CHANGESTAT_FN(x_edge_ages_mon){ 
  ZERO_ALL_CHANGESTATS();
  if(type == TICK) {
    CHANGE_STAT[0] = N_EDGES;
  }
}


C_CHANGESTAT_FN(c_edge_ages_mon){
  GET_AUX_STORAGE(StoreTimeAndLasttoggle, dur_inf);
  
  int age = ElapsedTimeToggle(tail,head,dur_inf,edgeflag);

  CHANGE_STAT[0] += edgeflag ? - age - 1 : age + 1;
}

S_CHANGESTAT_FN(s_edge_ages_mon){
  GET_AUX_STORAGE(StoreTimeAndLasttoggle, dur_inf);
  
  CHANGE_STAT[0] = 0;
  for (Edge k=1; k <= N_EDGES; k++){
    Vertex tail, head;
    FindithEdge(&tail, &head, k, nwp);
    int age = ElapsedTime(tail,head,dur_inf) + 1;
    CHANGE_STAT[0] += age;
  }
}

/*****************
 edgecov_ages

 Weighted sum of ages of all extant ties.

*****************/

X_CHANGESTAT_FN(x_edgecov_ages_mon){ 
  ZERO_ALL_CHANGESTATS();
  if(type == TICK) {
    int noffset = BIPARTITE, nrow;
    if(noffset > 0){
      nrow = noffset;
    }else{
      nrow = INPUT_PARAM[0];
    }
    
    // Sum of weights.
    for (Edge k=1; k <= N_EDGES; k++){
      Vertex tail, head;
      FindithEdge(&tail, &head, k, nwp);
      double val = INPUT_ATTRIB[(head-1-noffset)*nrow+(tail-1)];
      CHANGE_STAT[0] += val;
    }
  }
}

C_CHANGESTAT_FN(c_edgecov_ages_mon){
  GET_AUX_STORAGE(StoreTimeAndLasttoggle, dur_inf);
  
  int noffset = BIPARTITE, nrow;
  if(noffset > 0){
    nrow = noffset;
  }else{
    nrow = INPUT_PARAM[0];
  }

  double val = INPUT_ATTRIB[(head - 1 - noffset) * nrow + (tail - 1)];
  if(val!=0){
    int age = ElapsedTimeToggle(tail, head, dur_inf,edgeflag);

    CHANGE_STAT[0] += edgeflag ? - age*val - val : age*val + val;
  }
}

S_CHANGESTAT_FN(s_edgecov_ages_mon){
  GET_AUX_STORAGE(StoreTimeAndLasttoggle, dur_inf);
  
  int noffset = BIPARTITE, nrow;
  if(noffset > 0){
    nrow = noffset;
  }else{
    nrow = INPUT_PARAM[0];
  }

  CHANGE_STAT[0] = 0;
  for (Edge k=1; k <= N_EDGES; k++){
    Vertex tail, head;
    FindithEdge(&tail, &head, k, nwp);
    double val = INPUT_ATTRIB[(head-1-noffset)*nrow+(tail-1)];
    int age = ElapsedTime(tail,head,dur_inf) + 1;
    CHANGE_STAT[0] += age*val;
  }
}

/*****************
 mean_age

 Mean of (optionally log-) ages of all extant ties.

 The mean_ages of an empty network is defined to be emptyval.

 *****************/

I_CHANGESTAT_FN(i_mean_age_mon){
  ALLOC_STORAGE(2, double, sto);
  GET_AUX_STORAGE(StoreTimeAndLasttoggle, dur_inf);  
  int transform = INPUT_PARAM[1];
  
  for(Edge k=1; k <= N_EDGES; k++){
    Vertex tail, head;
    FindithEdge(&tail, &head, k, nwp);
    int et = ElapsedTime(tail,head,dur_inf);
    CSD_TRANSFORM_ET(et);
    sto[0] += ett1;
  }
}

X_CHANGESTAT_FN(x_mean_age_mon){
  ZERO_ALL_CHANGESTATS();
  GET_AUX_STORAGE(StoreTimeAndLasttoggle, dur_inf);  
  
  int transform = INPUT_PARAM[1];
  
  if(type == TICK) {
    GET_STORAGE(double, sto);
    
    if(transform == 0) {
      sto[0] += N_EDGES;
      CHANGE_STAT[0] = N_EDGES ? 1 : 0;
    } else {
      double oldval = sto[0];
      sto[0] = 0;
      for(Edge k=1; k <= N_EDGES; k++){
        Vertex tail, head;
        FindithEdge(&tail, &head, k, nwp);
        int et = ElapsedTime(tail,head,dur_inf) + 1;
        CSD_TRANSFORM_ET(et);
        sto[0] += ett1;
      }
      CHANGE_STAT[0] = N_EDGES ? (sto[0] - oldval)/N_EDGES : 0;
    }
  }
}

C_CHANGESTAT_FN(c_mean_age_mon){
  GET_STORAGE(double, sto);
  GET_AUX_STORAGE(StoreTimeAndLasttoggle, dur_inf);
    
  double s0 = sto[0], s1 = sto[0]; // Sum of age values of initial and final network.
  double zeroval = INPUT_PARAM[0]; // Empty network value.
  int transform = INPUT_PARAM[1]; // Transformation code.
  Edge e0, e1; // Number of edges in initial and final network.
  
  e0 = e1 = N_EDGES;
  
  int change = edgeflag ? -1 : 1;
  int et = ElapsedTimeToggle(tail,head,dur_inf,edgeflag);  
  CSD_TRANSFORM_ET(et);
  s1 += change*ett1;
  e1 += change;
  sto[1] = s1;

  CHANGE_STAT[0]=(e1==0?zeroval:s1/e1)-(e0==0?zeroval:s0/e0);
}

U_CHANGESTAT_FN(u_mean_age_mon){
  GET_STORAGE(double, sto);
  sto[0] = sto[1];
}

S_CHANGESTAT_FN(s_mean_age_mon){
  GET_AUX_STORAGE(StoreTimeAndLasttoggle, dur_inf);
  
  CHANGE_STAT[0] = 0;
  double zeroval = INPUT_PARAM[0];
  int transform = INPUT_PARAM[1]; // Transformation code.

  if(N_EDGES>0){
    for (Edge k=1; k <= N_EDGES; k++){
      Vertex tail, head;
      FindithEdge(&tail, &head, k, nwp);
      int et = ElapsedTime(tail,head,dur_inf);
      CSD_TRANSFORM_ET(et);
      CHANGE_STAT[0] += ett1;
    }
    
    CHANGE_STAT[0] /= N_EDGES;
  }else{
    CHANGE_STAT[0] = zeroval;
  }
}

/*****************
 edgecov_mean_age

 Weighted mean of ages of all extant ties.

 The edgecov_mean_ages of an empty network is defined to be emptyval.

 *****************/
I_CHANGESTAT_FN(i_edgecov_mean_age_mon) {
  ALLOC_STORAGE(4, double, sto);
  GET_AUX_STORAGE(StoreTimeAndLasttoggle, dur_inf);  
  int transform = INPUT_PARAM[1];
  
  int noffset = BIPARTITE, nrow;
  if(noffset > 0){
    nrow = noffset;
  }else{
    nrow = INPUT_PARAM[2];
  }
  
  for(Edge k=1; k <= N_EDGES; k++){
    Vertex tail, head;
    FindithEdge(&tail, &head, k, nwp);
    double val = INPUT_ATTRIB[(head - 1 - noffset) * nrow + (tail - 1)];   
    if(val!=0){   
      int et = ElapsedTime(tail,head,dur_inf);
      CSD_TRANSFORM_ET(et);
      sto[0] += ett1*val;
      sto[1] += val;
    }
  }
}

X_CHANGESTAT_FN(x_edgecov_mean_age_mon) {
  ZERO_ALL_CHANGESTATS();
  if(type == TICK) {
    GET_AUX_STORAGE(StoreTimeAndLasttoggle, dur_inf);  

    int noffset = BIPARTITE, nrow;
    if(noffset > 0){
      nrow = noffset;
    }else{
      nrow = INPUT_PARAM[2];
    }
  
    int transform = INPUT_PARAM[1];
    
    GET_STORAGE(double, sto);
    
    if(sto[1] != 0) {
      if(transform == 0) {
        sto[0] += sto[1];
        CHANGE_STAT[0] = 1;
      } else {
        double oldval = sto[0];
        sto[0] = 0;
        for(Edge k=1; k <= N_EDGES; k++) {
          Vertex tail, head;
          FindithEdge(&tail, &head, k, nwp);
          double val = INPUT_ATTRIB[(head - 1 - noffset) * nrow + (tail - 1)];   
          if(val!=0) {
            int et = ElapsedTime(tail,head,dur_inf) + 1;
            CSD_TRANSFORM_ET(et);
            sto[0] += ett1*val;
          }
        }
        CHANGE_STAT[0] = (sto[0] - oldval)/sto[1];
      }
    }
  }
}
  
C_CHANGESTAT_FN(c_edgecov_mean_age_mon){
  GET_AUX_STORAGE(StoreTimeAndLasttoggle, dur_inf);
  
  int noffset = BIPARTITE, nrow;
  if(noffset > 0){
    nrow = noffset;
  }else{
    nrow = INPUT_PARAM[2];
  }

  GET_STORAGE(double, sto);

  double s0 = sto[0], s1 = sto[0]; // Sum of age values of initial and final network.
  double zeroval = INPUT_PARAM[0]; // Empty network value.
  int transform = INPUT_PARAM[1]; // Transformation code.
  double e0 = sto[1], e1 = sto[1]; // Sum of edge weights in initial and final network.
  
  double val = INPUT_ATTRIB[(head - 1 - noffset) * nrow + (tail - 1)];   
  if(val!=0){
    if(edgeflag){
      int et = ElapsedTimeToggle(tail,head,dur_inf,edgeflag);
      CSD_TRANSFORM_ET(et);
      s1 -= ett1*val;
      e1 -= val;
    }else{
      int et = ElapsedTimeToggle(tail,head,dur_inf,edgeflag);
      CSD_TRANSFORM_ET(et);
      s1 += ett1*val;
      e1 += val;
    }
  }
  
  CHANGE_STAT[0]=(e1==0?zeroval:s1/e1)-(e0==0?zeroval:s0/e0);
  
  sto[2] = s1;
  sto[3] = e1;
}

U_CHANGESTAT_FN(u_edgecov_mean_age_mon){
  GET_STORAGE(double, sto);
  sto[0] = sto[2];
  sto[1] = sto[3];
}

S_CHANGESTAT_FN(s_edgecov_mean_age_mon){
  GET_AUX_STORAGE(StoreTimeAndLasttoggle, dur_inf);
  
  CHANGE_STAT[0] = 0;
  double zeroval = INPUT_PARAM[0], s=0, e=0;
  int transform = INPUT_PARAM[1]; // Transformation code.
  int noffset = BIPARTITE, nrow;
  if(noffset > 0){
    nrow = noffset;
  }else{
    nrow = INPUT_PARAM[2];
  }

  for (Edge k=1; k <= N_EDGES; k++){
    Vertex tail, head;
    FindithEdge(&tail, &head, k, nwp);
    double val = INPUT_ATTRIB[(head - 1 - noffset) * nrow + (tail - 1)];
    if(val!=0){
      int et = ElapsedTime(tail,head,dur_inf);    
      CSD_TRANSFORM_ET(et);
      s += ett1 * val;
      e += val;
    }
  }
   
  if(e!=0){
    CHANGE_STAT[0] = s/e;
  }else{
    CHANGE_STAT[0] = zeroval;
  }
}

/*****************
 degree_mean_age

 Mean of ages of all extant ties with a particular degree.

 The degree_mean_age of a network with no actors with degree of interest is defined to be emptyval.

 *****************/
I_CHANGESTAT_FN(i_degree_mean_age_mon){
  ALLOC_STORAGE(4, void *, sto);
  
  sto[0] = Calloc(N_CHANGE_STATS, double);
  sto[1] = Calloc(N_CHANGE_STATS, int);
  
  sto[2] = Calloc(N_CHANGE_STATS, double);
  sto[3] = Calloc(N_CHANGE_STATS, int);
  
  double *age = (double *)sto[0];
  int *count = (int *)sto[1];

  GET_AUX_STORAGE(StoreTimeAndLasttoggle, dur_inf);
    
  Vertex *id=IN_DEG, *od=OUT_DEG;
  int transform = INPUT_PARAM[1]; // Transformation code.
    
  for(unsigned int j = 0; j < N_CHANGE_STATS; j++){
    double s0 = 0;
    Edge e0 = 0;
  
    Vertex deg = INPUT_PARAM[j+2];
      
    for (Edge k=1; k <= N_EDGES; k++){
      Vertex tail, head;
      FindithEdge(&tail, &head, k, nwp);
        
      unsigned int w = (od[tail]+id[tail]==deg) + (od[head]+id[head]==deg);
        
      if(w){
        int et = ElapsedTime(tail,head,dur_inf);
        CSD_TRANSFORM_ET(et);
        s0 += ett1*w;
        e0+=w;
      }
    }
    
    age[j] = s0;
    count[j] = e0;
  }
}
 
X_CHANGESTAT_FN(x_degree_mean_age_mon){
  ZERO_ALL_CHANGESTATS();
  if(type == TICK) {
    GET_AUX_STORAGE(StoreTimeAndLasttoggle, dur_inf);
    GET_STORAGE(void *, sto);
    
    double *age = (double *)sto[0];
    int *count = (int *)sto[1];
   
    Vertex *id=IN_DEG, *od=OUT_DEG;
    double zeroval = INPUT_PARAM[0];
    int transform = INPUT_PARAM[1]; // Transformation code.

    for(unsigned int j = 0; j < N_CHANGE_STATS; j++) {
      double s0, s1;
      int e0;
      if(transform == 0) { // do it the fast way
        s0 = age[j];
        e0 = count[j];
        
        s1 = s0 + e0;
      } else { // transform == 1 and we need to do it the old way
        s0 = 0;
        s1 = 0;
        e0 = 0;
      
        Vertex deg = INPUT_PARAM[j+2];
        
        for (Edge k=1; k <= N_EDGES; k++){
          Vertex tail, head;
          FindithEdge(&tail, &head, k, nwp);
          
          unsigned int w = (od[tail]+id[tail]==deg) + (od[head]+id[head]==deg);
          
          if(w){
            int et = ElapsedTime(tail,head,dur_inf) + 1;
            CSD_TRANSFORM_ET(et);
            s0 += ett*w;
            s1 += ett1*w;
            e0+=w;
          }
        }      
      }
      
      CHANGE_STAT[j]=(e0==0?zeroval:s1/e0)-(e0==0?zeroval:s0/e0);
      
      age[j] = s1;      
    }
  }
}


C_CHANGESTAT_FN(c_degree_mean_age_mon){
  GET_AUX_STORAGE(StoreTimeAndLasttoggle, dur_inf);
  
  Vertex *id=IN_DEG, *od=OUT_DEG;
  double zeroval = INPUT_PARAM[0];
  int transform = INPUT_PARAM[1]; // Transformation code.
  
  GET_STORAGE(void *, sto);
    
  double *age = (double *)sto[0];
  int *count = (int *)sto[1];
  double *newage = (double *)sto[2];
  int *newcount = (int *)sto[3];

  
  
  for(unsigned int j = 0; j < N_CHANGE_STATS; j++){
    double s0 = age[j], s1 = age[j];
    Edge e0 = count[j], e1 = count[j];

    Vertex deg = INPUT_PARAM[j+2];
    
    int change = edgeflag ? -1 : +1;
    int taildiff = (od[tail]+id[tail] + change == deg)-(od[tail]+id[tail] == deg);
    int headdiff = (od[head]+id[head] + change == deg)-(od[head]+id[head] == deg);

    Edge e;
    Vertex head1, tail1;
    
    switch(taildiff){
      case -1: // tail was previously counted, but is no longer
        STEP_THROUGH_OUTEDGES(tail, e, head1){
          int et = ElapsedTime(tail,head1,dur_inf);
          CSD_TRANSFORM_ET(et);
          s1 -= ett1;
          e1--;
        }
        STEP_THROUGH_INEDGES(tail, e, head1){
          int et = ElapsedTime(head1,tail,dur_inf);
          CSD_TRANSFORM_ET(et);
          s1 -= ett1;
          e1--;
        }
        break;
      // We don't need to do anything special for the focus dyad here:
      // if it's formed, then it wasn't counted in the first place;
      // if it's dissolved, then it will have been subtracted off by the previous two loops.
      
      case +1: // tail was previously not counted, but is now
        STEP_THROUGH_OUTEDGES(tail, e, head1){
          int et = ElapsedTime(tail,head1,dur_inf);
          CSD_TRANSFORM_ET(et);
          s1 += ett1;
          e1++;
        }
        STEP_THROUGH_INEDGES(tail, e, head1){
          int et = ElapsedTime(head1,tail,dur_inf);
          CSD_TRANSFORM_ET(et);
          s1 += ett1;
          e1++;
        }
        // Here, we need to handle the focus dyad:
        if(change==+1){// if it's formed, add to s1
          int et = ElapsedTimeToggle(tail,head,dur_inf,edgeflag);
          CSD_TRANSFORM_ET(et);
          s1 += ett1;
          e1++;
        }else{// if it's dissolved, it had been counted in the previous two loops, and it should be subtracted
          int et = ElapsedTimeToggle(tail,head,dur_inf,edgeflag);
          CSD_TRANSFORM_ET(et);
          s1 -= ett1;
          e1--;
        }
        break;
    }

    switch(headdiff){
      case -1: // head was previously counted, but is no longer
        STEP_THROUGH_OUTEDGES(head, e, tail1){
          int et = ElapsedTime(head,tail1,dur_inf);
          CSD_TRANSFORM_ET(et);
          s1 -= ett1;
          e1--;
        }
        STEP_THROUGH_INEDGES(head, e, tail1){
          int et = ElapsedTime(tail1,head,dur_inf);
          CSD_TRANSFORM_ET(et);
          s1 -= ett1;
          e1--;
        }
        break;
      // We don't need to do anything special for the focus dyad here:
      // if it's formed, then it wasn't counted in the first place;
      // if it's dissolved, then it will have been subtracted off by the previous two loops.

      case +1: // head was previously not counted, but is now
        STEP_THROUGH_OUTEDGES(head, e, tail1){
          int et = ElapsedTime(head,tail1,dur_inf);
          CSD_TRANSFORM_ET(et);
          s1 += ett1;
          e1++;
        }
        STEP_THROUGH_INEDGES(head, e, tail1){
          int et = ElapsedTime(tail1,head,dur_inf);
          CSD_TRANSFORM_ET(et);
          s1 += ett1;
          e1++;
        }
        // Here, we need to handle the focus dyad:
        if(change==+1){// if it's formed, add to s1
          int et = ElapsedTimeToggle(tail,head,dur_inf,edgeflag);
          CSD_TRANSFORM_ET(et);
          s1 += ett1;
          e1++;
        }else{// if it's dissolved, it had been counted in the previous two loops, and it should be subtracted
          int et = ElapsedTimeToggle(tail,head,dur_inf,edgeflag);
          CSD_TRANSFORM_ET(et);
          s1 -= ett1;
          e1--;
        }
        break;
    }
  
    newage[j] = s1;
    newcount[j] = e1;
  
    CHANGE_STAT[j]=(e1==0?zeroval:s1/e1)-(e0==0?zeroval:s0/e0);
  }
}

U_CHANGESTAT_FN(u_degree_mean_age_mon){
  // FIXME: Do not assume that the dyad in this u_ call is the same as that in the last c_ call.
  GET_STORAGE(void *, sto);
    
  double *age = (double *)sto[0];
  int *count = (int *)sto[1];
  double *newage = (double *)sto[2];
  int *newcount = (int *)sto[3];

  memcpy(age, newage, N_CHANGE_STATS*sizeof(double));
  memcpy(count, newcount, N_CHANGE_STATS*sizeof(int));
}

F_CHANGESTAT_FN(f_degree_mean_age_mon){
  GET_STORAGE(void *, sto);

  Free(sto[3]);
  Free(sto[2]);
  Free(sto[1]);
  Free(sto[0]);  
}

S_CHANGESTAT_FN(s_degree_mean_age_mon){
  GET_AUX_STORAGE(StoreTimeAndLasttoggle, dur_inf);
  
  Vertex *id=IN_DEG, *od=OUT_DEG;
  double zeroval = INPUT_PARAM[0];
  int transform = INPUT_PARAM[1]; // Transformation code.

  ZERO_ALL_CHANGESTATS();

  for(unsigned int j = 0; j < N_CHANGE_STATS; j++){
    Vertex deg = INPUT_PARAM[j+2];
    Edge e=0;

    for (Edge k=1; k <= N_EDGES; k++){
      Vertex tail, head;
      FindithEdge(&tail, &head, k, nwp);
      
      unsigned int w = (od[tail]+id[tail]==deg ? 1:0) + (od[head]+id[head]==deg ? 1:0);
      
      if(w){
        int et = ElapsedTime(tail,head,dur_inf);
        CSD_TRANSFORM_ET(et);
        CHANGE_STAT[j] += ett1*w;
        e+=w;
      }
    }
    
    if(e) CHANGE_STAT[j] /= e;
    else CHANGE_STAT[j] = zeroval;
  }
}

/*****************
 degree_by_attr_mean_age

 Mean of ages of all extant ties with a particular degree, by actor attribute.

 The degree_by_attr_mean_age of a network with no actors with degree of interest is defined to be emptyval.

 *****************/
I_CHANGESTAT_FN(i_degree_by_attr_mean_age_mon){
  ALLOC_STORAGE(4, void *, sto);
  
  sto[0] = Calloc(N_CHANGE_STATS, double);
  sto[1] = Calloc(N_CHANGE_STATS, int);
  
  sto[2] = Calloc(N_CHANGE_STATS, double);
  sto[3] = Calloc(N_CHANGE_STATS, int);
  
  double *age = (double *)sto[0];
  int *count = (int *)sto[1];

  GET_AUX_STORAGE(StoreTimeAndLasttoggle, dur_inf);
    
  Vertex *id=IN_DEG, *od=OUT_DEG;
  int transform = INPUT_PARAM[1]; // Transformation code.
    
  for(unsigned int j = 0; j < N_CHANGE_STATS; j++){
    double s0 = 0;
    Edge e0 = 0;
  
    Vertex deg = INPUT_PARAM[2*j+2];
    int testattr = INPUT_PARAM[2*j+3];
      
    for (Edge k=1; k <= N_EDGES; k++){
      Vertex tail, head;
      FindithEdge(&tail, &head, k, nwp);
        
      Vertex taildeg = od[tail]+id[tail], headdeg = od[head]+id[head];
      int tailattr = INPUT_PARAM[2*N_CHANGE_STATS + tail + 1]; 
      int headattr = INPUT_PARAM[2*N_CHANGE_STATS + head + 1]; 
  
      unsigned int w = ((taildeg==deg && tailattr==testattr) ? 1 : 0) +
        ((headdeg==deg && headattr==testattr) ? 1 : 0);
        
      if(w){
        int et = ElapsedTime(tail,head,dur_inf);
        CSD_TRANSFORM_ET(et);
        s0 += ett1*w;
        e0+=w;
      }
    }
    
    age[j] = s0;
    count[j] = e0;
  }
}

X_CHANGESTAT_FN(x_degree_by_attr_mean_age_mon){
  ZERO_ALL_CHANGESTATS();
  if(type == TICK) {
    GET_AUX_STORAGE(StoreTimeAndLasttoggle, dur_inf);
    GET_STORAGE(void *, sto);
    
    double *age = (double *)sto[0];
    int *count = (int *)sto[1];
   
    Vertex *id=IN_DEG, *od=OUT_DEG;
    double zeroval = INPUT_PARAM[0];
    int transform = INPUT_PARAM[1]; // Transformation code.

    for(unsigned int j = 0; j < N_CHANGE_STATS; j++) {
      double s0, s1;
      int e0;
      if(transform == 0) { // do it the fast way
        s0 = age[j];
        e0 = count[j];
        
        s1 = s0 + e0;
      } else { // transform == 1 and we need to do it the old way
        s0 = 0;
        s1 = 0;
        e0 = 0;
      
        Vertex deg = INPUT_PARAM[2*j+2];
        int testattr = INPUT_PARAM[2*j+3];
        
        for (Edge k=1; k <= N_EDGES; k++){
          Vertex tail, head;
          FindithEdge(&tail, &head, k, nwp);
          
          Vertex taildeg = od[tail]+id[tail], headdeg = od[head]+id[head];
          int tailattr = INPUT_PARAM[2*N_CHANGE_STATS + tail + 1]; 
          int headattr = INPUT_PARAM[2*N_CHANGE_STATS + head + 1]; 
  
          unsigned int w = ((taildeg==deg && tailattr==testattr) ? 1 : 0) +
            ((headdeg==deg && headattr==testattr) ? 1 : 0);
          
          if(w){
            int et = ElapsedTime(tail,head,dur_inf) + 1;
            CSD_TRANSFORM_ET(et);
            s0 += ett*w;
            s1 += ett1*w;
            e0+=w;
          }
        }      
      }
      
      CHANGE_STAT[j]=(e0==0?zeroval:s1/e0)-(e0==0?zeroval:s0/e0);
      
      age[j] = s1;      
    }
  }
}

C_CHANGESTAT_FN(c_degree_by_attr_mean_age_mon){
  GET_AUX_STORAGE(StoreTimeAndLasttoggle, dur_inf);

  GET_STORAGE(void *, sto);
    
  double *age = (double *)sto[0];
  int *count = (int *)sto[1];
  double *newage = (double *)sto[2];
  int *newcount = (int *)sto[3];
    
  Vertex *id=IN_DEG, *od=OUT_DEG;
  double zeroval = INPUT_PARAM[0];
  int transform = INPUT_PARAM[1]; // Transformation code.
  
  for(unsigned int j = 0; j < N_CHANGE_STATS; j++){
    double s0 = age[j], s1 = age[j];
    Edge e0 = count[j], e1 = count[j];

    Vertex deg = INPUT_PARAM[2*j+2];
    int testattr = INPUT_PARAM[2*j+3];

    int tailattr = INPUT_PARAM[2*N_CHANGE_STATS + tail + 1]; 
    int headattr = INPUT_PARAM[2*N_CHANGE_STATS + head + 1]; 

    // If neither attribute matches, this toggle has no effect on the statistic.
    if(tailattr!=testattr && headattr!=testattr){
      newage[j] = age[j];
      newcount[j] = count[j];
      continue; 
    }

    int change = edgeflag ? -1 : +1;
    int taildiff = (od[tail]+id[tail] + change == deg)-(od[tail]+id[tail] == deg);
    int headdiff = (od[head]+id[head] + change == deg)-(od[head]+id[head] == deg);

    Edge e;
    Vertex head1, tail1;
    
    switch(taildiff * (tailattr==testattr)){ // If tailattr!=testattr, it'll look for case 0, i.e., do nothing.
      case -1: // tail was previously counted, but is no longer
        STEP_THROUGH_OUTEDGES(tail, e, head1){
          int et = ElapsedTime(tail,head1,dur_inf);
          CSD_TRANSFORM_ET(et);
          s1 -= ett1;
          e1--;
        }
        STEP_THROUGH_INEDGES(tail, e, head1){
          int et = ElapsedTime(head1,tail,dur_inf);
          CSD_TRANSFORM_ET(et);
          s1 -= ett1;
          e1--;
        }
        break;
      // We don't need to do anything special for the focus dyad here:
      // if it's formed, then it wasn't counted in the first place;
      // if it's dissolved, then it will have been subtracted off by the previous two loops.

      case +1: // tail was previously not counted, but is now
        STEP_THROUGH_OUTEDGES(tail, e, head1){
          int et = ElapsedTime(tail,head1,dur_inf);
          CSD_TRANSFORM_ET(et);
          s1 += ett1;
          e1++;
        }
        STEP_THROUGH_INEDGES(tail, e, head1){
          int et = ElapsedTime(head1,tail,dur_inf);
          CSD_TRANSFORM_ET(et);
          s1 += ett1;
          e1++;
        }
        // Here, we need to handle the focus dyad:
        if(change==+1){// if it's formed, add to s1
          int et = ElapsedTimeToggle(tail,head,dur_inf,edgeflag);
          CSD_TRANSFORM_ET(et);
          s1 += ett1;
          e1++;
        }else{// if it's dissolved, it had been counted in the previous two loops, and it should be subtracted
          int et = ElapsedTimeToggle(tail,head,dur_inf,edgeflag);
          CSD_TRANSFORM_ET(et);
          s1 -= ett1;
          e1--;
        }
        break;
    }

    switch(headdiff * (headattr==testattr)){ // If headattr!=testattr, it'll look for case 0, i.e., do nothing.
      case -1: // head was previously counted, but is no longer
        STEP_THROUGH_OUTEDGES(head, e, tail1){
          int et = ElapsedTime(head,tail1,dur_inf);
          CSD_TRANSFORM_ET(et);
          s1 -= ett1;
          e1--;
        }
        STEP_THROUGH_INEDGES(head, e, tail1){
          int et = ElapsedTime(tail1,head,dur_inf);
          CSD_TRANSFORM_ET(et);
          s1 -= ett1;
          e1--;
        }
        break;
      // We don't need to do anything special for the focus dyad here:
      // if it's formed, then it wasn't counted in the first place;
      // if it's dissolved, then it will have been subtracted off by the previous two loops.

      case +1: // head was previously not counted, but is now
        STEP_THROUGH_OUTEDGES(head, e, tail1){
          int et = ElapsedTime(head,tail1,dur_inf);
          CSD_TRANSFORM_ET(et);
          s1 += ett1;
          e1++;
        }
        STEP_THROUGH_INEDGES(head, e, tail1){
          int et = ElapsedTime(tail1,head,dur_inf);
          CSD_TRANSFORM_ET(et);
          s1 += ett1;
          e1++;
        }
        // Here, we need to handle the focus dyad:
        if(change==+1){// if it's formed, add to s1
          int et = ElapsedTimeToggle(tail,head,dur_inf,edgeflag);
          CSD_TRANSFORM_ET(et);
          s1 += ett1;
          e1++;
        }else{// if it's dissolved, it had been counted in the previous two loops, and it should be subtracted
          int et = ElapsedTimeToggle(tail,head,dur_inf,edgeflag);
          CSD_TRANSFORM_ET(et);
          s1 -= ett1;
          e1--;
        }
        break;
    }
  
    CHANGE_STAT[j]=(e1==0?zeroval:s1/e1)-(e0==0?zeroval:s0/e0);
    
    newcount[j] = e1;
    newage[j] = s1;
  }
}

U_CHANGESTAT_FN(u_degree_by_attr_mean_age_mon){
  // FIXME: Do not assume that the dyad in this u_ call is the same as that in the last c_ call.
  GET_STORAGE(void *, sto);
    
  double *age = (double *)sto[0];
  int *count = (int *)sto[1];
  double *newage = (double *)sto[2];
  int *newcount = (int *)sto[3];

  memcpy(age, newage, N_CHANGE_STATS*sizeof(double));
  memcpy(count, newcount, N_CHANGE_STATS*sizeof(int));
}

F_CHANGESTAT_FN(f_degree_by_attr_mean_age_mon){
  GET_STORAGE(void *, sto);

  Free(sto[3]);
  Free(sto[2]);
  Free(sto[1]);
  Free(sto[0]);  
}


S_CHANGESTAT_FN(s_degree_by_attr_mean_age_mon){
  GET_AUX_STORAGE(StoreTimeAndLasttoggle, dur_inf);
  
  Vertex *id=IN_DEG, *od=OUT_DEG;
  double zeroval = INPUT_PARAM[0];
  int transform = INPUT_PARAM[1]; // Transformation code.
  
  ZERO_ALL_CHANGESTATS(i);

  for(unsigned int j = 0; j < N_CHANGE_STATS; j++){
    Edge e=0;

    for (Edge k=1; k <= N_EDGES; k++){
      Vertex tail, head;
      FindithEdge(&tail, &head, k, nwp);
      Vertex taildeg = od[tail]+id[tail], headdeg = od[head]+id[head];
      int tailattr = INPUT_PARAM[2*N_CHANGE_STATS + tail + 1]; 
      int headattr = INPUT_PARAM[2*N_CHANGE_STATS + head + 1]; 
    
      Vertex deg = INPUT_PARAM[2*j+2];
      int testattr = INPUT_PARAM[2*j+3];

      unsigned int w = ((taildeg==deg && tailattr==testattr) ? 1 : 0) +
        ((headdeg==deg && headattr==testattr) ? 1 : 0);
      
      if(w){
        int et = ElapsedTime(tail,head,dur_inf);
        CSD_TRANSFORM_ET(et);
        CHANGE_STAT[j] += ett1*w;
        e+=w;
      }
    }
    
    if(e) CHANGE_STAT[j] /= e;
    else CHANGE_STAT[j] = zeroval;
  }
}


/*****************
 degrange_mean_age

 Mean of ages of all extant ties with a particular degree range.

 The degrange_mean_age of a network with no actors with degree of interest is defined to be emptyval.

 *****************/

// A macro indicating whether x is in [from,to)
#define FROM_TO(x, from, to) ((x)>=(from) && (x)<(to))

I_CHANGESTAT_FN(i_degrange_mean_age_mon){
  ALLOC_STORAGE(4, void *, sto);
  
  sto[0] = Calloc(N_CHANGE_STATS, double);
  sto[1] = Calloc(N_CHANGE_STATS, int);
  
  sto[2] = Calloc(N_CHANGE_STATS, double);
  sto[3] = Calloc(N_CHANGE_STATS, int);
  
  double *age = (double *)sto[0];
  int *count = (int *)sto[1];

  GET_AUX_STORAGE(StoreTimeAndLasttoggle, dur_inf);
    
  Vertex *id=IN_DEG, *od=OUT_DEG;
  int transform = INPUT_PARAM[1]; // Transformation code.
    
  for(unsigned int j = 0; j < N_CHANGE_STATS; j++){
    double s0 = 0;
    Edge e0 = 0;
  
    Vertex from = INPUT_PARAM[j*2+2], to = INPUT_PARAM[j*2+3];
      
    for (Edge k=1; k <= N_EDGES; k++){
      Vertex tail, head;
      FindithEdge(&tail, &head, k, nwp);
        
      unsigned int w = FROM_TO(od[tail]+id[tail],from,to) + FROM_TO(od[head]+id[head],from,to);
        
      if(w){
        int et = ElapsedTime(tail,head,dur_inf);
        CSD_TRANSFORM_ET(et);
        s0 += ett1*w;
        e0+=w;
      }
    }
    
    age[j] = s0;
    count[j] = e0;
  }
}

X_CHANGESTAT_FN(x_degrange_mean_age_mon){
  ZERO_ALL_CHANGESTATS();
  if(type == TICK) {
    GET_AUX_STORAGE(StoreTimeAndLasttoggle, dur_inf);
    GET_STORAGE(void *, sto);
    
    double *age = (double *)sto[0];
    int *count = (int *)sto[1];
   
    Vertex *id=IN_DEG, *od=OUT_DEG;
    double zeroval = INPUT_PARAM[0];
    int transform = INPUT_PARAM[1]; // Transformation code.

    for(unsigned int j = 0; j < N_CHANGE_STATS; j++) {
      double s0, s1;
      int e0;
      if(transform == 0) { // do it the fast way
        s0 = age[j];
        e0 = count[j];
        
        s1 = s0 + e0;
      } else { // transform == 1 and we need to do it the old way
        s0 = 0;
        s1 = 0;
        e0 = 0;
      
        Vertex from = INPUT_PARAM[j*2+2], to = INPUT_PARAM[j*2+3];        
        
        for (Edge k=1; k <= N_EDGES; k++){
          Vertex tail, head;
          FindithEdge(&tail, &head, k, nwp);
  
          unsigned int w = FROM_TO(od[tail]+id[tail],from,to) + FROM_TO(od[head]+id[head],from,to);
          
          if(w){
            int et = ElapsedTime(tail,head,dur_inf) + 1;
            CSD_TRANSFORM_ET(et);
            s0 += ett*w;
            s1 += ett1*w;
            e0+=w;
          }
        }      
      }
      
      CHANGE_STAT[j]=(e0==0?zeroval:s1/e0)-(e0==0?zeroval:s0/e0);
      
      age[j] = s1;      
    }
  }
}

C_CHANGESTAT_FN(c_degrange_mean_age_mon){
  GET_AUX_STORAGE(StoreTimeAndLasttoggle, dur_inf);
  
  Vertex *id=IN_DEG, *od=OUT_DEG;
  double zeroval = INPUT_PARAM[0];
  int transform = INPUT_PARAM[1]; // Transformation code.
  
  GET_STORAGE(void *, sto);
    
  double *age = (double *)sto[0];
  int *count = (int *)sto[1];
  double *newage = (double *)sto[2];
  int *newcount = (int *)sto[3];
  
  
  for(unsigned int j = 0; j < N_CHANGE_STATS; j++){
    double s0 = age[j], s1 = age[j];
    Edge e0 = count[j], e1 = count[j];

    Vertex from = INPUT_PARAM[j*2+2], to = INPUT_PARAM[j*2+3];

    int change = edgeflag ? -1 : +1;
    // In the degree range case, it's possible to gain or lose a tie without entering or exiting a given degree range.
    unsigned int tailin1 = FROM_TO(od[tail]+id[tail] + change, from, to),
      tailin0 = FROM_TO(od[tail]+id[tail], from, to),
      headin1 = FROM_TO(od[head]+id[head] + change, from, to),
      headin0 = FROM_TO(od[head]+id[head], from, to);
    
    Edge e;
    Vertex head1, tail1;

    if(tailin0 && !tailin1){ // tail was previously counted, but is no longer
      STEP_THROUGH_OUTEDGES(tail, e, head1){
        int et = ElapsedTime(tail,head1,dur_inf);
        CSD_TRANSFORM_ET(et);
        s1 -= ett1;
        e1--;
      }
      STEP_THROUGH_INEDGES(tail, e, head1){
        int et = ElapsedTime(head1,tail,dur_inf);
        CSD_TRANSFORM_ET(et);
        s1 -= ett1;
        e1--;
      }
      // We don't need to do anything special for the focus dyad here:
      // if it's formed, then it wasn't counted in the first place;
      // if it's dissolved, then it will have been subtracted off by the previous two loops.
    }else if(!tailin0 && tailin1){ // tail was previously not counted, but is now
      STEP_THROUGH_OUTEDGES(tail, e, head1){
        int et = ElapsedTime(tail,head1,dur_inf);
        CSD_TRANSFORM_ET(et);
        s1 += ett1;
        e1++;
      }
      STEP_THROUGH_INEDGES(tail, e, head1){
        int et = ElapsedTime(head1,tail,dur_inf);
        CSD_TRANSFORM_ET(et);
        s1 += ett1;
        e1++;
      }
      // Here, we need to handle the focus dyad:
      if(change==+1){// if it's formed, add to s1
        int et = ElapsedTimeToggle(tail,head,dur_inf,edgeflag);
        CSD_TRANSFORM_ET(et);
        s1 += ett1;
        e1++;
      }else{// if it's dissolved, it had been counted in the previous two loops, and it should be subtracted
        int et = ElapsedTimeToggle(tail,head,dur_inf,edgeflag);
        CSD_TRANSFORM_ET(et);
        s1 -= ett1;
        e1--;
      }
    }else if(tailin0 && tailin1){ // tail was counted both times, but we need to handle the focus dyad
      if(change==+1){// if it's formed, add to s1
        int et = ElapsedTimeToggle(tail,head,dur_inf,edgeflag);
        CSD_TRANSFORM_ET(et);
        s1 += ett1;
        e1++;
      }else{// if it's dissolved, it must be subtracted from s1
        int et = ElapsedTimeToggle(tail,head,dur_inf,edgeflag);
        CSD_TRANSFORM_ET(et);
        s1 -= ett1;
        e1--;
      }
    }
    // If !tailin0 && !tailin1, then it made no difference.
    
    if(headin0 && !headin1){ // head was previously counted, but is no longer
      STEP_THROUGH_OUTEDGES(head, e, tail1){
        int et = ElapsedTime(head,tail1,dur_inf);
        CSD_TRANSFORM_ET(et);
        s1 -= ett1;
        e1--;
      }
      STEP_THROUGH_INEDGES(head, e, tail1){
        int et = ElapsedTime(tail1,head,dur_inf);
        CSD_TRANSFORM_ET(et);
        s1 -= ett1;
        e1--;
      }
      // We don't need to do anything special for the focus dyad here:
      // if it's formed, then it wasn't counted in the first place;
      // if it's dissolved, then it will have been subtracted off by the previous two loops.
    }else if(!headin0 && headin1){ // head was previously not counted, but is now
      STEP_THROUGH_OUTEDGES(head, e, tail1){
        int et = ElapsedTime(head,tail1,dur_inf);
        CSD_TRANSFORM_ET(et);
        s1 += ett1;
        e1++;
      }
      STEP_THROUGH_INEDGES(head, e, tail1){
        int et = ElapsedTime(tail1,head,dur_inf);
        CSD_TRANSFORM_ET(et);
        s1 += ett1;
        e1++;
      }
      // Here, we need to handle the focus dyad:
      if(change==+1){// if it's formed, add to s1
        int et = ElapsedTimeToggle(tail,head,dur_inf,edgeflag);
        CSD_TRANSFORM_ET(et);
        s1 += ett1;
        e1++;
      }else{// if it's dissolved, it had been counted in the previous two loops, and it should be subtracted
        int et = ElapsedTimeToggle(tail,head,dur_inf,edgeflag);
        CSD_TRANSFORM_ET(et);
        s1 -= ett1;
        e1--;
      }
    }else if(headin0 && headin1){ // tail was counted both times, but we need to handle the focus dyad
      if(change==+1){// if it's formed, add to s1
        int et = ElapsedTimeToggle(tail,head,dur_inf,edgeflag);
        CSD_TRANSFORM_ET(et);
        s1 += ett1;
        e1++;
      }else{// if it's dissolved, it must be subtracted from s1
        int et = ElapsedTimeToggle(tail,head,dur_inf,edgeflag);
        CSD_TRANSFORM_ET(et);
        s1 -= ett1;
        e1--;
      }
    }
    // If !headin0 && !headin1, then it made no difference.
  
    CHANGE_STAT[j]=(e1==0?zeroval:s1/e1)-(e0==0?zeroval:s0/e0);
    
    newage[j] = s1;
    newcount[j] = e1;
  }
}

U_CHANGESTAT_FN(u_degrange_mean_age_mon){
  // FIXME: Do not assume that the dyad in this u_ call is the same as that in the last c_ call.
  GET_STORAGE(void *, sto);
    
  double *age = (double *)sto[0];
  int *count = (int *)sto[1];
  double *newage = (double *)sto[2];
  int *newcount = (int *)sto[3];

  memcpy(age, newage, N_CHANGE_STATS*sizeof(double));
  memcpy(count, newcount, N_CHANGE_STATS*sizeof(int));
}

F_CHANGESTAT_FN(f_degrange_mean_age_mon){
  GET_STORAGE(void *, sto);

  Free(sto[3]);
  Free(sto[2]);
  Free(sto[1]);
  Free(sto[0]);  
}


S_CHANGESTAT_FN(s_degrange_mean_age_mon){
  GET_AUX_STORAGE(StoreTimeAndLasttoggle, dur_inf);
  
  Vertex *id=IN_DEG, *od=OUT_DEG;
  double zeroval = INPUT_PARAM[0];
  int transform = INPUT_PARAM[1]; // Transformation code.
  
  ZERO_ALL_CHANGESTATS(i);

  for(unsigned int j = 0; j < N_CHANGE_STATS; j++){
    Vertex from = INPUT_PARAM[j*2+2], to = INPUT_PARAM[j*2+3];
    Edge e=0;

    for (Edge k=1; k <= N_EDGES; k++){
      Vertex tail, head;
      FindithEdge(&tail, &head, k, nwp);
      
      unsigned int w = FROM_TO(od[tail]+id[tail],from,to) + FROM_TO(od[head]+id[head],from,to);
      
      if(w){
        int et = ElapsedTime(tail,head,dur_inf);    
        CSD_TRANSFORM_ET(et);
        CHANGE_STAT[j] += ett1*w;
        e+=w;
      }
    }
    
    if(e) CHANGE_STAT[j] /= e;
    else CHANGE_STAT[j] = zeroval;
  }
}

/*****************
 degrange_by_attr_mean_age

 Mean of ages of all extant ties with a particular degree, by actor attribute.

 The degrange_by_attr_mean_age of a network with no actors with degree of interest is defined to be emptyval.

 *****************/

I_CHANGESTAT_FN(i_degrange_by_attr_mean_age_mon){
  ALLOC_STORAGE(4, void *, sto);
  
  sto[0] = Calloc(N_CHANGE_STATS, double);
  sto[1] = Calloc(N_CHANGE_STATS, int);
  
  sto[2] = Calloc(N_CHANGE_STATS, double);
  sto[3] = Calloc(N_CHANGE_STATS, int);
  
  double *age = (double *)sto[0];
  int *count = (int *)sto[1];

  GET_AUX_STORAGE(StoreTimeAndLasttoggle, dur_inf);
    
  Vertex *id=IN_DEG, *od=OUT_DEG;
  int transform = INPUT_PARAM[1]; // Transformation code.
    
  for(unsigned int j = 0; j < N_CHANGE_STATS; j++){
    double s0 = 0;
    Edge e0 = 0;
  
    Vertex from = INPUT_PARAM[3*j+2], to = INPUT_PARAM[3*j+3];
    int testattr = INPUT_PARAM[3*j+4];

      
    for (Edge k=1; k <= N_EDGES; k++){
      Vertex tail, head;
      FindithEdge(&tail, &head, k, nwp);
        
      Vertex taildeg = od[tail]+id[tail], headdeg = od[head]+id[head];
      int tailattr = INPUT_PARAM[3*N_CHANGE_STATS + tail + 1]; 
      int headattr = INPUT_PARAM[3*N_CHANGE_STATS + head + 1]; 
   
      unsigned int w = (FROM_TO(taildeg, from, to) && tailattr==testattr) +
        (FROM_TO(headdeg, from, to) && headattr==testattr);
        
      if(w){
        int et = ElapsedTime(tail,head,dur_inf);
        CSD_TRANSFORM_ET(et);
        s0 += ett1*w;
        e0+=w;
      }
    }
    
    age[j] = s0;
    count[j] = e0;
  }
}

X_CHANGESTAT_FN(x_degrange_by_attr_mean_age_mon){
  ZERO_ALL_CHANGESTATS();
  if(type == TICK) {
    GET_AUX_STORAGE(StoreTimeAndLasttoggle, dur_inf);
    GET_STORAGE(void *, sto);
    
    double *age = (double *)sto[0];
    int *count = (int *)sto[1];
   
    Vertex *id=IN_DEG, *od=OUT_DEG;
    double zeroval = INPUT_PARAM[0];
    int transform = INPUT_PARAM[1]; // Transformation code.

    for(unsigned int j = 0; j < N_CHANGE_STATS; j++) {
      double s0, s1;
      int e0;
      if(transform == 0) { // do it the fast way
        s0 = age[j];
        e0 = count[j];
        
        s1 = s0 + e0;
      } else { // transform == 1 and we need to do it the old way
        s0 = 0;
        s1 = 0;
        e0 = 0;
      
        Vertex from = INPUT_PARAM[3*j+2], to = INPUT_PARAM[3*j+3];
        int testattr = INPUT_PARAM[3*j+4];
        
        for (Edge k=1; k <= N_EDGES; k++){
          Vertex tail, head;
          FindithEdge(&tail, &head, k, nwp);
  
          Vertex taildeg = od[tail]+id[tail], headdeg = od[head]+id[head];
          int tailattr = INPUT_PARAM[3*N_CHANGE_STATS + tail + 1]; 
          int headattr = INPUT_PARAM[3*N_CHANGE_STATS + head + 1]; 
   
          unsigned int w = (FROM_TO(taildeg, from, to) && tailattr==testattr) +
            (FROM_TO(headdeg, from, to) && headattr==testattr);
          
          if(w){
            int et = ElapsedTime(tail,head,dur_inf) + 1;
            CSD_TRANSFORM_ET(et);
            s0 += ett*w;
            s1 += ett1*w;
            e0+=w;
          }
        }      
      }
      
      CHANGE_STAT[j]=(e0==0?zeroval:s1/e0)-(e0==0?zeroval:s0/e0);
      
      age[j] = s1;      
    }
  }
}

C_CHANGESTAT_FN(c_degrange_by_attr_mean_age_mon){
  GET_AUX_STORAGE(StoreTimeAndLasttoggle, dur_inf);
  
  Vertex *id=IN_DEG, *od=OUT_DEG;
  double zeroval = INPUT_PARAM[0];
  int transform = INPUT_PARAM[1]; // Transformation code.
  
  GET_STORAGE(void *, sto);
    
  double *age = (double *)sto[0];
  int *count = (int *)sto[1];
  double *newage = (double *)sto[2];
  int *newcount = (int *)sto[3];
  
  for(unsigned int j = 0; j < N_CHANGE_STATS; j++){
    double s0 = age[j], s1 = age[j];
    Edge e0 = count[j], e1 = count[j];

    Vertex from = INPUT_PARAM[3*j+2], to = INPUT_PARAM[3*j+3];
    int testattr = INPUT_PARAM[3*j+4];
    
    int tailattr = INPUT_PARAM[3*N_CHANGE_STATS + tail + 1]; 
    int headattr = INPUT_PARAM[3*N_CHANGE_STATS + head + 1]; 

    // If neither attribute matches, this toggle has no effect on the statistic.
    if(tailattr!=testattr && headattr!=testattr){
      newage[j] = age[j];
      newcount[j] = count[j];
      continue; 
    }

    int change = edgeflag ? -1 : +1;
    // In the degree range case, it's possible to gain or lose a tie without entering or exiting a given degree range.
    unsigned int tailin1 = FROM_TO(od[tail]+id[tail] + change, from, to),
      tailin0 = FROM_TO(od[tail]+id[tail], from, to),
      headin1 = FROM_TO(od[head]+id[head] + change, from, to),
      headin0 = FROM_TO(od[head]+id[head], from, to);

    Edge e;
    Vertex head1, tail1;
    
    if(tailattr==testattr){
      if(tailin0 && !tailin1){ // tail was previously counted, but is no longer
        STEP_THROUGH_OUTEDGES(tail, e, head1){
          int et = ElapsedTime(tail,head1,dur_inf);
          CSD_TRANSFORM_ET(et);
          s1 -= ett1;
          e1--;
        }
        STEP_THROUGH_INEDGES(tail, e, head1){
          int et = ElapsedTime(head1,tail,dur_inf);
          CSD_TRANSFORM_ET(et);
          s1 -= ett1;
          e1--;
        }
        // We don't need to do anything special for the focus dyad here:
        // if it's formed, then it wasn't counted in the first place;
        // if it's dissolved, then it will have been subtracted off by the previous two loops.
      }else if(!tailin0 && tailin1){ // tail was previously not counted, but is now
        STEP_THROUGH_OUTEDGES(tail, e, head1){
          int et = ElapsedTime(tail,head1,dur_inf);
          CSD_TRANSFORM_ET(et);
          s1 += ett1;
          e1++;
        }
        STEP_THROUGH_INEDGES(tail, e, head1){
          int et = ElapsedTime(head1,tail,dur_inf);
          CSD_TRANSFORM_ET(et);
          s1 += ett1;
          e1++;
        }
        // Here, we need to handle the focus dyad:
        if(change==+1){// if it's formed, add to s1
          int et = ElapsedTimeToggle(tail,head,dur_inf,edgeflag);
          CSD_TRANSFORM_ET(et);
          s1 += ett1;
          e1++;
        }else{// if it's dissolved, it had been counted in the previous two loops, and it should be subtracted
          int et = ElapsedTimeToggle(tail,head,dur_inf,edgeflag);
          CSD_TRANSFORM_ET(et);
          s1 -= ett1;
          e1--;
        }
      }else if(tailin0 && tailin1){ // tail was counted both times, but we need to handle the focus dyad
        if(change==+1){// if it's formed, add to s1
          int et = ElapsedTimeToggle(tail,head,dur_inf,edgeflag);
          CSD_TRANSFORM_ET(et);
          s1 += ett1;
          e1++;
        }else{// if it's dissolved, it must be subtracted from s1
          int et = ElapsedTimeToggle(tail,head,dur_inf,edgeflag);
          CSD_TRANSFORM_ET(et);
          s1 -= ett1;
          e1--;
        }
      }
      // If !tailin0 && !tailin1, then it made no difference.
    }
    
    if(headattr==testattr){
      if(headin0 && !headin1){ // head was previously counted, but is no longer
        STEP_THROUGH_OUTEDGES(head, e, tail1){
          int et = ElapsedTime(head,tail1,dur_inf);
          CSD_TRANSFORM_ET(et);
          s1 -= ett1;
          e1--;
        }
        STEP_THROUGH_INEDGES(head, e, tail1){
          int et = ElapsedTime(tail1,head,dur_inf);
          CSD_TRANSFORM_ET(et);
          s1 -= ett1;
          e1--;
        }
        // We don't need to do anything special for the focus dyad here:
        // if it's formed, then it wasn't counted in the first place;
        // if it's dissolved, then it will have been subtracted off by the previous two loops.
      }else if(!headin0 && headin1){ // head was previously not counted, but is now
        STEP_THROUGH_OUTEDGES(head, e, tail1){
          int et = ElapsedTime(head,tail1,dur_inf);
          CSD_TRANSFORM_ET(et);
          s1 += ett1;
          e1++;
        }
        STEP_THROUGH_INEDGES(head, e, tail1){
          int et = ElapsedTime(tail1,head,dur_inf);
          CSD_TRANSFORM_ET(et);
          s1 += ett1;
          e1++;
        }
        // Here, we need to handle the focus dyad:
        if(change==+1){// if it's formed, add to s1
          int et = ElapsedTimeToggle(tail,head,dur_inf,edgeflag);
          CSD_TRANSFORM_ET(et);
          s1 += ett1;
          e1++;
        }else{// if it's dissolved, it had been counted in the previous two loops, and it should be subtracted
          int et = ElapsedTimeToggle(tail,head,dur_inf,edgeflag);
          CSD_TRANSFORM_ET(et);
          s1 -= ett1;
          e1--;
        }
      }else if(headin0 && headin1){ // tail was counted both times, but we need to handle the focus dyad
        if(change==+1){// if it's formed, add to s1
          int et = ElapsedTimeToggle(tail,head,dur_inf,edgeflag);
          CSD_TRANSFORM_ET(et);
          s1 += ett1;
          e1++;
        }else{// if it's dissolved, it must be subtracted from s1
          int et = ElapsedTimeToggle(tail,head,dur_inf,edgeflag);
          CSD_TRANSFORM_ET(et);
          s1 -= ett1;
          e1--;
        }
      }
      // If !headin0 && !headin1, then it made no difference.
    }
  
    CHANGE_STAT[j]=(e1==0?zeroval:s1/e1)-(e0==0?zeroval:s0/e0);
    
    newage[j] = s1;
    newcount[j] = e1;
  }
}

U_CHANGESTAT_FN(u_degrange_by_attr_mean_age_mon){
  // FIXME: Do not assume that the dyad in this u_ call is the same as that in the last c_ call.
  GET_STORAGE(void *, sto);
    
  double *age = (double *)sto[0];
  int *count = (int *)sto[1];
  double *newage = (double *)sto[2];
  int *newcount = (int *)sto[3];

  memcpy(age, newage, N_CHANGE_STATS*sizeof(double));
  memcpy(count, newcount, N_CHANGE_STATS*sizeof(int));
}

F_CHANGESTAT_FN(f_degrange_by_attr_mean_age_mon){
  GET_STORAGE(void *, sto);

  Free(sto[3]);
  Free(sto[2]);
  Free(sto[1]);
  Free(sto[0]);  
}


S_CHANGESTAT_FN(s_degrange_by_attr_mean_age_mon){
  GET_AUX_STORAGE(StoreTimeAndLasttoggle, dur_inf);
  
  Vertex *id=IN_DEG, *od=OUT_DEG;
  double zeroval = INPUT_PARAM[0];
  int transform = INPUT_PARAM[1]; // Transformation code.
  ZERO_ALL_CHANGESTATS(i);

  for(unsigned int j = 0; j < N_CHANGE_STATS; j++){
    Edge e=0;

    for (Edge k=1; k <= N_EDGES; k++){
      Vertex tail, head;
      FindithEdge(&tail, &head, k, nwp);
      Vertex taildeg = od[tail]+id[tail], headdeg = od[head]+id[head];
      int tailattr = INPUT_PARAM[3*N_CHANGE_STATS + tail + 1]; 
      int headattr = INPUT_PARAM[3*N_CHANGE_STATS + head + 1]; 
    
      Vertex from = INPUT_PARAM[3*j+2], to = INPUT_PARAM[3*j+3];
      int testattr = INPUT_PARAM[3*j+4];

      unsigned int w = (FROM_TO(taildeg, from, to) && tailattr==testattr) +
        (FROM_TO(headdeg, from, to) && headattr==testattr);
      
      if(w){
        int et = ElapsedTime(tail,head,dur_inf);
        CSD_TRANSFORM_ET(et);
        CHANGE_STAT[j] += ett1*w;
        e+=w;
      }
    }
    
    if(e) CHANGE_STAT[j] /= e;
    else CHANGE_STAT[j] = zeroval;
  }
}

#undef FROM_TO
#undef CSD_TRANSFORM_ET
