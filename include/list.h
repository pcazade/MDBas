#ifndef LISTH_INCLUDED
#define LISTH_INCLUDED

#ifdef	__cplusplus
extern "C" {
#endif

    void allocate_heuristic(PARALLEL *parallel);

    void makelist(CTRL *ctrl,PARAM *param,PARALLEL *parallel,PBC *box,NEIGH *neigh,
                  CONSTRAINT constList[],BOND bond[],ANGLE angle[],DIHE dihe[],DIHE impr[],
		  double x[], double y[],double z[],double vx[], double vy[],double vz[],
		  int frozen[],int ***neighList,int **neighPair,int **neighList14,
                  int ***exclList,int **exclPair,int iBuffer[]);
    
    void heuristic_update(CTRL *ctrl,PARAM *param,PARALLEL *parallel,NEIGH *neigh,
		      double vx[],double vy[],double vz[],int iBuffer[]);

    void exclude_list(CTRL *ctrl,PARAM *param,PARALLEL *parallel,NEIGH *neigh,CONSTRAINT constList[],
                      BOND bond[],ANGLE angle[],DIHE dihe[],DIHE impr[],
                      int **neighList14,int ***exclList,int **exclPair);

    void init_verlet_list(PARAM *param,PARALLEL *parallel,PBC *box,NEIGH *neigh,double x[],double y[],double z[],
                          int frozen[],int ***neighList,int **neighPair,int **exclList,
                          int exclPair[]);

    void verlet_list(PARAM *param,PARALLEL *parallel,PBC *box,NEIGH *neigh,double x[],double y[],double z[],
                     int frozen[],int ***neighList,int neighPair[],int **exclList,
                     int exclPair[]);

    void init_link_cell_verlet_list(PARAM *param,PARALLEL *parallel,PBC *box,NEIGH *neigh,double x[],double y[],double z[],
                                    int frozen[],int ***neighList,int **neighPair,int **exclList,
                                    int exclPair[]);

    void link_cell_verlet_list(PARAM *param,PARALLEL *parallel,PBC *box,NEIGH *neigh,double x[],double y[],double z[],
                               int frozen[],int ***neighList,int neighPair[],int **exclList,
                               int exclPair[]);

// void fast_verlet_list(PARAM *param,PARALLEL *parallel,PBC *box,NEIGH *neigh,double x[],double y[],double z[],
// 		      int frozen[],int **neighList,int **neighPair,
// 		      int **exclList,int exclPair[]);

#ifdef	__cplusplus
}
#endif

#endif
