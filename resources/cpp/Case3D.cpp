#include "omp.h"
#include "stdio.h"
#include "time.h"
#include "math.h"
#include "memory.h"
#include "model.h"


//#include "stdafx.h"
int LoadDataToArrays();
//int LoadDensity(double r);
int SaveDataToFiles();
int SetActualPressure(int sz, double r, double **pPS, double *pPrs, int dir);
int Receive(int argc, char* argv[]);
int CalcParallelPart(double &TotTimer);
int RetArgMod(int d, double *parg, double *pmod);
int CompareRes(double r, double TotTimer, int Proc_Num);
int SetUVatTime(double r, double *pu, double *pv);
int SetTWatTime(double r, double *pt, double *pw);
double CompareArrays(double *p1, double *p2, double &Mx);
double QuadrInterpol(double r, double f0, double f1, double f2);
void SetBoundaryCond(int sz, double r, double **pRS, double **pRS2, double *pR, int val, int dir);
void StoreToArray(int str, int sz, double *pSrc, double *pDst);
void SetActualUV(int str, int sz, double *pU, double *pV, int dir);
void SetActualT(int str, int sz, double *pT, int dir);
void InitLocalArrays(int str, int sz, double **pPS, double **pWS, double **pRS, double **pRS2, int val, int dir);
// DEBUG
void SaveArray(char *name, double *pAr, int sz1, int sz2, int sz3, int ijk=1);


const char pVrtLvlFlName[VER_NMB][32]=
						{ "08000d", "09000d", "10000d", "11000d",
						  "12000d",	"13000d", "14000d", "15000d",
						  "16000d", "17000d", "18000d" };


//CALC

extern	double	*PP[SRS_NMB];
extern  double  WW_FOR_CALC[SRS_NMB][SZ_X1+1][SZ_X2+1][SZ_X3+1];
extern  double  tau;
extern  int		M_prm;



// Set as a command line
int Subdomains=2, // number of subdomains in each decoposion
	M_prm=1, // m-parameter of addit.-aver. splitting 
	TmKnotsPer12h=3600*12, // how many points will be on 12 hours period
	TmLimCalc=3600*1; //(in sec.) calc. will be done up to this time
	

double tau, Timer[3][3][5][4]={0.}; // [dir][val][dom][stg]

double	*TT[SRS_NMB], // for interpolating
		*PP[SRS_NMB], // for interpolating
		*UU[SRS_NMB],// for interpolating 
		*VV[SRS_NMB],// for interpolating
		*WW[SRS_NMB];// for interpolating

double	pD[(SZ_X3+1)*(SZ_X2+1)*(SZ_X1+1)]={0.}; // !!! is not used !!!
double	pUU[(SZ_X3+1)*(SZ_X2+1)*(SZ_X1+1)]={0.}; // keep the result
double	pVV[(SZ_X3+1)*(SZ_X2+1)*(SZ_X1+1)]={0.}; // keep the result
double	pTT[(SZ_X3+1)*(SZ_X2+1)*(SZ_X1+1)]={0.}; // keep the result
double	pTemp[(SZ_X3+1)*(SZ_X2+1)*(SZ_X1+1)]={0.}; // keep the result




int main(int argc, char* argv[])
{
 int Proc_Num=1;
 double TotTimer=0.; // timer for parallel part of task

 // P R E P A R E   P A R T
 if(Receive(argc, argv)!=0)	 return 1;
 LoadDataToArrays();
 //LoadDensity(0.5*TmLimCalc/CALC_TIME);
 //SaveArray("OUT/pU1", pUU, SZ_X3, SZ_X2, SZ_X1);

 // P A R A L L E L   P A R T
 Proc_Num=CalcParallelPart(TotTimer);
 // E N D   P A R T
 CompareRes((double)TmLimCalc/CALC_TIME, TotTimer, Proc_Num);
 
 SaveDataToFiles();
 //SaveDebug();

 return 0;
}


int CalcParallelPart(double &TotTimer)
{
 const int sdm_sz=(SZ_X1-1)/Subdomains+1; // lam-size of subdomains
 int NmbThreads;
 double *pOut[DIR][EVL_MT_VAL][MAX_SUBDMN];

 TotTimer=omp_get_wtime();
 omp_set_num_threads(Subdomains*DIR*EVL_MT_VAL);
 #pragma omp parallel
 {
  int dom, val, dir, proc, n, i, j, k, d, d1;
  unsigned int  sz_bite;
  double pTask1D[6*(SZ_X1+1)]={0.}, *pR1, *pR2, *pR,
		 *pF, *pC, *pV, *pU, *pW, *pP, *pMu, *pD,
		 *pPS[SRS_NMB], *pWS[SRS_NMB], *pRS[SRS_NMB],
		 *pRS2[SRS_NMB];

  double tick=0; 
  #pragma omp master
  NmbThreads=omp_get_num_threads(); 

  proc=omp_get_thread_num();
  dom=proc%Subdomains;
  val=proc/(Subdomains*DIR);
  dir=(proc%(Subdomains*DIR))/Subdomains; // 0 - lam, 1 - fi, 2 - z
    
  if(dir>0) sz_bite=(SZ_X3+1)*(SZ_X2+1)*(sdm_sz+1);
  else sz_bite=(SZ_X3+1)*(SZ_X1+1)*(sdm_sz+1);

  pR1=new double [sz_bite]; // overlapped domains
  pR2=new double [sz_bite];
  pF=new double [sz_bite];
  pC=new double [sz_bite];
  pU=new double [sz_bite];
  pV=new double [sz_bite];
  pW=new double [sz_bite];
  pP=new double [sz_bite];
  pD=new double [sz_bite];
  pMu=new double [sz_bite];
 
  pOut[dir][val][dom]=pR1;
  for(i=0; i<=(SRS_NMB-1); i++)
  {
   pPS[i]=new double [sz_bite];
   pWS[i]=new double [sz_bite];
   pRS[i]=new double [sz_bite];
   if(val<2) pRS2[i]=new double [sz_bite];	
  }
  // Init all arrays
 // memset(pMu, 0, sz_bite*sizeof(double));
  InitLocalArrays(dom*(sdm_sz-1), (dir>0)?sz_bite:sdm_sz, pPS, pWS, pRS, pRS2, val, dir);
/*
#pragma omp barrier

  #pragma omp master
  {
   for(i=0; i<=(SRS_NMB-1); i++)
   {
	delete []PP[i];
	delete []WW[i];
	delete []UU[i];
	delete []VV[i];
	delete []TT[i];
   }
  }
*/
  ///////////////// m - C Y C L E  ////////////////////
  for(n=M_prm; (n*tau)<=TmLimCalc; n+=M_prm)
  {
tick=omp_get_wtime();//1) підготовка до m-циклу

   SetActualUV(dom*(sdm_sz-1), (dir>0)?sz_bite:sdm_sz, pU, pV, dir);
   SetActualT(dom*(sdm_sz-1), (dir>0)?sz_bite:sdm_sz, pF, dir); // temporary pF == pT
   switch(val)
   {
    case 0:	memcpy(pR1, pU, sz_bite*sizeof(double));	break;
    case 1:	memcpy(pR1, pV, sz_bite*sizeof(double));	break;
    case 2:	memcpy(pR1, pF, sz_bite*sizeof(double));	break;
   }

   SetActualPressure(sdm_sz, (n-M_prm)*tau/CALC_TIME, pPS, pP, dir);
   SetActualDensity(sdm_sz, pP, pF /*pT*/, pD, dir);
   SetActualDiff(sdm_sz, pU, pV, pMu, pP, pD, dir);
   SetActualW(sdm_sz, (n-M_prm)*tau/CALC_TIME, pU, pV, pWS, pW, dir);
   SetFPart(sdm_sz, (n-M_prm)*tau/CALC_TIME, pU, pV, pW, pP, pPS, pD, pF, val, dir);//F1=F2=F3=tau*F
   SetActualAdvection(sdm_sz, pU, pV, pW, pMu, pC, dir);
   SetBoundaryCond(sdm_sz, n*tau/CALC_TIME, pRS, pRS2, pR2, val, dir);

Timer[dir][val][dom][0]+=(omp_get_wtime()-tick); 
tick=omp_get_wtime();//2) m-цикл

/// D E B U G /////////////
//if(val==0)
//{
/*
 if(dom==0) StoreToArray(0, (Subdomains>1)?sdm_sz:(sdm_sz+1), pP, pUU);
 else if(dom==(Subdomains-1)) StoreToArray(1+(sdm_sz-1)*dom, sdm_sz, &pP[(SZ_X3+1)*(SZ_X2+1)], pUU);
 else StoreToArray(1+(sdm_sz-1)*dom, sdm_sz-1, &pP[(SZ_X3+1)*(SZ_X2+1)], pUU);
*/
/*
 char buf[32];
 sprintf(buf, "OUT/pC_%d_%d", dir, dom);
 if(dir) SaveArray(buf, pC, SZ_X3, SZ_X2, sdm_sz);
 else SaveArray(buf, pC, SZ_X3, sdm_sz, SZ_X1, 0);
 */
//}
//////////////////////

  // The first stage of m-cycle
   switch(dir)
   {// lam, order access [j,k,i]
	case 0:		RunCalcLam(sdm_sz, pR1, pR2, pC, pMu, pD, pF, pTask1D);		break;
	// fi, order access [i,j,k]
	case 1:		RunCalcFi(sdm_sz, pR1, pR2, pC, pMu, pD, pF, pTask1D);		break;
	// z, order access [i,j,k]
	case 2:		RunCalcZ(sdm_sz, pR1, pR2, pC, pMu, pD, pF, pTask1D);		break;
   }

  if(M_prm%2)
  {// The solution is always in pR1 !
   pR=pR1;
   pR1=pR2;
   pR2=pR;
  }

Timer[dir][val][dom][1]+=(omp_get_wtime()-tick); 
  // The second stage of m-cycle (average)
#pragma omp barrier
tick=omp_get_wtime(); //3) усереднення

  if(dir==1)// fi
  {
   int sub, j1;
   // Average for each direction
   for(i=1; i<=(sdm_sz-1); i++)// lambda, without boundary
	for(j=dom*(sdm_sz-1); j<=(SZ_X2+dom*(sdm_sz-1)-1); j++)// fi, without boundary
	{
	 j1=j%SZ_X2;
	 if((j1==0)||(j1==(SZ_X2-1))) continue;
	 for(k=1; k<=(SZ_X3-1); k++)// z, without boundary
	  {
	   d=(i*(SZ_X2+1)+j1)*(SZ_X3+1)+k;
	   sub=(j1-1)/(sdm_sz-1);
	   d1=((j1-sub*(sdm_sz-1))*(SZ_X3+1)+k)*(SZ_X1+1)+(i+dom*(sdm_sz-1));	
	   pR1[d]+=(pOut[2][val][dom][d]+pOut[0][val][sub][d1]);
	   pR1[d]/=3.;
	  }
	}
   // Set boundary conditions from dir == 2 
   for(i=1; i<=(sdm_sz-1); i++)// lambda, without boundary
	for(j=1; j<=(SZ_X2-1); j++)// fi, without boundary
	 for(k=0; k<=SZ_X3; k+=SZ_X3)// z, only boundary
	  {
	   d=(i*(SZ_X2+1)+j)*(SZ_X3+1)+k;
	   pR1[d]=pOut[2][val][dom][d];
	  }
   //  Set boundary conditions from dir == 0, dom == 0
   if(dom==0)
	for(j=1; j<=(SZ_X2-1); j++)// fi, without boundary
	{
	 for(k=1; k<=(SZ_X3-1); k++)// z, without boundary
	  {
	   sub=(j-1)/(sdm_sz-1);
	   d1=((j-sub*(sdm_sz-1))*(SZ_X3+1)+k)*(SZ_X1+1)+(0+dom*(sdm_sz-1));
	   pR1[(0*(SZ_X2+1)+j)*(SZ_X3+1)+k]=pOut[0][val][sub][d1];
	  }
	}
   //  Set boundary conditions from dir == 0, dom == Subdomains-1
   if(dom==(Subdomains-1))
	for(j=(SZ_X2-1); j>=1; j--)// fi, without boundary
	{
	 for(k=1; k<=(SZ_X3-1); k++)// z, without boundary
	  {
	   sub=(j-1)/(sdm_sz-1);
	   d1=((j-sub*(sdm_sz-1))*(SZ_X3+1)+k)*(SZ_X1+1)+(sdm_sz+dom*(sdm_sz-1));
	   pR1[(sdm_sz*(SZ_X2+1)+j)*(SZ_X3+1)+k]=pOut[0][val][sub][d1];
	  }
	}
   }

Timer[dir][val][dom][2]+=(omp_get_wtime()-tick); 
#pragma omp barrier
tick=omp_get_wtime();//4) просторове об"єднання результатів 

 // Storing to global arrays 
 if(dir==1)
 {// There is data race; it is not divided for subdomains
  switch(val)
  {
   case 0: pR=pUU; break;
   case 1: pR=pVV; break;
   case 2: pR=pTT; break;
  }
  if(dom==0) StoreToArray(0, (Subdomains>1)?sdm_sz:(sdm_sz+1), pR1, pR);
  else if(dom==(Subdomains-1)) StoreToArray(1+(sdm_sz-1)*dom, sdm_sz, &pR1[(SZ_X3+1)*(SZ_X2+1)], pR);
  else StoreToArray(1+(sdm_sz-1)*dom, sdm_sz-1, &pR1[(SZ_X3+1)*(SZ_X2+1)], pR);
 }

Timer[dir][val][dom][3]+=(omp_get_wtime()-tick);
#pragma omp barrier 
 }
  delete []pR1;
  delete []pR2;
  delete []pF;
  delete []pC;
  delete []pU;
  delete []pV;
  delete []pW;
  delete []pP;
  delete []pD;
  delete []pMu;
  for(i=0; i<=(SRS_NMB-1); i++)
  {
   delete []pPS[i];
   delete []pWS[i];
   delete []pRS[i];
   if(val<2) delete []pRS2[i];
  }
 }

 TotTimer=omp_get_wtime()-TotTimer;
 return NmbThreads;
}



int LoadDataToArrays()
{
  char letter[ALL_MT_VAL-1]={'t', 'p', 'u', 'v', 'w'};
  int sz=(SZ_X1+1)*(SZ_X2+1)*(SZ_X3+1);
   
  for(int r=0; r<=(SRS_NMB-1); r++)
  {
   PP[r]= new double [sz];
   WW[r]= new double [sz];
   UU[r]= new double [sz];
   VV[r]= new double [sz];
   TT[r]= new double [sz];
  }

  omp_set_num_threads(ALL_MT_VAL-1);
  #pragma omp parallel
  {
   char buf[32], pFlName[SZ_X3+1][32];
   int i, j, k, n, proc;
   FILE *fl;

   // init input paths
   for(k=0; k<=SZ_X3; k++) sprintf(pFlName[k], "%s%s%s", FL_IN_DIR, pVrtLvlFlName[k], FL_EXT_IN);	
   proc=omp_get_thread_num();
   for(k=0; k<=SZ_X3; k++) pFlName[k][12]=letter[proc];
	for(n=0; n<=(SRS_NMB-1); n++)// t
	{
	 for(k=0; k<=SZ_X3; k++)// z
	 {
	  if(n>0)// change path
	  {
	   pFlName[k][5]=(char)('1'+n);// folder
	   sprintf(buf, "%.2u", 12*n);
	   pFlName[k][14]=buf[0];
	   pFlName[k][15]=buf[1];
	  }
	  fl=fopen(pFlName[k], "r");
	  for(j=0; j<=SZ_X2; j++)// fi
		for(i=0; i<=SZ_X1; i++)// lambda
	   { 
	    fgets(buf, 32, fl);
	    switch(proc)
	    {
	     case 0:	TT[n][(i*(SZ_X2+1)+j)*(SZ_X3+1)+k]=atof(buf);	break;
	     case 1:	PP[n][(i*(SZ_X2+1)+j)*(SZ_X3+1)+k]=atof(buf);	break;
	     case 2:	UU[n][(i*(SZ_X2+1)+j)*(SZ_X3+1)+k]=atof(buf);	break;
	     case 3:	VV[n][(i*(SZ_X2+1)+j)*(SZ_X3+1)+k]=atof(buf);	break;
	     case 4:	WW[n][(i*(SZ_X2+1)+j)*(SZ_X3+1)+k]=atof(buf);	break;
	    }
	   }
	  fclose(fl);
	 }
	}
  
   switch(proc)
   {
    case 0: memcpy(pTT, TT[0], sz*sizeof(double)); break;
	case 2: memcpy(pUU, UU[0], sz*sizeof(double)); break;
	case 3: memcpy(pVV, VV[0], sz*sizeof(double)); 	break;
	//case 4: memcpy(pWW, WW[0], (SZ_X3+1)*(SZ_X2+1)*(SZ_X1+1)*sizeof(double)); break;
   }
  }
 return 0;
}

/*
int LoadDensity(double r)
{
 char buf[32], pFlName[SZ_X3+1][32];
 int i, j, k, n;
 double t[SRS_NMB];
 FILE *fl[SRS_NMB];

 // init input paths
 for(k=0; k<=SZ_X3; k++)
    sprintf(pFlName[k], "%s%s%s", FL_IN_DIR, pVrtLvlFlName[k], FL_EXT_IN);	

 for(k=0; k<=SZ_X3; k++)// z
 {
  for(n=0; n<=(SRS_NMB-1); n++)// time series
  {
   if(n>0)
   {
    pFlName[k][5]=(char)('1'+n);// folder
    sprintf(buf, "%.2u", 12*n);
    pFlName[k][14]=buf[0];
    pFlName[k][15]=buf[1];
   }
   fl[n]=fopen(pFlName[k], "r");
  }

  for(j=0; j<=SZ_X2; j++)// fi
	for(i=0; i<=SZ_X1; i++)// lambda
   { 
	for(n=0; n<=(SRS_NMB-1); n++)
	{
     fgets(buf, 32, fl[n]);
     t[n]=atof(buf);	
	}
	pD[(i*(SZ_X2+1)+j)*(SZ_X3+1)+k]=QuadrInterpol(r, t[0], t[1], t[2]);
   }
  for(n=0; n<=(SRS_NMB-1); n++) fclose(fl[n]);
 }
 return 0;
}
*/

int SaveDataToFiles()
{
 char letter[ALL_MT_VAL]={'t', 'p', 'd', 'u', 'v', 'w'};
 omp_set_num_threads(ALL_MT_VAL);
 #pragma omp parallel
  {
   char buf[32], pFlName[SZ_X3+1][32];
   int i, j, k, proc;
   FILE *fl;

   // init input paths
   for(k=0; k<=SZ_X3; k++) sprintf(pFlName[k], "%s%s%s", FL_OUT_DIR, pVrtLvlFlName[k], FL_EXT_OUT);	
   proc=omp_get_thread_num();
   for(k=0; k<=SZ_X3; k++)// z
   {
    pFlName[k][9]=letter[proc];
    fl=fopen(pFlName[k], "w");
	for(j=0; j<=SZ_X2; j++)// fi
	  for(i=0; i<=SZ_X1; i++)// lambda
      {
       switch(proc)
	    {
	     case 0:	fprintf(fl, "%.2f\n", pTT[(i*(SZ_X2+1)+j)*(SZ_X3+1)+k]); break;// 
	     //case 1:	fprintf(fl, "%.2f\n", pPP[i][j][k]);		break;// out actual pres.
	     //case 2:	fprintf(fl, "%.4f\n", pD[(i*(SZ_X2+1)+j)*(SZ_X3+1)+k]);		break;// out actual density
	     case 3:	fprintf(fl, "%.2f\n", pUU[(i*(SZ_X2+1)+j)*(SZ_X3+1)+k]);	break;//
	     case 4:	fprintf(fl, "%.2f\n", pVV[(i*(SZ_X2+1)+j)*(SZ_X3+1)+k]);	break;//	
	     //case 5:	fprintf(fl, "%.4f\n", pW[(i*(SZ_X2+1)+j)*(SZ_X3+1)+k]);	break;//
	    }
      }
    fclose(fl);
  }
 }
 return 0;
}


void SaveArray(char *name, double *pAr, int szZ, int szF, int szL, int ijk)
{// default: ijk == 1  
 int i, j, k=5, m;
 FILE *fl;
 
 fl=fopen(name, "w");
 if(ijk)
  for(k=0; k<=szZ; k++)// z
   for(j=0; j<=szF; j++)// fi
	for(i=0; i<=szL; i++)// lambda
	 fprintf(fl, "%.6e\n", pAr[(i*(szF+1)+j)*(szZ+1)+k]);
 else //jki
  for(k=0; k<=szZ; k++)// z
   for(j=0; j<=szF; j++)// fi
    for(i=0; i<=szL; i++)// lambda
	 fprintf(fl, "%.6e\n", pAr[(j*(szZ+1)+k)*(szL+1)+i]);
 fclose(fl);
}


int SetActualPressure(int sz, double r, double **pPS, double *pPrs, int dir)
{
 int k, j, i, d;
   
 if(dir>0)  // lambda-domain decomposition
  for(i=0; i<=sz; i++)// lambda, overlapped domains
   for(j=0; j<=SZ_X2; j++)// fi
	 for(k=0; k<=SZ_X3; k++)// z
	 {
	  d=(i*(SZ_X2+1)+j)*(SZ_X3+1)+k;
	  pPrs[d]=QuadrInterpol(r, pPS[0][d], pPS[1][d], pPS[2][d]); 
	 }
 else
  for(j=0; j<=sz; j++)// fi
   for(k=0; k<=SZ_X3; k++)// z 
	 for(i=0; i<=SZ_X1; i++)// lam
	 {
	  d=(j*(SZ_X3+1)+k)*(SZ_X1+1)+i;
	  pPrs[d]=QuadrInterpol(r, pPS[0][d], pPS[1][d], pPS[2][d]); 
	  //pPrs[d]=pPS[0][d]-r*(0.5*pPS[2][d]-2.*pPS[1][d]+1.5*pPS[0][d]-r*(0.5*pPS[2][d]-pPS[1][d]+0.5*pPS[0][d]));
	 }

 return 0;
}


int RetArgMod(int d, double *parg, double *pmod)
{
 for(int s=0; s<=(SRS_NMB-1); s++)
 {
  pmod[s]=sqrt(UU[s][d]*UU[s][d]+VV[s][d]*VV[s][d]);
  parg[s]=atan2(VV[s][d], UU[s][d]);
  if(parg[s]==0.)
  {
   if((VV[s][d]>0.)&&(UU[s][d]==0.)) parg[s]=M_PI/2.;
   else if((VV[s][d]<0.)&&(UU[s][d]==0.)) parg[s]=-M_PI/2.;
   else if((VV[s][d]==0.)&&(UU[s][d]<0.)) parg[s]=M_PI;
  }
 }
 return 0;
}


inline double QuadrInterpol(double r, double f0, double f1, double f2)
{// r is from 0. to 2.
 return f0-r*(0.5*f2-2.*f1+1.5*f0-r*(0.5*f2-f1+0.5*f0));
}


int Receive(int argc, char* argv[])
{
 /*if(argc<5)
 {
  printf("\nToo few parameters in command line\n");
  return 1;
 }
 else if(argc>5)
 {
  printf("\nToo many parameters in command line\n");
  return 2;
 }
 Subdomains=atoi(argv[1]);
 M_prm=atoi(argv[2]);
 TmKnotsPer12h=atoi(argv[3]);
 TmLimCalc=atoi(argv[4]);

 tau=(double)CALC_TIME/TmKnotsPer12h; // time step*/

 Subdomains=2;
 M_prm=10;
 TmKnotsPer12h=10;
 TmLimCalc=100;

 tau=(double)CALC_TIME/TmKnotsPer12h; // time step*/

 return 0;
}


// distance between data (SZ_X1=SZ_X2=181, SZ_X3=10)
//					du			dv			dT				
//
// 0-12h	(L2)	32900		31000		1250   
//			(MAX)	55.3		55.4		13.8
//
// 0-24h	(L2)	35400		32000		2170
//			(MAX)	60.3		58.2		17.6
//
// 12-24h	(L2)	34400		31400		1090		
//			(MAX)	57.8		56.7		12.1		
//

int CompareRes(double r, double TotTimer, int Proc_Num) // Comparing interpol. and calc. values in L2 norm
{// r = 0 ,..., 2 
 int i, j, k, m;
 double perr[3]={0.}, perMax[3]={0.}, z;
 FILE *fl;
 
 SetUVatTime(r, pD, pTemp);
 perr[0]=MESH_ANG_STP*sqrt(MESH_VRT_STP*CompareArrays(pUU, pD, perMax[0]));
 perr[1]=MESH_ANG_STP*sqrt(MESH_VRT_STP*CompareArrays(pVV, pTemp, perMax[1]));
 
 SetTWatTime(r, pTemp, NULL);
 //SetActualW(1, SZ_X1-1);
 perr[2]=MESH_ANG_STP*sqrt(MESH_VRT_STP*CompareArrays(pTT, pTemp, perMax[2]));
 

 // Timer[dir][val][dom][stg]
 // Time for the first stage (avrg)
 z=0.;
 for(i=0; i<=(Subdomains-1); i++)// dom
  for(j=0; j<=2; j++) // val
   for(k=0; k<=2; k++) // dir
   {
    if(Timer[0][0][0][0]>Timer[k][j][i][0]) Timer[0][0][0][0]=Timer[k][j][i][0];
    if(z<Timer[k][j][i][0]) z=Timer[k][j][i][0];
   }
 Timer[0][0][0][0]=(z+Timer[0][0][0][0])/2.;
 // Time for the second stage (min)
 Timer[0][0][0][1]=1000000.;
 for(i=0; i<=(Subdomains-1); i++)// dom
  for(j=0; j<=2; j++) // val
   for(k=0; k<=2; k++) // dir
   {
    if(Timer[0][0][0][1]>Timer[k][j][i][1]) Timer[0][0][0][1]=Timer[k][j][i][1];
   }
 // Time for the thrird stage (min)
 Timer[0][0][0][2]=1000000.;
 for(i=0; i<=(Subdomains-1); i++)// dom
   for(j=0; j<=2; j++) // val
   {
    if(Timer[0][0][0][2]>Timer[1][j][i][2]) Timer[0][0][0][2]=Timer[1][j][i][2];
   }
 // Time for the fourth stage (max)
 Timer[0][0][0][3]=0.;
 for(i=0; i<=(Subdomains-1); i++)// dom
   for(j=0; j<=2; j++) // val
   {
    if(Timer[0][0][0][3]<Timer[1][j][i][3]) Timer[0][0][0][3]=Timer[1][j][i][3];
   }
 

 fl=fopen("res.txt", "a");
 fprintf(fl, "\n>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\n");
 //fprintf(fl, "Total time of paral. part is %.2f sec (Thread number was %d)\n", TotTimer, Proc_Num);
 fprintf(fl, "THREAD NUMBER: %d\t\t\tm=%d\tS=%d\n", Proc_Num, M_prm, Subdomains);
 fprintf(fl, "TIME: (%.2f+%.2f+%.2f+%.2f)\t\t%.2f sec\n", 
	  Timer[0][0][0][0], Timer[0][0][0][1], Timer[0][0][0][2], Timer[0][0][0][3],
	  Timer[0][0][0][0]+Timer[0][0][0][1]+Timer[0][0][0][2]+Timer[0][0][0][3]);
 fprintf(fl, "dt=%.4f\tSZ=%d*%d*%d\tTimeCalc=%d\nDiss=%.0f\n",
		tau, SZ_X1+1, SZ_X2+1, SZ_X3+1, TmLimCalc, COE_DISS);
 fprintf(fl, "ERROR\tL2\t\tMAX\nU\t%.2e\t%.2e\nV\t%.2e\t%.2e\nT\t%.2e\t%.2e\n",
		perr[0], perMax[0], perr[1], perMax[1], perr[2], perMax[2]);
 fclose(fl);

 printf("S=%d, m=%d: Time=%.2f\n", Subdomains, M_prm, 
		Timer[0][0][0][0]+Timer[0][0][0][1]+Timer[0][0][0][2]+Timer[0][0][0][3]);
 printf("ERROR\tL2\t\tMAX\nU\t%.2e\t%.2e\nV\t%.2e\t%.2e\nT\t%.2e\t%.2e\n",
		perr[0], perMax[0], perr[1], perMax[1], perr[2], perMax[2]);
 for(i=0; i<=(SRS_NMB-1); i++)
 {
  delete []PP[i];
  delete []WW[i];
  delete []UU[i];
  delete []VV[i];
  delete []TT[i];
 }

 return 0;  
}


int SetUVatTime(double r, double *pu, double *pv)
{
 int k, j, i, d;
 double a=0., m=1., pa[3], pm[3];
 
 for(i=1; i<=(SZ_X1-1); i++)// lambda, without boundary
   for(j=1; j<=(SZ_X2-1); j++)// fi, without boundary
     for(k=1; k<=(SZ_X3-1); k++)// z, without boundary
	 {
	  d=(i*(SZ_X2+1)+j)*(SZ_X3+1)+k;
	  RetArgMod(d, pa, pm);
	  a=QuadrInterpol(r, pa[0], pa[1], pa[2]);
	  m=QuadrInterpol(r, pm[0], pm[1], pm[2]);	
      pu[d]=m*cos(a);
	  pv[d]=m*sin(a);
	 }
 return 0;  
}


int SetTWatTime(double r, double *pt, double *pw)
{
 int k, j, i, d;
 
 for(i=1; i<=(SZ_X1-1); i++)// lambda, without boundary
   for(j=1; j<=(SZ_X2-1); j++)// fi, without boundary
	 for(k=1; k<=(SZ_X3-1); k++)// z, without boundary    
	 {
	  d=(i*(SZ_X2+1)+j)*(SZ_X3+1)+k;
	  pt[d]=QuadrInterpol(r, TT[0][d], TT[1][d], TT[2][d]);
	  //pw[d]=QuadrInterpol(r, WW[0][i][j][k], WW[1][i][j][k], WW[2][i][j][k]);
	 }
 return 0;  
}


double CompareArrays(double *p1, double *p2, double &Mx)
{// for L2-norm
 int i, j, k;
 double s=0., max=0., z;
 
 for(i=1; i<=(SZ_X1-1); i++)// lambda, without boundary
  for(j=1; j<=(SZ_X2-1); j++)// fi, without boundary
    for(k=1; k<=(SZ_X3-1); k++)// z, without boundary
    {
	 z=p1[(i*(SZ_X2+1)+j)*(SZ_X3+1)+k]-p2[(i*(SZ_X2+1)+j)*(SZ_X3+1)+k];
	 z*=z;
	 s+=z*z;
	 if(z>max) max=z;
	}
 Mx=sqrt(max);
 return s;
}

/*
void InitD(int str, int sz, double *pA, int dir)
{// dir: 0 - lam, 1 - fi, 2 - z
 if(dir>0) memcpy(pA, &pD[str*(SZ_X2+1)*(SZ_X3+1)], sz*sizeof(double)); // lambda-domain decomposition
 else // fi-domain decomposition
 {
  int i, j, k;
  for(j=str; j<=(str+sz); j++)// fi
   for(k=0; k<=SZ_X3; k++)// z
	for(i=0; i<=SZ_X1; i++)// lambda
	 pA[((j-str)*(SZ_X3+1)+k)*(SZ_X1+1)+i]=pD[(i*(SZ_X2+1)+j)*(SZ_X3+1)+k];	 
 }
}
*/

void SetBoundaryCond(int sz, double r, double **pRS, double **pRS2, double *pR, int val, int dir)
{
 int i, j, k, d;
 double a, m;

 switch(val)
 {
  case 0: // u
	  switch(dir)
	  {
	   case 0: // lam, [j k i]
	    for(j=0; j<=sz; j++)// fi, without boundary
		 for(k=0; k<=SZ_X3; k++)// z, 
		  for(i=0; i<=SZ_X1; i+=SZ_X1)// lambda, only boundary
		  {
		   d=(j*(SZ_X3+1)+k)*(SZ_X1+1)+i;
		   a=QuadrInterpol(r, pRS[0][d], pRS[1][d], pRS[2][d]);
		   m=QuadrInterpol(r, pRS2[0][d], pRS2[1][d], pRS2[2][d]);	
	       pR[d]=m*cos(a);
		  }
	   break;

	   case 1: // fi, [i j k]
	    for(i=0; i<=sz; i++)// lambda, without boundary
	     for(j=0; j<=SZ_X2; j+=SZ_X2)// fi, only boundary
		  for(k=0; k<=SZ_X3; k++)// z, 
		  {
		   d=(i*(SZ_X2+1)+j)*(SZ_X3+1)+k;
		   a=QuadrInterpol(r, pRS[0][d], pRS[1][d], pRS[2][d]);
		   m=QuadrInterpol(r, pRS2[0][d], pRS2[1][d], pRS2[2][d]);
	       pR[d]=m*cos(a);
		  }
	   break;

       case 2: // z, [i j k]
	    for(i=0; i<=sz; i++)// lambda,
	     for(j=0; j<=SZ_X2; j++)// fi, 
		  for(k=0; k<=SZ_X3; k+=SZ_X3)// z, only boundary
		  {
		   d=(i*(SZ_X2+1)+j)*(SZ_X3+1)+k;
		   a=QuadrInterpol(r, pRS[0][d], pRS[1][d], pRS[2][d]);
		   m=QuadrInterpol(r, pRS2[0][d], pRS2[1][d], pRS2[2][d]);
	       pR[d]=m*cos(a);
		  }
	   break;
	  }
  break;

  case 1: // v
	  switch(dir)
	  {
	   case 0: // lam, [j k i]
	    for(j=0; j<=sz; j++)// fi, without boundary
		 for(k=0; k<=SZ_X3; k++)// z, 
		  for(i=0; i<=SZ_X1; i+=SZ_X1)// lambda, only boundary
		  {
		   d=(j*(SZ_X3+1)+k)*(SZ_X1+1)+i;
		   a=QuadrInterpol(r, pRS[0][d], pRS[1][d], pRS[2][d]);
		   m=QuadrInterpol(r, pRS2[0][d], pRS2[1][d], pRS2[2][d]);	
	       pR[d]=m*sin(a);
		  }
	   break;

	   case 1: // fi, [i j k]
	    for(i=0; i<=sz; i++)// lambda, without boundary
	     for(j=0; j<=SZ_X2; j+=SZ_X2)// fi, only boundary
		  for(k=0; k<=SZ_X3; k++)// z, 
		  {
		   d=(i*(SZ_X2+1)+j)*(SZ_X3+1)+k;
		   a=QuadrInterpol(r, pRS[0][d], pRS[1][d], pRS[2][d]);
		   m=QuadrInterpol(r, pRS2[0][d], pRS2[1][d], pRS2[2][d]);
	       pR[d]=m*sin(a);
		  }
	   break;

	   case 2: // z, [i j k]
	    for(i=0; i<=sz; i++)// lambda,
	     for(j=0; j<=SZ_X2; j++)// fi, 
		  for(k=0; k<=SZ_X3; k+=SZ_X3)// z, only boundary
		  {
		   d=(i*(SZ_X2+1)+j)*(SZ_X3+1)+k;
		   a=QuadrInterpol(r, pRS[0][d], pRS[1][d], pRS[2][d]);
		   m=QuadrInterpol(r, pRS2[0][d], pRS2[1][d], pRS2[2][d]);
	       pR[d]=m*sin(a);
		  }
	   break;
	  }
  break;

  case 2: // T
	  switch(dir)
	  {
	   case 0: // lam, [j k i]
	    for(j=0; j<=sz; j++)// fi, without boundary
		 for(k=0; k<=SZ_X3; k++)// z, 
		  for(i=0; i<=SZ_X1; i+=SZ_X1)// lambda, only boundary
		  {
		   d=(j*(SZ_X3+1)+k)*(SZ_X1+1)+i;
	       pR[d]=QuadrInterpol(r, pRS[0][d], pRS[1][d], pRS[2][d]);
		  }
	   break;

	   case 1: // fi, [i j k]
	    for(i=0; i<=sz; i++)// lambda, without boundary
	     for(j=0; j<=SZ_X2; j+=SZ_X2)// fi, only boundary
		  for(k=0; k<=SZ_X3; k++)// z,
		  {
		   d=(i*(SZ_X2+1)+j)*(SZ_X3+1)+k;
	       pR[d]=QuadrInterpol(r, pRS[0][d], pRS[1][d], pRS[2][d]);
		  }
	   break;

	   case 2: // z, [i j k]
	    for(i=0; i<=sz; i++)// lambda, without boundary
	     for(j=0; j<=SZ_X2; j++)// fi, only boundary
		  for(k=0; k<=SZ_X3; k+=SZ_X3)// z, 
		  {
	       d=(i*(SZ_X2+1)+j)*(SZ_X3+1)+k;
	       pR[d]=QuadrInterpol(r, pRS[0][d], pRS[1][d], pRS[2][d]);
		  }
	   break;
	  }
  break;
 }
}


void StoreToArray(int str, int sz, double *pSrc, double *pDst)
{
 memcpy(&pDst[str*(SZ_X3+1)*(SZ_X2+1)], pSrc, sz*(SZ_X3+1)*(SZ_X2+1)*sizeof(double));
}


void SetActualUV(int str, int sz, double *pU, double *pV, int dir)
{// dir: 0 - lam, 1 - fi, 2 - z
 if(dir>0)
 {
  memcpy(pU, &pUU[str*(SZ_X2+1)*(SZ_X3+1)], sz*sizeof(double)); // lambda-domain decomposition
  memcpy(pV, &pVV[str*(SZ_X2+1)*(SZ_X3+1)], sz*sizeof(double)); // lambda-domain decomposition
 }
 else // fi-domain decomposition
 {
  int i, j, k;
  for(j=str; j<=(str+sz); j++)// fi
   for(k=0; k<=SZ_X3; k++)// z
	for(i=0; i<=SZ_X1; i++)// lambda
	{
	 pU[((j-str)*(SZ_X3+1)+k)*(SZ_X1+1)+i]=pUU[(i*(SZ_X2+1)+j)*(SZ_X3+1)+k];
	 pV[((j-str)*(SZ_X3+1)+k)*(SZ_X1+1)+i]=pVV[(i*(SZ_X2+1)+j)*(SZ_X3+1)+k];
	}
 }
}


void SetActualT(int str, int sz, double *pT, int dir)
{// dir: 0 - lam, 1 - fi, 2 - z
 if(dir>0) memcpy(pT, &pTT[str*(SZ_X2+1)*(SZ_X3+1)], sz*sizeof(double)); // lambda-domain decomposition
 else // fi-domain decomposition
 {
  int i, j, k;
  for(j=str; j<=(str+sz); j++)// fi
   for(k=0; k<=SZ_X3; k++)// z
	for(i=0; i<=SZ_X1; i++)// lambda
	 pT[((j-str)*(SZ_X3+1)+k)*(SZ_X1+1)+i]=pTT[(i*(SZ_X2+1)+j)*(SZ_X3+1)+k];
 }
}


int SetActualW(int sz, double r, double *pU, double *pV, double **pWS, double *pW, int dir)
{
 int k, j, i, d;
 double a, b, L0, L1;

 if(dir>0) // lambda-domain decomposition
 {
  for(i=1; i<=(sz-1); i++)// lambda, without boundary
   for(j=1; j<=(SZ_X2-1); j++)// fi, without boundary 
   {
    d=(i*(SZ_X2+1)+j)*(SZ_X3+1)+0;
	pW[d]=QuadrInterpol(r, pWS[0][d], pWS[1][d], pWS[2][d]);
   }
  
  for(i=1; i<=(sz-1); i++)// lambda, without boundary 
   for(j=1; j<=(SZ_X2-1); j++)// fi, without boundary
   {
    a=cos(j*MESH_ANG_STP);
    b=tan(j*MESH_ANG_STP);

    L0=(pU[((i+1)*(SZ_X2+1)+j)*(SZ_X3+1)+0]-pU[((i-1)*(SZ_X2+1)+j)*(SZ_X3+1)+0])/(2.*a*MESH_ANG_STP)
		  +(pV[(i*(SZ_X2+1)+j+1)*(SZ_X3+1)+0]-pV[(i*(SZ_X2+1)+j-1)*(SZ_X3+1)+0])/(2.*MESH_ANG_STP)
		  -pV[(i*(SZ_X2+1)+j)*(SZ_X3+1)+0]*b;
    for(k=1; k<=SZ_X3; k++)// z
    {
     d=(i*(SZ_X2+1)+j)*(SZ_X3+1)+k;
    
	 L1=(pU[((i+1)*(SZ_X2+1)+j)*(SZ_X3+1)+k]-pU[((i-1)*(SZ_X2+1)+j)*(SZ_X3+1)+k])/(2.*a*MESH_ANG_STP)
		  +(pV[(i*(SZ_X2+1)+j+1)*(SZ_X3+1)+k]-pV[(i*(SZ_X2+1)+j-1)*(SZ_X3+1)+k])/(2.*MESH_ANG_STP)
		  -pV[d]*b;

	 pW[d]=pW[d-1]-MESH_VRT_STP*(L0+L1)/(2.*ERH_RDS);
	 L0=L1;
    }
   }
 }
 else // fi-domain decomposition
 {
  for(j=1; j<=(sz-1); j++)// fi, without boundary
   for(i=1; i<=(SZ_X1-1); i++)// lan, without boundary
   {
    d=(j*(SZ_X3+1)+0)*(SZ_X1+1)+i;
	pW[d]=QuadrInterpol(r, pWS[0][d], pWS[1][d], pWS[2][d]);
   }

  for(j=1; j<=(sz-1); j++)// fi, without boundary
  {
   a=cos(j*MESH_ANG_STP);
   b=tan(j*MESH_ANG_STP);
   for(i=1; i<=(SZ_X1-1); i++)// lambda, without boundary
   {
    L0=(pU[(j*(SZ_X3+1)+0)*(SZ_X1+1)+i+1]-pU[(j*(SZ_X3+1)+0)*(SZ_X1+1)+i-1])/(2.*a*MESH_ANG_STP)
		  +(pV[((j+1)*(SZ_X3+1)+0)*(SZ_X1+1)+i]-pV[((j-1)*(SZ_X3+1)+0)*(SZ_X1+1)+i])/(2.*MESH_ANG_STP)
		  -pV[(j*(SZ_X3+1)+0)*(SZ_X1+1)+i]*b;

    for(k=1; k<=SZ_X3; k++)// z
    {
     d=(j*(SZ_X3+1)+k)*(SZ_X1+1)+i;
	 L1=(pU[d+1]-pU[d-1])/(2.*a*MESH_ANG_STP)
		+(pV[((j+1)*(SZ_X3+1)+k)*(SZ_X1+1)+i]-pV[((j-1)*(SZ_X3+1)+k)*(SZ_X1+1)+i])/(2.*MESH_ANG_STP)
		-pV[d]*b;

	 pW[d]=pW[(j*(SZ_X3+1)+k-1)*(SZ_X1+1)+i]-MESH_VRT_STP*(L0+L1)/(2.*ERH_RDS);
	 L0=L1;
    }
   }
  }
 }
 return 0;
}


void InitLocalArrays(int str, int sz, double **pPS, double **pWS, double **pRS, double **pRS2, int val, int dir)
{// dir: 0 - lam, 1 - fi, 2 - z
 int i, j, k, d;

 if(dir>0)
 {
  for(d=0; d<=(SRS_NMB-1); d++)
  {
   memcpy(pPS[d], &PP[d][str*(SZ_X2+1)*(SZ_X3+1)], sz*sizeof(double)); // lambda-domain decomposition
   memcpy(pWS[d], &WW[d][str*(SZ_X2+1)*(SZ_X3+1)], sz*sizeof(double)); // lambda-domain decomposition
  }
 }
 else // fi-domain decomposition
 {
  for(d=0; d<=(SRS_NMB-1); d++)
   for(j=str; j<=(str+sz); j++)// fi
    for(k=0; k<=SZ_X3; k++)// z
	 for(i=0; i<=SZ_X1; i++)// lambda
	 {
	  pPS[d][((j-str)*(SZ_X3+1)+k)*(SZ_X1+1)+i]=PP[d][(i*(SZ_X2+1)+j)*(SZ_X3+1)+k];
	  pWS[d][((j-str)*(SZ_X3+1)+k)*(SZ_X1+1)+i]=WW[d][(i*(SZ_X2+1)+j)*(SZ_X3+1)+k];
	 }
  }
 

 if(val==2)// T
 {
  if(dir>0)
  {
   for(d=0; d<=(SRS_NMB-1); d++)
     memcpy(pRS[d], &TT[d][str*(SZ_X2+1)*(SZ_X3+1)], sz*sizeof(double)); // lambda-domain decomposition
  }
  else
  {
   for(d=0; d<=(SRS_NMB-1); d++)
    for(j=str; j<=(str+sz); j++)// fi
     for(k=0; k<=SZ_X3; k++)// z
	  for(i=0; i<=SZ_X1; i++)// lambda
	  {
	   pRS[d][((j-str)*(SZ_X3+1)+k)*(SZ_X1+1)+i]=TT[d][(i*(SZ_X2+1)+j)*(SZ_X3+1)+k];
	  }
  }
 }
 else// U, V
 {
  double mod[SRS_NMB], arg[SRS_NMB];
  if(dir>0)
  {
   sz/=((SZ_X3+1)*(SZ_X2+1));
   sz--;
   for(i=str; i<=(str+sz); i++)//  lambda
    for(j=0; j<=SZ_X2; j++)// fi
     for(k=0; k<=SZ_X3; k++)// z
	 {
	  RetArgMod((i*(SZ_X2+1)+j)*(SZ_X3+1)+k, arg, mod);
	  for(d=0; d<=(SRS_NMB-1); d++)
	  {
	   pRS[d][((i-str)*(SZ_X2+1)+j)*(SZ_X3+1)+k]=arg[d];
	   pRS2[d][((i-str)*(SZ_X2+1)+j)*(SZ_X3+1)+k]=mod[d];
	  }
	 }
  }
  else
  {
   for(j=str; j<=(str+sz); j++)// fi
    for(k=0; k<=SZ_X3; k++)// z
	 for(i=0; i<=SZ_X1; i++)// lambda
	 {
	  RetArgMod((i*(SZ_X2+1)+j)*(SZ_X3+1)+k, arg, mod);
	  for(d=0; d<=(SRS_NMB-1); d++)
	  {
	   pRS[d][((j-str)*(SZ_X3+1)+k)*(SZ_X1+1)+i]=arg[d];
	   pRS2[d][((j-str)*(SZ_X3+1)+k)*(SZ_X1+1)+i]=mod[d];
	  }
	 }
  }
 }
}







// -------------------------------------------------------------------------------------------------
// CALC


int Task1D_Minus(double dt, double dx, int beg, int end, double *p1D)
{// beg - the first point to calc., end -- the last one
 int i;
 double a, b, h=2.*dx*dx;

 for(i=end; i>=beg; i--)
 {
  a=dt*(p1D[T1D_C+i]*p1D[T1D_D+i]*dx+(p1D[T1D_MU+i]+p1D[T1D_MU+i-1]))/(h*p1D[T1D_D+i]);
  b=dt*(p1D[T1D_C+i]*p1D[T1D_D+i]*dx-(p1D[T1D_MU+i]+p1D[T1D_MU+i+1]))/(h*p1D[T1D_D+i]);
  p1D[T1D_NEW+i]=(p1D[T1D_F+i]+p1D[T1D_OLD+i]+
	  a*(p1D[T1D_OLD+i-1]-p1D[T1D_OLD+i])-b*p1D[T1D_NEW+i+1])/(1.-b);
 }
 return 0;
}


int Task1D_Plus(double dt, double dx, int beg, int end, double *p1D)
{// beg - the first point to calc., end -- the last one
 int i;
 double a, b, h=2.*dx*dx;

 for(i=beg; i<=end; i++)
 {
  a=dt*(p1D[T1D_C+i]*p1D[T1D_D+i]*dx-(p1D[T1D_MU+i+1]+p1D[T1D_MU+i]))/(h*p1D[T1D_D+i]);
  b=dt*(p1D[T1D_C+i]*p1D[T1D_D+i]*dx+(p1D[T1D_MU+i-1]+p1D[T1D_MU+i]))/(h*p1D[T1D_D+i]);
  p1D[T1D_NEW+i]=(p1D[T1D_F+i]+p1D[T1D_OLD+i]+
	  a*(p1D[T1D_OLD+i]-p1D[T1D_OLD+i+1])+b*p1D[T1D_NEW+i-1])/(1.+b);
 }
 return 0;
}


int Calc1DTask(double dt, double dx, double *p1D, int sz)// solving 1D task
{
 int i, beg=1;
 double a_p, a_m, b_p, b_m, d, h=2.*dx*dx,
		zn;

 for(i=2; i<=(sz-1); i++)// The main part of domain
 {
  zn=p1D[T1D_C+i-1]*p1D[T1D_C+i];
  if(zn<0.)
  {
   if(p1D[T1D_C+i-1]<0.)// change from "-" (i-1) to "+" (i)
   {// calc.  p1D[T1D_NEW+i-1] and p1D[T1D_NEW+i]  
    a_p=dt*(p1D[T1D_C+i]*p1D[T1D_D+i]*dx-(p1D[T1D_MU+i+1]+p1D[T1D_MU+i]))/(h*p1D[T1D_D+i]);
	a_m=dt*(p1D[T1D_C+i-1]*p1D[T1D_D+i-1]*dx+(p1D[T1D_MU+i-1]+p1D[T1D_MU+i-2]))/(h*p1D[T1D_D+i-1]);
	b_p=dt*(p1D[T1D_C+i]*p1D[T1D_D+i]*dx+(p1D[T1D_MU+i-1]+p1D[T1D_MU+i]))/(h*p1D[T1D_D+i]);
	b_m=dt*(p1D[T1D_C+i-1]*p1D[T1D_D+i-1]*dx-(p1D[T1D_MU+i-1]+p1D[T1D_MU+i]))/(h*p1D[T1D_D+i-1]);
	d=1.+b_p-b_m;
	a_m=p1D[T1D_F+i-1]+a_m*p1D[T1D_OLD+i-2]+(1.-a_m)*p1D[T1D_OLD+i-1];
	a_p=p1D[T1D_F+i]-a_p*p1D[T1D_OLD+i+1]+(1.+a_p)*p1D[T1D_OLD+i];
	p1D[T1D_NEW+i-1]=((1.+b_p)*a_m-b_m*a_p)/d;
	p1D[T1D_NEW+i]=(b_p*a_m+(1.-b_m)*a_p)/d;
	if((i-1)>beg) Task1D_Minus(dt, dx, beg, i-2, p1D);
	beg=++i;
   }
   else if(p1D[T1D_C+i-1]>0.)
   {
	Task1D_Plus(dt, dx, beg, i-1, p1D);
    beg=i;	
   }
  }
  else if(zn==0.)// p1D[T1D_C+i] == 0
  {
   if(p1D[T1D_C+i-1]<0.)  p1D[T1D_C+i]=-1.e-9;
   else p1D[T1D_C+i]=1.e-9;
  }
 }
 // The end-part of domain
 if(beg<=(sz-1))
 {
  if((p1D[T1D_C+sz-1]+p1D[T1D_C+beg])<0.) Task1D_Minus(dt, dx, beg, sz-1, p1D);
  else Task1D_Plus(dt, dx, beg, sz-1, p1D);
 }
 return 0;
}


int RunCalcLam(int sz, double *p1, double *p2, double *pC, double *pMu, 
			   double *pD, double *pF, double *pT1D) // fi-domain decomposition
{
 int i, j, k, m, d;
 double *pOld=p1, *pNew=p2, *pTemp;
 
 for(m=1; m<=M_prm; m++)// m-cycle
 {
  for(j=1; j<=(sz-1); j++)// fi, without boundary
   for(k=1; k<=(SZ_X3-1); k++)// z, without boundary
   {
	for(i=0; i<=SZ_X1; i++)// lambda, without boundary
    {
	 d=(j*(SZ_X3+1)+k)*(SZ_X1+1)+i;
	 pT1D[T1D_OLD+i]=pOld[d];
	 pT1D[T1D_C+i]=pC[d];
	 pT1D[T1D_MU+i]=pMu[d];
	 pT1D[T1D_D+i]=pD[d];
	 pT1D[T1D_F+i]=pF[d];
	}
	pT1D[T1D_NEW+0]=p2[(j*(SZ_X3+1)+k)*(SZ_X1+1)+0]; // bound. cond. do not change
	pT1D[T1D_NEW+SZ_X1]=p2[(j*(SZ_X3+1)+k)*(SZ_X1+1)+SZ_X1]; // bound. cond. do not change
	Calc1DTask(3.*tau, MESH_ANG_STP, pT1D, SZ_X1); // 3*tau
	for(i=0; i<=SZ_X1; i++) pNew[(j*(SZ_X3+1)+k)*(SZ_X1+1)+i]=pT1D[T1D_NEW+i];
   }
  pTemp=pOld;
  pOld=pNew;
  pNew=pTemp;
 }
 return 0;
}


int RunCalcFi(int sz, double *p1, double *p2, double *pC, double *pMu, 
			  double *pD, double *pF, double *pT1D) // lambda-domain decomposition
{
 int i, j, k, m, d;
 double *pOld=p1, *pNew=p2, *pTemp;
 
 for(m=1; m<=M_prm; m++)// m-cycle
 {
  for(i=1; i<=(sz-1); i++)// lambda, without boundary
   for(k=1; k<=(SZ_X3-1); k++)// z, without boundary
   {/// !!!
	for(j=0; j<=SZ_X2; j++)// fi, without boundary
	{
	 d=(i*(SZ_X2+1)+j)*(SZ_X3+1)+k;
	 pT1D[T1D_OLD+j]=pOld[d];
	 pT1D[T1D_C+j]=pC[d];
	 pT1D[T1D_MU+j]=pMu[d];
	 pT1D[T1D_D+j]=pD[d];
	 pT1D[T1D_F+j]=pF[d];
	}
	pT1D[T1D_NEW+0]=p2[(i*(SZ_X2+1)+0)*(SZ_X3+1)+k]; // bound. cond. do not change
	pT1D[T1D_NEW+SZ_X2]=p2[(i*(SZ_X2+1)+SZ_X2)*(SZ_X3+1)+k]; // bound. cond. do not change
	Calc1DTask(3.*tau, MESH_ANG_STP, pT1D, SZ_X2); // 3*tau
	for(j=0; j<=SZ_X2; j++) pNew[(i*(SZ_X2+1)+j)*(SZ_X3+1)+k]=pT1D[T1D_NEW+j];
  }
  pTemp=pOld;
  pOld=pNew;
  pNew=pTemp;
 }
 return 0;
}


int RunCalcZ(int sz, double *p1, double *p2, double *pC, double *pMu, 
			 double *pD, double *pF, double *pT1D) // lambda-domain decomposition
{
 int i, j, k, m, d;
 double *pOld=p1, *pNew=p2, *pTemp;
 
 for(m=1; m<=M_prm; m++)// m-cycle
 {
  for(i=1; i<=(sz-1); i++)// lambda, without boundary
   for(j=1; j<=(SZ_X2-1); j++)// fi, without boundary 
   {
	for(k=0; k<=SZ_X3; k++)// z, without boundary
    {
	 d=(i*(SZ_X2+1)+j)*(SZ_X3+1)+k;
	 pT1D[T1D_OLD+k]=pOld[d];
	 pT1D[T1D_C+k]=pC[d];
	 pT1D[T1D_MU+k]=pMu[d];
	 pT1D[T1D_D+k]=pD[d];
	 pT1D[T1D_F+k]=pF[d];
	}
	pT1D[T1D_NEW+0]=p2[(i*(SZ_X2+1)+j)*(SZ_X3+1)+0]; // bound. cond. do not change
	pT1D[T1D_NEW+SZ_X3]=p2[(i*(SZ_X2+1)+j)*(SZ_X3+1)+SZ_X3]; // bound. cond. do not change
	Calc1DTask(3.*tau, MESH_VRT_STP, pT1D, SZ_X3); // 3*tau
	for(k=0; k<=SZ_X3; k++) pNew[(i*(SZ_X2+1)+j)*(SZ_X3+1)+k]=pT1D[T1D_NEW+k];	
   }
  pTemp=pOld;
  pOld=pNew;
  pNew=pTemp;
 }
return 0;
}


void SetActualAdvection(int sz, double *pU, double *pV, double *pW, double *pMu, double *pC, int dir)
{
 int i, j, k, d;
 double a;
 
 switch(dir)
 {
  case 0:// lam, !!! j k i
	for(j=1; j<=(sz-1); j++)// fi 
	{
	 a=ERH_RDS*cos(j*MESH_ANG_STP);
	 for(k=1; k<=(SZ_X3-1); k++)// z, without boundary
	  for(i=1; i<=(SZ_X1-1); i++)// lambda, without boundary
	  {
	   d=(j*(SZ_X3+1)+k)*(SZ_X1+1)+i;
	   pC[d]=pU[d]/a;
	  }
	}
  break;
 
  case 1: // fi, !!! i j k 
	for(i=1; i<=(sz-1); i++)// lambda, without boundary
	 for(j=1; j<=(SZ_X2-1); j++)// fi, without boundary
	  for(k=1; k<=(SZ_X3-1); k++)// z, without boundary
	  {
	   d=(i*(SZ_X2+1)+j)*(SZ_X3+1)+k;
	   pC[d]=pV[d]/ERH_RDS;
	  }
  break;

  case 2: // z, !!! i j k 
	for(i=1; i<=(sz-1); i++)// lambda
	 for(j=1; j<=(SZ_X2-1); j++)// fi, without boundary
	  for(k=1; k<=(SZ_X3-1); k++)// z, without boundary
	  {
	   d=(i*(SZ_X2+1)+j)*(SZ_X3+1)+k;
	   pC[d]=pW[d];
	   //pC[d]=pW[d]-2.*pMu[d]/ERH_RDS;
	  }
  break;
 }
}


void SetActualDensity(int sz, double *pP, double *pT, double *pD, int dir)
{
 int i, j, k, d;

 if(dir>0)
  for(i=1; i<=(sz-1); i++)// lambda, without boundary
   for(j=1; j<=(SZ_X2-1); j++)// fi, without boundary
	for(k=1; k<=(SZ_X3-1); k++)// z, without boundary
	{
	 d=(i*(SZ_X2+1)+j)*(SZ_X3+1)+k;
	 pD[d]=pP[d]/(GAS_R*pT[d]);
	}
 else
  for(j=1; j<=(sz-1); j++)// fi 
   for(k=1; k<=(SZ_X3-1); k++)// z, without boundary
	for(i=1; i<=(SZ_X1-1); i++)// lambda, without boundary
	{
	 d=(j*(SZ_X3+1)+k)*(SZ_X1+1)+i;
	 pD[d]=pP[d]/(GAS_R*pT[d]);
	} 
}


void SetActualDiff(int sz, double *pU, double *pV, double *pMu, double *pP, double *pD, int dir)
{
 int i, j, k, d;
 double d1, d2, a, b, c, q;

  switch(dir)
 {
  case 0:// lam, !!! j k i
	q=0.08*COE_DISS*MESH_ANG_STP*MESH_ANG_STP;
	for(j=1; j<=(sz-1); j++)// fi
	{
	 a=ERH_RDS*cos(j*MESH_ANG_STP);
	 b=sin(j*MESH_ANG_STP);
	 c=(a*a+ERH_RDS_2)*q;
	 for(k=1; k<=(SZ_X3-1); k++)// z, without boundary
	  for(i=1; i<=(SZ_X1-1); i++)// lam, without boundary
	  {
	   d=(j*(SZ_X3+1)+k)*(SZ_X1+1)+i;
	   d1=(pU[d+1]-pU[d-1])/(2.*a*MESH_ANG_STP)-pV[d]*b/a
		  -(pV[((j+1)*(SZ_X3+1)+k)*(SZ_X1+1)+i]-pV[((j-1)*(SZ_X3+1)+k)*(SZ_X1+1)+i])/(2.*ERH_RDS*MESH_ANG_STP);
	   d2=(pV[d+1]-pV[d-1])/(2.*a*MESH_ANG_STP)+pU[d]*b/a
		  +(pU[((j+1)*(SZ_X3+1)+k)*(SZ_X1+1)+i]-pU[((j-1)*(SZ_X3+1)+k)*(SZ_X1+1)+i])/(2.*ERH_RDS*MESH_ANG_STP);
	   pMu[d]=pD[d]*(COE_DISS*COE_K0+c*sqrt(d1*d1+d2*d2))/(a*a);
	   if(i==1) pMu[(j*(SZ_X3+1)+k)*(SZ_X1+1)+0]=pMu[d];
	   else if(i==(SZ_X1-1)) pMu[(j*(SZ_X3+1)+k)*(SZ_X1+1)+SZ_X1]=pMu[d];
	  }
    }
  break;
 
  case 1: // fi, !!! i j k 
	q=0.08*COE_DISS*MESH_ANG_STP*MESH_ANG_STP;
	for(i=1; i<=(sz-1); i++)// lambda
	 for(j=1; j<=(SZ_X2-1); j++)// fi, without boundary
	 {
	  a=ERH_RDS*cos(j*MESH_ANG_STP);
	  b=sin(j*MESH_ANG_STP);
	  c=(a*a+ERH_RDS_2)*q;
	  for(k=1; k<=(SZ_X3-1); k++)// z, without boundary
	  {
	   d=(i*(SZ_X2+1)+j)*(SZ_X3+1)+k;
	   d1=(pU[((i+1)*(SZ_X2+1)+j)*(SZ_X3+1)+k]-pU[((i-1)*(SZ_X2+1)+j)*(SZ_X3+1)+k])/(2.*a*MESH_ANG_STP)
		  -(pV[(i*(SZ_X2+1)+j+1)*(SZ_X3+1)+k]-pV[(i*(SZ_X2+1)+j-1)*(SZ_X3+1)+k])/(2.*ERH_RDS*MESH_ANG_STP)
		  -pV[d]*b/a;
	   d2=(pV[((i+1)*(SZ_X2+1)+j)*(SZ_X3+1)+k]-pV[((i-1)*(SZ_X2+1)+j)*(SZ_X3+1)+k])/(2.*a*MESH_ANG_STP)
		  +(pU[(i*(SZ_X2+1)+j+1)*(SZ_X3+1)+k]-pU[(i*(SZ_X2+1)+j-1)*(SZ_X3+1)+k])/(2.*ERH_RDS*MESH_ANG_STP)
		  +pU[d]*b/a;
	   pMu[d]=pD[d]*(COE_DISS*COE_K0+c*sqrt(d1*d1+d2*d2))/ERH_RDS_2;
	   if(j==1) pMu[(i*(SZ_X2+1)+0)*(SZ_X3+1)+k]=pMu[d];
	   else if(j==(SZ_X2-1)) pMu[(i*(SZ_X2+1)+SZ_X2)*(SZ_X3+1)+k]=pMu[d];
	  }
    }
  break;

  case 2: // z, !!! i j k
    a=600.;//0.4*1500.;
	q=a*a/(2.*MESH_VRT_STP);
	for(i=1; i<=(sz-1); i++)// lambda, without boundary
	 for(j=1; j<=(SZ_X2-1); j++)// fi, without boundary
      for(k=1; k<=(SZ_X3-1); k++)// z, without boundary
      {
	   d=(i*(SZ_X2+1)+j)*(SZ_X3+1)+k;
	   d1=pU[d+1]-pU[d-1];
	   d2=pV[d+1]-pV[d-1];
	   pMu[d]=pD[d]*q*sqrt(d1*d1+d2*d2);
	   if(k==1) pMu[d-1]=pMu[d];
	   else if(k==(SZ_X3-1)) pMu[d+1]=pMu[d];
	 }
  break;
 } 
}


void SetFPart(int sz, double r, double *pU, double *pV, double *pW, double *pP,
			  double **pPS, double *pD, double *pF, int val, int dir)// F1=F2=F3=tau*F
{
 int i, j, k, d;
 double a, b;
 
 switch(val)
 {
  case 0: // u
	  if(dir>0)// fi, z, [i j k]
	  {
	   for(i=1; i<=(sz-1); i++)// lambda, without boundary
		for(j=1; j<=(SZ_X2-1); j++)// fi, without boundary
		{
		 a=cos(j*MESH_ANG_STP);
		 b=sin(j*MESH_ANG_STP);
		 for(k=1; k<=(SZ_X3-1); k++)// z, without boundary
		 {
		  d=(i*(SZ_X2+1)+j)*(SZ_X3+1)+k;
		  pF[d]=tau*(-(pP[((i+1)*(SZ_X2+1)+j)*(SZ_X3+1)+k]-pP[((i-1)*(SZ_X2+1)+j)*(SZ_X3+1)+k])/
		  (2.*MESH_ANG_STP*a*ERH_RDS*pD[d])+
		  (pV[d]*b-pW[d]*a)*(2.*ERH_OMEGA+pU[d]/(a*ERH_RDS)));
		 }
	   }
	  }
	  else// lam, [j k i]
	  {
	   for(j=1; j<=(sz-1); j++)// fi, without boundary
	   {
	    a=cos(j*MESH_ANG_STP);
		b=sin(j*MESH_ANG_STP);
		for(k=1; k<=(SZ_X3-1); k++)// z, without boundary
	     for(i=1; i<=(SZ_X1-1); i++)// lambda, without boundary
		 {
		  d=(j*(SZ_X3+1)+k)*(SZ_X1+1)+i;
		  pF[d]=tau*(-(pP[d+1]-pP[d-1])/(2.*MESH_ANG_STP*a*ERH_RDS*pD[d])+
		  (pV[d]*b-pW[d]*a)*(2.*ERH_OMEGA+pU[d]/(a*ERH_RDS)));
		 }
	   }
	  }
  break;

  case 1:// v
	  if(dir>0)// fi, z, [i j k]
	  {
	   for(i=1; i<=(sz-1); i++)// lambda, without boundary
		for(j=1; j<=(SZ_X2-1); j++)// fi, without boundary
		{
		 a=ERH_RDS*cos(j*MESH_ANG_STP);
		 b=sin(j*MESH_ANG_STP);
		 for(k=1; k<=(SZ_X3-1); k++)// z, without boundary
		 {
		  d=(i*(SZ_X2+1)+j)*(SZ_X3+1)+k;
		  pF[d]=-tau*((pP[(i*(SZ_X2+1)+j+1)*(SZ_X3+1)+k]-pP[(i*(SZ_X2+1)+j-1)*(SZ_X3+1)+k])/
		  (2.*MESH_ANG_STP*ERH_RDS*pD[d])+pU[d]*b*(2.*ERH_OMEGA+pU[d]/a
		  +pV[d]*pW[d]/ERH_RDS));
		}
	   }
	  }
	  else// lam, [j k i]
	  {
	   for(j=1; j<=(sz-1); j++)// fi, without boundary
	   {
		a=ERH_RDS*cos(j*MESH_ANG_STP);
		b=sin(j*MESH_ANG_STP);
		for(k=1; k<=(SZ_X3-1); k++)// z, without boundary
		 for(i=1; i<=(SZ_X1-1); i++)// lam, without boundary
		 {
		  d=(j*(SZ_X3+1)+k)*(SZ_X1+1)+i;
		  pF[d]=-tau*((pP[((j+1)*(SZ_X3+1)+k)*(SZ_X1+1)+i]-pP[((j-1)*(SZ_X3+1)+k)*(SZ_X1+1)+i])/
		  (2.*MESH_ANG_STP*ERH_RDS*pD[d])+pU[d]*b*(2.*ERH_OMEGA+pU[d]/a
		  +pV[d]*pW[d]/ERH_RDS));
		 }
	   }
	  }
  break;

  case 2:// T
	 if(dir>0)// fi, z, [i j k]
	 {
	  for(i=1; i<=(sz-1); i++)// lambda, without boundary
	   for(j=1; j<=(SZ_X2-1); j++)// fi, without boundary
	   {
	    //a=ERH_RDS*cos(j*MESH_ANG_STP);
		for(k=1; k<=(SZ_X3-1); k++)// z, without boundary
		{
		 d=(i*(SZ_X2+1)+j)*(SZ_X3+1)+k;
		 pF[d]=0.;
		 /*
		 if(r<1.) b=(pPS[1][(i*(SZ_X2+1)+j)*(SZ_X3+1)+k]-pPS[0][(i*(SZ_X2+1)+j)*(SZ_X3+1)+k])/(3600.*12.);
		 else b=(pPS[2][(i*(SZ_X2+1)+j)*(SZ_X3+1)+k]-pPS[1][(i*(SZ_X2+1)+j)*(SZ_X3+1)+k])/(3600.*12.);
		 pF[d]=tau*(b+
		 pU[d]*(pP[((i+1)*(SZ_X2+1)+j)*(SZ_X3+1)+k]-pP[((i-1)*(SZ_X2+1)+j)*(SZ_X3+1)+k])/(2.*MESH_ANG_STP*a) 
		 +pV[d]*(pP[(i*(SZ_X2+1)+j+1)*(SZ_X3+1)+k]-pP[(i*(SZ_X2+1)+j-1)*(SZ_X3+1)+k])/(2.*MESH_ANG_STP*ERH_RDS)
		 +pW[d]*(pP[d+1]-pP[d-1])/(2.*MESH_VRT_STP))
		 /(GAS_Cp*pD[d]);
		 */
		}
	   }
	 }
	 else// lam, [j k i]
	 {
	  for(j=1; j<=(sz-1); j++)// fi, without boundary
	  {
	   //a=ERH_RDS*cos(j*MESH_ANG_STP);
	   for(k=1; k<=(SZ_X3-1); k++)// z, without boundary
	    for(i=1; i<=(SZ_X1-1); i++)// lambda, without boundary
		{
		 d=(j*(SZ_X3+1)+k)*(SZ_X1+1)+i;
		 pF[d]=0.;
		 /*
		 if(r<1.) b=(pPS[1][(i*(SZ_X2+1)+j)*(SZ_X3+1)+k]-pPS[0][(i*(SZ_X2+1)+j)*(SZ_X3+1)+k])/(3600.*12.);
		 else b=(pPS[2][(i*(SZ_X2+1)+j)*(SZ_X3+1)+k]-pPS[1][(i*(SZ_X2+1)+j)*(SZ_X3+1)+k])/(3600.*12.);
		 pF[d]=tau*(b+
		 pU[d]*(pP[d+1]-pP[d-1])/(2.*MESH_ANG_STP*a) 
		 +pV[d]*(pP[((j+1)*(SZ_X3+1)+k)*(SZ_X1+1)+i]-pP[((j-1)*(SZ_X3+1)+k)*(SZ_X1+1)+i])/(2.*MESH_ANG_STP*ERH_RDS)
		 +pW[d]*(pP[(j*(SZ_X3+1)+k+1)*(SZ_X1+1)+i]-pP[(j*(SZ_X3+1)+k-1)*(SZ_X1+1)+i])/(2.*MESH_VRT_STP))
		 /(GAS_Cp*pD[d]);
		 */
		}
	   }
	 }
  break;
 }
}

