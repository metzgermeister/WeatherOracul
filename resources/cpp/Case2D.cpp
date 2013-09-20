// Case2D.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"

const double	T=24.*3600.;  // max time interval (in seconds) -- do not change
const double	T1=1.*3600.;  // calculate time interval (in seconds)
const double	dt=100.; // seconds
const double	dh=1.5*M_PI/180.; // radians
const double	S=3.0; // parameter of model dissipation


char  flNames[9][9]={"01Hs.inp", "01Us.inp", "01Vs.inp",
					"02Hs.inp", "02Us.inp", "02Vs.inp",
					"03Hs.inp", "03Us.inp", "03Vs.inp"};
double pU[SZ_FI*SZ_LAM], pV[SZ_FI*SZ_LAM],
	   pU1[SZ_FI*SZ_LAM], pV1[SZ_FI*SZ_LAM],
	   pU2[SZ_FI*SZ_LAM], pV2[SZ_FI*SZ_LAM],
	   pF1[SZ_FI*SZ_LAM], pF2[SZ_FI*SZ_LAM],
	   pMu1[SZ_FI*SZ_LAM], pMu2[SZ_FI*SZ_LAM], 
	   pC[SZ_FI*SZ_LAM],
	   pH1[SZ_FI*SZ_LAM], pH2[SZ_FI*SZ_LAM], 
	   pH3[SZ_FI*SZ_LAM], pH[SZ_FI*SZ_LAM],
	   pFiBndCnd[2][3][2][SZ_LAM], pLamBndCnd[2][3][2][SZ_FI]; 



int LoadDataFromFiles();
int SplitSum(double *pUU, double *pUU1, double *pUU2);
int SetBoundaryCnd(double s, double *pU1, double *pU2, double *pV1, double *pV2);
int SplitLam(double *pR, double *pCC, double *pMMu, double *pFF1, double *pR1);
int MinusX1(int i1, int i2, int j, double *pR, double *pCC,
			double *pMMu, double *pFF1, double *pR1);
double QuadrInterpol(double r, double f0, double f1, double f2);
int SplitFi(double *pR, double *pCC, double *pMMu, double *pFF2, double *pR2);
int MinusX2(int j1, int j2, int i, double *pR, double *pCC,
			double *pMMu, double *pFF2, double *pR2);
int SetMu(double *pU, double *pV, double *pMu1, double *pMu2);
double fRo(double z);
double fG(double z, double fi);
int SetFU(double *pUU, double *pVV, double *pHH, double *pFF1, double *pFF2);
int SetFV(double *pUU, double *pVV, double *pHH, double *pFF1, double *pFF2);
int SetActualH(double r, double *pHH);
int LoadOrginalValues(double r, double *pUU1, double *pVV1);
int CalcNorma(double r);
int SaveToFile(char ch, double *pAr);





int _tmain(int argc, _TCHAR* argv[])
{
 int i;

 LoadDataFromFiles();
 for(i=0; (i*dt)<=(T1-dt); i++)
 {
  SetMu(pU, pV, pMu1, pMu2);
  SetBoundaryCnd(2.*(i+1)*dt/T, pU1, pU2, pV1, pV2);
  SetActualH(2.*i*dt/T, pH); 

  // Calc u
  SetFU(pU, pV, pH, pF1, pF2);
  memcpy(pC, pU, SZ_LAM*SZ_FI*sizeof(double));
  SplitLam(pU, pC, pMu1, pF1, pU1);
  memcpy(pC, pV, SZ_LAM*SZ_FI*sizeof(double));
  SplitFi(pU, pC, pMu2, pF2, pU2);

  // Calc v
  SetFV(pU, pV, pH, pF1, pF2);
  memcpy(pC, pU, SZ_LAM*SZ_FI*sizeof(double));
  SplitLam(pV, pC, pMu1, pF1, pV1);
  memcpy(pC, pV, SZ_LAM*SZ_FI*sizeof(double));
  SplitFi(pV, pC, pMu2, pF2, pV2);

  // Final
  SplitSum(pU, pU1, pU2);
  SplitSum(pV, pV1, pV2);
 }

 CalcNorma(2.*(i-1)*dt/T);
 SaveToFile('u', pU);
 SaveToFile('v', pV);
 printf("The end");
 getch();

 return 0;
}



int SplitLam(double *pR, double *pCC, double *pMMu, double *pFF1, double *pR1)
{
 char flg1, flg2;
 int i, j, beg;
 double p, q, pm, qm, g, z, h;

   for(j=1; j<=(SZ_FI-2); j++)
   {
	flg1=0x01;
	z=dt/(2*dh*dh*ERH_RDS_2*cos(j*dh)*cos(j*dh));
	h=dh*ERH_RDS*cos(j*dh);
    for(i=1; i<=(SZ_LAM-2); i++)
	{
	 if(pCC[j*SZ_LAM+i]>0.) flg2=0x01;
	 else if(pCC[j*SZ_LAM+i]==0.) flg2=flg1;
	 else flg2=0x00;

	 if(flg2)
	 {
	  p=z*(h*pCC[j*SZ_LAM+i]+pMMu[j*SZ_LAM+i]+pMMu[j*SZ_LAM+i-1]);
	  q=z*(h*pCC[j*SZ_LAM+i]-pMMu[j*SZ_LAM+i]-pMMu[j*SZ_LAM+i+1]);

	  if(flg1==0x00)
	  {// from "-" (i-1) to "+" (i)
	   pm=z*(-h*pCC[j*SZ_LAM+i-1]+pMMu[j*SZ_LAM+i-1]+pMMu[j*SZ_LAM+i]);
	   qm=z*(h*pCC[j*SZ_LAM+i-1]+pMMu[j*SZ_LAM+i-1]+pMMu[j*SZ_LAM+i-2]);

	   g=1.+p+pm;
	   q=pR[j*SZ_LAM+i]+pFF1[j*SZ_LAM+i]+q*(pR[j*SZ_LAM+i]-pR[j*SZ_LAM+i+1]);
	   qm=pR[j*SZ_LAM+i-1]+pFF1[j*SZ_LAM+i-1]+qm*(pR[j*SZ_LAM+i-2]-pR[j*SZ_LAM+i-1]);
	   p=pm*q+p*qm;

	   pR1[j*SZ_LAM+i-1]=(qm+p)/g;
	   pR1[j*SZ_LAM+i]=(q+p)/g;

	   MinusX1(beg, i-2, j, pR, pCC, pMMu, pFF1, pR1);
	  }
	  else // "+", "+" (pCC[i]>0)
	  {
	   pR1[j*SZ_LAM+i]=pR1[j*SZ_LAM+i-1]+(pR[j*SZ_LAM+i]-pR1[j*SZ_LAM+i-1]
			+q*(pR[j*SZ_LAM+i]-pR[j*SZ_LAM+i+1])+pFF1[j*SZ_LAM+i])/(1.+p);
	  }
	 }
	 else if(flg1==0x01) beg=i; // from "+" (i-1) to "-" (i)

	 flg1=flg2;
	}

	if(flg1==0x00) MinusX1(beg, SZ_LAM-2, j, pR, pCC, pMMu, pFF1, pR1);
   }
 return 0;
}



int MinusX1(int i1, int i2, int j, double *pR, double *pCC,
			double *pMMu, double *pFF1, double *pR1)
{
 int i;
 double p, q, z=dt/(2.*dh*dh*ERH_RDS_2*cos(j*dh)*cos(j*dh)),
		h=dh*ERH_RDS*cos(j*dh);

 for(i=i2; i>=i1; i--)
 {
  p=z*(-h*pCC[j*SZ_LAM+i]+pMMu[j*SZ_LAM+i]+pMMu[j*SZ_LAM+i+1]);
  q=z*(h*pCC[j*SZ_LAM+i]+pMMu[j*SZ_LAM+i]+pMMu[j*SZ_LAM+i-1]);

  pR1[j*SZ_LAM+i]=pR1[j*SZ_LAM+i+1]+(pR[j*SZ_LAM+i]-pR1[j*SZ_LAM+i+1]
			+q*(pR[j*SZ_LAM+i-1]-pR[j*SZ_LAM+i])+pFF1[j*SZ_LAM+i])/(1.+p);
 }
 return 0;
}



int SplitFi(double *pR, double *pCC, double *pMMu, double *pFF2, double *pR2)
{
 char flg1, flg2;
 int i, j, beg;
 double p, q, pm, qm, g, z=dt/(2*dh*dh*ERH_RDS_2),
		h=dh*ERH_RDS;

   for(i=1; i<=(SZ_LAM-2); i++)
   {
	flg1=0x01;
	for(j=1; j<=(SZ_FI-2); j++)
	{
	 if(pCC[j*SZ_LAM+i]>0.) flg2=0x01;
	 else if(pCC[j*SZ_LAM+i]==0.) flg2=flg1;
	 else flg2=0x00;

	 if(flg2)
	 {
	  p=z*(h*pCC[j*SZ_LAM+i]+pMMu[j*SZ_LAM+i]+pMMu[(j-1)*SZ_LAM+i]);
	  q=z*(h*pCC[j*SZ_LAM+i]-pMMu[j*SZ_LAM+i]-pMMu[(j+1)*SZ_LAM+i]);

	  if(flg1==0x00)
	  {// from "-" (j-1) to "+" (j)
	   pm=z*(-h*pCC[(j-1)*SZ_LAM+i]+pMMu[(j-1)*SZ_LAM+i]+pMMu[j*SZ_LAM+i]);
	   qm=z*(h*pCC[(j-1)*SZ_LAM+i]+pMMu[(j-1)*SZ_LAM+i]+pMMu[(j-2)*SZ_LAM+i]);

	   g=1.+p+pm;
	   q=pR[j*SZ_LAM+i]+pFF2[j*SZ_LAM+i]+q*(pR[j*SZ_LAM+i]-pR[(j+1)*SZ_LAM+i]);
	   qm=pR[(j-1)*SZ_LAM+i]+pFF2[(j-1)*SZ_LAM+i]+qm*(pR[(j-2)*SZ_LAM+i]-pR[(j-1)*SZ_LAM+i]);
	   p=pm*q+p*qm;

	   pR2[(j-1)*SZ_LAM+i]=(qm+p)/g;
	   pR2[j*SZ_LAM+i]=(q+p)/g;

	   MinusX2(beg, j-2, i, pR, pCC, pMMu, pFF2, pR2);
	  }
	  else // "+", "+" (pCC[j]>0)
	  {
	   pR2[j*SZ_LAM+i]=pR2[(j-1)*SZ_LAM+i]+(pR[j*SZ_LAM+i]-pR2[(j-1)*SZ_LAM+i]
			+q*(pR[j*SZ_LAM+i]-pR[(j+1)*SZ_LAM+i])+pFF2[j*SZ_LAM+i])/(1.+p);
	  }
	 }
	 else if(flg1==0x01) beg=j; // from "+" (j-1) to "-" (j)

	 flg1=flg2;
	}

	if(flg1==0x00) MinusX2(beg, SZ_FI-2, i, pR, pCC, pMMu, pFF2, pR2);
   }
 return 0;
}



int MinusX2(int j1, int j2, int i, double *pR, double *pCC,
			double *pMMu, double *pFF2, double *pR2)
{
 int j;
 double p, q, z=dt/(2.*dh*dh*ERH_RDS_2),
		h=dh*ERH_RDS;

 for(j=j2; j>=j1; j--)
 {
  p=z*(-h*pCC[j*SZ_LAM+i]+pMMu[j*SZ_LAM+i]+pMMu[(j+1)*SZ_LAM+i]);
  q=z*(h*pCC[j*SZ_LAM+i]+pMMu[j*SZ_LAM+i]+pMMu[(j-1)*SZ_LAM+i]);

  pR2[j*SZ_LAM+i]=pR2[(j+1)*SZ_LAM+i]+(pR[j*SZ_LAM+i]-pR2[(j+1)*SZ_LAM+i]
			+q*(pR[(j-1)*SZ_LAM+i]-pR[j*SZ_LAM+i])+pFF2[j*SZ_LAM+i])/(1.+p);
 }
 return 0;
}


int SetFU(double *pUU, double *pVV, double *pHH, double *pFF1, double *pFF2)
{
 int i, j;
 double cs, sn, g;

 for(j=1; j<=(SZ_FI-2); j++)
 {
  cs=ERH_RDS*cos(j*dh);
  sn=dt*sin(j*dh);
  g=dt*fG(SRF_ALT, j*dh)/(2.*cs*dh);
  for(i=1; i<=(SZ_LAM-2); i++)
  {
   pFF1[j*SZ_LAM+i]=-g*(pHH[j*SZ_LAM+i+1]-pHH[j*SZ_LAM+i-1]);
   pFF2[j*SZ_LAM+i]=sn*pVV[j*SZ_LAM+i]*(pUU[j*SZ_LAM+i]/cs+
						2.*ERH_OMEGA);
  }
 }
 return 0;
}



int SetFV(double *pUU, double *pVV, double *pHH, double *pFF1, double *pFF2)
{
 int i, j;
 double cs, sn, g;

 for(j=1; j<=(SZ_FI-2); j++)
 {
  cs=ERH_RDS*cos(j*dh);
  sn=dt*sin(j*dh);
  g=dt*fG(SRF_ALT, j*dh)/(2.*dh*ERH_RDS);
  for(i=1; i<=(SZ_LAM-2); i++)
  {
   pFF1[j*SZ_LAM+i]=-sn*pUU[j*SZ_LAM+i]*(pUU[j*SZ_LAM+i]/cs+
						2.*ERH_OMEGA);
   pFF2[j*SZ_LAM+i]=-g*(pHH[(j+1)*SZ_LAM+i]-pHH[(j-1)*SZ_LAM+i]);
  }
 }
 return 0;
}


int SetActualH(double r, double *pHH)
{
 int i;
 
 for(i=0; i<=(SZ_FI*SZ_LAM-1); i++)
	pHH[i]=QuadrInterpol(r, pH1[i], pH2[i], pH3[i]);

 return 0;
}



int SetBoundaryCnd(double r, double *pUU1, double *pUU2, double *pVV1, double *pVV2)
{
 int i, j;
 // lam
 for(i=0; i<=(SZ_FI-1); i++)
  for(j=0; j<=(SZ_LAM-1); j+=(SZ_LAM-1))
    pUU1[i*SZ_LAM+j]=QuadrInterpol(r, pLamBndCnd[0][0][j/(SZ_LAM-1)][i],
		pLamBndCnd[0][1][j/(SZ_LAM-1)][i], pLamBndCnd[0][2][j/(SZ_LAM-1)][i]);
 for(i=0; i<=(SZ_FI-1); i++)
  for(j=0; j<=(SZ_LAM-1); j+=(SZ_LAM-1))
    pVV1[i*SZ_LAM+j]=QuadrInterpol(r, pLamBndCnd[1][0][j/(SZ_LAM-1)][i],
		pLamBndCnd[1][1][j/(SZ_LAM-1)][i], pLamBndCnd[1][2][j/(SZ_LAM-1)][i]);
 // fi
 for(j=0; j<=(SZ_LAM-1); j++)
  for(i=0; i<=(SZ_FI-1); i+=(SZ_FI-1))
	pUU2[i*SZ_LAM+j]=QuadrInterpol(r, pFiBndCnd[0][0][i/(SZ_FI-1)][j],
		pFiBndCnd[0][1][i/(SZ_FI-1)][j], pFiBndCnd[0][2][i/(SZ_FI-1)][j]);
 for(j=0; j<=(SZ_LAM-1); j++)
  for(i=0; i<=(SZ_FI-1); i+=(SZ_FI-1))
	pVV2[i*SZ_LAM+j]=QuadrInterpol(r, pFiBndCnd[1][0][i/(SZ_FI-1)][j],
		pFiBndCnd[1][1][i/(SZ_FI-1)][j], pFiBndCnd[1][2][i/(SZ_FI-1)][j]);
 return 0;
}


int SetMu(double *pUU, double *pVV, double *pM1, double *pM2)
{
 int i, j;
 double dd, dc, tg, cs, ss=S/fRo(SRF_ALT), 
		k0=50000.;

 for(j=1; j<=(SZ_FI-2); j++)
 {
  tg=tan(j*dh);
  cs=dh*cos(j*dh);
  for(i=1; i<=(SZ_LAM-2); i++)
  {
   dc=(pUU[j*SZ_LAM+i+1]-pUU[j*SZ_LAM+i-1])/(2.*cs)-
	   (pVV[(j+1)*SZ_LAM+i]-pVV[(j-1)*SZ_LAM+i])/(2.*dh)-
	   pVV[j*SZ_LAM+i]*tg;
   dd=(pVV[j*SZ_LAM+i+1]-pVV[j*SZ_LAM+i-1])/(2.*cs)+
	   (pUU[(j+1)*SZ_LAM+i]-pUU[(j-1)*SZ_LAM+i])/(2.*dh)+
	   pUU[j*SZ_LAM+i]*tg;
   pM2[j*SZ_LAM+i]=ss*(k0+0.08*ERH_RDS*(cs*cs+dh*dh)*sqrt(dd*dd+dc*dc));
   pM1[j*SZ_LAM+i]=pM2[j*SZ_LAM+i];
  }
 }
 return 0;
}



int CalcNorma(double r)
{
 int i, j;
 double maxU=0., L2U=0., maxV=0., L2V=0., z;

 LoadOrginalValues(r, pU2, pV2);
 // For u
 for(j=1; j<=(SZ_FI-2); j++)
   for(i=1; i<=(SZ_LAM-2); i++)
   {
    z=pU2[j*SZ_LAM+i]-pU[j*SZ_LAM+i];
	z*=z;
	L2U+=z;
	if(z>maxU) maxU=z;
   }
  maxU=sqrt(maxU);
  L2U=dh*sqrt(L2U);

 // For v
 for(j=1; j<=(SZ_FI-2); j++)
   for(i=1; i<=(SZ_LAM-2); i++)
   {
    z=pV2[j*SZ_LAM+i]-pV[j*SZ_LAM+i];
	z*=z;
	L2V+=z;
	if(z>maxV) maxV=z;
   }
  maxV=sqrt(maxV);
  L2V=dh*sqrt(L2V);

 printf("maxU=%9.2e L2U=%9.2e\nmaxV=%9.2e L2V=%9.2e\n",
		maxU, L2U, maxV, L2V); 
 
 return 0;
}


int LoadOrginalValues(double r, double *pUU1, double *pVV1)
{
 char buf[64];
 int i, j, k;
 double *pArOut[2]={pUU1, pVV1}, z[3];
 FILE *pFile[3];

 for(i=0; i<=1; i++)
 {
  for(k=0; k<=2; k++)
  {
   sprintf(buf, "IN\\%s", flNames[3*k+1+i]); 
   pFile[k]=fopen(buf, "r");
  }

  for(j=0; j<=(SZ_FI*SZ_LAM-1); j++)
  {
   for(k=0; k<=2; k++)
   {
    fgets(buf, 64, pFile[k]);
    z[k]=atof(buf);
   }
   pArOut[i][j]=QuadrInterpol(r, z[0], z[1], z[2]);
  }

  for(k=0; k<=2; k++) fclose(pFile[k]);
 }
 
 return 0;
}


double QuadrInterpol(double r, double f0, double f1, double f2)
{// f(0)=f0,  f(1)=f1, f(2)=f2 
 return f0-r*(0.5*f2-2.*f1+1.5*f0-r*(0.5*f2-f1+0.5*f0));
}


int SplitSum(double *pUU, double *pUU1, double *pUU2)
{
 int i, j;

 for(i=1; i<=(SZ_FI-2); i++)
  for(j=1; j<=(SZ_LAM-2); j++)
    pUU[i*SZ_LAM+j]=pUU1[i*SZ_LAM+j]+pUU2[i*SZ_LAM+j]-pUU[i*SZ_LAM+j];

 // Boundary condition
 for(j=0; j<=(SZ_LAM-1); j+=(SZ_LAM-1))
   for(i=0; i<=(SZ_FI-1); i++)
     pUU[i*SZ_LAM+j]=pUU1[i*SZ_LAM+j];

 for(i=0; i<=(SZ_FI-1); i+=(SZ_FI-1))
	for(j=0; j<=(SZ_LAM-1); j++)
     pUU[i*SZ_LAM+j]=pUU2[i*SZ_LAM+j];
 
 return 0;
}



int SaveToFile(char ch, double *pAr)
{
 char buf[64];
 int i;
 FILE *pFile;
 
 sprintf(buf, "OUT\\%c.out", ch); 
 pFile=fopen(buf, "w");
 for(i=0; i<=(SZ_FI*SZ_LAM-1); i++)
 {
  sprintf(buf, "%.2f\n", pAr[i]);
  fputs(buf, pFile);
 }
 fclose(pFile);
 
 return 0;
}


int LoadDataFromFiles()
{
 char buf[64];
 int i, j, k, m[2];
 double *pHUV[3]={pU, pV, pH3};
 FILE *pFile;

 // Set initial condition for u, v
 for(i=0; i<=1; i++)
 {
  sprintf(buf, "IN\\%s", flNames[1+i]); 
  pFile=fopen(buf, "r");
  for(j=0; j<=(SZ_FI*SZ_LAM-1); j++)
  {
   fgets(buf, 64, pFile);
   pHUV[i][j]=atof(buf);	
  }
  fclose(pFile);
 }

 // Set H
 pHUV[0]=pH1;
 pHUV[1]=pH2;
 for(i=0; i<=2; i++)
 {
  sprintf(buf, "IN\\%s", flNames[3*i]); 
  pFile=fopen(buf, "r");
  for(j=0; j<=(SZ_FI*SZ_LAM-1); j++)
  {
   fgets(buf, 64, pFile);
   pHUV[i][j]=atof(buf);	
  }
  fclose(pFile);
 }

 // For boundary condition
 // [u, v] [time: 1, 2, 3] [boundary: 0, SZ] [knots]
 for(i=0; i<=1; i++) // u, v
  for(j=0; j<=2; j++) // time 1, 2, 3
  {
   sprintf(buf, "IN\\%s", flNames[3*j+i+1]); 
   pFile=fopen(buf, "r");
   for(k=0; k<=(SZ_FI*SZ_LAM-1); k++)
   {
    m[0]=k/SZ_LAM;
	m[1]=k%SZ_LAM;
	fgets(buf, 64, pFile);	
	if((m[0])==0) pFiBndCnd[i][j][0][m[1]]=atof(buf);
	else if((m[0])==(SZ_FI-1)) pFiBndCnd[i][j][1][m[1]]=atof(buf);
	if((m[1])==0) pLamBndCnd[i][j][0][m[0]]=atof(buf);
	else if((m[1])==(SZ_LAM-1)) pLamBndCnd[i][j][1][m[0]]=atof(buf);
   }
   fclose(pFile);
  }
 
 return 0;
}


double fRo(double z)
{
 return 1.22*exp(-RO_RATE*z);
}


double fG(double z, double fi)
{
 double cs=cos(2.*fi);
 return ERH_G-0.025928*cs+(6.9e-5)*cs*cs-(3.086e-6)*z;
}