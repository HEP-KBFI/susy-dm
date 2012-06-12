/*======  Spectrum calculator  ========= 
   Choose RGE from the list below. SuSpect is included 
   in micrOMEGAs, to use another code define the path 
   to the corresponding package in lib/Makefile
=====================================*/ 
#define RGE suspect
     /* choose 'suspect','isajet','softSusy','spheno'*/

/*=========   SUSY scenario  ==========
  One can define SUGRA, AMSB, EWSB (for low scale input). 
  By default the program reads SLHA data file 
=======================================*/
#define SUGRA 
/* #define AMSB  */
/* #define EWSB  */

/*====== Modules ===============
   Keys to switch on 
   various modules of micrOMEGAs
================================*/

#define MASSES_INFO      
      /* Display information about SUSY and Higgs masses 
      */
#define CONSTRAINTS     
      /* Display  deltarho, B_>sgamma, Bs->mumu, gmuon and
         check LEP mass limits 
      */ 
#define OMEGA            
      /* Calculate relic density and display contribution of
         individual channels 
      */
//#define INDIRECT_DETECTION  
      /* Compute spectra of gamma/positron/neutrinos
         for DM annihilation; calculate <sigma*v> and
         integrate gamma signal over DM galactic squared
         density for given line of sight.  
      */

#define DIRECT_DETECTION
      
/*#define RESET_FORMFACTORS*/
      /* Modify default nucleus form factors, 
         DM velocity distribution,
         A-dependence of Fermi-dencity
      */     
//#define CDM_NUCLEON 
      /* Calculate amplitudes and cross-sections for 
         CDM-mucleon collisions 
      */  
/* #define TEST_Direct_Detection */
      /* 
        Compare analytical formula for DD against micrOMEGAS calculation.
        As well as compare tree level and box improved approaches.
       */      
#define CDM_NUCLEUS
      /* Calculate number of events for 1kg*day 
         and recoil energy distibution for various nuclei
      */
//#define DECAYS
      /* Calculate decay widths and branchings  */      

#define CROSS_SECTIONS
      /* Calculate cross sections of reactions specified by the user */

/*===== end of Modules  ======*/

/*===== Options ========*/
/* #define SHOWPLOTS */
     /* Display  graphical plots on the screen */ 

//#define DEBUG

/*===== End of DEFINE  settings ===== */

#define PB_IN_CM2 pow(10,-36)

#include"../sources/micromegas.h"
#include"../sources/micromegas_aux.h"
#include"lib/pmodel.h"

//These are needed to call suspect
#include"lib/suspect_path.h"
#include<sys/wait.h>
#include<unistd.h>
#define FIN  "suspect2_lha.in"
#define FOUT "suspect2_lha.out"

#define SUGRAMODEL_(A) A ## SUGRA
#define SUGRAMODEL(A) SUGRAMODEL_(A)

#define AMSBMODEL_(A) A ## AMSB
#define AMSBMODEL(A) AMSBMODEL_(A)

#define EWSBMODEL_(A) A ## EwsbMSSM
#define EWSBMODEL(A) EWSBMODEL_(A)

#define PRINTRGE_(A) printf(" Spectrum calculator is %s\n", #A)
#define PRINTRGE(A)  PRINTRGE_(A)


int runSuspect() {
  char buff[2000];
  int err;

  if(!access(FIN, R_OK)) unlink(FIN);

  sprintf(buff,"%s/suspect.exe",SUSPECT);
  if(access( buff,X_OK))
  { printf("Executable \n %s\n is not found. Program stops.\n",buff);
    exit(13);
  }  
  
  err=System(buff);   
  if(err>=0)  err=slhaRead(FOUT,0); else cleanSLHAdata();   
  if(delFiles){unlink(FIN);unlink(FOUT);unlink("suspect2.out"); }
  return err;
}

typedef struct{
	double u;
	double l;
	char name[20];
} bounds;

bounds* load_bounds(FILE* file, int* n_bounds) {
	bounds* b = NULL;
	float u=0;
	float l=0;
	char str[10];
	
	int count = 0;
	while(fscanf(file, "%[^=]=%f~%f\n", str, &l, &u) == 3) {
		count++;
		bounds* new_b = realloc(b, count * sizeof(bounds));
		if(new_b != NULL) {
			b = new_b;
			b[count-1].u = u;	
			b[count-1].l = l;
			strcpy(b[count-1].name, str);
		} else {
			fprintf(stderr, "Error allocating memory\n");
			return NULL;
		}	
	}
	*n_bounds = count;
	return b;
}

int is_in_bounds(bounds* b, int n_bounds, const char* name, double val) {
	int i = 0;
	for(i = 0; i < n_bounds; i++) {
		if(strcmp(b[i].name, name) == 0) {
			if (val >= b[i].l && val <= b[i].u) {
				return 1;
			}
			else {
				return 0;
			}
		}
	}
	printf("Warning: %s bounds not specified\n", name);
	return -1;
}

void check_value(bounds* b, int n_bounds, const char* name, const double value) {
	if (is_in_bounds(b, n_bounds, name, value) == 0) {
		printf("%s outside limits: %.2E\n", name, value);
		#ifndef DEBUG
		exit(1);
		#endif
	}
}

int main(int argc,char** argv)
{ 
	//Outputs
	double mu = 0;
	double m_a = 0;
	
	bounds* b = NULL;
	int n_bounds = 0;
	FILE* f = fopen("bounds.txt", "r");
	if(f == NULL) {
		fprintf(stderr, "Error opening file bounds.txt for bounds\n");
		exit(1);
	}
	b = load_bounds(f, &n_bounds);
	if(b == NULL) {
		fprintf(stderr, "Error loading bounds\n");
		exit(1);
	}
	fclose(f);


	int err;
	char cdmName[10];
	int spin2, charge3,cdim;

	delFiles=0; /* switch to save/delete RGE input/output */
	ForceUG=0;  /* to Force Unitary Gauge assign 1 */

	double m0,mhf,a0,tb;
	double gMG1, gMG2, gMG3,  gAl, gAt, gAb,  sgn, gMHu,  gMHd,
		gMl2, gMl3, gMr2, gMr3, gMq2, gMq3, gMu2, gMu3, gMd2, gMd3;
        
	printf("This program accepts input from the standard input in the following format:\n");
	printf("m0 mhf a0 tb Mtp MbMb alfSMZ sgn\n");
	printf("example: .\\main < sample.txt\n");
 
	printf("\n\n\n========= mSUGRA scenario =====\n");
	PRINTRGE(RGE);
	double Mtp,MbMb,alfSMZ;
	
	int input_result = 0;
	#ifndef DEBUG	
	//input_result = scanf("%lf %lf %lf %lf %lf %lf %lf %lf", &m0, &mhf, &a0, &tb, &Mtp, &MbMb, &alfSMZ, &sgn);
  input_result = scanf("%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf",
                      &tb, &Mtp, &MbMb, &alfSMZ, &sgn,
                      &gMG1, &gMG2, &gMG3,
                      &gAl, &gAt, &gAb,
                      &gMHu, &gMHd,
                      &gMl2, &gMl3, &gMr2 ,&gMr3,
                      &gMq2, &gMq3, &gMu2, &gMu3, &gMd2, &gMd3); //23 parameters in total	

  #endif
	
	//#ifdef DEBUG
	//input_result = sscanf("120 500 -350 10 173.1 160.0 0.1 1.0", "%lf %lf %lf %lf %lf %lf %lf %lf", &m0, &mhf, &a0, &tb, &Mtp, &MbMb, &alfSMZ, &sgn);
	//#endif
	
	if(input_result != 23) {
		fprintf(stderr, "Format error in input\n");
		exit(1);
	}
  else {
	//printf("m0 = %.5f\n", m0);
	//printf("mhf = %.5f\n", mhf);
	//printf("a0 = %.5f\n", a0);
	//printf("tb = %.5f\n", tb);
	//printf("Mtp = %.5f\n", Mtp);
	//printf("MbMb = %.5f\n", MbMb);
	//printf("alfSMZ = %.5f\n", alfSMZ);
	//printf("sgn = %.1f\n", sgn);
    assignValW("Mtp", Mtp);
    assignValW("MbMb", MbMb);
    assignValW("alfSMZ", alfSMZ);
  }

	///*==== simulation of mSUGRA =====*/
	//gMG1=mhf, gMG2=mhf,gMG3=mhf;
	//gAl=a0,   gAt=a0,  gAb=a0;  gMHu=m0,  gMHd=m0;
	//gMl2=m0,  gMl3=m0, gMr2=m0, gMr3=m0;
	//gMq2=m0,  gMq3=m0, gMu2=m0, gMd2=m0, gMu3=m0, gMd3=m0;
	
	err= suspectSUGRA(tb, gMG1, gMG2, gMG3,  gAl,  gAt, gAb,  sgn, gMHu, gMHd, gMl2, gMl3, gMr2, gMr3, gMq2,  gMq3, gMu2, gMu3, gMd2, gMd3); 
//  runSuspect();
	
	int nw;
	printf("Warnings from spectrum calculator:\n");
	nw=slhaWarnings(stdout);
	if(nw==0) printf(" .....none\n");
	else exit(1);
	  
	if(err) exit(1);
	err=sortOddParticles(cdmName);
	if(err) { printf("Can't calculate %s\n",cdmName); return 1;}
	
	qNumbers(cdmName,&spin2, &charge3, &cdim);
	printf("\nDark matter candidate is '%s' with spin=%d/2  mass=%.2E\n", cdmName, spin2, Mcdm); 
	  
	
	if(charge3) { printf("Dark Matter has electric charge %d/3\n",charge3); exit(1);}
	if(cdim!=1) { printf("Dark Matter is a color particle\n"); exit(1);}
	if(strcmp(cdmName,"~o1")) printf(" ~o1 is not CDM\n"); 
	                            else o1Contents(stdout);
	
	double higgs_mass = 0;
	double higgs_width = 0;
    
    int ret=slhaWrite("slha.in");
    if(ret!=0) {
        printf("Error writing SLHA file: ret code: %d", ret);
        printf("Quitting");
        return -1;
    }
    slhaRead("slha.out",0);

//	#ifdef DEBUG
//	int i = 0;
//	for(i = 0;i < nModelParticles; i++) {
//		printf("%s = %.5f\n", ModelPrtcls[i].name, findValW(ModelPrtcls[i].mass));	
//	}
//	for(i = 0;i < nModelVars; i++) {
//		printf("%s = %.5f\n", varNames[i], varValues[i]);	
//	}
//	#endif
		
	#ifdef MASSES_INFO
	{
		printf("\n=== MASSES OF HIGGS AND SUSY PARTICLES: ===\n");
		
		int i;
		txtList LL;
		int dim;
		for(i = 0;i < nModelParticles; i++) {
			if(strcmp(ModelPrtcls[i].name, "h") == 0) {
				higgs_mass = findValW(ModelPrtcls[i].mass);
	       			higgs_width = pWidth(ModelPrtcls[i].name, &LL, &dim);
			}	
		}
		printHiggs(stdout);
		printMasses(stdout,1);
		if (is_in_bounds(b, n_bounds, "higgs", higgs_mass) == 0) {
			printf("Higgs mass outside limits: %.5f\n", higgs_mass);
			exit(1);
		}
	}
	#endif

	//Displays the bounds	
	int bound_c = 0;
	printf("\n===Physical constraints for precision data===\n");
	for(bound_c = 0; bound_c < n_bounds; bound_c++) {
		printf("%.2E <= %s <= %.2E\n", b[bound_c].l, b[bound_c].name, b[bound_c].u);
	}

	double v_deltarho = deltarho();
	double v_gmuon  = gmuon();
	double v_bsgnlo = bsgnlo(0); //In micromegas 2.4.5 bsgnlo returns a pointer to a datastructure containing ...(?) 
	double v_bsmumu = bsmumu();
	double v_btaunu = btaunu();
	#ifdef CONSTRAINTS
	{
		printf("\n\n==== Physical Constraints: =====\n"); 
		
		printf("deltarho=%.2E\n", v_deltarho);
		printf("gmuon=%.2E\n", v_gmuon);
		printf("bsgnlo=%.2E\n", v_bsgnlo);
		printf("bsmumu=%.2E\n", v_bsmumu);
		printf("btaunu=%.2E\n", v_btaunu);

		check_value(b, n_bounds, "deltarho", v_deltarho);
		check_value(b, n_bounds, "gmuon", v_gmuon);	
		check_value(b, n_bounds, "bsgnlo", v_bsgnlo);	
		check_value(b, n_bounds, "bsmumu", v_bsmumu);	
		check_value(b, n_bounds, "btaunu", v_btaunu);
		
		if(masslimits()==0) printf("MassLimits OK\n");
		else exit(1);
	}
	#endif

	double Omega,Xf;   
	#ifdef OMEGA
	{
		int fast=1;
		double Beps=1.E-5, cut=0.01;
		printf("\n==== Calculation of relic density =====\n");  
		Omega=darkOmega(&Xf,fast,Beps);
		printf("Xf=%.2e Omega=%.2e\n",Xf,Omega);
		printChannels(Xf,cut,Beps,1,stdout);
		if (is_in_bounds(b, n_bounds, "Omega", Omega) == 0) {
			printf("Omega outside limits: %.5f\n", Omega);
			exit(1);
		}
	}
	#endif
	
	//Find SUSY parameters mu and M_A(MH3)
	int i = 0;	
	for(i = 0;i < nModelVars; i++) {
		if(strcmp(varNames[i], "mu") == 0) {
			mu = varValues[i];
		}
		else if(strcmp(varNames[i], "MH3") == 0) {
			m_a = varValues[i];
		}
	}
  
	double nEvents;
	#ifdef CDM_NUCLEUS
	{
		double dNdE[300];
		
		printf("\n======== Direct Detection ========\n");    
		
		nEvents=nucleusRecoil(Maxwell,131,Z_Xe,J_Xe131,S00Xe131,S01Xe131,S11Xe131,FeScLoop,dNdE);
		printf("131Xe: Total number of events=%.2E /day/kg\n",nEvents);
		printf("Number of events in 10 - 50 KeV region=%.2E /day/kg\n", cutRecoilResult(dNdE,10,50));                                   
	}
	#endif 

	#ifdef DECAYS
	{  
	  txtList L;
	   int dim;
	   double width,br;
	   char * pname;
	
	   pname = "h";
	    width=pWidth(pname,&L,&dim);
	    printf("%s->%d*x :   total width=%E \n and Branchings:\n",pname,dim,width);
	    printTxtList(L,stdout);
	
	   pname = "l";
	    width=pWidth(pname,&L,&dim);
	    printf("%s->%d*x :   total width=%E \n and Branchings:\n",pname,dim,width);
	    printTxtList(L,stdout);
	    printf("Br(e,Ne,nl)= %E\n",findBr(L,"e,Ne,nl"));
	
	   pname = "~o2";
	    width=pWidth(pname,&L,&dim);
	    printf("%s->%d*x :   total width=%E \n and Branchings:\n",pname,dim,width);
	    printTxtList(L,stdout);
	    
	   pname = "~g";
	    width=pWidth(pname,&L,&dim);
	    printf("%s->%d*x :   total width=%E \n and Branchings:\n",pname,dim,width);
	    printTxtList(L,stdout);
	    
	    
	}
	#endif

//	double xs[4];
	double xs_p;
	double amplitude_p;
	
	#ifdef CROSS_SECTIONS
	{
//		double Pcm=500, cosmin=-0.99, cosmax=0.99, cs;
//		numout* cc;
//		printf("\n====== Calculation of cross section ====\n");  
//	
//		printf(" e^+, e^- annihilation\n");
//		Pcm=500.;
//		Helicity[0]=0.5;    /* helicity : spin projection on direction of motion   */    
//		Helicity[1]=-0.5;   /* helicities ={ 0.5, -0.5} corresponds to vector state */
//		printf("Process e,E->2*x at Pcm=%.3E GeV\n",Pcm);
//		cc=newProcess("e%,E%->2*x","eE_2x");
//		if(cc) {
//			int ntot,l;
//			char * name[4];
//	 		procInfo1(cc,&ntot,NULL,NULL);
//			for(l=1;l<=ntot; l++) {
//				int err;
//				double cs;
//				char txt[100];
//				procInfo2(cc,l,name,NULL);
//				sprintf(txt,"%3s,%3s -> %3s %3s  ",name[0],name[1],name[2],name[3]);
//				cs= cs22(cc,l,Pcm,cosmin,cosmax,&err);
//				if(err) printf("%-20.20s    Error\n",txt);
//				else if(cs*PB_IN_CM2>pow(10,-37)) printf("%-20.20s  %.6E [cm^2]\n",txt,cs*PB_IN_CM2);
//
//				#ifdef DEBUG	
//				printf("name[2]=%s\tname[3]=%s\tcs = %.6E\n", name[2], name[3], cs*PB_IN_CM2);	
//				#endif	
//
//				if(strcmp(name[2], "e") == 0 && strcmp(name[3], "E") == 0) {
//					xs[0] = cs*PB_IN_CM2;
//				} else if(strcmp(name[2], "A") == 0 && strcmp(name[3], "A") == 0) {
//					xs[1] = cs*PB_IN_CM2;
//				} else if(strcmp(name[2], "A") == 0 && strcmp(name[3], "Z") == 0) {
//					xs[2] = cs*PB_IN_CM2;
//				} else if(strcmp(name[2], "W+") == 0 && strcmp(name[3], "W-") == 0) {
//					xs[3] = cs*PB_IN_CM2;
//				}
//			}
//		}

		printf("\n====Calculations for the scattering cross-section and amplitude====\n");
		double pAsi[2];
		double nAsi[2];
		double pAsd[2];
		double nAsd[2];

		//Rest energy of the proton in the natural system of units [Gev]
		double nuc = 0.939;

		//The mass of the LOP is Mcdm, the mass of dark matter
		double kfac = 4.0/(M_PI*2.57e27) * pow(nuc*Mcdm/(nuc+Mcdm),2);
	
		nucleonAmplitudes(NULL, pAsi, nAsi, pAsd, nAsd);
		amplitude_p = pAsi[0]*pAsi[0];

		#ifdef DEBUG
		printf("pAsi[0] = %.3E\t pAsi[1] = %.3E\n", pAsi[0], pAsi[1]);
		printf("nAsi[0] = %.3E\t nAsi[1] = %.3E\n", nAsi[0], nAsi[1]);
		printf("pAsd[0] = %.3E\t pAsd[1] = %.3E\n", pAsd[0], pAsd[1]);
		printf("nAsd[0] = %.3E\t nAsd[1] = %.3E\n", nAsd[0], nAsd[1]);
		#endif
		
		xs_p = kfac*amplitude_p;
		printf("scattering cross-section for proton xs = %.6E\n", xs_p);
		printf("scattering amplitude for proton amp = %.6E\n", amplitude_p);
	}
	#endif
	printf("\nThis program prints its output in the following format to stdout\n");
	printf("m0 mhf a0 tb Mtp MbMb alfSMZ sgn");
	printf(" gmuon bsgamma deltarho bsmumu btaunu");
	printf(" Mcdm Omega Xf Zn11 Zn12 Zn13 Zn14 Zt11 Zt12 At Ab Al mu MG1 MG2 MNE1 MNE2 MNE3 MNE4 M_A(MH3) MSt1 MSt2 MSb1 MSb2 MSeL MSmL MSl1 MSG MSuL MC1 MC2 Xe131 higgs_mass higgs_width proton_xs proton_amplitude DMname\n");

	printf("RESULT: ");
	printf("%.6f %.6f %.6f %.6f %.6f %.6f %.6f %.0f ", m0, mhf, a0, tb, Mtp, MbMb, alfSMZ, sgn);
	printf("%.3e %.3e %.3e %.3e %.3e ", v_gmuon, v_bsgnlo, v_deltarho, v_bsmumu, v_btaunu);

	printf("%.6f %.6f %.6f %.6f %.6f %.6f %.6f %.6f %.6f %.6f %.6f %.6f %.6f %.6f %.6f %.6f %.6f %.6f %.6f %.6f %.6f %.6f %.6f %.6f %.6f %.6f %.6f %.6f %.6f %.6f %.6f %.6E %.6f %.3E %.6E %.6E %s\n",
		Mcdm, Omega, Xf, findValW("Zn11"), findValW("Zn12"), findValW("Zn13"), findValW("Zn14"), findValW("Zt11"), findValW("Zt12"), findValW("At"), findValW("Ab"), findValW("Al"), mu, findValW("MG1"), findValW("MG2"), findValW("MNE1"), findValW("MNE2"), findValW("MNE3"), findValW("MNE4"), m_a, findValW("MSt1"), findValW("MSt2"), findValW("MSb1"), findValW("MSb2"), findValW("MSeL"), findValW("MSmL"), findValW("MSl1"),findValW("MSG"),findValW("MSuL"),findValW("MC1"), findValW("MC2"),nEvents, higgs_mass, higgs_width, xs_p, amplitude_p, cdmName); 
   
	return 0;
}
