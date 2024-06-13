/************************************************************************** 
*                         NSSS.c
* 
* Computes constant neutron star sequences with constant baryonic mass
* 
**************************************************************************/
//
#include <stdio.h>
#include <math.h>
#include <stddef.h>
#include <stdlib.h>
#include <string.h> 

#include "consts.h"
#include "struct.h"

#include "nrutil.h"
#include "equil.h"
#include "equil_util.h"
#include "findmodel.h"
#include "surface.h"
#include "stableorbit.h"
#include "interpol.h"


/* Main; where it all starts and ends */

int main(int argc, char **argv)     /* Number of command line arguments, Command line arguments */
{ NeutronStar star;
  EOS eos;
  NeutronStar starDM;
  EOS eosDM;  
  int i, ierr;
  double
    e_min,e_minDM,  e_max, e_max_mass, ratio_r = 1.0,
    e_center=10e15,  e_centerDM=10e15,                   /* central en. density */
    B,                            /* Quark Bag Constant */
    K=3.0,                        /* Second parameter in "quark" eos */
    spin_freq=100,                  /* Spin Frequency */
    Gamma_P=0.0;                      /* Gamma for polytropic EOS */  
                
  int j;

  int a = 0, numseq=2;
  int spin_lim = 0;
  float e_c[4], M_0[4];
  float M0, Mtot,Mstat, Rstat, Radius, freq,freqK,energy_value, temp_energy, ej,Volp;
  float maxmass, maxradius;   // Mass and radius of the maximum mass neutron star for an EOS
  float T, W;
  long int angmom;
  //double Kfreq, Kfreq_j;

  //andreas
  double lumi;
  double u_lumi;
  double poten[4];
  int call;
  double ratio_ch=0.005;
  double energy_min=0.2;//10^{15}$ g/$cm^3
  double e_ch=0.05;//10^{15}$ g/$cm^3
  FILE *out;

  //andreas

  FILE *fpointer;

  char eos_file[80] = "no EOS file specified";   /* EOS file name */
  char eos_fileDM[80] = "no EOS file specified";   /* EOS file name */
  char eos_type[80] = "tab";                     /* EOS type (poly or tab) */
  char eos_typeDM[80] = "DM";                     /* EOS type (poly or tab) */
  char data_dir[80] = "junk";                    /* Data output directory */
  char data_dirDM[80] = "junkDM";                    /* Data output directory */
  char filename[100] = "Many Fermion";


  //andreas

  /* READ IN THE COMMAND LINE OPTIONS */
  for(i=1;i<argc;i++) 
    if(argv[i][0]=='-'){
      switch(argv[i][1]){

      case 'b':
	/* IF A TABULATED EOS WAS CHOSEN, CHOOSE THE
	   NAME OF THE FILE */
	sscanf(argv[i+1],"%s",eos_file);
	break;
	
      case 'd':
	/* IF A TABULATED EOS WAS CHOSEN, CHOOSE THE
	   NAME OF THE FILE */
	sscanf(argv[i+1],"%s",eos_fileDM);
	break;
	
      case 'e':
	/* CHOOSE THE CENTRAL ENERGY DENSITY OF THE 
	   NEUTRON STAR (IN g/cm^3) */
	sscanf(argv[i+1],"%lf",&e_min);
	if(strcmp(eos_type,"poly")!=0)
	  e_min *= C*C*KSCALE;
	break;
	
      case 'c':
	/* CHOOSE THE CENTRAL ENERGY DENSITY OF THE 
	   NEUTRON STAR (IN g/cm^3) */
	sscanf(argv[i+1],"%lf",&e_minDM);
	if(strcmp(eos_typeDM,"poly")!=0)
	  e_minDM *= C*C*KSCALE;
	break;
      case 'r':
	/* IF A TABULATED EOS WAS CHOSEN, CHOOSE THE
	   NAME OF THE FILE */
	sscanf(argv[i+1],"%lf",&ratio_r);
	break;
      }
    }

  //strncat(filename, eos_file, 65);
  strncat(filename, ".txt", 65);
  //printf("%s\n", filename);
  fpointer = fopen(filename, "a");


  //strncat(filename, eos_file, 65);
  //strncat(filename, "_table.txt", 65);
  //printf("%s\n", filename);
  out = fopen(filename, "a"); 
  /* PRINT THE HEADER */
  if(strcmp(eos_type,"tab")==0)
    printf("%s,  MDIVxSDIV=%dx%d\n",eos_file,MDIV,SDIV);
  if(strcmp(eos_type,"quark")==0)
    printf("Quark star with B=%f, MDIVxSDIV=%dx%d\n",B/1.602e33/KSCALE,MDIV,SDIV);

  /* SetUpStar loads in the eos and sets up the grid */
  /* Source code for SetUpStar can be found in findmodel.c */





  //printf("The star infrastructure has been set up! \n");

  e_center = e_min;
  e_centerDM = e_minDM;
  //temp_energy = e_center;

  // Computing the star with the maximum mass and its corresponding radius
    //Computing the non-rotating spherical neutron star
      ierr = SetUpStar(eos_file, eos_type, eos_fileDM, eos_typeDM, data_dir,data_dirDM, Gamma_P, B, K,&eos, &star,&eosDM, &starDM);

  /*
  e_center=1;
   
   //while(e_center>0.45){
      e_centerDM=0.00185;
   while(e_centerDM>0.00001){
  */
   ratio_r=1.00;
   ierr = MakeSphere(&eos, &star, e_center,&eosDM, &starDM, e_centerDM); 
      /*if(star.r_eDM==-1){
        continue;
      }*/
  printf("Rotating NSs\n" );
     printf(" \n" );
        printf("DM_file  e_c e_cDM MassTot MassBM MassDM RadiusBM RadiusDM Rrat_S_BM RratioBM RratioDM  Omega OmegaDM  Omega_K \n" );
while(ratio_r>=0.7){


      rns(ratio_r, e_center, &eos, &star, e_centerDM, &eosDM, &starDM);  
 
   //if(ratio_r==1.00) printf(" %f \n",star.MassDM/MSUN/(star.Mass/MSUN+star.MassDM/MSUN));
  if(!((isnan((star.Mass))&&isnan((star.MassDM)))||(isnan((star.OmegaDM))&&isnan((star.Omega))))){
  //if(star.R_eDM<star.R_e) exit(0);
  /*call = Surface(&eos,&star,&poten[0],&poten[1],&poten[2],&poten[3]);
    for(i=1; i<= MDIV; i++){
        printf("%.5f %.5f %.5f %.5f %.5f %.5f %.5f\n", acos((i-1.0)/(MDIV-1.0)),star.r_surf[i],star.gravity_surf[i],star.Mass/MSUN+star.MassDM/MSUN, star.R_e*1e-5,ratio_r, star.Omega/(2.0*PI));
    }*/
      printf("%s %g %g %.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f %.3f %.3f %.5f\n",
       eos_fileDM, star.e_center, star.e_centerDM, star.Mass/MSUN+star.MassDM/MSUN, star.Mass/MSUN, star.MassDM/MSUN, star.R_e*1e-5, star.R_eDM*1e-5, star.Ratio_sch,ratio_r, star.r_ratioDM, star.Omega/(2.0*PI),star.OmegaDM/(2.0*PI), star.Omega_K/(2.0*PI));
           
      fprintf(fpointer,"%s %g %g %.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f %.3f %.3f %.5f\n",
       eos_fileDM, star.e_center, star.e_centerDM, star.Mass/MSUN+star.MassDM/MSUN, star.Mass/MSUN, star.MassDM/MSUN, star.R_e*1e-5, star.R_eDM*1e-5, star.Ratio_sch,ratio_r, star.r_ratioDM, star.Omega/(2.0*PI),star.OmegaDM/(2.0*PI), star.Omega_K/(2.0*PI));
     
  //if(star.MassDM/(star.Mass+star.MassDM)<0.005) exit(0);

 }
 /*fprintf(fpointer, "%g %g  %.5f  %.5f  %.5f %.5f %.3f %.5f %.3f %.5f\n", 
              star.e_center, star.e_centerDM, star.Mass/MSUN+star.MassDM/MSUN, star.Mass/MSUN, star.MassDM/MSUN, star.R_e*1e-5, star.R_eDM*1e-5, ratio_r, star.Omega/(2.0*PI), star.Omega_K/(2.0*PI));
*/      

  ratio_r-=0.1;       

   }
  /*   e_centerDM-=0.00001;
   }*/
     //e_center-=0.05;
  //}


/*
 
  // Computing one neutron star
  if(1){

    
    ratio_r=1.00;
    double MaxM;
    //Computing the non-rotating spherical neutron star
    ierr = SetUpStar(eos_file, eos_type, data_dir, Gamma_P, B, K,&eos, &star);
    ierr = MakeSphere(&eos, &star, e_center); 
    rns(1.00, e_center, &eos, &star); 
    MaxM = star.Mass/MSUN;
    

    
    while(e_center>energy_min){

      //Computing the non-rotating spherical neutron star
      ierr = SetUpStar(eos_file, eos_type, data_dir, Gamma_P, B, K,&eos, &star);
      ierr = MakeSphere(&eos, &star, e_center); 
      rns(1.00, e_center, &eos, &star); 
      Mstat = star.Mass/MSUN;
      Rstat = star.R_e*1e-5;   

      // Computing the one star with given value of ratio_r

      while(star.Omega/(2.0*PI)<star.Omega_K/(2.0*PI)){


	rns(ratio_r, e_center, &eos, &star); 


	//andreas
	call = Surface(&eos,&star,&poten[0],&poten[1],&poten[2],&poten[3]);


	W = star.Mp + T - star.Mass;

	Mtot=star.Mass/MSUN;
	M0=star.Mass_0/MSUN;
	Radius=star.R_e*1e-5;
	freq=star.Omega/(2.0*PI);
	freqK=star.Omega_K/(2.0*PI);
	angmom=star.ang_mom;
	Volp=star.Vp;


          printf("%g %.5f  %.5f  %.5f %.5f %.3f %.5f %.3f %.5f\n",
            star.e_center, star.Mass/MSUN, star.Mass_0/MSUN, Mstat, star.R_e*1e-5, ratio_r, Rstat, star.Omega/(2.0*PI), star.Omega_K/(2.0*PI));
           fprintf(fpointer, "%7g %7g %7g %6g %8g %4g %6g %7g %6g  %g  %g  %g %g\n", 
             star.e_center, star.Mass/MSUN, star.Mass_0/MSUN, Mstat, star.R_e*1e-5, ratio_r, Rstat, star.Omega/(2.0*PI), star.Omega_K/(2.0*PI), star.ang_mom, T/MSUN, W/MSUN,star.r_surf[MDIV]/star.r_surf[1]);
             //fprintf(fpointer, "%7g %7g %7g %6g %f %8g %4g %6g %7g %6g  %g  %g  %g %f %f %f\n", 
             //star.e_center, star.Mass/MSUN, star.Mass_0/MSUN, Mstat, 0.0, star.R_e*1e-5, ratio_r, Rstat, star.Omega/(2.0*PI), star.Omega_K/(2.0*PI), star.ang_mom, T,W, 0.0,0.0,0.0 );
	if(ratio_r==1.00){
             fprintf(out, "%s %g %s %g %s %g %s \n","&",star.e_center,"&",star.Mass/MSUN,"&",star.R_e*1e-5,"\\\\");
	}

	ratio_r-=ratio_ch; 
      }
      //break;
      ratio_r=1.00;
      e_center-=e_ch;
    }

    //andreas
  }
*/

  fclose(out);//andreas
  fclose(fpointer);
  return 0;
}









