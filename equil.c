/*****************************************************************************
*	equil.c
*
*		The code in this file is a set of procedures written by
*	Nikolaos Stergioulas. These are the procedures used to integrate
*	the field equations for a rapidly rotating neutron star.
*
* 	The most important procedures are:
*	
*	make_grid:	Create the MDIV x SDIV grid. 
*			MDIV = number of divisions of variable mu=cos theta
*			SDIV = number of divisions of radial variable s
*	load_eos:	Load the equation of state file
*	make_center:	Calculate the central pressure and enthalpy
*	sphere:		Compute the metric of a spherical star 
*	TOV:		Integrates the Tolman-Oppenheimer-Volkoff
*			equations for spherically symmetric star
*	spin:		Integrates the equations for a rapidly rotating
*			neutron star with oblateness = r_ratio = 
*				radius of pole/radius of equator
*	mass_radius:	Calculates the gravitational mass and equatorial
*			radius of the rotating star, along with other
*			equilibrium quantities. 
*
******************************************************************************/



#include <stdio.h>
#include <string.h> 
#include <math.h>
#include "equil_util.h"
#include "consts.h"
#include "nrutil.h"
#include "equil.h"
#include <stdbool.h>

#define C 2.9979e10                  /* speed of light in vacuum */
#define G 6.6732e-8                  /* gravitational constant */ 
#define KAPPA 1.346790806509621e+13  /* square of length scale = 1e-15*C*C/G */
#define KSCALE 1.112668301525780e-36 /* KAPPA*G/(C*C*C*C) */  
#define MSUN 1.987e33                /* Mass of Sun */
#define PI 3.1415926535  

/*******************************************************************/
/* Create computational grid.                                      */
/* Points in the mu-direction are stored in the array mu[i].       */
/* Points in the s-direction are stored in the array s_gp[j].      */
/*******************************************************************/

double max(double v1,double v2){
     if(v1>=v2){
     return v1;
     }
     else{ 
     return v2;
     }
}


void make_grid(double s_gp[SDIV+1], 
               double mu[MDIV+1])                        
{ 
  int m, s;                         /* counters */
    
      for(s=1;s<=SDIV;s++) 
         s_gp[s] = SMAX*(s-1.0)/(SDIV-1.0);

	/* s_gp[1] = 0.0     corresponds to the center of the star
	   s_gp[SDIV] = SMAX corresponds to infinity */

	/* SMAX is defined in the file consts.h */

      for(m=1;m<=MDIV;m++) 
         mu[m] = (m-1.0)/(MDIV-1.0);

	/* mu[1] = 0.0    corresponds to the plane of the equator 
	   mu[MDIV] = 1.0 corresponds to the axis of symmetry */

	/* s_gp[0] and mu[0] are not used by the program */

}



/*************************************************************************/
/* Load EOS file.                                                        */ 
/*************************************************************************/
void load_eos( char eos_file[], 
               double log_e_tab[2001], 
               double log_p_tab[2001], 
               double log_h_tab[2001],
               double log_n0_tab[2001], 
               int *n_tab)
{
 int i;                    /* counter */

 double p,                 /* pressure */
        rho,               /* density */
        h,                 /* enthalpy */
        n0;                /* number density */    
        //g;                 /* Gamma */

 FILE *f_eos;              /* pointer to eos_file */
  

    /* OPEN FILE TO READ */

    if((f_eos=fopen(eos_file,"r")) == NULL ) {    
       printf("cannot open file:  %s\n",eos_file); 
       exit(0);
    }

 
    /* READ NUMBER OF TABULATED POINTS */

    fscanf(f_eos,"%d\n",n_tab);


    /* READ EOS, H, N0 AND MAKE THEM DIMENSIONLESS */
 
    for(i=1;i<=(*n_tab);i++) {  
      /*fscanf(f_eos,"%lf %lf %lf %lf %lf\n",&rho,&p,&h,&n0,&g) ; */
       fscanf(f_eos,"%lf %lf %lf %lf\n",&rho,&p,&h,&n0) ;
       log_e_tab[i]=log10(rho*C*C*KSCALE);     /* multiply by C^2 to get */ 
       log_p_tab[i]=log10(p*KSCALE);           /* energy density. */
       log_h_tab[i]=log10(h/(C*C));        
       log_n0_tab[i]=log10(n0);
       /*Gamma_tab[i]=g;*/
    }
}



/*******************************************************************/
double e_of_rho0(double rho0, double Gamma_P)
{
 return(pow(rho0,Gamma_P)/(Gamma_P-1.0)+rho0);
}
   

/*C*/
/*******************************************************************/
double e_at_p(double pp, 
              double pp_surface,
              double log_e_tab[2001], 
              double log_p_tab[2001],
              int    n_tab, 
              int    *n_nearest_pt,
              char eos_type[],
              double Gamma_P)
{
 if((strcmp(eos_type,"tab")==0) || (strcmp(eos_type,"DM")==0)){
   if(pp<pp_surface){
     return 0;
   }else{
     return pow(10.0,interp(log_p_tab,log_e_tab,n_tab,log10(pp), n_nearest_pt));
   }
 }
}

/*C*/
/*******************************************************************/
double p_at_e(double ee,
              double ee_surface, 
              double log_p_tab[2001], 
              double log_e_tab[2001],
              int    n_tab, 
              int    *n_nearest_pt)
{ 
 if(ee<ee_surface){
  return 0;
 }else{
 return pow(10.0,interp(log_e_tab,log_p_tab,n_tab,log10(ee), n_nearest_pt));
 }
} 

/*C*/
/*******************************************************************/
double p_at_h(double hh,
              double hh_surface, 
              double log_p_tab[2001], 
              double log_h_tab[2001],
              int    n_tab, 
              int    *n_nearest_pt)
{

 if(hh<hh_surface){
  return 0;
 }else{
 return pow(10.0,interp(log_h_tab,log_p_tab,n_tab,log10(hh), n_nearest_pt));
 }
}

/*C*/
/*******************************************************************/
double h_at_p(double pp, 
              double pp_surface,
              double log_h_tab[2001], 
              double log_p_tab[2001],
              int    n_tab, 
              int    *n_nearest_pt)
{

 if(pp<pp_surface){
  return 0;
 }else{
 return pow(10.0,interp(log_p_tab,log_h_tab,n_tab,log10(pp), n_nearest_pt));
 }
}

/*C*/
/*******************************************************************/
double n0_at_e(double ee, 
	       double ee_surface,
               double log_n0_tab[2001], 
               double log_e_tab[2001],
               int    n_tab, 
               int    *n_nearest_pt)
{
 if(ee<ee_surface){
  return 0;
 }else{
  return pow(10.0,interp(log_e_tab,log_n0_tab,n_tab,log10(ee), n_nearest_pt));
 }
}
 
/*C*/
/***************************************************************/
void make_center(
	       char eos_file[], 
               double log_e_tab[2001], 
               double log_p_tab[2001], 
               double log_h_tab[2001],
               double log_n0_tab[2001], 
               int n_tab,                 
	       char eos_type[],
	       double Gamma_P, 
	       double e_center,
	       double *p_center, 
	       double *h_center,
	       double e_surface,
	       double p_surface)

{

 int n_nearest;

 double rho0_center;

 n_nearest=n_tab/2; 


   (*p_center) = p_at_e( e_center, e_surface, log_p_tab, log_e_tab, n_tab, &n_nearest);
   (*h_center) = h_at_p( (*p_center), p_surface, log_h_tab, log_p_tab, n_tab, &n_nearest);
 



}

/*C*/
/***********************************************************************/
/* Computes the gravitational mass, equatorial radius, angular momentum
 *	of the star
 * 	and the velocity of co- and counter-rotating particles      
 *	with respect to a ZAMO                                         */
/***********************************************************************/
void mass_radius(
		 double s_gp[SDIV+1],
		 double mu[MDIV+1],
		 double log_e_tab[2001], 
		 double log_p_tab[2001], 
		 double log_h_tab[2001],
		 double log_n0_tab[2001], 
		 int n_tab,                 
		 char eos_type[],
		 double Gamma_P, 
		 double log_e_tabDM[2001], 
		 double log_p_tabDM[2001], 
		 double log_h_tabDM[2001],
		 double log_n0_tabDM[2001], 
		 int n_tabDM,                 
		 char eos_typeDM[],
		 double Gamma_PDM, 
		 double **rho,
		 double **gama,
		 double **alpha,
		 double **omega,
		 double **energy,
		 double **pressure,
		 double **enthalpy,
		 double **velocity_sq,
		 double **energyDM,
		 double **pressureDM,
		 double **enthalpyDM,
		 double **velocity_sqDM,
                 double r_ratio,
                 double r_ratioDM,
                 double *Ratio_sch,
		 double e_surface,
		 double e_surfaceDM,
                 double r_e,
                 double r_eDM,
                 double Omega,
                 double *Mass, 
		 double *Mass_0,
		 double *ang_mom,
                 double *R_e,
                 double *MassDM, 
		 double *Mass_0DM,
		 double *ang_momDM,
                 double *R_eDM,
		 double *v_plus,
		 double *v_minus,
		 double *Omega_K,
		 double *Vp,
     double *Mp)

{


 int s,
     m,
     n_nearest,
     n_nearestDM;

int index;

 
 double   
   **rho_0, /*rest mass density*/
   **velocity,
   gama_equator,              /* gama at equator */
   rho_equator,               /* rho at equator */
   omega_equator,             /* omega at equator */
   s1,
   s_1,r_eDM_old,
   d_gama_s,
   d_rho_s,
   d_omega_s,
   sqrt_v,
   D_m[SDIV+1],               /* int. quantity for M */
   D_m_0[SDIV+1],             /* int. quantity for M_0 */ 
   D_J[SDIV+1],               /* int. quantity for J */
   D_mDM[SDIV+1],               /* int. quantity for M */
   D_m_0DM[SDIV+1],             /* int. quantity for M_0 */ 
   D_JDM[SDIV+1],               /* int. quantity for J */
   s_e,                 
   d_o_e[SDIV+1],
   d_g_e[SDIV+1],
   d_r_e[SDIV+1],
   d_v_e[SDIV+1],
   doe,
   dge, 
   dre,
   dve,
   vek, 
   ratio_old,ratio_oldDM,    
   gama_mu_0[SDIV+1],                   
   rho_mu_0[SDIV+1], 
   gama_mu_1[SDIV+1],                   
   rho_mu_1[SDIV+1],                      
   omega_mu_0[SDIV+1],
   J,
   JDM,
   r_p,
   s_p,
   D_vp[SDIV+1],
   D_mp[SDIV+1],
   Rv,r_eq;        
   r_eDM_old=r_eDM;
   /* Circumferential radius */
   bool cond=(r_e>=r_eDM);

   ratio_old=r_ratio;
   ratio_oldDM=r_ratioDM;

   r_ratio = ((cond)?ratio_old:ratio_oldDM);
   r_ratioDM=((!cond)?ratio_old:ratio_oldDM);


   r_eq=((cond)?r_e:r_eDM);
   r_eDM=((!cond)?r_e:r_eDM);

   r_p= r_ratio*r_eq;                              /* radius at pole */
   s_p= r_p/(r_p+r_eq);                            /* s-coordinate at pole */
   s_e=0.5;

   rho_0 = dmatrix(1,SDIV,1,MDIV);
   velocity = dmatrix(1,SDIV,1,MDIV);

   for(s=1;s<=SDIV;s++) {               
      gama_mu_0[s]=gama[s][1];                   
      rho_mu_0[s]=rho[s][1];                                                    
      gama_mu_1[s]=gama[s][MDIV];                   
      rho_mu_1[s]=rho[s][MDIV];        
   }

   n_nearest= SDIV/2;
   gama_equator=interp(s_gp,gama_mu_0,SDIV,s_e, &n_nearest);  
   rho_equator=interp(s_gp,rho_mu_0,SDIV,s_e, &n_nearest);   

/* Circumferential radius */

 /*(*R_e) = sqrt(KAPPA)*r_eq*exp((gama_equator-rho_equator)/2.0);
 (*R_eDM) = sqrt(KAPPA)*r_eDM*exp((interp(s_gp,gama_mu_0,SDIV,(r_eDM/r_eq)/(1.+(r_eDM/r_eq)), &n_nearest) -interp(s_gp,rho_mu_0,SDIV,(r_eDM/r_eq)/(1.+(r_eDM/r_eq)), &n_nearest))/2.0);*/
 if(cond){
 

 (*R_e) = sqrt(KAPPA)*r_e*exp((gama_equator-rho_equator)/2.0);
 (*R_eDM) = sqrt(KAPPA)*r_eDM_old*exp((interp(s_gp,gama_mu_0,SDIV,(r_eDM_old/r_e)/(1.+(r_eDM_old/r_e)), &n_nearest) -interp(s_gp,rho_mu_0,SDIV,(r_eDM_old/r_e)/(1.+(r_eDM_old/r_e)), &n_nearest))/2.0);
 double R_p=sqrt(KAPPA)*r_p*exp((interp(s_gp,gama_mu_1,SDIV,s_p, &n_nearest)-interp(s_gp,rho_mu_1,SDIV,s_p, &n_nearest))/2.0);
 (*Ratio_sch)=R_p/(*R_e);
 }else{
 (*R_eDM) = sqrt(KAPPA)*r_eDM_old*exp((gama_equator-rho_equator)/2.0);
 (*R_e) = sqrt(KAPPA)*r_e*exp((interp(s_gp,gama_mu_0,SDIV,(r_e)/(r_eDM_old+(r_e)), &n_nearest) -interp(s_gp,rho_mu_0,SDIV,(r_e)/(r_eDM_old+(r_e)), &n_nearest))/2.0); 
 
  double R_p=sqrt(KAPPA)*ratio_old*r_e*exp((interp(s_gp,gama_mu_1,SDIV,(ratio_old*r_e/r_eDM_old)/(1.+(ratio_old*r_e/r_eDM_old)), &n_nearest)-interp(s_gp,rho_mu_1,SDIV,(ratio_old*r_e/r_eDM_old)/(1.+(ratio_old*r_e/r_eDM_old)), &n_nearest))/2.0);
 (*Ratio_sch)=R_p/(*R_e);
 }

 /* Masses and angular momentum */
 
   (*Mass) = 0.0;              /* initialize */
   (*Mass_0) = 0.0;
   (*MassDM) = 0.0;              /* initialize */
   (*Mass_0DM) = 0.0;
   (*Vp) = 0.0;
   (*Mp) = 0.0;
   J=0.0;
   JDM=0.0;
   Rv=0.0;
   /* CALCULATE THE REST MASS DENSITY */
 if((strcmp(eos_type,"tab")==0) || (strcmp(eos_type,"DM")==0)) {
   n_nearest=n_tab/2;
   for(s=1;s<=SDIV;s++)
      for(m=1;m<=MDIV;m++) {
           if(energy[s][m]>e_surface)
             rho_0[s][m]=n0_at_e(energy[s][m],e_surface, log_n0_tab, log_e_tab, n_tab,
                                             &n_nearest)*MB*KSCALE*SQ(C);
           else
             rho_0[s][m]=0.0;
      }  
 }
 double s_eDM;
   for(s=1;s<=SDIV;s++) {
    D_m[s]=0.0;           /* initialize */
    D_mDM[s]=0.0;           /* initialize */
    D_m_0[s]=0.0;
    D_J[s]=0.0;
    D_JDM[s]=0.0;
    D_vp[s]=0.0;
    D_mp[s]=0.0;
    
    for(m=1;m<=MDIV-2;m+=2) {

    if(1){
     D_m[s] += (1.0/(3.0*(MDIV-1)))*( exp(2.0*alpha[s][m]+gama[s][m])*
              (((energy[s][m]+pressure[s][m])/(1.0-velocity_sq[s][m]))*
              (1.0+velocity_sq[s][m]+(2.0*s_gp[s]*sqrt(velocity_sq[s][m])/
              (1.0-s_gp[s]))*sqrt(1.0-mu[m]*mu[m])*r_eq*omega[s][m]*
              exp(-rho[s][m])) + 2.0*pressure[s][m])

            + 4.0*exp(2.0*alpha[s][m+1]+gama[s][m+1])*
              (((energy[s][m+1]+pressure[s][m+1])/(1.0-velocity_sq[s][m+1]))*
              (1.0+velocity_sq[s][m+1]+(2.0*s_gp[s]*sqrt(velocity_sq[s][m+1])/
              (1.0-s_gp[s]))*sqrt(1.0-mu[m+1]*mu[m+1])*r_eq*omega[s][m+1]*
              exp(-rho[s][m+1])) + 2.0*pressure[s][m+1]) 

            + exp(2.0*alpha[s][m+2]+gama[s][m+2])*
              (((energy[s][m+2]+pressure[s][m+2])/(1.0-velocity_sq[s][m+2]))*
              (1.0+velocity_sq[s][m+2]+(2.0*s_gp[s]*sqrt(velocity_sq[s][m+2])/
              (1.0-s_gp[s]))*sqrt(1.0-mu[m+2]*mu[m+2])*r_eq*omega[s][m+2]*
              exp(-rho[s][m+2])) + 2.0*pressure[s][m+2]));


    D_mDM[s] += (1.0/(3.0*(MDIV-1)))*( exp(2.0*alpha[s][m]+gama[s][m])*
              (((energyDM[s][m]+pressureDM[s][m])/(1.0-velocity_sqDM[s][m]))*
              (1.0+velocity_sqDM[s][m]+(2.0*s_gp[s]*sqrt(velocity_sqDM[s][m])/
              (1.0-s_gp[s]))*sqrt(1.0-mu[m]*mu[m])*r_eq*omega[s][m]*
              exp(-rho[s][m])) + 2.0*pressureDM[s][m])

            + 4.0*exp(2.0*alpha[s][m+1]+gama[s][m+1])*
              (((energyDM[s][m+1]+pressureDM[s][m+1])/(1.0-velocity_sqDM[s][m+1]))*
              (1.0+velocity_sqDM[s][m+1]+(2.0*s_gp[s]*sqrt(velocity_sqDM[s][m+1])/
              (1.0-s_gp[s]))*sqrt(1.0-mu[m+1]*mu[m+1])*r_eq*omega[s][m+1]*
              exp(-rho[s][m+1])) + 2.0*pressureDM[s][m+1]) 

            + exp(2.0*alpha[s][m+2]+gama[s][m+2])*
              (((energyDM[s][m+2]+pressureDM[s][m+2])/(1.0-velocity_sqDM[s][m+2]))*
              (1.0+velocity_sqDM[s][m+2]+(2.0*s_gp[s]*sqrt(velocity_sqDM[s][m+2])/
              (1.0-s_gp[s]))*sqrt(1.0-mu[m+2]*mu[m+2])*r_eq*omega[s][m+2]*
              exp(-rho[s][m+2])) + 2.0*pressureDM[s][m+2]));
              //printf("%d %d %f %f %f\n",s,m, energyDM[s][m+2]+pressureDM[s][m+2], energyDM[s][m+1]+pressureDM[s][m+1],energyDM[s][m]+pressureDM[s][m]); 
      }else{

     D_mDM[s] += (1.0/(3.0*(MDIV-1)))*( exp(2.0*alpha[s][m]+gama[s][m])*
              (((energyDM[s][m]+pressureDM[s][m])/(1.0-velocity_sq[s][m]))*
              (1.0+velocity_sq[s][m]+(2.0*s_gp[s]*sqrt(velocity_sq[s][m])/
              (1.0-s_gp[s]))*sqrt(1.0-mu[m]*mu[m])*r_eq*omega[s][m]*
              exp(-rho[s][m])) + 2.0*pressureDM[s][m])

            + 4.0*exp(2.0*alpha[s][m+1]+gama[s][m+1])*
              (((energyDM[s][m+1]+pressureDM[s][m+1])/(1.0-velocity_sq[s][m+1]))*
              (1.0+velocity_sq[s][m+1]+(2.0*s_gp[s]*sqrt(velocity_sq[s][m+1])/
              (1.0-s_gp[s]))*sqrt(1.0-mu[m+1]*mu[m+1])*r_eq*omega[s][m+1]*
              exp(-rho[s][m+1])) + 2.0*pressureDM[s][m+1]) 

            + exp(2.0*alpha[s][m+2]+gama[s][m+2])*
              (((energyDM[s][m+2]+pressureDM[s][m+2])/(1.0-velocity_sq[s][m+2]))*
              (1.0+velocity_sq[s][m+2]+(2.0*s_gp[s]*sqrt(velocity_sq[s][m+2])/
              (1.0-s_gp[s]))*sqrt(1.0-mu[m+2]*mu[m+2])*r_eq*omega[s][m+2]*
              exp(-rho[s][m+2])) + 2.0*pressureDM[s][m+2]));


    D_m[s] += (1.0/(3.0*(MDIV-1)))*( exp(2.0*alpha[s][m]+gama[s][m])*
              (((energy[s][m]+pressure[s][m])/(1.0-velocity_sqDM[s][m]))*
              (1.0+velocity_sqDM[s][m]+(2.0*s_gp[s]*sqrt(velocity_sqDM[s][m])/
              (1.0-s_gp[s]))*sqrt(1.0-mu[m]*mu[m])*r_eq*omega[s][m]*
              exp(-rho[s][m])) + 2.0*pressure[s][m])

            + 4.0*exp(2.0*alpha[s][m+1]+gama[s][m+1])*
              (((energy[s][m+1]+pressure[s][m+1])/(1.0-velocity_sqDM[s][m+1]))*
              (1.0+velocity_sqDM[s][m+1]+(2.0*s_gp[s]*sqrt(velocity_sqDM[s][m+1])/
              (1.0-s_gp[s]))*sqrt(1.0-mu[m+1]*mu[m+1])*r_eq*omega[s][m+1]*
              exp(-rho[s][m+1])) + 2.0*pressure[s][m+1]) 

            + exp(2.0*alpha[s][m+2]+gama[s][m+2])*
              (((energy[s][m+2]+pressure[s][m+2])/(1.0-velocity_sqDM[s][m+2]))*
              (1.0+velocity_sqDM[s][m+2]+(2.0*s_gp[s]*sqrt(velocity_sqDM[s][m+2])/
              (1.0-s_gp[s]))*sqrt(1.0-mu[m+2]*mu[m+2])*r_eq*omega[s][m+2]*
              exp(-rho[s][m+2])) + 2.0*pressure[s][m+2])); 
      
      
      }



     D_m_0[s] += (1.0/(3.0*(MDIV-1)))*( exp(2.0*alpha[s][m]+(gama[s][m]
              -rho[s][m])/2.0)*rho_0[s][m]/sqrt(1.0-velocity_sq[s][m])

             + 4.0* exp(2.0*alpha[s][m+1]+(gama[s][m+1]
             -rho[s][m+1])/2.0)*rho_0[s][m+1]/sqrt(1.0-velocity_sq[s][m+1])
         
             + exp(2.0*alpha[s][m+2]+(gama[s][m+2]
             -rho[s][m+2])/2.0)*rho_0[s][m+2]/sqrt(1.0-velocity_sq[s][m+2]));

     D_J[s] += (1.0/(3.0*(MDIV-1)))*( sqrt(1.0-mu[m]*mu[m])*
              exp(2.0*alpha[s][m]+gama[s][m]-rho[s][m])*(energy[s][m]
              +pressure[s][m])*sqrt(velocity_sq[s][m])/(1.0-velocity_sq[s][m])
  
              +4.0*sqrt(1.0-mu[m+1]*mu[m+1])*
              exp(2.0*alpha[s][m+1]+gama[s][m+1]-rho[s][m+1])*(energy[s][m+1]
              +pressure[s][m+1])*sqrt(velocity_sq[s][m+1])/
              (1.0-velocity_sq[s][m+1])

              + sqrt(1.0-mu[m+2]*mu[m+2])*
              exp(2.0*alpha[s][m+2]+gama[s][m+2]-rho[s][m+2])*(energy[s][m+2]
              +pressure[s][m+2])*sqrt(velocity_sq[s][m+2])/
              (1.0-velocity_sq[s][m+2]));


     D_JDM[s] += (1.0/(3.0*(MDIV-1)))*( sqrt(1.0-mu[m]*mu[m])*
              exp(2.0*alpha[s][m]+gama[s][m]-rho[s][m])*(energyDM[s][m]
              +pressureDM[s][m])*sqrt(velocity_sqDM[s][m])/(1.0-velocity_sqDM[s][m])
  
              +4.0*sqrt(1.0-mu[m+1]*mu[m+1])*
              exp(2.0*alpha[s][m+1]+gama[s][m+1]-rho[s][m+1])*(energyDM[s][m+1]
              +pressureDM[s][m+1])*sqrt(velocity_sqDM[s][m+1])/
              (1.0-velocity_sqDM[s][m+1])

              + sqrt(1.0-mu[m+2]*mu[m+2])*
              exp(2.0*alpha[s][m+2]+gama[s][m+2]-rho[s][m+2])*(energyDM[s][m+2]
              +pressureDM[s][m+2])*sqrt(velocity_sqDM[s][m+2])/
              (1.0-velocity_sqDM[s][m+2]));

     D_mp[s] += (1.0/(3.0*(MDIV-1)))*( exp(2.0*alpha[s][m]+(gama[s][m]
              -rho[s][m])/2.0)*((energy[s][m]+pressure[s][m])/sqrt(1.0-velocity_sq[s][m]))

             + 4.0* exp(2.0*alpha[s][m+1]+(gama[s][m+1]
             -rho[s][m+1])/2.0)*((energy[s][m+1]+pressure[s][m+1])/sqrt(1.0-velocity_sq[s][m+1]))
         
             + exp(2.0*alpha[s][m+2]+(gama[s][m+2]
             -rho[s][m+2])/2.0)*((energy[s][m+2]+pressure[s][m+2])/sqrt(1.0-velocity_sq[s][m+2])));

     if(energy[s][m+1]+pressure[s][m+1]==0){
      D_vp[s]+=0.0;
     }else{
 
     D_vp[s] += (1.0/(3.0*(MDIV-1)))*( exp(2.0*alpha[s][m]+(gama[s][m]
              -rho[s][m])/2.0)*(1.0/sqrt(1.0-velocity_sq[s][m]))

             + 4.0* exp(2.0*alpha[s][m+1]+(gama[s][m+1]
             -rho[s][m+1])/2.0)*(1.0/sqrt(1.0-velocity_sq[s][m+1]))
         
             + exp(2.0*alpha[s][m+2]+(gama[s][m+2]
             -rho[s][m+2])/2.0)*(1.0/sqrt(1.0-velocity_sq[s][m+2])));
             Rv= s_gp[s]*r_eq*sqrt(KAPPA)/((1.-s_gp[s])*100000.);
    }    

    }
   }
   index = 1;
    for(s=1;s<=SDIV-2;s+=2) { 
        //printf("%f %f %f ",log(velocity_sq[s][0]),energy[s][0],omega[s][0]);
     (*Mass) += (SMAX/(3.0*(SDIV-1)))*(pow(sqrt(s_gp[s])/(1.0-s_gp[s]),4.0)*
          D_m[s]+4.0*pow(sqrt(s_gp[s+1])/(1.0-s_gp[s+1]),4.0)*D_m[s+1]
          +pow(sqrt(s_gp[s+2])/(1.0-s_gp[s+2]),4.0)*D_m[s+2]);

     (*MassDM) += (SMAX/(3.0*(SDIV-1)))*(pow(sqrt(s_gp[s])/(1.0-s_gp[s]),4.0)*
          D_mDM[s]+4.0*pow(sqrt(s_gp[s+1])/(1.0-s_gp[s+1]),4.0)*D_mDM[s+1]
          +pow(sqrt(s_gp[s+2])/(1.0-s_gp[s+2]),4.0)*D_mDM[s+2]);
          
     //printf("%d \t %g %g %g\n", index, (SMAX/(3.0*(SDIV-1)))*(pow(sqrt(s_gp[s])/(1.0-s_gp[s]),4.0)*
      //    D_mDM[s],4.0*pow(sqrt(s_gp[s+1])/(1.0-s_gp[s+1]),4.0)*D_mDM[s+1],pow(sqrt(s_gp[s+2])/(1.0-s_gp[s+2]),4.0)*D_mDM[s+2])); 

     //printf("%d \t %g\n", index, (*Mass));      /***************************/

     (*Mass_0) += (SMAX/(3.0*(SDIV-1)))*(pow(sqrt(s_gp[s])/(1.0-s_gp[s]),4.0)*
          D_m_0[s]+4.0*pow(sqrt(s_gp[s+1])/(1.0-s_gp[s+1]),4.0)*D_m_0[s+1]
          +pow(sqrt(s_gp[s+2])/(1.0-s_gp[s+2]),4.0)*D_m_0[s+2]);
 
     J += (SMAX/(3.0*(SDIV-1)))*((pow(s_gp[s],3.0)/pow(1.0-s_gp[s],5.0))*
          D_J[s]+ 4.0*(pow(s_gp[s+1],3.0)/pow(1.0-s_gp[s+1],5.0))*
          D_J[s+1] + (pow(s_gp[s+2],3.0)/pow(1.0-s_gp[s+2],5.0))*
          D_J[s+2]);

     JDM += (SMAX/(3.0*(SDIV-1)))*((pow(s_gp[s],3.0)/pow(1.0-s_gp[s],5.0))*
          D_JDM[s]+ 4.0*(pow(s_gp[s+1],3.0)/pow(1.0-s_gp[s+1],5.0))*
          D_JDM[s+1] + (pow(s_gp[s+2],3.0)/pow(1.0-s_gp[s+2],5.0))*
          D_JDM[s+2]);

     (*Mp) += (SMAX/(3.0*(SDIV-1)))*(pow(sqrt(s_gp[s])/(1.0-s_gp[s]),4.0)*
          D_mp[s]+4.0*pow(sqrt(s_gp[s+1])/(1.0-s_gp[s+1]),4.0)*D_mp[s+1]
          +pow(sqrt(s_gp[s+2])/(1.0-s_gp[s+2]),4.0)*D_mp[s+2]);

     (*Vp) += (SMAX/(3.0*(SDIV-1)))*(pow(sqrt(s_gp[s])/(1.0-s_gp[s]),4.0)*
          D_vp[s]+4.0*pow(sqrt(s_gp[s+1])/(1.0-s_gp[s+1]),4.0)*D_vp[s+1]
          +pow(sqrt(s_gp[s+2])/(1.0-s_gp[s+2]),4.0)*D_vp[s+2]);
  
     index++;
    }
   
    if((strcmp(eos_type,"tab")==0) || (strcmp(eos_type,"DM")==0)) {
      (*Mass) *= 4.*PI*sqrt(KAPPA)*C*C*pow(r_eq,3.0)/G;
      (*MassDM) *= 4.*PI*sqrt(KAPPA)*C*C*pow(r_eq,3.0)/G;
      (*Mass_0) *= 4.*PI*sqrt(KAPPA)*C*C*pow(r_eq,3.0)/G;
      (*Vp) *= 4.*PI*sqrt(KAPPA)*sqrt(KAPPA)*sqrt(KAPPA)/3.0;
      (*Mp) *= 4.*PI*sqrt(KAPPA)*C*C*pow(r_eq,3.0)/G;
  
    }
     if(isnan((*Mass))&&isnan((*MassDM))){
       //printf("The system is unstable\n"); 
       //exit(0);
       return;
     }
    if(r_ratio==1.0) 
         J=0.0; 
    else {    
          if((strcmp(eos_type,"tab")==0) || (strcmp(eos_type,"DM")==0)) 
              J *= 4.0*PI*KAPPA*C*C*C*pow(r_eq,4.0)/G;
    }

    (*ang_mom) = J;

    if(r_ratioDM==1.0) 
         JDM=0.0; 
    else {    
          if((strcmp(eos_type,"tab")==0) || (strcmp(eos_type,"DM")==0)) 
              JDM *= 4.0*PI*KAPPA*C*C*C*pow(r_eq,4.0)/G;
    }

    (*ang_momDM) = JDM;

    //printf(" J = %g \n", J);


  /* Compute the velocities of co-rotating and counter-rotating particles
	with respect to a ZAMO 	*/

  for(s=1+(SDIV-1)/2;s<=SDIV;s++) {
    s1= s_gp[s]*(1.0-s_gp[s]);
    s_1=1.0-s_gp[s];
        
    d_gama_s=deriv_s(gama,s,1);
    d_rho_s=deriv_s(rho,s,1);
    d_omega_s=deriv_s(omega,s,1);

    sqrt_v= exp(-2.0*rho[s][1])*r_eq*r_eq*pow(s_gp[s],4.0)*pow(d_omega_s,2.0) 
            + 2.0*s1*(d_gama_s+d_rho_s)+s1*s1*(d_gama_s*d_gama_s-d_rho_s*d_rho_s);

    if(sqrt_v>0.0) sqrt_v= sqrt(sqrt_v);
     else {
      sqrt_v=0.0;
     }

    v_plus[s]=(exp(-rho[s][1])*r_eq*s_gp[s]*s_gp[s]*d_omega_s + sqrt_v)/
              (2.0+s1*(d_gama_s-d_rho_s));

    v_minus[s]=(exp(-rho[s][1])*r_eq*s_gp[s]*s_gp[s]*d_omega_s - sqrt_v)/
               (2.0+s1*(d_gama_s-d_rho_s));
  }


/* Kepler angular velocity */

   for(s=1;s<=SDIV;s++) { 
     d_o_e[s]=deriv_s(omega,s,1);
     d_g_e[s]=deriv_s(gama,s,1);
     d_r_e[s]=deriv_s(rho,s,1);
     d_v_e[s]=deriv_s(velocity,s,1);
     /* Value of omega on the equatorial plane*/
     omega_mu_0[s] = omega[s][1];
   }

   n_nearest=SDIV/2; 
   doe=interp(s_gp,d_o_e,SDIV,s_e, &n_nearest);
   dge=interp(s_gp,d_g_e,SDIV,s_e, &n_nearest);
   dre=interp(s_gp,d_r_e,SDIV,s_e, &n_nearest);
   dve=interp(s_gp,d_v_e,SDIV,s_e, &n_nearest);

  vek=(doe/(8.0+dge-dre))*r_eq*exp(-rho_equator) + sqrt(((dge+dre)/(8.0+dge
        -dre)) + pow((doe/(8.0+dge-dre))*r_eq*exp(-rho_equator),2.0));


  //if (r_ratio ==1.0)
  //  omega_equator = 0.0;
  //else
    omega_equator = interp(s_gp,omega_mu_0,SDIV,s_e, &n_nearest);




   (*Omega_K) = (C/sqrt(KAPPA))*(omega_equator+vek*exp(rho_equator)/r_eq);


   free_dmatrix(velocity,1,SDIV,1,MDIV);
   free_dmatrix(rho_0,1,SDIV,1,MDIV);

              //printf("%.5f\n",Rv);


}

/*C*/
/**************************************************************************/
double dm_dr_is(double r_is, 
                double r, 
                double m, 
                double p, 
                double e_center, 
                double p_surface,
                double log_e_tab[SDIV+1],
                double log_p_tab[SDIV+1],
                int    n_tab,
                int    *n_nearest_pt,
                char eos_type[],
                double Gamma_P)
{
 double dmdr,
        e_d;

 if(p<p_surface) 
    e_d=0.0;
 else  
    e_d = e_at_p(p,p_surface, log_e_tab, log_p_tab, n_tab, n_nearest_pt, eos_type, 
                                                                  Gamma_P);
 
 if(r_is<RMIN) 
    dmdr=4.0*PI*e_center*r*r*(1.0+4.0*PI*e_center*r*r/3.0);
 else
    dmdr=4.0*PI*e_d*r*r*r*sqrt(1.0-2.0*m/r)/r_is;
 
return dmdr;
}
 
/*C*/
/**************************************************************************/
double dp_dr_is(double r_is, 
                double r, 
                double m, 
                double p,
                double p_other, 
                double e_center, 
                double p_surface,
                double log_e_tab[SDIV+1],
                double log_p_tab[SDIV+1],
                int    n_tab,
                int    *n_nearest_pt,
                char eos_type[],
                double Gamma_P)
{ double dpdr,
         e_d; 

  if(p<p_surface){
   e_d=0.0;
   p=0.0; 
  }
  else        
   e_d=e_at_p(p,p_surface, log_e_tab, log_p_tab, n_tab, n_nearest_pt, eos_type, 
                                                                  Gamma_P);
  
  if(r_is<RMIN) dpdr = -4.0*PI*(e_center+p)*(e_center+3.0*p)*r*(1.0
                     +4.0*e_center*r*r/3.0)/3.0;

  else 
   dpdr = -(e_d+p)*(m+4.0*PI*r*r*r*(p+p_other))/(r*r_is*sqrt(1.0-2.0*m/r));

 return dpdr;
}


/**************************************************************************/
double dr_dr_is(double r_is, double r, double m)
{
 double drdris;

 if(r_is<RMIN) drdris=1.0;
  else
   drdris=(r/r_is)*sqrt(1.0-2.0*m/r);

 return drdris;
}

/*C*/
/************************************************************************/
void TOV(
	       int    i_check, 
               char   eos_type[],
               double e_center,
               double p_center,
               double p_surface,
               double e_surface,
               double Gamma_P, 
               double log_e_tab[2001],
               double log_p_tab[2001],
               double log_n0_tab[2001],
               double log_h_tab[2001],
               int    n_tab,
               double r_is_gp[RDIV+1], 
               double lambda_gp[RDIV+1], 
               char   eos_typeDM[],
               double e_centerDM,
               double p_centerDM,
               double p_surfaceDM,
               double e_surfaceDM,
               double Gamma_PDM, 
               double log_e_tabDM[2001],
               double log_p_tabDM[2001],
               double log_n0_tabDM[2001],
               double log_h_tabDM[2001],
               int    n_tabDM,
               double nu_gp[RDIV+1], 
               double *r_is_final, 
               double *r_final, 
               double *m_final,
               double *rDM_is_final, 
               double *rDM_final, 
               double *mDM_final)
{
  int i=2,
      n_nearest,n_nearestDM;

  double r,                           /* radius */
         r_is,                        /* isotropic radial coordinate */
         r_is_est,                    /* estimate on final isotr. radius */ 
         r_is_check,                  /*                      */    
         dr_is_save,                  /* r_is saving interval */  
         rho0,
         e_d,                         /* density */
         p,                           /* pressure */
         h,                           /* stepsize during integration */
         m,                           /* mass   */
         rho0DM,
         e_dDM,                         /* density */
         pDM,                           /* pressure */
         mDM,                           /* mass   */
         nu_s,
         hh,
         a1,a2,a3,a4,b1,b2,b3,b4,     /* coeff. in Runge-Kutta equations */
         c1,c2,c3,c4,
         a1DM,a2DM,a3DM,a4DM,b1DM,b2DM,b3DM,b4DM,     /* coeff. in Runge-Kutta equations */
         c1DM,c2DM,c3DM,c4DM,
         k_rescale, 
         r_gp[RDIV+1],
         m_gp[RDIV+1],
         e_d_gp[RDIV+1],
         p_d_gp[RDIV+1],
         m_gpDM[RDIV+1],
         e_d_gpDM[RDIV+1],
         p_d_gpDM[RDIV+1];   


    if(i_check==1) {
      if((strcmp(eos_type,"tab")==0) || (strcmp(eos_type,"DM")==0))
        r_is_est=1.5e6/sqrt(KAPPA);


      h=r_is_est/10000.0;     
    }
    else {
          r_is_est= max(*rDM_is_final,*r_is_final);
          h=r_is_est/100000.0;   
      	  dr_is_save = max(*rDM_is_final,*r_is_final)/RDIV;
    	  r_is_check = dr_is_save;
	 }

    r_is=0.0;                            /* initial isotropic radius */
    r=0.0;                               /* initial radius */
    m=0.0;                               /* initial mass */ 
    p=p_center;                          /* initial pressure */ 
    mDM=0.0;
    pDM=p_centerDM;
    double RBM,RDM,RBMis,RDMis=0.0;
    n_nearest = n_tab/2;
    n_nearestDM = n_tabDM/2;
    r_is_gp[1]=0.0;
    r_gp[1]=0.0;
    m_gp[1]=0.0;
    m_gpDM[1]=0.0;
    lambda_gp[1]=0.0;
    e_d_gp[1] = e_at_p(p_center,p_surface, log_e_tab, log_p_tab, n_tab, &n_nearest, eos_type, Gamma_P);; 
    e_d_gpDM[1] = e_at_p(pDM,p_surfaceDM, log_e_tabDM, log_p_tabDM, n_tabDM, &n_nearestDM, eos_typeDM, Gamma_PDM);
    
    p_d_gp[1] =p;
    
    p_d_gpDM[1] = pDM;
    
    



    while ( (p>=p_surface) || (pDM>=p_surfaceDM) ) { 
 
      e_d = e_at_p(p,p_surface, log_e_tab, log_p_tab, n_tab, &n_nearest, eos_type, 
                                                                    Gamma_P);
      e_dDM = e_at_p(pDM,p_surfaceDM, log_e_tabDM, log_p_tabDM, n_tabDM, &n_nearestDM, eos_typeDM, 
                                                                    Gamma_PDM);
     if((i_check==3) && (r_is>r_is_check) && (i<=RDIV)) {
      r_is_gp[i]=r_is;
      r_gp[i]=r;
      if((p>=p_surface)){
        m_gp[i]=m;
        e_d_gp[i]=e_d;
        p_d_gp[i]=p;
      }else{
        m_gp[i]=0.0;
        e_d_gp[i]=0.0;
        p_d_gp[i]=0.0;    
      } 
      if((pDM>=p_surfaceDM)){
        m_gpDM[i]=mDM;
        e_d_gpDM[i]=e_dDM; 
        p_d_gpDM[i]=pDM; 
      }else{
        m_gpDM[i]=0.0;
        e_d_gpDM[i]=0.0;       
        p_d_gpDM[i]=0.0; 
      }
      i++;   
      r_is_check += dr_is_save;
     }    
       

     if((p>=p_surface)){
      RBM=r;
      RBMis=r_is;
     (*r_is_final)=r_is;
     (*r_final)=r;
     }
     if((pDM>=p_surfaceDM)){
      RDM=r;
      RDMis=r_is;
      (*rDM_is_final)=r_is;
      (*rDM_final)=r;
     }
     (*m_final)=m;
     (*mDM_final)=mDM;

 
     a1=dr_dr_is(r_is,r,m+mDM);
     
     if(p>p_surface){  
     b1=dm_dr_is(r_is,r,m+mDM,p, e_center, p_surface, log_e_tab, log_p_tab, n_tab,
                                                &n_nearest, eos_type, Gamma_P);
     c1=dp_dr_is(r_is,r,m+mDM,p,pDM, e_center, p_surface, log_e_tab, log_p_tab, n_tab,
                                                &n_nearest, eos_type, Gamma_P);
     }

     if(pDM>p_surfaceDM){  
     b1DM=dm_dr_is(r_is,r,m+mDM,pDM, e_centerDM, p_surfaceDM, log_e_tabDM, log_p_tabDM, n_tabDM,
                                                &n_nearestDM, eos_typeDM, Gamma_PDM);
     c1DM=dp_dr_is(r_is,r,m+mDM,pDM,p, e_centerDM, p_surfaceDM, log_e_tabDM, log_p_tabDM, n_tabDM,
                                                &n_nearestDM, eos_typeDM, Gamma_PDM);
     }

     
     a2=dr_dr_is(r_is+h/2.0, r+h*a1/2.0, m+h*b1/2.0+mDM+h*b1DM/2.0);

     if(p+h*b1/2.0>p_surface){                           
     b2=dm_dr_is(r_is+h/2.0, r+h*a1/2.0, m+h*b1/2.0+mDM+h*b1DM/2.0, p+h*c1/2.0, e_center, 
                          p_surface, log_e_tab, log_p_tab, n_tab,&n_nearest, 
                          eos_type, Gamma_P);

     c2=dp_dr_is(r_is+h/2.0, r+h*a1/2.0, m+h*b1/2.0+mDM+h*b1DM/2.0, p+h*c1/2.0, pDM+h*c1DM/2.0, e_center, 
                          p_surface, log_e_tab, log_p_tab, n_tab,&n_nearest,  
                          eos_type, Gamma_P);
     }
     if(pDM+h*b1DM/2.0>p_surfaceDM){                           
     b2DM=dm_dr_is(r_is+h/2.0, r+h*a1/2.0, m+h*b1/2.0+mDM+h*b1DM/2.0, pDM+h*c1DM/2.0, e_centerDM, 
                          p_surfaceDM, log_e_tabDM, log_p_tabDM, n_tabDM,&n_nearestDM, 
                          eos_typeDM, Gamma_PDM);

     c2DM=dp_dr_is(r_is+h/2.0, r+h*a1/2.0, m+h*b1/2.0+mDM+h*b1DM/2.0, pDM+h*c1DM/2.0, p+h*c1/2.0, e_centerDM, p_surfaceDM, log_e_tabDM, log_p_tabDM, n_tabDM,&n_nearestDM,eos_typeDM, Gamma_PDM);
     }                  
                          

     a3=dr_dr_is(r_is+h/2.0, r+h*a2/2.0, m+h*b2/2.0+ mDM+h*b2DM/2.0);
     
     if(p+h*c2/2.0>p_surface){                           
     b3=dm_dr_is(r_is+h/2.0, r+h*a2/2.0, m+h*b2/2.0+ mDM+h*b2DM/2.0, p+h*c2/2.0, e_center, 
                          p_surface, log_e_tab, log_p_tab, n_tab,&n_nearest, 
                          eos_type, Gamma_P);

     c3=dp_dr_is(r_is+h/2.0, r+h*a2/2.0, m+h*b2/2.0+ mDM+h*b2DM/2.0, p+h*c2/2.0, pDM+h*c2DM/2.0, e_center, 
                          p_surface, log_e_tab, log_p_tab, n_tab,&n_nearest, 
                          eos_type, Gamma_P);
     }
     
     
     if(pDM+h*c2DM/2.0>p_surfaceDM){                           
     b3DM=dm_dr_is(r_is+h/2.0, r+h*a2/2.0, m+h*b2/2.0+ mDM+h*b2DM/2.0, pDM+h*c2DM/2.0, e_centerDM, 
                          p_surfaceDM, log_e_tabDM, log_p_tabDM, n_tabDM,&n_nearestDM, 
                          eos_typeDM, Gamma_PDM);

     c3DM=dp_dr_is(r_is+h/2.0, r+h*a2/2.0, m+h*b2/2.0+ mDM+h*b2DM/2.0, pDM+h*c2DM/2.0, p+h*c2/2.0, e_centerDM, 
                          p_surfaceDM, log_e_tabDM, log_p_tabDM, n_tabDM,&n_nearestDM, 
                          eos_typeDM, Gamma_PDM);
     }

     a4=dr_dr_is(r_is+h, r+h*a3, m+h*b3+mDM+h*b3DM);

     if(p+h*c3>p_surface){                           
     b4=dm_dr_is(r_is+h, r+h*a3, m+h*b3+mDM+h*b3DM, p+h*c3, e_center, p_surface, 
                                log_e_tab, log_p_tab, n_tab,&n_nearest, 
                                eos_type, Gamma_P);

     c4=dp_dr_is(r_is+h, r+h*a3, m+h*b3+mDM+h*b3DM, p+h*c3,pDM+h*c3DM, e_center, p_surface, 
                                log_e_tab, log_p_tab, n_tab,&n_nearest, 
                                eos_type, Gamma_P);
     }     
     if(pDM+h*c3DM>p_surfaceDM){                           
     b4DM=dm_dr_is(r_is+h, r+h*a3, m+h*b3+mDM+h*b3DM, pDM+h*c3DM, e_centerDM, p_surfaceDM, 
                                log_e_tabDM, log_p_tabDM, n_tabDM,&n_nearestDM, 
                                eos_typeDM, Gamma_PDM);

     c4DM=dp_dr_is(r_is+h, r+h*a3, m+h*b3+mDM+h*b3DM, pDM+h*c3DM,p+h*c3, e_centerDM, p_surfaceDM, 
                                log_e_tabDM, log_p_tabDM, n_tabDM,&n_nearestDM, 
                                eos_typeDM, Gamma_PDM);
     }

     
          
     r += (h/6.0)*(a1+2.*a2+2.*a3+a4);
     if((p+(h/6.0)*(c1+2.*c2+2.*c3+c4))>p_surface){
       p += (h/6.0)*(c1+2.*c2+2.*c3+c4);
       m += (h/6.0)*(b1+2.*b2+2.*b3+b4);
     }else{
       p=0.0;
     }

     if((pDM+(h/6.0)*(c1DM+2.*c2DM+2.*c3DM+c4DM))>p_surfaceDM){
       pDM += (h/6.0)*(c1DM+2.*c2DM+2.*c3DM+c4DM);
       mDM += (h/6.0)*(b1DM+2.*b2DM+2.*b3DM+b4DM);
     }else{
       pDM=0.0;
     }       

     
     r_is += h;

    }
    (*r_is_final)=RBMis;
    (*r_final)=RBM;
    (*rDM_is_final)=RDMis;
    (*rDM_final)=RDM;
    r_is_gp[RDIV]=max(RBMis,RDMis);
    r_gp[RDIV]=max(RBM,RDM);
    m_gp[RDIV]=(*m_final);
    m_gpDM[RDIV]=(*mDM_final);

/* Rescale r_is and compute lambda */

   if(i_check==3) {
      k_rescale=0.5*(r_gp[RDIV]/r_is_gp[RDIV])*(1.0-(*m_final+*mDM_final)/r_gp[RDIV]+
                sqrt(1.0-2.0*(*m_final+*mDM_final)/r_gp[RDIV]));
 
      (*r_is_final) *= k_rescale;
      (*rDM_is_final) *= k_rescale;
       
      nu_s = log((1.0-(*m_final+*mDM_final)/(2.0*r_is_gp[RDIV]*k_rescale))/(1.0+(*m_final+*mDM_final)/
               (2.0*r_is_gp[RDIV]*k_rescale)));

      for(i=1;i<=RDIV;i++) {
         r_is_gp[i] *= k_rescale;
 
         if(i==1){ 
           lambda_gp[1]= log(1.0/k_rescale);
  

         }
           else {lambda_gp[i]=log(r_gp[i]/r_is_gp[i]); }

         if(((e_d_gp[i]+e_d_gpDM[i])<(e_surfaceDM+e_surface))){
 
           hh=0.0;
          
          /*else if((e_d_gp[i]>e_surface)&&(e_d_gpDM[i]<e_surfaceDM)){
              p=p_at_e(e_d_gp[i],e_surface, log_p_tab, log_e_tab, n_tab, &n_nearest);
	      hh=h_at_p(p,p_surface, log_h_tab, log_p_tab, n_tab, &n_nearest);
          }
          else if((e_d_gp[i]<e_surface)&&(e_d_gpDM[i]>e_surfaceDM)){
              pDM=p_at_e(e_d_gpDM[i],e_surfaceDM, log_p_tabDM, log_e_tabDM, n_tabDM, &n_nearestDM);
	      hh=h_at_p(pDM,p_surfaceDM, log_h_tabDM, log_p_tabDM, n_tabDM, &n_nearestDM);
          }*/
          }else {  
          
                  p=p_at_e(e_d_gp[i],e_surface, log_p_tab, log_e_tab, n_tab, &n_nearest);
                  double hB=exp(h_at_p(p,p_surface, log_h_tab, log_p_tab, n_tab, &n_nearest));
                  double eB=e_at_p(p,p_surface, log_e_tab, log_p_tab, n_tab, &n_nearest, eos_type,Gamma_P);
                  double rhoB=(eB+p)/hB;
                     
                  pDM=p_at_e(e_d_gpDM[i],e_surfaceDM, log_p_tabDM, log_e_tabDM, n_tabDM, &n_nearestDM);
                  double hx=exp(h_at_p(pDM,p_surfaceDM, log_h_tabDM, log_p_tabDM, n_tabDM, &n_nearestDM));
                  double ex=e_at_p(pDM,p_surfaceDM, log_e_tabDM, log_p_tabDM, n_tabDM, &n_nearestDM, eos_typeDM,Gamma_PDM);
                  double rhox=(ex+pDM)/hx;               
                  
                 // hh=((hB*rhoB+hx*rhox)/(rhoB+rhox));
 

int Num=i;
 int size=RDIV-Num;

  double sum=0.0,          /* final sum */
         h0,h1,hph,hdh,hmh;


  /*if((RDIV-Num-1)%2==0){
*/
     if(p_d_gp[i]==0.0){
       sum=h_at_p(pDM,p_surfaceDM, log_h_tabDM, log_p_tabDM, n_tabDM, &n_nearestDM);
     }else if(p_d_gpDM[i]==0.0){
       sum=h_at_p(p,p_surface, log_h_tab, log_p_tab, n_tab, &n_nearest);
     }else{
      for(int j=RDIV;j>=Num+2;j-=2) {
    
      h0=-(p_d_gpDM[j]+p_d_gp[j]-p_d_gpDM[j-1]-p_d_gp[j-1]);
      h1=-(p_d_gpDM[j-1]+p_d_gp[j-1]-p_d_gpDM[j-2]-p_d_gp[j-2]);
      hph=h1+h0;
      hdh=h1/h0;
      hmh=h1*h0;
      //if((dx!=0)&&(dx1!=0)){
       sum+= (hph)/6.*((2.-hdh)/(p_d_gpDM[j]+p_d_gp[j]+e_d_gpDM[j]+e_d_gp[j])+hph*hph/(p_d_gpDM[j-1]+p_d_gp[j-1]+e_d_gpDM[j-1]+e_d_gp[j-1])/(hmh)+(2.-1./hdh)/(p_d_gpDM[j-2]+p_d_gp[j-2]+e_d_gpDM[j-2]+e_d_gp[j-2]));


       //}
           //  printf("%lf %lf %lf %lf\n",(2-dx1/dx)*f(p_d_gpDM[j]+p_d_gp[j],e_d_gpDM[j]+e_d_gp[j]),(dx+dx1)*(dx+dx1)*f(p_d_gpDM[j+1]+p_d_gp[j+1],e_d_gpDM[j+1]+e_d_gp[j+1])/(dx*dx1),(2-dx/dx1)*f(p_d_gpDM[j+2]+p_d_gp[j+2],e_d_gpDM[j+2]+e_d_gp[j+2]),p_d_gpDM[j]);
       
      }
      
if((RDIV-Num)%2==1){
      h1=(p_d_gpDM[Num]+p_d_gp[Num]-p_d_gpDM[Num+1]-p_d_gp[Num+1]);
      h0=(p_d_gpDM[Num+1]+p_d_gp[Num+1]-p_d_gpDM[Num+2]-p_d_gp[Num+2]);

     sum += (2.*h1*h1+3.*h0*h1)/(6.*(h0+h1))/(p_d_gpDM[Num]+p_d_gp[Num]+e_d_gpDM[Num]+e_d_gp[Num]);
     sum += (h1*h1+3.*h1*h0)/(p_d_gpDM[Num+1]+p_d_gp[Num+1]+e_d_gpDM[Num+1]+e_d_gp[Num+1])/(6.*h0);
     sum -= h1*h1*h1/(p_d_gpDM[Num+2]+p_d_gp[Num+2]+e_d_gpDM[Num+2]+e_d_gp[Num+2])/(6.*h0*(h0+h1));
}      
      }
      hh=sum;

    /*}else{
      for(int j=Num;j<=RDIV-2;j+=1) {

      dx=(p_d_gpDM[j]+p_d_gp[j]-p_d_gpDM[j+1]-p_d_gp[j+1]);
      dx1=(p_d_gpDM[j+1]+p_d_gp[j+1]-p_d_gpDM[j+2]-p_d_gp[j+2]);
       sum+= (dx+dx1)/6.*((2-dx1/dx)*f(p_d_gpDM[j]+p_d_gp[j],e_d_gpDM[j]+e_d_gp[j])+(dx+dx1)*(dx+dx1)*f(p_d_gpDM[j+1]+p_d_gp[j+1],e_d_gpDM[j+1]+e_d_gp[j+1])/(dx*dx1)+(2-dx/dx1)*f(p_d_gpDM[j+2]+p_d_gp[j+2],e_d_gpDM[j+2]+e_d_gp[j+2]));
       
      }   
      dx=(p_d_gpDM[RDIV-1]+p_d_gp[RDIV-1]-p_d_gpDM[RDIV]-p_d_gp[RDIV]);
      dx1=(p_d_gpDM[RDIV-2]+p_d_gp[RDIV-2]-p_d_gpDM[RDIV-1]-p_d_gp[RDIV-1]);

     sum+=  (2*dx*dx+3*dx*dx1)/(6*(dx1+dx2))*f(p_d_gpDM[RDIV]+p_d_gp[RDIV],e_d_gpDM[RDIV]+e_d_gp[RDIV])+(dx*dx+3*dx*dx1)/(6*(dx2))*f(p_d_gpDM[RDIV-1]+p_d_gp[RDIV-1],e_d_gpDM[RDIV-1]+e_d_gp[RDIV-1])-(dx*dx*dx)/(6*dx2*(dx1+dx))*f(p_d_gpDM[RDIV-2]+p_d_gp[RDIV-2],e_d_gpDM[RDIV-2]+e_d_gp[RDIV-2]);
    }
*/
  
  
  /*
  for(int j=Num;j<=RDIV-2;j+=2){
     dx=p_d_gpDM[j]+p_d_gp[j]-p_d_gpDM[j+1]-p_d_gp[j+1];
     dx1=p_d_gpDM[j+1]+p_d_gp[j+1]-p_d_gpDM[j+2]-p_d_gp[j+2];
     //prjntf("%lf %lf \n",sum,( f( p_d_gp[j]+p_d_gpDM[j],e_d_gp[j]+e_d_gpDM[j])  )*(dx/3.0));
     sum =sum+ (f( p_d_gp[j]+p_d_gpDM[j],e_d_gp[j]+e_d_gpDM[j])*dx  +4*f(p_d_gp[j+1]+p_d_gpDM[j+1],e_d_gp[j+1]+e_d_gpDM[j+1])*(dx+dx1)/2+f(p_d_gp[j+2]+p_d_gpDM[j+2],e_d_gp[j+2]+e_d_gpDM[j+2])*dx1  )/3.0;
   }*/





              }

         nu_gp[i]=nu_s-hh;
                //printf("%f %f\n",hh, h_at_p(p,p_surface, log_h_tab, log_p_tab, n_tab, &n_nearest));
               //printf("%lf %lf\n",h_integ_tab(p_d_gp,e_d_gp,p_d_gpDM,e_d_gpDM,log_p_tab[1],i),p_d_gpDM[i]); 
      }
      nu_gp[RDIV]=nu_s;


                   

   }
  printf("RK4 for non-rotating NSs\n" );
  printf("RBM RDM Mtot MBM MDM MDM_fract\n" );
    printf( "%7g %7g %7g %7g %7g %7g\n", 
     RBM*sqrt(KAPPA)/100000.,RDM*sqrt(KAPPA)/100000.,(*m_final+*mDM_final)*sqrt(KAPPA)*C*C/(G*MSUN),(*m_final)*sqrt(KAPPA)*C*C/(G*MSUN),(*mDM_final)*sqrt(KAPPA)*C*C/(G*MSUN),(*mDM_final)/(*m_final+*mDM_final));
}
//
/*C*/
/*************************************************************************/
void sphere(double s_gp[SDIV+1], 
	    double log_e_tab[2001], 
	    double log_p_tab[2001], 
	    double log_h_tab[2001],
	    double log_n0_tab[2001], 
	    int n_tab,                 
	    char eos_type[],
	    double Gamma_P, 
	    double e_center,
	    double p_center, 
	    double h_center,
	    double p_surface,
	    double e_surface,
	    double log_e_tabDM[2001], 
	    double log_p_tabDM[2001], 
	    double log_h_tabDM[2001],
	    double log_n0_tabDM[2001], 
	    int n_tabDM,                 
	    char eos_typeDM[],
	    double Gamma_PDM, 
	    double e_centerDM,
	    double p_centerDM, 
	    double h_centerDM,
	    double p_surfaceDM,
	    double e_surfaceDM,
	    double **rho,
	    double **gama,
	    double **alpha,
	    double **omega,
	    double *r_e,
	    double *r_eDM)

{
 int s,
     m,
     n_nearest;

 double r_is_s,
        r_is_final,
        r_final, 
        m_final,
        rDM_is_final,
        rDM_final, 
        mDM_final,
        lambda_s,
        nu_s,
        r_is_gp[RDIV+1],
        lambda_gp[RDIV+1],
        nu_gp[RDIV+1],
        gama_mu_0[SDIV+1],
        rho_mu_0[SDIV+1],
        gama_eq,
        rho_eq,
        R_is,
        s_e=0.5;

 /* The function TOV integrates the TOV equations. The function
	can be found in the file equil.c */

 TOV(1, eos_type, e_center, p_center, p_surface, e_surface, Gamma_P,
              log_e_tab, log_p_tab,log_n0_tab, log_h_tab, n_tab, r_is_gp, lambda_gp, 
               eos_typeDM, e_centerDM, p_centerDM, p_surfaceDM, e_surfaceDM, Gamma_PDM,
              log_e_tabDM, log_p_tabDM,log_n0_tabDM, log_h_tabDM, n_tabDM,
              nu_gp, &r_is_final, &r_final, &m_final,&rDM_is_final, &rDM_final, &mDM_final);

 TOV(2, eos_type, e_center, p_center, p_surface, e_surface, Gamma_P,
              log_e_tab, log_p_tab,log_n0_tab, log_h_tab, n_tab, r_is_gp, lambda_gp, 
               eos_typeDM, e_centerDM, p_centerDM, p_surfaceDM, e_surfaceDM, Gamma_PDM,
              log_e_tabDM, log_p_tabDM,log_n0_tabDM, log_h_tabDM, n_tabDM,
              nu_gp, &r_is_final, &r_final, &m_final,&rDM_is_final, &rDM_final , &mDM_final);
              
 TOV(3, eos_type, e_center, p_center, p_surface, e_surface, Gamma_P,
              log_e_tab, log_p_tab,log_n0_tab, log_h_tab, n_tab, r_is_gp, lambda_gp, 
               eos_typeDM, e_centerDM, p_centerDM, p_surfaceDM, e_surfaceDM, Gamma_PDM,
              log_e_tabDM, log_p_tabDM,log_n0_tabDM, log_h_tabDM, n_tabDM,
              nu_gp, &r_is_final, &r_final, &m_final,&rDM_is_final, &rDM_final , &mDM_final);


 R_is=r_is_gp[RDIV];
 n_nearest=RDIV/2;
 for(s=1;s<=SDIV;s++) {
    r_is_s=R_is*(s_gp[s]/(1.0-s_gp[s]));

    if(r_is_s<R_is) {
      lambda_s=interp(r_is_gp,lambda_gp,RDIV,r_is_s,&n_nearest);
      nu_s=interp(r_is_gp,nu_gp,RDIV,r_is_s,&n_nearest);

    }
    else {

      lambda_s=2.0*log(1.0+(m_final+mDM_final)/(2.0*r_is_s));
      nu_s=log((1.0-(m_final+mDM_final)/(2.0*r_is_s))/(1.0+(m_final+mDM_final)/(2.0*r_is_s)));
    }

    gama[s][1]=nu_s+lambda_s;
    rho[s][1]=nu_s-lambda_s;

    for(m=1;m<=MDIV;m++) {
        gama[s][m]=gama[s][1];        
        rho[s][m]=rho[s][1];
        alpha[s][m]=(gama[s][1]-rho[s][1])/2.0;
        omega[s][m]=0.0; 
    }
 
    gama_mu_0[s]=gama[s][1];                   /* gama at \mu=0 */
    rho_mu_0[s]=rho[s][1];                     /* rho at \mu=0 */

 }

   n_nearest=SDIV/2;
   gama_eq = interp(s_gp,gama_mu_0,SDIV,s_e,&n_nearest); /* gama at equator */
   rho_eq = interp(s_gp,rho_mu_0,SDIV,s_e,&n_nearest);   /* rho at equator */
    
    

    
   if((r_final)>(rDM_final)){

     (*r_e)= r_final*exp(0.5*(rho_eq-gama_eq)); 
     (*r_eDM)= rDM_final*exp(interp(s_gp,rho_mu_0,SDIV,rDM_final/(rDM_final+r_final),&n_nearest)-interp(s_gp,gama_mu_0,SDIV,rDM_final/(rDM_final+r_final),&n_nearest)); 
   }else{
      (*r_e)= r_final*exp(interp(s_gp,rho_mu_0,SDIV,r_final/(rDM_final+r_final),&n_nearest)-interp(s_gp,gama_mu_0,SDIV,r_final/(rDM_final+r_final),&n_nearest)); 
     (*r_eDM)= rDM_final*exp(0.5*(rho_eq-gama_eq)); 

   }
   //printf("%f %f %f \n",*r_e,*r_eDM, (gama_eq+rho_eq-gama[1][1]-rho[1][1])/2.);
   //printf("%f %f %f \n",s_gp[s],r_is_gp[RDIV-1]-r_is_gp[RDIV-2],r_is_gp[RDIV-2]-r_is_gp[RDIV-3]);
   //printf("%f %f \n",interp(s_gp,gama_mu_0,SDIV,0.5,&n_nearest),log((1.0-(m_final+mDM_final)/(2.0*max(r_is_final,rDM_is_final)))/(1.0+(m_final+mDM_final)/(2*max(r_is_final,rDM_is_final))))+2.0*log(1.0+(m_final+mDM_final*max(r_is_final,rDM_is_final))/(2.0)));
 /*if(mDM_final/(m_final+mDM_final)>0.1){
   (*r_eDM)=-1;
 }*/

}


/*C*/
/*************************************************************************/
/* Main iteration cycle for computation of the rotating star's metric    */
/*************************************************************************/
void spin(double s_gp[SDIV+1],
	  double mu[MDIV+1],
	  double log_e_tab[2001], 
	  double log_p_tab[2001], 
	  double log_h_tab[2001],
	  double log_n0_tab[2001], 
	  int n_tab,                 
	  char eos_type[],
	  double Gamma_P, 
	  double h_center,
	  double enthalpy_min,
	  double log_e_tabDM[2001], 
	  double log_p_tabDM[2001], 
	  double log_h_tabDM[2001],
	  double log_n0_tabDM[2001], 
	  int n_tabDM,                 
	  char eos_typeDM[],
	  double Gamma_PDM, 
	  double h_centerDM,
	  double enthalpy_minDM,
	  double **rho,
	  double **gama,
	  double **alpha,
	  double **omega,
	  double **energy,
	  double **pressure,
	  double **enthalpy,
	  double **velocity_sq,
	  double **energyDM,
	  double **pressureDM,
	  double **enthalpyDM,
	  double **velocity_sqDM,
	  int    a_check, 
	  double accuracy,
	  double cf,
	  double r_ratio,
	  double r_ratioDM,
	  double *r_e_new,
	  double *rDM_e_new,
	  double *Omega,
  	  double *OmegaDM)

 {
 int m,                      /* counter */
     s,                      /* counter */
     n,                      /* counter */
     k,                      /* counter */
     n_of_it=0,              /* number of iterations */
     n_nearest,
     n_nearestDM,
     print_dif = 0,
     i,
     j;

double **D2_rho,
  **D2_gama,
  **D2_omega;

float  ***f_rho,
       ***f_gama;

 
double   sum_rho=0.0,         /* intermediate sum in eqn for rho */
	 sum_gama=0.0,        /* intermediate sum in eqn for gama */
	 sum_omega=0.0,       /* intermediate sum in eqn for omega */
         r_e_old,             /* equatorial radius in previus cycle */
   	 dif=1.0,dif2=1.0,             /* difference | r_e_old - r_e | */
         d_gama_s,            /* derivative of gama w.r.t. s */
         d_gama_m,            /* derivative of gama w.r.t. m */
         d_rho_s,             /* derivative of rho w.r.t. s */
         d_rho_m,             /* derivative of rho w.r.t. m */
         d_omega_s,           /* derivative of omega w.r.t. s */
         d_omega_m,           /* derivative of omega w.r.t. m */
         d_gama_ss,           /* 2nd derivative of gama w.r.t. s */
         d_gama_mm,           /* 2nd derivative of gama w.r.t. m */
         d_gama_sm,           /* derivative of gama w.r.t. m and s */
         temp1,                /* temporary term in da_dm */ 
         temp2, 
         temp3,
         temp4,
         temp5,
         temp6,
         temp7,
         temp8,
         m1,                  
         s1,
         s2,
         ea,eaDM,
         rsm,
         gsm,
         omsm,
         esm,
         psm,
         esmDM,
         psmDM,
         v2sm,
         v2smDM,
         mum,
         sgp,
         s_1,
         e_gsm,e_gsmDM,
         e_rsm,e_rsmDM, 
         rho0sm,
         term_in_Omega_h,
         term_in_Omega_hDM,
         r_p,
         s_p,
         gama_pole_h,                  /* gama^hat at pole */  
         gama_pole_hDM,                  /* gama^hat at pole */  
         gama_center_h,                /* gama^hat at center */
         gama_equator_h,               /* gama^hat at equator */
         gama_equator_hDM,               /* gama^hat at equator */
         rho_pole_h,                   /* rho^hat at pole */ 
         rho_pole_hDM,                   /* rho^hat at pole */ 
         rho_center_h,                 /* rho^hat at center */
         rho_equator_h,                /* rho^hat at equator */ 
         rho_equator_hDM,                /* rho^hat at equator */ 
         omega_equator_h,              /* omega^hat at equator */         
         omega_equator_hDM,              /* omega^hat at equator */         
         gama_mu_1[SDIV+1],            /* gama at \mu=1 */
         gama_mu_0[SDIV+1],            /* gama at \mu=0 */
         rho_mu_1[SDIV+1],             /* rho at \mu=1 */
         rho_mu_0[SDIV+1],             /* rho at \mu=0 */
         omega_mu_0[SDIV+1],           /* omega at \mu=0 */
       **da_dm,
       **dgds,
       **dgdm,
       **D1_rho,
       **D1_gama,
       **D1_omega,
       **S_gama,
       **S_rho,
       **S_omega,
       **f2n,
       **P_2n,   
       **P1_2n_1,
         Omega_h,
         Omega_hDM,
         sin_theta[MDIV+1],
         theta[MDIV+1],
         sk,
         sj,
         sk1,
         sj1,
         r_e,
         dum_e,
         k0,k1,k2,k3,
         kDM0,kDM1,kDM2,kDM3,
         l0,l1,l2,l3,
         lDM0,lDM1,lDM2,lDM3,r_eq;

         double s_e=0.5,
         s_eDM_new=0.5,
         s_pDM_new=0.5,
         s_e_new=0.0,
         r_eDM, r_eDM_old,
         ratio_old, ratio_oldDM;
        

       
    f2n = dmatrix(1,LMAX+1,1,SDIV);
    f_rho = f3tensor(1,SDIV,1,LMAX+1,1,SDIV);
    f_gama = f3tensor(1,SDIV,1,LMAX+1,1,SDIV);
 
    P_2n = dmatrix(1,MDIV,1,LMAX+1);   
    P1_2n_1 = dmatrix(1,MDIV,1,LMAX+1);


    for(n=0;n<=LMAX;n++) 
       for(i=2;i<=SDIV;i++) f2n[n+1][i] = pow((1.0-s_gp[i])/s_gp[i],2.0*n);

    if(SMAX!=1.0) {

     for(j=2;j<=SDIV;j++)
        for(n=1;n<=LMAX;n++)
           for(k=2;k<=SDIV;k++) {
                 sk=s_gp[k];
                 sj=s_gp[j];
                 sk1=1.0-sk;
                 sj1=1.0-sj;

                 if(k<j) {   
                          f_rho[j][n+1][k] = f2n[n+1][j]*sj1/(sj*
                                  f2n[n+1][k]*sk1*sk1);
                          f_gama[j][n+1][k] = f2n[n+1][j]/(f2n[n+1][k]*sk*sk1);
	                 }else {     
                          f_rho[j][n+1][k] = f2n[n+1][k]/(f2n[n+1][j]*sk*sk1);
                          f_gama[j][n+1][k] = f2n[n+1][k]*sj1*sj1*sk/(sj*sj
                                            *f2n[n+1][j]*sk1*sk1*sk1);
                 }
	    }
     j=1;
 
       n=0; 
       for(k=2;k<=SDIV;k++) {
          sk=s_gp[k];
          f_rho[j][n+1][k]=1.0/(sk*(1.0-sk));
       }

       n=1;
       for(k=2;k<=SDIV;k++) {
          sk=s_gp[k];
          sk1=1.0-sk;         
          f_rho[j][n+1][k]=0.0;
          f_gama[j][n+1][k]=1.0/(sk*sk1);
       }

       for(n=2;n<=LMAX;n++)
          for(k=1;k<=SDIV;k++) {
             f_rho[j][n+1][k]=0.0;
             f_gama[j][n+1][k]=0.0;
          }


     k=1;

       n=0;
       for(j=1;j<=SDIV;j++)
          f_rho[j][n+1][k]=0.0;

       for(j=1;j<=SDIV;j++)
          for(n=1;n<=LMAX;n++) {
             f_rho[j][n+1][k]=0.0;
             f_gama[j][n+1][k]=0.0;
          }


     n=0;
     for(j=2;j<=SDIV;j++)
        for(k=2;k<=SDIV;k++) {
               sk=s_gp[k];
               sj=s_gp[j];
               sk1=1.0-sk;
               sj1=1.0-sj;

               if(k<j) 
                 f_rho[j][n+1][k] = sj1/(sj*sk1*sk1);
               else     
                 f_rho[j][n+1][k] = 1.0/(sk*sk1);
      }

   }
   else{      
        for(j=2;j<=SDIV-1;j++)
           for(n=1;n<=LMAX;n++)
              for(k=2;k<=SDIV-1;k++) {
                 sk=s_gp[k];
                 sj=s_gp[j];
                 sk1=1.0-sk;
                 sj1=1.0-sj;

                 if(k<j) {   
                          f_rho[j][n+1][k] = f2n[n+1][j]*sj1/(sj*
                                           f2n[n+1][k]*sk1*sk1);
                          f_gama[j][n+1][k] = f2n[n+1][j]/(f2n[n+1][k]*sk*sk1);
                 }else {     
                          f_rho[j][n+1][k] = f2n[n+1][k]/(f2n[n+1][j]*sk*sk1);

                          f_gama[j][n+1][k] = f2n[n+1][k]*sj1*sj1*sk/(sj*sj
                                            *f2n[n+1][j]*sk1*sk1*sk1);
                 }
	      }
   
        j=1;
 
          n=0; 
          for(k=2;k<=SDIV-1;k++) {
             sk=s_gp[k];
             f_rho[j][n+1][k]=1.0/(sk*(1.0-sk));
          }

          n=1;
          for(k=2;k<=SDIV-1;k++) {
             sk=s_gp[k];
             sk1=1.0-sk;         
             f_rho[j][n+1][k]=0.0;
             f_gama[j][n+1][k]=1.0/(sk*sk1);
          }

          for(n=2;n<=LMAX;n++)
             for(k=1;k<=SDIV-1;k++) {
                f_rho[j][n+1][k]=0.0;
                f_gama[j][n+1][k]=0.0;
             }

        k=1;
 
          n=0;
          for(j=1;j<=SDIV-1;j++)
             f_rho[j][n+1][k]=0.0;

          for(j=1;j<=SDIV-1;j++)
             for(n=1;n<=LMAX;n++) {
                f_rho[j][n+1][k]=0.0;
                f_gama[j][n+1][k]=0.0;
             }
 
 
        n=0;
          for(j=2;j<=SDIV-1;j++)
             for(k=2;k<=SDIV-1;k++) {
                sk=s_gp[k];
                sj=s_gp[j];
                sk1=1.0-sk;
                sj1=1.0-sj;

                if(k<j) 
                  f_rho[j][n+1][k] = sj1/(sj*sk1*sk1);
                else     
                  f_rho[j][n+1][k] = 1.0/(sk*sk1);
             }
 
        j=SDIV;
          for(n=1;n<=LMAX;n++)
             for(k=1;k<=SDIV;k++) {
                f_rho[j][n+1][k] = 0.0;
                f_gama[j][n+1][k] = 0.0;
             }

        k=SDIV;
          for(j=1;j<=SDIV;j++)
              for(n=1;n<=LMAX;n++) {
                 f_rho[j][n+1][k] = 0.0;
                 f_gama[j][n+1][k] = 0.0;
              }
   }

  n=0;
   for(i=1;i<=MDIV;i++)
      P_2n[i][n+1]=legendre(2*n,mu[i]);

   for(i=1;i<=MDIV;i++)
     for(n=1;n<=LMAX;n++) {
      P_2n[i][n+1]=legendre(2*n,mu[i]);
      P1_2n_1[i][n+1] = plgndr(2*n-1 ,1,mu[i]);
    }

  free_dmatrix(f2n,1,LMAX+1,1,SDIV);


  for(m=1;m<=MDIV;m++) { 
     sin_theta[m] = sqrt(1.0-mu[m]*mu[m]);  
     theta[m] = asin(sin_theta[m]);
  }

  
        bool cond=((*r_e_new)>=(*rDM_e_new));
        
    n_nearest = n_tab/2;
    n_nearestDM = n_tabDM/2;
      ratio_old=r_ratio;
      ratio_oldDM=r_ratioDM;

      r_ratio = ((cond)?ratio_old:ratio_oldDM);
      r_ratioDM=((!cond)?ratio_old:ratio_oldDM);

      r_eq = ((cond)?(*r_e_new):(*rDM_e_new));
      r_eDM=(!(cond)?(*r_e_new):(*rDM_e_new));


  while(((dif> accuracy)||(dif2> accuracy)) || n_of_it<2) { 

      if(print_dif!=0)
         printf("%4.3e\n",dif);

 
      /* Rescale potentials and construct arrays with the potentials along
       | the equatorial and polar directions.
      */        

      for(s=1;s<=SDIV;s++) {
                 //printf("%f %f %f %f %d\n",enthalpy[s][1],h_center,enthalpyDM[s][1],h_centerDM,s);
         for(m=1;m<=MDIV;m++) {
            rho[s][m] /= SQ(r_eq);
            gama[s][m] /= SQ(r_eq); 
            alpha[s][m] /= SQ(r_eq);
            omega[s][m] *= r_eq;
         }
         rho_mu_0[s]=rho[s][1];     
         gama_mu_0[s]=gama[s][1];   
         omega_mu_0[s]=omega[s][1]; 
         rho_mu_1[s]=rho[s][MDIV];  
         gama_mu_1[s]=gama[s][MDIV];
      }
 
      /* Compute new r_e. */ 

      r_e_old=r_eq;
      r_eDM_old=r_eDM;
      r_p=r_ratio*r_eq;                          
      s_p=r_p/(r_p+r_eq); 
      s_pDM_new=r_ratioDM*r_eDM/(r_ratioDM*r_eDM+r_eq);            
      s_eDM_new=r_eDM/(r_eDM+r_eq);            
  
      n_nearest= SDIV/2;
      gama_pole_h=interp(s_gp,gama_mu_1,SDIV,s_p,&n_nearest); 
      gama_pole_hDM=interp(s_gp,gama_mu_1,SDIV,s_pDM_new,&n_nearest); 
      gama_equator_h=interp(s_gp,gama_mu_0,SDIV,s_e,&n_nearest);
      gama_equator_hDM=interp(s_gp,gama_mu_0,SDIV,s_eDM_new,&n_nearest);
      gama_center_h=gama[1][1];                    
  
      rho_pole_h=interp(s_gp,rho_mu_1,SDIV,s_p,&n_nearest);   
      rho_pole_hDM=interp(s_gp,rho_mu_1,SDIV,s_pDM_new,&n_nearest);   
      rho_equator_h=interp(s_gp,rho_mu_0,SDIV,s_e,&n_nearest);
      rho_equator_hDM=interp(s_gp,rho_mu_0,SDIV,s_eDM_new,&n_nearest);
      rho_center_h=rho[1][1];                      
    n_nearest = n_tab/2;
    n_nearestDM = n_tabDM/2;


      //r_eq=sqrt(2*( log((exp(h_center)*rho_center +exp(h_centerDM)*rho_centerDM)/(rho_center+rho_centerDM)) )/(gama_pole_h+rho_pole_h-gama_center_h-rho_center_h));
       if((cond)){
       r_eq=sqrt(2.0*( h_center)/(gama_pole_h+rho_pole_h-gama_center_h-rho_center_h));

       r_eDM=sqrt(SQ(r_eDM_old)*2.0*( h_centerDM)/((gama_pole_hDM+rho_pole_hDM-gama_center_h-rho_center_h)*SQ(r_e_old)));

       if(h_centerDM<=enthalpy_minDM){
        r_eDM=0.;
        r_eDM_old=1.;
       }
       if(h_center<=enthalpy_min){
        r_eq=0.;
        r_e_old=1.;
       }
       }else{
       r_eq=sqrt(2.0*( h_centerDM)/(gama_pole_h+rho_pole_h-gama_center_h-rho_center_h));

       r_eDM=sqrt(SQ(r_eDM_old)*2.0*( h_center)/((gama_pole_hDM+rho_pole_hDM-gama_center_h-rho_center_h)*SQ(r_e_old)));

       if(h_centerDM<=enthalpy_minDM){
        r_eq=0.;
        r_e_old=1.;
       }
       if(h_center<=enthalpy_min){
        r_eDM=0.;
        r_eDM_old=1.;
       }
       }       

              //printf("%f %f %f %f \n",r_e_old,r_eDM_old,(sqrt(2*( h_center)/(gama_pole_h+rho_pole_h-gama_center_h-rho_center_h))/r_e_old)/((sqrt(2*( h_centerDM)/((gama_pole_hDM+rho_pole_hDM-gama_center_h-rho_center_h)*SQ(r_e_old))))),0.04);


      /* Compute angular velocity Omega. */
      
      if(cond){
      
      if(r_ratio==1.0) {
        Omega_h=0.0;
        omega_equator_h=0.0;
      } 
      else {
            omega_equator_h=interp(s_gp,omega_mu_0,SDIV,s_e, &n_nearest);
            term_in_Omega_h=1.0-exp(SQ(r_eq)*(gama_pole_h+rho_pole_h
                                             -gama_equator_h-rho_equator_h));
            if(term_in_Omega_h>=0.0) 
               Omega_h = omega_equator_h + exp(SQ(r_eq)*rho_equator_h)
                                            *sqrt(term_in_Omega_h);
            else {
                Omega_h=0.0;
	    }
      }

       if(r_ratioDM==1.0) {
        Omega_hDM=0.0;
        omega_equator_hDM=0.0;        
      } 
      else {
            omega_equator_hDM=interp(s_gp,omega_mu_0,SDIV,s_eDM_new, &n_nearestDM);
            term_in_Omega_hDM=1.0-exp(SQ(r_eq)*(gama_pole_hDM+rho_pole_hDM-gama_equator_hDM-rho_equator_hDM));
            if(term_in_Omega_hDM>=0.0) {
               Omega_hDM = omega_equator_hDM + exp(SQ(r_eq)*rho_equator_hDM) *sqrt(term_in_Omega_hDM);
                                            
            }
            else {
                Omega_hDM=0.0;
	    }
      }
      }else{
      
      if(r_ratio==1.0) {
        Omega_hDM=0.0;
        omega_equator_hDM=0.0;
      } 
      else {
            omega_equator_hDM=interp(s_gp,omega_mu_0,SDIV,s_e, &n_nearest);
            term_in_Omega_hDM=1.0-exp(SQ(r_eq)*(gama_pole_h+rho_pole_h
                                             -gama_equator_h-rho_equator_h));
            if(term_in_Omega_hDM>=0.0) 
               Omega_hDM = omega_equator_hDM + exp(SQ(r_eq)*rho_equator_h)
                                            *sqrt(term_in_Omega_hDM);
            else {
                Omega_hDM=0.0;
	    }
      }

       if(r_ratioDM==1.0) {
        Omega_h=0.0;
        omega_equator_h=0.0;        
      } 
      else {
            omega_equator_h=interp(s_gp,omega_mu_0,SDIV,s_eDM_new, &n_nearest);
            term_in_Omega_h=1.0-exp(SQ(r_eq)*(gama_pole_hDM+rho_pole_hDM-gama_equator_hDM-rho_equator_hDM));
            if(term_in_Omega_h>=0.0) {
               Omega_h = omega_equator_h + exp(SQ(r_eq)*rho_equator_hDM) *sqrt(term_in_Omega_h);
                                            
            }
            else {
                Omega_h=0.0;
	    }
      }      
      
      
      }

      /* Compute velocity, energy density and pressure. */
      
      n_nearest=n_tab/2; 
      n_nearestDM=n_tabDM/2; 
      
      for(s=1;s<=SDIV;s++) {
         sgp=s_gp[s];

         for(m=1;m<=MDIV;m++) {
            rsm=rho[s][m];
            
            if((r_ratio==1.0)&&(r_ratioDM==1.0)) 
                velocity_sq[s][m]=0.0;
            else 
            velocity_sq[s][m]=SQ((Omega_h-omega[s][m])*(sgp/(1.0-sgp))
                                  *sin_theta[m]*exp(-rsm*SQ(r_eq)));

            if(velocity_sq[s][m]>=1.0) 
              velocity_sq[s][m]=0.0;

            if((r_ratio==1.0)&&(r_ratioDM==1.0)) 
                velocity_sqDM[s][m]=0.0;
            else               
            velocity_sqDM[s][m]=SQ((Omega_hDM-omega[s][m])*(sgp/(1.0-sgp))
                                  *sin_theta[m]*exp(-rsm*SQ(r_eq)));

            if(velocity_sqDM[s][m]>=1.0) 
              velocity_sqDM[s][m]=0.0;
              
              
        if(cond){
        

             enthalpy[s][m]=enthalpy_min + 0.5*(SQ(r_eq)*(gama_pole_h+rho_pole_h
                           -gama[s][m]-rsm)-log(1.0-velocity_sq[s][m]));
  
            if((enthalpy[s][m]<=enthalpy_min) || (sgp>s_e)) {
                  pressure[s][m]=0.0;
                  energy[s][m]=0.0; 
	    }
            else {

                     pressure[s][m]=p_at_h(enthalpy[s][m],pow(10.0,log_h_tab[1]), log_p_tab, 
                                           log_h_tab, n_tab, &n_nearest);
                     energy[s][m]=e_at_p(pressure[s][m],pow(10.0,log_p_tab[1]), log_e_tab, 
                                      log_p_tab, n_tab, &n_nearest, eos_type,
                                       Gamma_P);
       //printf("%f %f %f %f %f %f %f %f \n",r_eq*s_gp[s]/(1-s_gp[s]),dif,dif2,r_eDM,rho_pole_hDM,gama_center_h,rho_center_h);
	           }
	           
	           
	    enthalpyDM[s][m]=enthalpy_minDM + 0.5*(SQ(r_eq)*(gama_pole_hDM+rho_pole_hDM
                           -gama[s][m]-rsm)-log(1.0-velocity_sqDM[s][m]));
  
            if((enthalpyDM[s][m]<=enthalpy_minDM) || (sgp>s_eDM_new)) {

                  pressureDM[s][m]=0.0;
                  energyDM[s][m]=0.0; 
	    }
	    
            else { 
            

                     pressureDM[s][m]=p_at_h(enthalpyDM[s][m],pow(10.0,log_h_tabDM[1]), log_p_tabDM, 
                                           log_h_tabDM, n_tabDM, &n_nearestDM);
                     energyDM[s][m]=e_at_p(pressureDM[s][m],pow(10.0,log_p_tabDM[1]), log_e_tabDM, 
                                      log_p_tabDM, n_tabDM, &n_nearestDM, eos_typeDM,
                                       Gamma_PDM);
	           }
           }else{
           

             enthalpyDM[s][m]=enthalpy_minDM + 0.5*(SQ(r_eq)*(gama_pole_h+rho_pole_h
                           -gama[s][m]-rsm)-log(1.0-velocity_sqDM[s][m]));
  
            if((enthalpyDM[s][m]<=enthalpy_minDM) || (sgp>s_e)) {
                  pressureDM[s][m]=0.0;
                  energyDM[s][m]=0.0; 
	    }
            else {

                     pressureDM[s][m]=p_at_h(enthalpyDM[s][m],pow(10.0,log_h_tabDM[1]), log_p_tabDM, 
                                           log_h_tabDM, n_tabDM, &n_nearestDM);
                     energyDM[s][m]=e_at_p(pressureDM[s][m],pow(10.0,log_p_tabDM[1]), log_e_tabDM, 
                                      log_p_tabDM, n_tabDM, &n_nearestDM, eos_typeDM,
                                       Gamma_PDM);

       //printf("%f %f %f %f %f %f %f %f \n",r_eq*s_gp[s]/(1-s_gp[s]),dif,dif2,r_eDM,rho_pole_hDM,gama_center_h,rho_center_h);
	           }
	           
	           
	    enthalpy[s][m]=enthalpy_min + 0.5*(SQ(r_eq)*(gama_pole_hDM+rho_pole_hDM
                           -gama[s][m]-rsm)-log(1.0-velocity_sq[s][m]));
  
            if((enthalpy[s][m]<=enthalpy_min) || (sgp>s_eDM_new)) {

                  pressure[s][m]=0.0;
                  energy[s][m]=0.0; 
	    }
	    
            else { 
            

                     pressure[s][m]=p_at_h(enthalpy[s][m],pow(10.0,log_h_tab[1]), log_p_tab, 
                                           log_h_tab, n_tab, &n_nearest);
                     energy[s][m]=e_at_p(pressure[s][m],pow(10.0,log_p_tab[1]), log_e_tab, 
                                      log_p_tab, n_tab, &n_nearest, eos_type,
                                       Gamma_P);
	           }
           }
           


            rho[s][m] *= SQ(r_eq);
            gama[s][m] *= SQ(r_eq);
            alpha[s][m] *= SQ(r_eq);

	 }

      }
      /* Compute metric potentials */
     // r_eq=sqrt(2*( h_center-enthalpy_min)/(gama_pole_h+rho_pole_h-gama_center_h-rho_center_h));
      S_gama = dmatrix(1,SDIV,1,MDIV);
      S_rho = dmatrix(1,SDIV,1,MDIV);
      S_omega = dmatrix(1,SDIV,1,MDIV);
      //printf("%5.4e \n",accuracy);//s_gp[SDIV]*r_e*sqrt(KAPPA)/((1.-s_gp[SDIV])*100000.));
      
      for(s=1;s<=SDIV;s++)
         for(m=1;m<=MDIV;m++) {
            rsm=rho[s][m];
            gsm=gama[s][m];
            omsm=omega[s][m];
            esm=energy[s][m];
            psm=pressure[s][m];
            esmDM=energyDM[s][m];
            psmDM=pressureDM[s][m];
            e_gsm=exp(0.5*gsm);
            e_rsm=exp(-rsm);
            e_gsmDM=exp(0.5*gsm*SQ(r_e_old*r_eDM/(r_eDM_old*r_eq)));
            e_rsmDM=exp(-rsm*SQ(r_e_old*r_eDM/(r_eDM_old*r_eq)));
            eaDM=16.0*PI*exp(2.0*alpha[s][m]*SQ(r_e_old*r_eDM/(r_eDM_old*r_eq)))*SQ(r_eq);
            v2sm=velocity_sq[s][m];
            v2smDM=velocity_sqDM[s][m];
            mum=mu[m];            
            m1=1.0-SQ(mum);
            sgp=s_gp[s];
            s_1=1.0-sgp;
            s1=sgp*s_1;
            s2=SQ(sgp/s_1);  

            ea=16.0*PI*exp(2.0*alpha[s][m])*SQ(r_eq);
 
            if(s==1) {
              d_gama_s=0.0;
              d_gama_m=0.0;
              d_rho_s=0.0;
              d_rho_m=0.0;
              d_omega_s=0.0;
              d_omega_m=0.0;
            }else{
                 d_gama_s=deriv_s(gama,s,m);
                 d_gama_m=deriv_m(gama,s,m);
                 d_rho_s=deriv_s(rho,s,m);
                 d_rho_m=deriv_m(rho,s,m);
                 d_omega_s=deriv_s(omega,s,m);
                 d_omega_m=deriv_m(omega,s,m);
	     }
          if((1)){
            S_rho[s][m] = e_gsm*(0.5*ea*(esm + psm)*s2*(1.0+v2sm)/(1.0-v2sm) + 0.5*ea*(esmDM + psmDM)*s2*(1.0+v2smDM)/(1.0-v2smDM)
  
                          + s2*m1*SQ(e_rsm)*(SQ(s1*d_omega_s) 
                       
                          + m1*SQ(d_omega_m))
                         
                          + s1*d_gama_s - mum*d_gama_m + 0.5*rsm*(ea*(psm+psmDM)*s2  
 
                          - s1*d_gama_s*(0.5*s1*d_gama_s+1.0) 
 
                          - d_gama_m*(0.5*m1*d_gama_m-mum)));

            S_gama[s][m] = e_gsm*(ea*(psm+psmDM)*s2 + 0.5*gsm*(ea*(psm+psmDM)*s2 - 0.5*SQ(s1

                           *d_gama_s) - 0.5*m1*SQ(d_gama_m)));

            S_omega[s][m]=e_gsm*e_rsm*( -ea*(Omega_hDM-omsm)*((esmDM+psmDM))

                          *s2/(1.0-v2smDM)  -ea*(Omega_h-omsm)*(esm+psm) *s2/(1.0-v2sm) 
                          
                          + omsm*( -0.5*ea*(((1.0+v2sm)*(esm) 
                           
                          + 2.0*v2sm*(psm))/(1.0-v2sm))*s2 
                          
                          -0.5*ea*(((1.0+v2smDM)*(esmDM) 
                           
                          + 2.0*v2smDM*(psmDM))/(1.0-v2smDM))*s2 

                          - s1*(2.0*d_rho_s+0.5*d_gama_s)

                          + mum*(2.0*d_rho_m+0.5*d_gama_m) + 0.25*SQ(s1)*(4.0

                          *SQ(d_rho_s)-SQ(d_gama_s)) + 0.25*m1*(4.0*SQ(d_rho_m)

                          - SQ(d_gama_m)) - m1*SQ(e_rsm)*(SQ(SQ(sgp)*d_omega_s)

                          + s2*m1*SQ(d_omega_m))));
             }else{
            S_rho[s][m] = e_gsm*(0.5*ea*(esmDM + psmDM)*s2*(1.0+v2sm)/(1.0-v2sm) + 0.5*ea*(esm + psm)*s2*(1.0+v2smDM)/(1.0-v2smDM)
  
                          + s2*m1*SQ(e_rsm)*(SQ(s1*d_omega_s) 
                       
                          + m1*SQ(d_omega_m))
                         
                          + s1*d_gama_s - mum*d_gama_m + 0.5*rsm*(ea*(psmDM+psm)*s2  
 
                          - s1*d_gama_s*(0.5*s1*d_gama_s+1.0) 
 
                          - d_gama_m*(0.5*m1*d_gama_m-mum)));

            S_gama[s][m] = e_gsm*(ea*(psm+psmDM)*s2 + 0.5*gsm*(ea*(psm+psmDM)*s2 - 0.5*SQ(s1

                           *d_gama_s) - 0.5*m1*SQ(d_gama_m)));

            S_omega[s][m]=e_gsm*e_rsm*( -ea*(Omega_hDM-omsm)*((esm+psm))

                          *s2/(1.0-v2smDM)  -ea*(Omega_h-omsm)*(esmDM+psmDM) *s2/(1.0-v2sm) 
                          
                          + omsm*( -0.5*ea*(((1.0+v2sm)*(esmDM) 
                           
                          + 2.0*v2sm*(psmDM))/(1.0-v2sm))*s2 
                          
                          -0.5*ea*(((1.0+v2smDM)*(esm) 
                           
                          + 2.0*v2smDM*(psm))/(1.0-v2smDM))*s2 

                          - s1*(2.0*d_rho_s+0.5*d_gama_s)

                          + mum*(2.0*d_rho_m+0.5*d_gama_m) + 0.25*SQ(s1)*(4.0

                          *SQ(d_rho_s)-SQ(d_gama_s)) + 0.25*m1*(4.0*SQ(d_rho_m)

                          - SQ(d_gama_m)) - m1*SQ(e_rsm)*(SQ(SQ(sgp)*d_omega_s)

                          + s2*m1*SQ(d_omega_m))));       
             
             
             
             }



	 }



      /* ANGULAR INTEGRATION */
   
      D1_rho = dmatrix(1,LMAX+1,1,SDIV);
      D1_gama = dmatrix(1,LMAX+1,1,SDIV);
      D1_omega = dmatrix(1,LMAX+1,1,SDIV);

      n=0;
      for(k=1;k<=SDIV;k++) {      

         for(m=1;m<=MDIV-2;m+=2) {
               sum_rho += (DM/3.0)*(P_2n[m][n+1]*S_rho[k][m]
                          + 4.0*P_2n[m+1][n+1]*S_rho[k][m+1] 
                          + P_2n[m+2][n+1]*S_rho[k][m+2]);
	 }

         D1_rho[n+1][k]=sum_rho;
         D1_gama[n+1][k]=0.0;
         D1_omega[n+1][k]=0.0;
         sum_rho=0.0;

      }

      for(n=1;n<=LMAX;n++)
         for(k=1;k<=SDIV;k++) {      
            for(m=1;m<=MDIV-2;m+=2) {

               sum_rho += (DM/3.0)*(P_2n[m][n+1]*S_rho[k][m]
                          + 4.0*P_2n[m+1][n+1]*S_rho[k][m+1] 
                          + P_2n[m+2][n+1]*S_rho[k][m+2]);
                       
               sum_gama += (DM/3.0)*(sin((2.0*n-1.0)*theta[m])*S_gama[k][m]
                           +4.0*sin((2.0*n-1.0)*theta[m+1])*S_gama[k][m+1]
                           +sin((2.0*n-1.0)*theta[m+2])*S_gama[k][m+2]);
  
               sum_omega += (DM/3.0)*(sin_theta[m]*P1_2n_1[m][n+1]*S_omega[k][m]
                            +4.0*sin_theta[m+1]*P1_2n_1[m+1][n+1]*S_omega[k][m+1]
                            +sin_theta[m+2]*P1_2n_1[m+2][n+1]*S_omega[k][m+2]);
	    }
            D1_rho[n+1][k]=sum_rho;
            D1_gama[n+1][k]=sum_gama;
            D1_omega[n+1][k]=sum_omega;
            sum_rho=0.0;
            sum_gama=0.0;
            sum_omega=0.0;
	}


      free_dmatrix(S_gama,1,SDIV,1,MDIV);
      free_dmatrix(S_rho,1,SDIV,1,MDIV);
      free_dmatrix(S_omega,1,SDIV,1,MDIV);



      /* RADIAL INTEGRATION */

      D2_rho = dmatrix(1,SDIV,1,LMAX+1);
      D2_gama = dmatrix(1,SDIV,1,LMAX+1);
      D2_omega = dmatrix(1,SDIV,1,LMAX+1);



      n=0;
      for(s=1;s<=SDIV;s++) {
            for(k=1;k<=SDIV-2;k+=2) { 
               sum_rho += (DS/3.0)*( f_rho[s][n+1][k]*D1_rho[n+1][k] 
                          + 4.0*f_rho[s][n+1][k+1]*D1_rho[n+1][k+1]
                          + f_rho[s][n+1][k+2]*D1_rho[n+1][k+2]);
 	    }
	    D2_rho[s][n+1]=sum_rho;
	    D2_gama[s][n+1]=0.0;
	    D2_omega[s][n+1]=0.0;
            sum_rho=0.0;
	 }

 
      for(s=1;s<=SDIV;s++)
         for(n=1;n<=LMAX;n++) {
            for(k=1;k<=SDIV-2;k+=2) { 
               sum_rho += (DS/3.0)*( f_rho[s][n+1][k]*D1_rho[n+1][k] 
                          + 4.0*f_rho[s][n+1][k+1]*D1_rho[n+1][k+1]
                          + f_rho[s][n+1][k+2]*D1_rho[n+1][k+2]);
 
               sum_gama += (DS/3.0)*( f_gama[s][n+1][k]*D1_gama[n+1][k] 
                           + 4.0*f_gama[s][n+1][k+1]*D1_gama[n+1][k+1]
                           + f_gama[s][n+1][k+2]*D1_gama[n+1][k+2]);
     
               if(k<s && k+2<=s) 
                 sum_omega += (DS/3.0)*( f_rho[s][n+1][k]*D1_omega[n+1][k] 
                              + 4.0*f_rho[s][n+1][k+1]*D1_omega[n+1][k+1]
                              + f_rho[s][n+1][k+2]*D1_omega[n+1][k+2]);
               else {
                 if(k>=s) 
                   sum_omega += (DS/3.0)*( f_gama[s][n+1][k]*D1_omega[n+1][k] 
                                + 4.0*f_gama[s][n+1][k+1]*D1_omega[n+1][k+1]
                                + f_gama[s][n+1][k+2]*D1_omega[n+1][k+2]);
                 else
                   sum_omega += (DS/3.0)*( f_rho[s][n+1][k]*D1_omega[n+1][k] 
                                + 4.0*f_rho[s][n+1][k+1]*D1_omega[n+1][k+1]
                                + f_gama[s][n+1][k+2]*D1_omega[n+1][k+2]);
               }
	    }
	    D2_rho[s][n+1]=sum_rho;
	    D2_gama[s][n+1]=sum_gama;
	    D2_omega[s][n+1]=sum_omega;
            sum_rho=0.0;
            sum_gama=0.0;
            sum_omega=0.0;
	 }
 
      free_dmatrix(D1_rho,1,LMAX+1,1,SDIV);
      free_dmatrix(D1_gama,1,LMAX+1,1,SDIV);
      free_dmatrix(D1_omega,1,LMAX+1,1,SDIV);


      /* SUMMATION OF COEFFICIENTS */

      for(s=1;s<=SDIV;s++) 
         for(m=1;m<=MDIV;m++) {

            gsm=gama[s][m];
            rsm=rho[s][m];
            omsm=omega[s][m];             
            e_gsm=exp(-0.5*gsm);
            e_rsm=exp(rsm);
            temp1=sin_theta[m];

            sum_rho += -e_gsm*P_2n[m][0+1]*D2_rho[s][0+1]; 

            for(n=1;n<=LMAX;n++) {

               sum_rho += -e_gsm*P_2n[m][n+1]*D2_rho[s][n+1]; 

               if(m==MDIV) {             
                 sum_omega += 0.5*e_rsm*e_gsm*D2_omega[s][n+1]; 
                 sum_gama += -(2.0/PI)*e_gsm*D2_gama[s][n+1];   
	       }
               else { 
                     sum_omega += -e_rsm*e_gsm*(P1_2n_1[m][n+1]/(2.0*n
                                  *(2.0*n-1.0)*temp1))*D2_omega[s][n+1];
  
                     sum_gama += -(2.0/PI)*e_gsm*(sin((2.0*n-1.0)*theta[m])
                                 /((2.0*n-1.0)*temp1))*D2_gama[s][n+1];   
	       }
	    }

            rho[s][m]=rsm + cf*(sum_rho-rsm);
            gama[s][m]=gsm + cf*(sum_gama-gsm);
            omega[s][m]=omsm + cf*(sum_omega-omsm);

            sum_omega=0.0;
            sum_rho=0.0;
            sum_gama=0.0; 
	  }

      free_dmatrix(D2_rho,1,SDIV,1,LMAX+1);
      free_dmatrix(D2_gama,1,SDIV,1,LMAX+1);
      free_dmatrix(D2_omega,1,SDIV,1,LMAX+1);


      /* CHECK FOR DIVERGENCE */

      if(fabs(omega[2][1])>100.0 || fabs(rho[2][1])>100.0 
         || fabs(gama[2][1])>300.0) {
         a_check=200; 
         break;
      }


      /* TREAT SPHERICAL CASE */
      
      if((r_ratio==1.0)&&(r_ratioDM==1.0)) {
        for(s=1;s<=SDIV;s++)
           for(m=1;m<=MDIV;m++) {
              rho[s][m]=rho[s][1];
              gama[s][m]=gama[s][1];
              omega[s][m]=0.0;          
	   }
      }
      

      /* TREAT INFINITY WHEN SMAX=1.0 */

      if(SMAX==1.0) {
         for(m=1;m<=MDIV;m++) {
            rho[SDIV][m]=0.0;
            gama[SDIV][m]=0.0;
            omega[SDIV][m]=0.0;
	 }
      } 

      
      /* COMPUTE FIRST ORDER DERIVATIVES OF GAMA */ 


      da_dm = dmatrix(1,SDIV,1,MDIV);
      dgds = dmatrix(1,SDIV,1,MDIV);
      dgdm = dmatrix(1,SDIV,1,MDIV); 
 
      for(s=1;s<=SDIV;s++)
         for(m=1;m<=MDIV;m++) {
            dgds[s][m]=deriv_s(gama,s,m);
            dgdm[s][m]=deriv_m(gama,s,m);
	 }



      /* ALPHA */
 
      if((r_ratio==1.0)&&(r_ratioDM==1.0)) {
        for(s=1;s<=SDIV;s++)
           for(m=1;m<=MDIV;m++)
              da_dm[s][m]=0.0; 
      } 
      else {
            for(s=2;s<=SDIV;s++)
               for(m=1;m<=MDIV;m++) {

                  da_dm[1][m]=0.0; 
       
                  sgp=s_gp[s];
                  s1=sgp*(1.0-sgp);
                  mum=mu[m]; 
                  m1=1.0-SQ(mum);
          
                  d_gama_s=dgds[s][m];
                  d_gama_m=dgdm[s][m];
                  d_rho_s=deriv_s(rho,s,m);
                  d_rho_m=deriv_m(rho,s,m);
                  d_omega_s=deriv_s(omega,s,m);
                  d_omega_m=deriv_m(omega,s,m);
                  d_gama_ss=s1*deriv_s(dgds,s,m)+(1.0-2.0*sgp)
                                               *d_gama_s;
                  d_gama_mm=m1*deriv_m(dgdm,s,m)-2.0*mum*d_gama_m;  
                  d_gama_sm=deriv_sm(gama,s,m);

           temp1=2.0*SQ(sgp)*(sgp/(1.0-sgp))*m1*d_omega_s*d_omega_m

                *(1.0+s1*d_gama_s) - (SQ(SQ(sgp)*d_omega_s) - 
 
                SQ(sgp*d_omega_m/(1.0-sgp))*m1)*(-mum+m1*d_gama_m); 
  
           temp2=1.0/(m1 *SQ(1.0+s1*d_gama_s) + SQ(-mum+m1*d_gama_m));

           temp3=s1*d_gama_ss + SQ(s1*d_gama_s);
  
           temp4=d_gama_m*(-mum+m1*d_gama_m);
   
           temp5=(SQ(s1*(d_rho_s+d_gama_s)) - m1*SQ(d_rho_m+d_gama_m))

                 *(-mum+m1*d_gama_m);

           temp6=s1*m1*(0.5*(d_rho_s+d_gama_s)* (d_rho_m+d_gama_m) 
  
                + d_gama_sm + d_gama_s*d_gama_m)*(1.0 + s1*d_gama_s); 

           temp7=s1*mum*d_gama_s*(1.0+s1*d_gama_s);

           temp8=m1*exp(-2.0*rho[s][m]);
 
          da_dm[s][m] = -0.5*(d_rho_m+d_gama_m) - temp2*(0.5*(temp3 - 

            d_gama_mm - temp4)*(-mum+m1*d_gama_m) + 0.25*temp5 

            - temp6 +temp7 + 0.25*temp8*temp1);	 
       }
   }

      for(s=1;s<=SDIV;s++) {
         alpha[s][1]=0.0;
         for(m=1;m<=MDIV-1;m++) 
            alpha[s][m+1]=alpha[s][m]+0.5*DM*(da_dm[s][m+1]+
                          da_dm[s][m]);
      } 
 

   free_dmatrix(da_dm,1,SDIV,1,MDIV);
   free_dmatrix(dgds,1,SDIV,1,MDIV);
   free_dmatrix(dgdm,1,SDIV,1,MDIV);


      for(s=1;s<=SDIV;s++)          
         for(m=1;m<=MDIV;m++) {   

            alpha[s][m] += -alpha[s][MDIV]+0.5*(gama[s][MDIV]-rho[s][MDIV]);


            if(alpha[s][m]>=300.0) {
              a_check=200; 
              break;
            }
            omega[s][m] /= r_eq;
         } 

      if(SMAX==1.0) {
         for(m=1;m<=MDIV;m++)      
            alpha[SDIV][m] = 0.0;
      }

      if(a_check==200)
        break;


      if(r_eq!=0){
        dif=fabs(r_e_old-r_eq)/r_eq;
      }else{
        dif=0;
      }      
      
      if(r_eDM!=0){
        dif2=fabs(r_eDM_old-r_eDM)/r_eDM;
      }else{
        dif2=0;
      }


      n_of_it++;

 }   /* end while */
//printf("N of it= %d\n", n_of_it);
//printf("%g %g %g %g %g\n",velocity_sq[1][1],alpha[1][1], gama[1][1], rho[1][1], omega[1][1]);

     /* COMPUTE OMEGA */  
 
    /* UPDATE r_e_new */
    if(cond){
      (*r_e_new) = r_eq;
      (*rDM_e_new) = r_eDM;
      (*Omega) = Omega_h*C/(r_eq*sqrt(KAPPA));
      (*OmegaDM) = Omega_hDM*C/(r_eq*sqrt(KAPPA));
    }else{
      (*r_e_new) = r_eDM;
      (*rDM_e_new) = r_eq;
      (*Omega) = Omega_h*C/(r_eDM*sqrt(KAPPA));
      (*OmegaDM) = Omega_hDM*C/(r_eq*sqrt(KAPPA));
    }


//printf("%g %g %g %g\n",log(Omega_h), log(Omega_hDM), log(r_ratio),  log(r_ratioDM));
    free_f3tensor(f_rho, 1,SDIV,1,LMAX+1,1,SDIV);
    free_f3tensor(f_gama,1,SDIV,1,LMAX+1,1,SDIV);
    free_dmatrix(P_2n,   1,MDIV,1,LMAX+1);   
    free_dmatrix(P1_2n_1,1,MDIV,1,LMAX+1);  
}



