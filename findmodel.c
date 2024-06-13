/************************************************************************** 
 *                    findmodel.c
 *   These routines are used to find special stellar models based on
 *         spin, mass, etc.
 *   Includes print and printpoly, MakeSphere, Surface,  Kepler, rns, 
 *         SetUpStar, and Set Spin
 *
 *************************************************************************/

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

#include "surface.h"

#include "findmodel.h"

/***************************************************************************
 * Printing Routine for tabulated stars. 
 ***************************************************************************/
									    
void print(double r_ratio,
	   double e_center, double Mass, double Mass_0, double R_e,
	   double Omega, double Omega_K, double J
	   )
{ double I_45;
  
  if( Omega == 0.0) I_45 = 0.0;
  else I_45 = J/(Omega*1.0e45);

  //printf("Someone is trying to print\n");

  printf(
	 "%5.4f \t%4.3f \t%4.3f \t%4.3f \t%4.3f \t%4.1f \t%5.1f \t%4.2f \t%4.3f \n",
	 r_ratio,
	 e_center,
	 Mass/MSUN,
	 Mass_0/MSUN,
	 R_e/1.0e5,
	 Omega/(2.0*PI),
	 Omega_K/(2.0*PI),
	 I_45,
	 ( C*J/(G*Mass*Mass)))
	 ;
}

/****************************************************************************
 * MakeSphere finds the spherical model with the same central energy density
 *       as the desired rotating star.
 ***************************************************************************/

int MakeSphere(EOS *eos,
	       NeutronStar *star,
	       double e_center,
	       EOSDM *eosDM,
	       NeutronStar *starDM,
	       double e_centerDM
)

{ 
  int  ierr; 
  double p_surface,p_surfaceDM;

  //  printf("Entering MakeSphere\n");
  /*
  if(strcmp(eos->eos_type,"poly") !=0) {
    star->e_surface=7.8*C*C*KSCALE;
    p_surface=1.01e8*KSCALE;
    star->enthalpy_min=1.0/(C*C);

    if (strcmp(eos->eos_type,"quark") == 0)
      star->e_surface=(eos->B)*(eos->K)*4.0/3.0;

  }
  else{
    star->e_surface=0.0;
    p_surface=0.0;
    star->enthalpy_min=0.0;

  }
  
  star->e_surface=pow(10,eos->log_e_tab[1]);
  p_surface=pow(10,eos->log_p_tab[1]);
  star->enthalpy_min=pow(10,eos->log_h_tab[1]);
    */
  star->e_surface=pow(10.,eos->log_e_tab[1]);
  p_surface=pow(10.,eos->log_p_tab[1]);
  star->enthalpy_min=pow(10.,eos->log_h_tab[1]);

    
  star->r_ratio = 1.00;
  star->Omega = 0.0;
  star->e_center = e_center;
 
 
  star->e_surfaceDM=pow(10.,eosDM->log_e_tabDM[1]);
  p_surfaceDM=pow(10.,eosDM->log_p_tabDM[1]);
  star->enthalpy_minDM=pow(10.,eosDM->log_h_tabDM[1]);
   
  star->e_centerDM = e_centerDM;

  /* CALCULATE THE PRESSURE AND ENTHALPY AT THE CENTRE OF THE STAR*/

  make_center(eos->eos_file, eos->log_e_tab, eos->log_p_tab, 
	      eos->log_h_tab, eos->log_n0_tab, eos->n_tab,eos->eos_type, 
	      eos->Gamma_P, 
	      star->e_center, &star->p_center, &star->h_center,star->e_surface,p_surface);
	      
  make_center(eosDM->eos_fileDM, eosDM->log_e_tabDM, eosDM->log_p_tabDM, 
	      eosDM->log_h_tabDM, eosDM->log_n0_tabDM, eosDM->n_tabDM,eosDM->eos_typeDM, 
	      eosDM->Gamma_PDM, 
	      star->e_centerDM, &star->p_centerDM, &star->h_centerDM,star->e_surfaceDM,p_surfaceDM);

  //printf("Entered make_center\n");


  /* COMPUTE A SPHERICAL STAR AS A FIRST GUESS FOR THE ROTATING STAR */

  sphere( star->metric.s_gp, eos->log_e_tab, eos->log_p_tab, 
	  eos->log_h_tab, eos->log_n0_tab, eos->n_tab, eos->eos_type, 
	  eos->Gamma_P, 
	  star->e_center, star->p_center, star->h_center, 
	  p_surface, star->e_surface, eosDM->log_e_tabDM, eosDM->log_p_tabDM, 
	  eosDM->log_h_tabDM, eosDM->log_n0_tabDM, eosDM->n_tabDM, eosDM->eos_typeDM, 
	  eosDM->Gamma_PDM, 
	  star->e_centerDM, star->p_centerDM, star->h_centerDM, 
	  p_surfaceDM, star->e_surfaceDM,
	  star->metric.rho, star->metric.gama, 
	  star->metric.alpha, star->metric.omega, 
	  &star->r_e,&star->r_eDM);


  //printf("Finished Sphere\n");

  /* THE PROCEDURE SPIN() WILL COMPUTE THE METRIC OF A STAR WITH
     GIVEN OBLATENESS. THE OBLATENESS IS SPECIFIED BY GIVING 
     THE RATIO OF THE LENGTH OF THE AXIS CONNECTING THE CENTRE OF THE STAR 
     TO ONE OF THE POLES TO THE RADIUS OF THE STAR'S EQUATOR. 
     THIS RATIO IS NAMED r_ratio.
     WHEN r_ratio = 1.0, THE STAR IS SPHERICAL */

//  ierr = rns( 1.0,  e_center, eos, star, e_centerDM, eosDM, starDM);

  //printf("Finished rns\n");

  return 0;
}

/****************************************************************************
 * Kepler finds the maximum spin rate for a star with the specified
 *        central energy density. It returns the kepler frequency.
 ***************************************************************************/

double Kepler( EOS *eos, NeutronStar *star,
	       double e_center, EOS *eosDM, NeutronStar *starDM,
	       double e_centerDM)
{

  int j, ierr;

  double r_ratio=1.00,
    dr=0.1,
    diff_Omega, old_diff_Omega;


  float ans, fh,fl,fm,fnew,sroot,xh,xl,xm,xnew,xacc=1e-4;

  

  /* THIS LOOP STARTS WITH A NON-ROTATING STAR AND INCREASES
     THE STAR'S OBLATENESS (BY DECREASING R_RATIO) AND 
     THEN CALCULATES THE STAR'S ANGULAR VELOCITY. ONCE THE
     COMPUTED VALUE OF ANGULAR VELOCITY IS LARGER THAN 
     THE ANGULAR VELOCITY OF A PARTICLE ORBITING THE STAR
     AT THE EQUATOR, (Omega_K), THE LOOP STOPS */
    
  diff_Omega = star->Omega_K - star->Omega;
  old_diff_Omega = diff_Omega;
  
  while( diff_Omega >0){
    /* Find the interval of r_ratio where the star has the
       correct angular velocity	*/
    r_ratio -= dr;
	
    /* Compute the star with the specified value of r_ratio	*/

    ierr = rns(r_ratio, e_center, eos, star, e_centerDM, eosDM, starDM);

    old_diff_Omega = diff_Omega;
    diff_Omega = star->Omega_K - star->Omega;
   
  } 

  /* The correct star lies between r_ratio and r_ratio + dr */
  xl = r_ratio;
  xh = r_ratio + dr;
  fl = diff_Omega;
  fh = old_diff_Omega;

  /* Use Ridder's method to find the correct star (Taken from 
     Numerical Recipes)	*/


  if ((fl > 0.0 && fh < 0.0) || (fl < 0.0 && fh > 0.0)) {
    ans=-1.11e30;
    for (j=1;j<=60;j++) {
      xm=0.5*(xl+xh);
      r_ratio = xm;

      ierr = rns( r_ratio, e_center, eos, star, e_centerDM, eosDM, starDM);

      fm= star->Omega_K - star->Omega;

      sroot=sqrt(fm*fm-fl*fh);
      if (sroot == 0.0) {
	r_ratio = ans;
	break;
      }
       
      xnew=xm+(xm-xl)*((fl >= fh ? 1.0 : -1.0)*fm/sroot);
      if (fabs(xnew-ans) <= xacc) {
	r_ratio = ans;
	break;
      }
      ans=xnew;
      r_ratio = ans;

      ierr = rns(r_ratio, e_center, eos, star, e_centerDM, eosDM, starDM);

      fnew =  star->Omega_K - star->Omega;

      if (fnew == 0.0){
	r_ratio = ans;
	break;
      }
       
      if (SIGN(fm,fnew) != fm) {
	xl=xm;
	fl=fm;
	xh=ans;
	fh=fnew;
      } else if (SIGN(fl,fnew) != fl) {
	xh=ans;
	fh=fnew;
      } else if (SIGN(fh,fnew) != fh) {
	xl=ans;
	fl=fnew;
      } else nrerror("never get here.");
      if (fabs(xh-xl) <= xacc){
	r_ratio = ans;
	break;
      }
    }
  }
  else {
    if (fh == 0.0){
      r_ratio +=dr;
    }
    nrerror("root must be bracketed in zriddr. (Kepler)");
  }
 
  /* THE RIDDER ZERO-FINDING ROUTINE HAS FOUND THE VALUE
     OF R_RATIO WHICH GIVES THE DESIRED STAR. */



  return (star->Omega/(2.0*PI));
}

/****************************************************************************
 * rns finds the rotating stellar model from the central energy density
 *        and the axis ratio.
 ***************************************************************************/

int rns(double r_ratio, 
	double e_center,
	EOS *eos,
	NeutronStar *star,
	double e_centerDM,
	EOSDM *eosDM,
	NeutronStar *starDM)
{
  int a_check = 0;
  int s,m;
  double cf=1.0;
  double accuracy=1e-5;
  
 

  //printf("Entering rns with r_ratio=%g\n",r_ratio);

  if(strcmp(eos->eos_type,"quark")==0) cf = 0.6;


  star->r_ratio = r_ratio;

  star->e_center = e_center;
  star->e_centerDM = e_centerDM;


  star->r_ratioDM = 1;
  double count=1;
  double dummy,dummyold= 1;
  double step=0.05;
  
  /*while((star->r_ratio!=1)&&((e_center>pow(10,eos->log_e_tab[1]))&&(e_centerDM>pow(10,eosDM->log_e_tabDM[1])))){
  starDM=star;
  spin(starDM->metric.s_gp, starDM->metric.mu, 
       eos->log_e_tab, eos->log_p_tab, eos->log_h_tab, eos->log_n0_tab, 
       eos->n_tab, eos->eos_type, eos->Gamma_P, 
       starDM->h_center, starDM->enthalpy_min,
       eosDM->log_e_tabDM, eosDM->log_p_tabDM, eosDM->log_h_tabDM, eosDM->log_n0_tabDM, 
       eosDM->n_tabDM, eosDM->eos_typeDM, eosDM->Gamma_PDM, 
       starDM->h_centerDM, starDM->enthalpy_minDM,
       starDM->metric.rho, starDM->metric.gama, starDM->metric.alpha, starDM->metric.omega, 
       starDM->energy, starDM->pressure, 
       starDM->enthalpy, starDM->velocity_sq,
       starDM->energyDM, starDM->pressureDM, 
       starDM->enthalpyDM, starDM->velocity_sqDM,
       a_check, accuracy, cf,
       starDM->r_ratio, starDM->r_ratioDM, &starDM->r_e, &starDM->r_eDM, &starDM->Omega, &starDM->OmegaDM);
   //printf("%f %f %f %f \n", starDM->r_ratioDM,starDM->OmegaDM,starDM->r_ratio,starDM->Omega);
   
     if(isnan((starDM->OmegaDM))&&isnan((starDM->Omega))){
       //printf("The system is unstable\n"); 
       //exit(0);
       return;
     }
   
      if(fabs(starDM->OmegaDM-starDM->Omega)/starDM->Omega<0.0005){
         star->r_ratioDM = starDM->r_ratioDM;
         break;
      }
       if(starDM->OmegaDM<starDM->Omega){
             dummy=1;
             if(dummy==dummyold){
                star->r_ratioDM -= step/count;
             }
             else{
                count++;
                star->r_ratioDM -= step/count;
             }
       }else{
             dummy=2;
             if(dummy==dummyold){
                star->r_ratioDM += step/count;
             }
             else{
                count++;
                star->r_ratioDM += step/count;
             }
       }
       
       dummyold=dummy;
      //printf("%f \n",fabs(starDM->OmegaDM-star->Omega)/starDM->Omega);
  }*/

  spin(star->metric.s_gp, star->metric.mu, 
       eos->log_e_tab, eos->log_p_tab, eos->log_h_tab, eos->log_n0_tab, 
       eos->n_tab, eos->eos_type, eos->Gamma_P, 
       star->h_center, star->enthalpy_min,
       eosDM->log_e_tabDM, eosDM->log_p_tabDM, eosDM->log_h_tabDM, eosDM->log_n0_tabDM, 
       eosDM->n_tabDM, eosDM->eos_typeDM, eosDM->Gamma_PDM, 
       star->h_centerDM, star->enthalpy_minDM,
       star->metric.rho, star->metric.gama, star->metric.alpha, star->metric.omega, 
       star->energy, star->pressure, 
       star->enthalpy, star->velocity_sq,
       star->energyDM, star->pressureDM, 
       star->enthalpyDM, star->velocity_sqDM,
       a_check, accuracy, cf,
       star->r_ratio, star->r_ratioDM, &star->r_e, &star->r_eDM, &star->Omega, &star->OmegaDM);
  

 
  mass_radius( star->metric.s_gp, star->metric.mu, 
	       eos->log_e_tab, eos->log_p_tab, eos->log_h_tab, eos->log_n0_tab, 
	       eos->n_tab, eos->eos_type, eos->Gamma_P, 
               eosDM->log_e_tabDM, eosDM->log_p_tabDM, eosDM->log_h_tabDM, eosDM->log_n0_tabDM, 
	       eosDM->n_tabDM, eosDM->eos_typeDM, eosDM->Gamma_PDM, 
	       star->metric.rho, star->metric.gama, 
	       star->metric.alpha, star->metric.omega, 
	       star->energy, star->pressure, 
	       star->enthalpy, star->velocity_sq,
	       star->energyDM, star->pressureDM, 
	       star->enthalpyDM, star->velocity_sqDM,
	       star->r_ratio,star->r_ratioDM,&star->Ratio_sch, star->e_surface, star->e_surfaceDM, star->r_e, star->r_eDM, star->Omega,
	       &star->Mass, &star->Mass_0, &star->ang_mom, &star->R_e, 
	       &star->MassDM, &star->Mass_0DM, &star->ang_momDM, &star->R_eDM, 
	       star->v_plus, star->v_minus, 
	       &star->Omega_K, &star->Vp, &star->Mp);

  for (m=1;m<=MDIV;m++)
    for (s=1;s<=SDIV;s++){
  
      star->metric.alpha_s[s][m] = deriv_s(star->metric.alpha,s,m);
      star->metric.alpha_mu[s][m] = deriv_m(star->metric.alpha,s,m);
      star->metric.alpha_ms[s][m] = deriv_sm(star->metric.alpha,s,m);

      star->metric.rho_s[s][m] = deriv_s(star->metric.rho,s,m);
      star->metric.rho_mu[s][m] = deriv_m(star->metric.rho,s,m);
      star->metric.rho_ms[s][m] = deriv_sm(star->metric.rho,s,m);

      star->metric.gama_s[s][m] = deriv_s(star->metric.gama,s,m);
      star->metric.gama_mu[s][m] = deriv_m(star->metric.gama,s,m);
      star->metric.gama_ms[s][m] = deriv_sm(star->metric.gama,s,m);

      star->metric.omega_s[s][m] = deriv_s(star->metric.omega,s,m);
      star->metric.omega_mu[s][m] = deriv_m(star->metric.omega,s,m);
      star->metric.omega_ms[s][m] = deriv_sm(star->metric.omega,s,m);

    }


  //printf("Finished rns with r_ratio=%g spin=%g \n",r_ratio, star->Omega/(2.0*PI));

  return 0;
}

/****************************************************************************
 * SetUpStar loads in the eos and sets up grid.
 ***************************************************************************/

int SetUpStar( char eos_file[80],
	       char eos_type[10],
	       char eos_fileDM[80],
	       char eos_typeDM[80],
	       char data_dir[80],
	       char data_dirDM[80],
	       double Gamma_P,
	       double B,
	       double K,
	       EOS *eos,
	       NeutronStar *star,
	       EOSDM *eosDM,
	       NeutronStar *starDM)
{

  double
   **rho,                          /* potential \rho */ 
   **gama,                         /* potential \gamma */ 
   **omega,                        /* potential \omega */ 
   **alpha,                        /* potential \alpha */ 
   **rho_s,
   **gama_s,
   **alpha_s,
   **omega_s,
   **rho_mu,
   **gama_mu,
   **alpha_mu,
   **omega_mu,
   **rho_ms,
   **gama_ms,
   **alpha_ms,
   **omega_ms,
   **energy,                       /* energy density \epsilon */
   **pressure,                     /* pressure */ 
   **enthalpy,                     /* enthalpy */
   **velocity_sq,                  /* square of velocity */ 
   **energyDM,                       /* energy density \epsilon */
   **pressureDM,                     /* pressure */ 
   **enthalpyDM,                     /* enthalpy */
   **velocity_sqDM,                  /* square of velocity */ 
   *s_surf,
   *r_is_surf,
   *r_surf,
    *gravity_surf,
   *v_plus,			/* vel. of co-rot. particle wrt ZAMO */
   *v_minus,			/* vel. of counter-rot. ... */
   *dpot_dr_dr;

  sprintf(eos->eos_file,"%s",eos_file);
  sprintf(eos->eos_type,"%s",eos_type);
  sprintf(eos->data_dir,"%s",data_dir);


  sprintf(eosDM->eos_fileDM,"%s",eos_fileDM);
  sprintf(eosDM->eos_typeDM,"%s",eos_typeDM);
  sprintf(eosDM->data_dirDM,"%s",data_dirDM);


  eos->Gamma_P = Gamma_P;
  eos->B = B;
  eos->K = K;

  /* LOAD TABULATED EOS */ 
  if(strcmp(eos_type,"tab")==0) 
    load_eos( eos->eos_file, eos->log_e_tab, eos->log_p_tab, 
	      eos->log_h_tab, eos->log_n0_tab, &eos->n_tab );
  if(strcmp(eos_typeDM,"DM")==0) 
    load_eos( eosDM->eos_fileDM, eosDM->log_e_tabDM, eosDM->log_p_tabDM, 
	      eosDM->log_h_tabDM, eosDM->log_n0_tabDM, &eosDM->n_tabDM );

  /* SET UP GRID */
  make_grid(star->metric.s_gp, star->metric.mu);



  /* ALLLOCATE MEMORY */

  rho = dmatrix(1,SDIV,1,MDIV);
  gama = dmatrix(1,SDIV,1,MDIV);
  alpha = dmatrix(1,SDIV,1,MDIV);
  omega = dmatrix(1,SDIV,1,MDIV);

  rho_s = dmatrix(1,SDIV,1,MDIV);
  gama_s = dmatrix(1,SDIV,1,MDIV);
  alpha_s = dmatrix(1,SDIV,1,MDIV);
  omega_s = dmatrix(1,SDIV,1,MDIV);

  rho_mu = dmatrix(1,SDIV,1,MDIV);
  gama_mu = dmatrix(1,SDIV,1,MDIV);
  alpha_mu = dmatrix(1,SDIV,1,MDIV);
  omega_mu = dmatrix(1,SDIV,1,MDIV);

  rho_ms = dmatrix(1,SDIV,1,MDIV);
  gama_ms = dmatrix(1,SDIV,1,MDIV);
  alpha_ms = dmatrix(1,SDIV,1,MDIV);
  omega_ms = dmatrix(1,SDIV,1,MDIV);

  energy = dmatrix(1,SDIV,1,MDIV);
  pressure = dmatrix(1,SDIV,1,MDIV);
  enthalpy = dmatrix(1,SDIV,1,MDIV);
  velocity_sq = dmatrix(1,SDIV,1,MDIV);

  s_surf = dvector(1,MDIV);
  r_is_surf = dvector(1,MDIV);
  r_surf = dvector(1,MDIV);
  gravity_surf = dvector(1,MDIV);

  v_plus = dvector(1,SDIV);
  v_minus = dvector(1,SDIV);
  dpot_dr_dr = dvector(1,SDIV);

  energyDM = dmatrix(1,SDIV,1,MDIV);
  pressureDM = dmatrix(1,SDIV,1,MDIV);
  enthalpyDM = dmatrix(1,SDIV,1,MDIV);
  velocity_sqDM = dmatrix(1,SDIV,1,MDIV);

/*  s_surfDM = dvector(1,MDIV);
  r_is_surfDM = dvector(1,MDIV);
  r_surfDM = dvector(1,MDIV);
  gravity_surfDM = dvector(1,MDIV);

  v_plusDM = dvector(1,SDIV);
  v_minusDM = dvector(1,SDIV);
  dpot_dr_drDM = dvector(1,SDIV);
*/
  /* PUT REST OF MODEL INFORMATION INTO star */
  star->metric.alpha = alpha;
  star->metric.gama = gama;
  star->metric.rho = rho;
  star->metric.omega = omega;
  star->metric.alpha_s = alpha_s;
  star->metric.gama_s = gama_s;
  star->metric.rho_s = rho_s;
  star->metric.omega_s = omega_s;
  star->metric.alpha_mu = alpha_mu;
  star->metric.gama_mu = gama_mu;
  star->metric.rho_mu = rho_mu;
  star->metric.omega_mu = omega_mu;
  star->metric.alpha_ms = alpha_ms;
  star->metric.gama_ms = gama_ms;
  star->metric.rho_ms = rho_ms;
  star->metric.omega_ms = omega_ms;

  star->energy = energy;
  star->pressure = pressure;
  star->enthalpy = enthalpy;
  star->velocity_sq = velocity_sq;

  star->energyDM = energyDM;
  star->pressureDM = pressureDM;
  star->enthalpyDM = enthalpyDM;
  star->velocity_sqDM = velocity_sqDM;

  star->v_plus = v_plus;
  star->v_minus = v_minus;
  star->dpot_dr_dr = dpot_dr_dr;

  star->s_surf = s_surf;
  star->r_is_surf = r_is_surf;
  star->r_surf = r_surf;
  star->gravity_surf = gravity_surf;


  return 0;

}

/****************************************************************************
 *  SetSpin finds a star with the specified central energy and the
 *          specified spin rate.
 ***************************************************************************/

int SetSpin( EOS *eos, NeutronStar *star,
	     double e_center,EOS *eosDM, NeutronStar *starDM,
	     double e_centerDM, double spinfreq)
{

  int j, ierr;

  double r_ratio=1.0,
    dr=0.02,
    diff_Omega, old_diff_Omega;


  float ans, fh,fl,fm,fnew,sroot,xh,xl,xm,xnew,xacc=1e-5;

  double newomega;

  //printf("entering SetSpin spin=%lf\n",spinfreq);

  newomega = spinfreq*2.0*PI;
  
 

  diff_Omega = newomega - star->Omega;
  old_diff_Omega = diff_Omega;
  
  /*printf("r_ratio=%g diff_Omega=%6.5e\n",r_ratio, diff_Omega);
  print(r_ratio,e_center, star->Mass, star->Mass_0, 
	  star->R_e, star->Omega, star->Omega_K, star->ang_mom);
  */


  while( diff_Omega >0){
    /* Find the interval of r_ratio where the star has the
       correct angular velocity	*/
    r_ratio -= dr;
	
    /* Compute the star with the specified value of r_ratio	*/

    ierr = rns(r_ratio, e_center, eos, star, e_centerDM, eosDM, starDM);
    /*
    printf("r_ratio=%g diff_Omega=%6.5e\n",r_ratio, diff_Omega);

    print(r_ratio,e_center, star->Mass, star->Mass_0, 
    star->R_e, star->Omega, star->Omega_K, star->ang_mom);*/
 
    old_diff_Omega = diff_Omega;
    diff_Omega = newomega - star->Omega;
   
  } 



  /* The correct star lies between r_ratio and r_ratio + dr */
  xl = r_ratio;
  xh = r_ratio + dr;
  fl = diff_Omega;
  fh = old_diff_Omega;

  /* Use Ridder's method to find the correct star (Taken from 
     Numerical Recipes)	*/


  if ((fl > 0.0 && fh < 0.0) || (fl < 0.0 && fh > 0.0)) {
    ans=-1.11e30;
    for (j=1;j<=60;j++) {
      xm=0.5*(xl+xh);
      r_ratio = xm;

      ierr = rns( r_ratio, e_center, eos, star, e_centerDM, eosDM, starDM);
      fm= newomega - star->Omega;

      sroot=sqrt(fm*fm-fl*fh);
      if (sroot == 0.0) {
	r_ratio = ans;
	break;
      }
       
      xnew=xm+(xm-xl)*((fl >= fh ? 1.0 : -1.0)*fm/sroot);
      if (fabs(xnew-ans) <= xacc) {
	r_ratio = ans;
	ierr = rns(r_ratio, e_center, eos, star, e_centerDM, eosDM, starDM);
	/*
	print(r_ratio,e_center, star->Mass, star->Mass_0, 
	star->R_e, star->Omega, star->Omega_K, star->ang_mom);*/

	break;
      }

      ans=xnew;
      r_ratio = ans;

      ierr = rns(r_ratio, e_center, eos, star, e_centerDM, eosDM, starDM);
      /*
      print(r_ratio,e_center, star->Mass, star->Mass_0, 
      star->R_e, star->Omega, star->Omega_K, star->ang_mom);*/

      fnew =  newomega - star->Omega;

      if (fnew == 0.0){
	r_ratio = ans;
	break;
      }
       
      if (SIGN(fm,fnew) != fm) {
	xl=xm;
	fl=fm;
	xh=ans;
	fh=fnew;
      } else if (SIGN(fl,fnew) != fl) {
	xh=ans;
	fh=fnew;
      } else if (SIGN(fh,fnew) != fh) {
	xl=ans;
	fl=fnew;
      } else nrerror("never get here.");
      if (fabs(xh-xl) <= xacc){
	r_ratio = ans;
	break;
      }
    }
  }
  else {
    if (fh == 0.0){
      r_ratio +=dr;
    }
    nrerror("root must be bracketed in zriddr. (SetSpin)");
  }
 
  /* THE RIDDER ZERO-FINDING ROUTINE HAS FOUND THE VALUE
     OF R_RATIO WHICH GIVES THE DESIRED STAR. */

  // printf("Finished Set Spin\n");

  return 0;
}

/****************************************************************************
 *  SetJ finds a star with the specified moment of inertia and the
 *          central energy density
 ***************************************************************************/

int SetJ( EOS *eos, NeutronStar *star,
	     double e_center,EOS *eosDM, NeutronStar *starDM,
	     double e_centerDM, double J)
{ /*J is given in cgs units*/
  int j, ierr;
  double r_ratio=1.0,
    dr=0.01,
    diff=0, old_diff=0;
   
  double ans, fh,fl,fm,fnew,sroot,xh,xl,xm,xnew,xacc=1e-4;

    printf("\t\t\t\t\t\tin set J ep = %6.5e, J = %6.5e, ang_mom = %6.5e, r_ratio = %6.5e\n", e_center, J, star->ang_mom, r_ratio);  

  diff = J - star->ang_mom;
  old_diff = diff;

  while( diff >0){

  printf("\t\t\t\t\t\tlooping in set J ep = %6.5e, J = %6.5e, ang_mom = %6.5e, r_ratio = %6.5e\n", e_center, J, star->ang_mom, r_ratio);  
    /* Find the interval of r_ratio where the star has the
       correct angular velocity	*/
    r_ratio -= dr;
	
    /* Compute the star with the specified value of r_ratio	*/

    ierr = rns(r_ratio, e_center, eos, star, e_centerDM, eosDM, starDM);

    old_diff = diff;
    diff = J - star->ang_mom;
 
  } 

  /* The correct star lies between r_ratio and r_ratio + dr */
  xl = r_ratio;
  xh = r_ratio + dr;
  fl = diff;
  fh = old_diff;

  /* Use Ridder's method to find the correct star (Taken from 
     Numerical Recipes)	*/

  if ((fl > 0.0 && fh < 0.0) || (fl < 0.0 && fh > 0.0)) {
    ans=-1.11e30;
    for (j=1;j<=60;j++) {

      xm=0.5*(xl+xh);
      r_ratio = xm;

      ierr = rns( r_ratio, e_center, eos, star, e_centerDM, eosDM, starDM);

      fm= J - star->ang_mom;
      sroot=sqrt(fm*fm-fl*fh);
      if (sroot == 0.0) {
	r_ratio = ans;
	break;
      }
       
      xnew=xm+(xm-xl)*((fl >= fh ? 1.0 : -1.0)*fm/sroot);
      if (fabs(xnew-ans) <= xacc) {
	r_ratio = ans;
	break;
      }
      ans=xnew;
      r_ratio = ans;

      ierr = rns(r_ratio, e_center, eos, star, e_centerDM, eosDM, starDM);

      fnew =  J - star->ang_mom;
      if (fnew == 0.0){
	r_ratio = ans;
	break;
      }
       
      if (SIGN(fm,fnew) != fm) {
	xl=xm;
	fl=fm;
	xh=ans;
	fh=fnew;
      } else if (SIGN(fl,fnew) != fl) {
	xh=ans;
	fh=fnew;
      } else if (SIGN(fh,fnew) != fh) {
	xl=ans;
	fl=fnew;
      } else nrerror("never get here.");
      if (fabs(xh-xl) <= xacc){
	r_ratio = ans;
	break;
      }
    }
  }
  else {
    if (fh == 0.0){
      r_ratio +=dr;
    }
    nrerror("root must be bracketed in zriddr.");
  }
 
  /* THE RIDDER ZERO-FINDING ROUTINE HAS FOUND THE VALUE
     OF R_RATIO WHICH GIVES THE DESIRED STAR. */

  ierr = rns(r_ratio, e_center, eos, star, e_centerDM, eosDM, starDM);
  printf("\t\t\t\t\t\tlooping in set J ep = %6.5e, J = %6.5e, ang_mom = %6.5e, r_ratio = %6.5e\n", 
	 e_center, J, star->ang_mom, r_ratio);  

  return 0;
}


/*  MaxMass finds the Maximum Mass non-rotating star. */

int MaxMass(double e_min, double e_max, EOS *eos, NeutronStar *star,  EOS *eosDM, NeutronStar *starDM, char name[80])
{

  int i, n=40;
  int ierr=0;
  int imax;
  int flag=0;

  double evec[100], mvec[100];
  double e_center,e_centerDM, e_step;
  double maxmass, xx;

  double ea, eb, ec, ma, mb, mc;
  double e0,e1,e2,e3, m1,m2;
  double gold=0.61803399;
  double tol=1e-4;

  FILE *output, *franges;
  franges = fopen("E_range.txt", "a");

  printf("Entering maxmass\n");

  output = fopen(name,"w");

  e_step = (e_max-e_min)/(1.0*n);

  /* Bracket the Maximum */

  maxmass = 0.0;
  imax = 0;

  fprintf(output,"#e_c   Mass   Mass_0   Radius    M/R \n");
  printf("#e_c   Mass   Mass_0   Radius    M/R \n");

  for (i=0;i<=n;i++){

    e_center = e_min + i*e_step;
    ierr = MakeSphere( eos, star, e_center, eosDM, starDM, e_centerDM);

    printf("MadeSphere\n");

    evec[i] = e_center;
    mvec[i] = star->Mass/MSUN;

    if (mvec[i] > maxmass){
      maxmass = mvec[i];
      imax = i;
    }

    if (mvec[i] < maxmass) flag=1;
    printf("i=%d  e=%g  Mass=%g\n",
	   i, evec[i],mvec[i]);

    //    fprintf(output,"%g %g \n",star->R_e*1e-5,star->Mass/MSUN);


    xx =  star->Mass*G/star->R_e*pow(C,-2.);
    

   fprintf(output,"%g %g %g %g %g \n",
	    e_center, star->Mass/MSUN, star->R_e*1e-5, xx,
	   e_center*4.0*PI*pow(star->R_e,3.)/(3.*star->Mass)*1e15);


  }

 

  if (flag==1){
    /* Store the bracketed values: 
       a: values before maximum
       b: values before maximum
       c: values after maximum
    */
    ea = evec[imax-1];
    ma = mvec[imax-1];

    eb = evec[imax];
    mb = mvec[imax];
    
    ec = evec[imax+1];
    mc = mvec[imax+1];

    printf("Maximum mass bracketed at: i=%d e=%g Mass=%g \n",
	   imax, evec[imax], mvec[imax]);

    /* Now use Golden Ratio method to find the maximum. */

    e0=ea;
    e3=ec;

    if (fabs(ec-eb) > fabs(eb-ea)){
      e1=eb;
      m1=mb;
      e2=eb + (1.0-gold)*(ec-eb);
      ierr = MakeSphere( eos, star, e2,eosDM, starDM, e2);
      m2 = star->Mass/MSUN;
    }
    else{
      e2=eb;
      m2=mb;
      e1=eb - (1.0-gold)*(eb-ea);
      ierr = MakeSphere( eos, star, e1, eosDM, starDM, e1);
      m1 = star->Mass/MSUN;
    }

    //printf("e1=%g m1=%g\n",e1,m1);
    //printf("e2=%g m2=%g\n",e2,m2);


    while ( fabs(e3-e0) > tol*(e1+e2)){

      if (m2>m1){
	e0=e1;
	e1=e2;
	e2=gold*e1+(1.0-gold)*e3;

	m1=m2;
	ierr = MakeSphere( eos, star, e2,eosDM, starDM, e2);
	m2 = star->Mass/MSUN;

	//printf("e1=%g m1=%g\n",e1,m1);
	//printf("e2=%g m2=%g\n",e2,m2);
      }
      else{
	e3=e2;
	e2=e1;
	e1=gold*e2+(1.0-gold)*e0;

	m2=m1;
	ierr = MakeSphere( eos, star, e1,eosDM, starDM, e1);
	m1 = star->Mass/MSUN;

	//	printf("e1=%g m1=%g\n",e1,m1);
	//printf("e2=%g m2=%g\n",e2,m2);

      }
    } /* end while loop */

    /* Maximum mass found */

    if (m1>m2)
      e_max = e1;
    else
      e_max = e2;

    ierr = MakeSphere( eos, star, e_max,eosDM, starDM, e_max);

    xx =  star->Mass*G/star->R_e*pow(C,-2.);
    
  
 
  fprintf(output,"%g %g %g %g %g \n",
	    e_max, star->Mass/MSUN, star->R_e*1e-5, xx,
	   e_max*4.0*PI*pow(star->R_e,3.)/(3.*star->Mass)*1e15);




    printf("Maximum Mass at: e=%g e15    m=%g Msun\n",e_max, star->Mass/MSUN);
    printf("Baryon Mass = %g Msun\n",star->Mass_0/MSUN);
    printf("Radius = %g  zeta=M/R=%g \n",star->R_e*1e-5, G*star->Mass/(star->R_e*C*C));

    fprintf(franges, "%g %g %g %g \n", e_max, star->Mass/MSUN, star->Mass_0/MSUN, star->R_e*1e-5);


  }
  else
    {
      printf("Maximum mass not bracketed!! Use a larger maximum mass.\n");
      ierr = 1;
      fprintf(franges, "Modify range of energy \n");
    }



  fclose(output);
  fclose(franges);

  return ierr;

}




/*  SetMassRatio computes a star with a specified value of r_ratio
    and a specified value of Rest Mass. */

int SetMassRatio(double Ratio,
                double Mass,
                double e_center,
                EOS *eos,
                NeutronStar *star,
                double e_centerDM,
                EOS *eosDM,
                NeutronStar *starDM,
                int printflag)
{

  int j, n;
  int ierr=0;

 double evec[100], mvec[100];
 double diff, olddiff, err, tol, sign;
 double estep;
 double olde, newe, oldm, newm;

  estep = 0.05;
  err = 1.0;
  tol = 1e-2;
  j=0;
  sign=1.0;

  printf("\nBeginning of SetMassRatio  Mass = %g  r_ratio=%g \n",Mass,Ratio);

  while ((err > tol) & (j<1)) {

    /* First Bracket the root */

    n=1;
    diff = sign*1.0;
    olddiff = diff;

    newe=0;
    newm=0;

    while ( olddiff*diff > 0  ){

      olde = newe;
      oldm = newm;

      ierr = MakeSphere( eos, star, e_center,eosDM, starDM, e_centerDM);
      if (Ratio<0.8) ierr = rns(0.8, e_center, eos, star, e_centerDM, eosDM, starDM);
      if (Ratio<0.7) ierr = rns(0.7, e_center, eos, star, e_centerDM, eosDM, starDM);
      if (Ratio<0.6) ierr = rns(0.6, e_center, eos, star, e_centerDM, eosDM, starDM);
      ierr = rns(Ratio, e_center, eos, star, e_centerDM, eosDM, starDM);

      printf("e_center = %g  mass_0=%g \n",
	     e_center, star->Mass_0/MSUN);

      newe=e_center;
      newm=star->Mass_0/MSUN;

      olddiff = diff;
      diff = Mass - star->Mass_0/MSUN;

      if ( j==0 && n==1 && diff < 0 ){
        sign = -1.0;
        estep *= sign;
        olddiff *= sign;
      }

      e_center += estep;
      n++;
    }

    /* root is bracketed by olde and newe */

    evec[1] = olde;
    mvec[1] = oldm;
    evec[2] = newe;
    mvec[2] = newm;

    n=2;

    e_center = printpolint( mvec, evec, n, Mass, &err,printflag);

    ierr = MakeSphere( eos, star, e_center,eos, starDM, e_centerDM);
    if (Ratio<0.8) ierr = rns(0.8, e_center, eos, star, e_centerDM, eosDM, starDM);
    if (Ratio<0.7) ierr = rns(0.7, e_center, eos, star, e_centerDM, eosDM, starDM);
    if (Ratio<0.6) ierr = rns(0.6, e_center, eos, star, e_centerDM, eosDM, starDM);
    rns(Ratio, e_center, eos, star, e_centerDM, eosDM, starDM);

    if (printflag)
      if(strcmp(eos->eos_type,"poly")!=0)
        print(Ratio,e_center, star->Mass, star->Mass_0,
            star->R_e, star->Omega, star->Omega_K, star->ang_mom);
      

    err = fabs(star->Mass_0/MSUN - Mass)/Mass;
    printf("error = %6.5e \n",err);

    if (err > tol){
      //if (n>2)
      //e_center = evec[n-2];
      //else e_center = evec[1];
      estep *= 0.1;
      j++;
      if (printflag)     printf("Trying new level of refinement: estep=%g\n",estep);
    }


  }


  if (star->Omega > star->Omega_K) ierr = 1;
  else{
    if (j>0) ierr = 2;
    if ( fabs(star->Omega_K-star->Omega)/star->Omega_K  <= 1e-3) ierr=-1;
  }

  // printf("ierr=%d\n",ierr);

  return ierr;

}









double printpolint(double *xp, double *yp, int order, double xb,
                   double *err, int printflag){
  int i, m, ns=1, j;
  double tmp, diff, den, dnum, cnum, yb;
  double *c, *d;

  if(printflag) printf("printpolint: xb=%g\n",xb);

  c = (double *) malloc((order+1)*sizeof(double));
  d = (double *) malloc((order+1)*sizeof(double));
  diff = fabs(xb-xp[1]);
  for (i=1; i<=order; i++){
    if (printflag) printf("x[%d]=%g   y[%d]=%g \n",
                          i,xp[i],i,yp[i]);
    if ( (tmp=fabs(xb-xp[i])) < diff ){
      ns=i;
      diff = tmp;
    }
    c[i] = yp[i];
    d[i] = yp[i];
  }

  yb = yp[ns--];
  for (m=1; m<order; m++){
    for (i=1; i<=order-m; i++){
      cnum = xp[i]-xb;
      dnum = xp[i+m]-xb;
      if ( (den=cnum-dnum) == 0.0 ){
        /*Two values of xp are equal:no polynomial passes through the points.*/
        printf("error in printpolint: xp[%d]==xp[%d]\n", i, i+m);
        printf("xb=%g\n",xb);
        for(j=1;j<=order;j++){
          printf("xp[%d]=%g    yp[%d]=%g\n",j,xp[j],j,yp[j]);
        }
        exit(1);
      }
      tmp = (c[i+1]-d[i])/den;
      c[i] = cnum*tmp;
      d[i] = dnum*tmp;
    }
    *err = 2*ns<order-m ? c[ns+1] : d[ns--];
    yb += *err;
  }
  free(c);
  free(d);
  if (printflag) printf("interpolated answer: yb=%g  err=%g\n",
                        yb,*err);
  return yb;
}
