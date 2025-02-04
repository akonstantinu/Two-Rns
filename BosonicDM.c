/***************************************************************************/
/*                             HnG.C                                       */
/*                                                                         */
/* Compile: cc -O HnG.c -lm -o HnG                                         */
/*                                                                         */
/* Usage: HnG  e_vs_P_filename   number_of_points_in_eos_file              */
/*                                                                         */
/* Output: number_of_points_in_eos_file                                    */
/*         energy_density/c^2 Pressure enthalpy number_density             */
/*         (cgs units)                                                     */
/*                                                                         */
/* Nikolaos Stergioulas,   Jan. 1997                                       */ 
/***************************************************************************/
   
#include <stdio.h>
#include <string.h>
#include <math.h>

 
#define PDIV 16003                /* divisions used in Simpson's integration */
#define PDIV_1 (PDIV-1)
#define NDIV 400                  /* number of output points */
#define G 6.6732e-8               /* gravitational constant */ 
#define KAPPA (1.0e-15*C*C/G)     /* scaling factor */
#define KSCALE (KAPPA*G/(C*C*C*C))
#define C 2.99792458e10
#define  hbar 1.05457182*pow(10,-27)
/*#define MB 1.659e-24 */
#define MB 1.66e-24               /* baryon mass in cgs */
#define MEVFM3 1.60217657e33      /* 1 Mev/fm^3 in dynes/cm^2 */
#define MEVC2FM3 1.7827e12      /* 1 Mev/fm^3 in g/cm^3 */
#define  mx 1.5 //100 MeV 3,2.5,1.5
#define lambda  0.1*3.14159265358979//0.1,5,0.2
 
 
static int imaxarg1, imaxarg2;
#define IMAX(a,b) (imaxarg1=(a),imaxarg2=(b),(imaxarg1)>(imaxarg2) ?\
   (imaxarg1):(imaxarg2))

static int iminarg1, iminarg2;
#define IMIN(a,b) (iminarg1=(a),iminarg2=(b),(iminarg1)<(iminarg2) ?\
   (iminarg1):(iminarg2))


int n_tab,                    /* number of points in EOS file */
    n_nearest=1;

double error,                 /* error in ratint (not used) */
       e_min,                 /* minimum energy density */ 
       e_max,                 /* maximum energy density */
       log_p_min,             /* minimum pressure */
       log_e_tab[501],        /* log of energy density points in EOS file */
       log_p_tab[501],        /* log of pressure points in EOS file */
       x[PDIV+1],             /* log of pressure points used in integration */
       energy_d[NDIV+1],      /* density output */
       log_p,                 /* pressure at density output */ 
       h[NDIV+1];             /* enthalpy output */ 

/***************************************************************************/
void hunt(double xx[], int n, double x, int *jlo)
{
	int jm,jhi,inc,ascnd;

	ascnd=(xx[n] > xx[1]);
	if (*jlo <= 0 || *jlo > n) {
		*jlo=0;
		jhi=n+1;
	} else {
		inc=1;
		if (x >= xx[*jlo] == ascnd) {
			if (*jlo == n) return;
			jhi=(*jlo)+1;
			while (x >= xx[jhi] == ascnd) {
				*jlo=jhi;
				inc += inc;
				jhi=(*jlo)+inc;
				if (jhi > n) {
					jhi=n+1;
					break;
				}
			}
		} else {
			if (*jlo == 1) {
				*jlo=0;
				return;
			}
			jhi=(*jlo);
			*jlo -= 1;
			while (x < xx[*jlo] == ascnd) {
				jhi=(*jlo);
				inc += inc;
				*jlo=jhi-inc;
				if (*jlo < 1) {
					*jlo=0;
					break;
				}
			}
		}
	}
	while (jhi-(*jlo) != 1) {
		jm=(jhi+(*jlo)) >> 1;
		if (x > xx[jm] == ascnd)
			*jlo=jm;
		else
			jhi=jm;
	}
}

/*************************************************************************/
/* Driver for the interpolation routine. First we find the tab. point    */
/* nearest to xb, then we interpolate using three points around xb.      */  
/*************************************************************************/
double interp(double xp[], double yp[], int np ,double xb)
{ 
 int k,        /* offset */
     m=4;      /* degree of interpolation */ 
 
 double y;

 hunt(xp,np,xb,&n_nearest);

 k=IMIN(IMAX(n_nearest-(m-1)/2,1),np+1-m);

 if( xb==xp[k] ||  xb==xp[k+1] || xb==xp[k+2] || xb==xp[k+3]) xb += 1e-16;

 y= (xb-xp[k+1])*(xb-xp[k+2])*(xb-xp[k+3])*yp[k]/
        ((xp[k]-xp[k+1])*(xp[k]-xp[k+2])*(xp[k]-xp[k+3]))
 
    +(xb-xp[k])*(xb-xp[k+2])*(xb-xp[k+3])*yp[k+1]/
       ((xp[k+1]-xp[k])*(xp[k+1]-xp[k+2])*(xp[k+1]-xp[k+3]))
 
    +(xb-xp[k])*(xb-xp[k+1])*(xb-xp[k+3])*yp[k+2]/
       ((xp[k+2]-xp[k])*(xp[k+2]-xp[k+1])*(xp[k+2]-xp[k+3]))
 
    +(xb-xp[k])*(xb-xp[k+1])*(xb-xp[k+2])*yp[k+3]/
       ((xp[k+3]-xp[k])*(xp[k+3]-xp[k+1])*(xp[k+3]-xp[k+2]));

 return (y);
}

/*************************************************************************/
/* Load EOS file.                                                        */
/* Format:  rho  p  #number                                              */ 
/*************************************************************************/

void load_eos( char *eos_file)
{
 int i; 
 double rho,                  /* mass density (input)*/
        pr;                   /* pressure */
      
 double junk1,junk2,junk3,junk4,junk5;

 FILE *fe;

 if((fe=fopen(eos_file,"r")) == NULL ) {
       printf("cannot open file\n");
       exit(0);
 }

 for (i=1;i<=n_tab;i++) {
   fscanf(fe,"%lf %lf %lf %lf \n",&rho,&pr,&junk1,&junk2) ;
       double dens=4.3*mx*mx*mx*mx/(lambda);//mx^4/3lambda
       log_e_tab[i]=log10(rho*C*C*KSCALE);
       log_p_tab[i]=log10(dens*MEVFM3/3.*pow(pow(1.+rho/(dens*MEVC2FM3),0.5)-1.,2.)*KSCALE);

       // printf("i=%d rho = %4.3e loge=%4.3e rho=%4.3e p = %4.3e logp=%4.3e\n", i, rho, log_e_tab[i], pow(10.0,log_e_tab[i])/(KSCALE*C*C),
       //     pr,log_p_tab[i]);

 }


}
/*************************************************************************/
/* Returns the function to be integrated at pressure p_val.              */
/*************************************************************************/
double f(double x_val)
{
 double log_e;                 /* temp. storage for log of energy density */
          
   log_e=interp(log_p_tab,log_e_tab,n_tab,x_val);
   return (1.0/(1.0+pow(10.0,log_e-x_val)));
}

/**************************************************************************/
/* Calculate the energy_d points where the enthalpy will be tabulated.    */ 
/**************************************************************************/
void make_e_grid(void)
{
 int i;

 double a;                    /* multiplication factor */

 a=pow(e_max/e_min, 1.0/(NDIV-1));

 for(i=1;i<=NDIV;i++) energy_d[i]=pow(a,1.0*(i-1))*e_min;
} 

/**************************************************************************/
/* Use Simpson's rule to integrate.                                       */   
/**************************************************************************/
double h_integ_tab(double x_val)
{ 
 register int i;
 
  double sum=0.0,          /* final sum */
         dx;

  /* Create x-grid for integration */
  
  for(i=1;i<=PDIV;i++) x[i]=log_p_min+(x_val-log_p_min)*(i-1.0)/PDIV_1; 

  dx=(x_val-log_p_min)/PDIV_1;    /* step in x-integration */

  for(i=1;i<=PDIV-2;i+=2) 
     sum += (dx/3.0)*(f(x[i])+4*f(x[i+1])+f(x[i+2]));

  return (sum);
}

/*******************************************************************/
/* Returns the derivative w.r.t. x of an array f[501].          */ 
/*******************************************************************/
double deriv_e(double logp[501], double loge[501], int i, int n_tab)
{
double d_temp;   
 switch(i) { 
            case 1    : d_temp=(logp[i+1]-logp[i])/(loge[i+1]-loge[i]);
                        break; 
      
            default     : d_temp=(logp[i]-logp[i-1])/(loge[i]-loge[i-1]);
                          break; 
 }
 
 return d_temp;
}

/**************************************************************************/
main(int argc, char *argv[])
{
 int i;                     /* counter */
 
 double dlogp_dloge[501],
        Gamma_ad[501],
        e_i, 
        p_i,
        n_0,
        mu_x_i;


 n_tab=atoi(argv[2]); 
 load_eos(argv[1]);       
 
 i=1;
 //printf("i=%d  loge=%4.3e e=%4.3e logp=%4.3e\n", i, log_e_tab[i], pow(10.0,log_e_tab[i])/(KSCALE*C*C),
 //	log_p_tab[i]);

 e_min=pow(10.0,log_e_tab[1]);
 e_max=pow(10.0,log_e_tab[n_tab]);

 log_p_min=interp(log_e_tab,log_p_tab,n_tab,log10(e_min));
  
 make_e_grid();

 //printf("i=%d  loge=%4.3e e=%4.3e logp=%4.3e\n", i, log_e_tab[i], pow(10.0,log_e_tab[i])/(KSCALE*C*C),
 //	log_p_tab[i]);

 for(i=1;i<=n_tab;i++) {

   if(i==1) h[i]=1.0/(C*C);   
       else {
       h[i]=log(10.0)*h_integ_tab(log_p_tab[i]);
     }
 }


 // for(i=1;i<=n_tab;i++)  dlogp_dloge[i]=deriv_e(log_p_tab,log_e_tab,i,n_tab);
 //
 // for(i=1;i<=n_tab;i++) {
 // Gamma_ad[i]=(pow(10.0,log_e_tab[i])+pow(10.0,log_p_tab[i]))
 //	         /pow(10.0,log_e_tab[i])*dlogp_dloge[i];
 //
 //}


printf("%d\n",n_tab);

 for(i=1;i<=n_tab;i++) {
 
   e_i=pow(10.0,log_e_tab[i])/(KSCALE*C*C);
   p_i=pow(10.0,log_p_tab[i])/KSCALE; 
   double dens=4.3*mx*mx*mx*mx/(lambda);
   double dummy = dens/3.*pow(pow(1.+e_i/(dens*MEVC2FM3),0.5)-1.,2.);

   mu_x_i=sqrt(sqrt(dummy*4./(3.*dens))+1.)*mx;
   n_0=mu_x_i*(mu_x_i*mu_x_i/(mx*mx)-1.)*3*4.3*mx*mx*mx/lambda*MEVC2FM3/(1.79 * pow(10,-27)*mx*100);
   //  printf("i=%d rho = %4.3e loge=%4.3e rho=%4.3e p = %4.3e logp=%4.3e\n", i, e_i, log_e_tab[i], pow(10.0,log_e_tab[i])/(KSCALE*C*C),
   //	      p_i,log_p_tab[i]);

   //n_0= (e_i+p_i/(C*C))*exp(-h[i])/MB;
   //printf("%6.5e %6.5e %16.15e %16.15e \n",e_i/MEVC2FM3, dummy,-log(mu_x_i/(mx*100))/*h[i]*/*C*C, n_0);
   printf("%6.5e %6.5e %16.15e %16.15e \n",e_i, p_i,log(mu_x_i/(mx))/*h[i]*/*C*C, n_0);

  }
}
