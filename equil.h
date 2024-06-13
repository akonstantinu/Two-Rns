double f(double ptot,double etot);

double h_integ_tab(double p1[RDIV+1],double e1[RDIV+1], double p2[RDIV+1],double e2[RDIV+1],int Num);

double max(double v1,double v2);

double min(double v1,double v2);

double dN_dr_is(double r_is, double r,double m, double p, double pDM, double p_surface, double p_surfaceDM,double e_center);

void make_grid( double s_gp[SDIV+1], 
                double mu[MDIV+1]);                        

void load_eos( char eos_file[], 
               double log_e_tab[2001], 
               double log_p_tab[2001], 
               double log_h_tab[2001],
               double log_n0_tab[2001], 
               int *n_tab);

double e_of_rho0(double rho0, double Gamma_P);

double e_at_p(double pp, 
              double pp_surface,
              double log_e_tab[2001], 
              double log_p_tab[2001],
              int    n_tab, 
              int    *n_nearest_pt,
              char eos_type[],
              double Gamma_P);

double p_at_e(double ee, 
              double ee_surface,
              double log_p_tab[2001], 
              double log_e_tab[2001],
              int    n_tab, 
              int    *n_nearest_pt);

double p_at_h(double hh, 
              double hh_surface,
              double log_p_tab[2001], 
              double log_h_tab[2001],
              int    n_tab, 
              int    *n_nearest_pt);

double h_at_p(double pp, 
              double pp_surface,
              double log_h_tab[2001], 
              double log_p_tab[2001],
              int    n_tab, 
              int    *n_nearest_pt);

double n0_at_e(double ee, 
               double ee_surface,
               double log_n0_tab[2001], 
               double log_e_tab[2001],
               int    n_tab, 
               int    *n_nearest_pt);


double dpds(double p,
            double p_surface,
	    double e,
	    double rho_s,
	    double gamma_s,
            double s,
            double v,
            double Omega,
            double omega,
            double omega_s);

double dpdm(double p,
            double p_surface,
	    double e,
	    double rho_m,
	    double gamma_m,
            double m,
            double v,
            double Omega,
            double omega,
            double omega_m);

double y_at_x(double x,
	      double x1,
	      double x2,
	      double y1,
	      double y2);





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
	       double p_surface);




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
     double *Mp);

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
                char   eos_type[],
                double Gamma_P);

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
                double Gamma_P);

double dr_dr_is(double r_is, double r, double m);

void TOV(int    i_check, 
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
               double *mDM_final);

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
	    double *r_eDM);


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
	  double *r_eDM_new,
	  double *Omega,
	  double *OmegaDM) ;
 
