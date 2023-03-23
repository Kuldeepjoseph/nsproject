#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_poly.h>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_multiroots.h>
#include <time.h>


//============== CONSTANTS ==========
#define pi 3.14159265359
#define hc1 197.32
#define rmp 938.0   // avg mass of nucleon in nucleus
#define rmla 1116.0   // mass of lambda hyperon
#define rmsi 1193.0 // mass of sigma hyperon
#define rmca 1313.0   // mass of cascade hyperon 
#define rmsig  550.0   // mass of sigma meson
#define rmome  783.0   // mass of omega meson
#define rmrho  770.0   // mass of rho meson
#define MIN(a,b) (((a)<(b))?(a):(b))

//============================================================
//  Program to calculate coupling constants for RMF parameters 
//  includes calc of isovector couplings 
//  Debarati CHATTERJEE
//============================================================
/*
  Compile and link with:
gcc calccoup_isovec.c -o calccoupisovec -Wall -lgsl -lgslcblas -lm
./calccoupisovec calccoupisovec_test.out calccoup_isovec.out
*/

//======= DEFINITIONS ============= 
 float pi2  = pi*pi ;
 float pi3  = pi*pi*pi ;
 float hc2 = hc1*hc1 ;
 float hc3 = hc1*hc1*hc1 ;


// ============== STRUCTURES & FUNCTION DEF======================
// Structure for Input empirical data
struct Empirical {
   float rho0;
   float lasat0;
   float ksat0 ;
   float jsym0;
   float lsym0;
   float ksym0;
   float effm;
//   float del;
//   float rho;
};

// Structure for Parameters to det EoS
struct Couplings {
    double gsigm;
    double gomeg;
    double bsig;
    double csig;
//    double grho;
//    double grwn;
//    double gw2n;
}eos_couplings ;

// Structure for EOS Output
struct EOS {
    double sigmo;
    double endens;
    double pres;
    double comp;
//    double jsym;
//    double lsym;
//    struct Couplings eos_couplings;
};

// ==================================================================
// function to calculate Couplings using satdata (of struct Empirical)
// and returns Output (of struct Couplings)
   struct Couplings Calc_Couplings( struct Empirical satdata );

// ==================================================================
// function to calculate EoS using satdata (of struct Empirical)
// and returns EOS Output (of struct EOS)
   struct EOS Calc_EOS( struct Empirical satdata, struct Couplings eos_couplings );
  
//=====================================================================   
// FUNCTIONS FOR MULTIROOT SOLVER

// Structure for parameters of multiroot solver
struct rparams
  {
   double p1;
   double p2;
   double p3;
   double p4;
   double p5;
   double p6;
   double p7;
  };
   
// Define Functions to solve
  int any_f(const gsl_vector * x, void *params,
              gsl_vector * f)
{
  double p1 = ((struct rparams *) params)->p1;
  double p2= ((struct rparams *) params)->p2;
  double p3 = ((struct rparams *) params)->p3;
  double p4 = ((struct rparams *) params)->p4;
  double p5 = ((struct rparams *) params)->p5;
  double p6 = ((struct rparams *) params)->p6;
  double p7 = ((struct rparams *) params)->p7;
//   printf("params p1,p2,p3,p7=%f %f %f %f\n",p1,p2,p3,p7) ;

// x1,x2,x3,x4 = couplings gsig, gomeg,bsig,csig
  const double x1 = gsl_vector_get (x, 0);
  const double x2 = gsl_vector_get (x, 1);
  const double x3 = gsl_vector_get (x, 2);
  const double x4 = gsl_vector_get (x, 3);
 // printf("x1,x2,x3,x4 =%f %f %f %f \n",x1,x2,x3,x4);

// Call Function to calc EoS parameters for SNM
// given satdata and couplings
 struct Couplings eos_couplings;
     eos_couplings.gsigm = x1 ;
     eos_couplings.gomeg = x2 ;
     eos_couplings.bsig = x3 ;
     eos_couplings.csig = x4 ;
//  printf("eos_couplings gsigm,gomeg,bsig,csig=%f %f %f %f\n",eos_couplings.gsigm,eos_couplings.gomeg,eos_couplings.bsig,eos_couplings.csig); 
 struct EOS eos_snm ;
 struct Empirical satdata;
   satdata.rho0 = p1 ;
   satdata.lasat0 = p2 ;
   satdata.ksat0 = p3 ;
   satdata.jsym0 = p4 ;
   satdata.lsym0 = p5;
   satdata.ksym0 = p6;
   satdata.effm = p7 ;
//   printf("satdata rho0,lasat0,ksat0,effm0=%f %f %f %f\n",satdata.rho0,satdata.lasat0,satdata.ksat0,satdata.effm) ;
     eos_snm = Calc_EOS( satdata, eos_couplings );
//     printf("eos sigmo,endens,pres,comp=%f %f %f %f\n",eos_snm.sigmo,eos_snm.endens,eos_snm.pres,eos_snm.comp);

// FUNCTIONS DEFINED HERE!!!
// printf("rmp,p2,p3,p7 =%f %f %f %f\n",rmp,p2,p3,p7);
  const double y1 = p7*rmp - rmp + x1*eos_snm.sigmo ; // effm
  const double y2 = eos_snm.endens - p2 ; // binding energy
  const double y3 = eos_snm.pres ;  // Pres 
  const double y4 = eos_snm.comp - p3;  // comp 
//  printf("y1,y2,y3,y4=%f %f %f %f\n",y1,y2,y3,y4);

  gsl_vector_set (f, 0, y1);
  gsl_vector_set (f, 1, y2);
  gsl_vector_set (f, 2, y3);
  gsl_vector_set (f, 3, y4);

  return GSL_SUCCESS;
}

//double min (double a, double b);

//==================== MAIN ===========================


int main(int argn, char** argv)
{
   FILE *fin,*fout_test,*fout_plot ;
   struct Empirical satdata0 ;
//   struct Couplings eos_couplings ;
   float rho0in,lasat0in,ksat0in,jsym0in,lsym0in,ksym0in,effmin ;
   struct EOS eos_snm ;
   double gsigsol,gomegsol,bsigsol,csigsol ;

//    printf("pi,pi2,rmp,hc1=%f %f %f %f\n",pi,pi2,rmp,hc1); 

// read input parameters from input file
 fin=fopen("calccoup_isovecGM1.in","r");
 fscanf(fin,"%f %f %f",&rho0in,&lasat0in,&ksat0in);
 fscanf(fin,"%f %f %f",&jsym0in,&lsym0in,&ksym0in);
 fscanf(fin,"%f",&effmin);

// open files to write
   fout_test=fopen(argv[1],"w");
//   fout_plot=fopen(argv[2],"a+");
   fout_plot=fopen(argv[2],"w");

// Calc of coupling constants

//  empirical data at saturation 
   satdata0.rho0 = rho0in ;
   satdata0.lasat0 = lasat0in ;
   satdata0.ksat0 = ksat0in ;
   satdata0.jsym0 = jsym0in ;
   satdata0.lsym0 = lsym0in;
   satdata0.ksym0 = ksym0in;
   satdata0.effm = effmin ;

 //  printf("satdata rho0,lasat0,ksat0,effm0=%f %f %f %f\n",satdata0.rho0,satdata0.lasat0,satdata0.ksat0,satdata0.effm) ;
   
// Call function to calculate EoS Couplings for given satdata
   eos_couplings = Calc_Couplings( satdata0 );
//      struct Couplings eos_couplings = {8.,10.,0.005,-0.002};
//      struct Couplings eos_couplings = {8.782,8.712,0.00865,-0.00242}; // TEST
      eos_snm = Calc_EOS( satdata0, eos_couplings) ; 
       printf(" sigmo=%f, endens=%f, pres=%f, comp=%f\n",eos_snm.sigmo,eos_snm.endens,eos_snm.pres,eos_snm.comp);

     gsigsol = eos_couplings.gsigm ;
     gomegsol = eos_couplings.gomeg ;
     bsigsol = eos_couplings.bsig ;
     csigsol = eos_couplings.csig ;
     printf("gsigsol,gomegsol,bsigsol,csigsol=%f %f %f %f\n",gsigsol,gomegsol,bsigsol,csigsol);

// Conversions 
         double gsms2 = pow(gsigsol,2)/rmsig/rmsig ; // in MeV^-2
         double gwmw2 = pow(gomegsol,2)/rmome/rmome ;
       //  b1 = bsig*rmp*sqrt(gsms2)*gsms2*rmsig**3
       //  c1 = csig*gsms2*gsms2*rmsig**4
         double b1 = bsigsol*rmp*pow(gsigsol,3) ;
         double c1 = csigsol*pow(gsigsol,4) ;
	 printf("gsms2,gwmw2,b1,c1=%f %f %f %f\n",gsms2,gwmw2,b1,c1);

         double gsms2f = gsms2*hc2 ;  // in fm^2
         double gwmw2f = gwmw2*hc2 ;  // in fm^2
         double bgle = b1/(rmp*sqrt(gsms2)*gsms2*pow(rmsig,3)) ;
         double cgle = c1/(gsms2*gsms2*pow(rmsig,4)) ;
	 printf("gsms2f,gwmw2f,bgle,cgle=%f %f %f %f\n",gsms2f,gwmw2f,bgle,cgle);


// CALC OF ISOVECTOR PARAMETERS     
// ASYMMETRY PARAMETER ==============
          double rhob=satdata0.rho0*hc3 ; // in MeV^3
	  double rhop = rhob/2. ;
          double rmn = satdata0.effm*rmp ;
          double pfp = pow(1.5*pi2*rhob,1./3.) ;
          double pfp2 = pfp*pfp ;
          double rmn2 = rmn*rmn ;
	  double rmup = sqrt(pfp2 + rmn2) ; // in MeV
	  double rmup2 = rmup*rmup ;
//	  printf("rhop/hc3,rmup=%g %g \n",rhop/hc3,rmup);

// Calc of Lambda_omega parameter =====				     
         double j0 = pow(pfp,2.)/6./rmup ;      // in MeV
	 double j1 = satdata0.jsym0 - j0 ;  // in MeV
//	 printf("j0,j1=%f %f\n",j0,j1) ;

         double mstgs2 = 1./gsms2 + 2.*bsigsol*rmp*gsigsol*eos_snm.sigmo
	                          + 3.*csigsol*pow(gsigsol*eos_snm.sigmo,2.) ;	 
	 double dnsdrmn = ( pfp/rmup*(rmup2 + 2.*rmn2)
		  - 3.*rmn2*log((pfp + rmup)/rmn) )/pi2 ; 
	 double drmndrho = - rmn/rmup*1./(mstgs2 + dnsdrmn) ;
	 double l0 = j0*( 1.+rmn2/rmup2*(1.-3.*rhob/rmn*drmndrho ) );
	 double l1 = satdata0.lsym0 - l0 ; 

// Calc of isovec couplings 	 
// older expr for grmr2
     //    double grmr2 = hc2*( satdata0.jsym0 - pow(pfp,2)/6./sqrt(pfp2 + rmn2) )
     //    *12.*pi2/pow(pfp,3) ;   // in fm^2
     //    double grmr2 = 1./( rhob/8./j1 ) ; // in MeV^-2
         double omeg0=sqrt(gwmw2)*2.*rhop/rmome ;

///         double lamomeg = 0.; // TEST
         double lamomeg = (1.-l1/3./j1)/(32.*gomegsol*omeg0*j1)/gwmw2 ;  
	 printf("l0,l1,lamomeg=%f %f %f\n",l0,l1,lamomeg) ;


	 double grmr2 = 1./( rhob/8./j1 - 2.*lamomeg*pow(gomegsol*omeg0,2.) ) ; // in fm^2
	 double grmr2f = hc2*grmr2 ;        // in fm^2
         double grho = sqrt(grmr2f)*rmrho/hc1 ;  // dimless
	 printf("grmr2f,grho=%f %f\n",grmr2f,grho) ;


//==================== PRINT OUTPUT ====================================
//   fout_test=fopen(argv[1],"w");
//   fout_plot=fopen(argv[2],"a+");

// TEST FILE

    fprintf(fout_test,"rho0=%f,lasat0= %f,ksat0= %f\n",rho0in,lasat0in,ksat0in);
    fprintf(fout_test,"jsym0=%f,lsym0= %f,ksym0= %f\n",jsym0in,lsym0in,ksym0in);
    fprintf(fout_test,"effmin=%f\n",effmin);

        fprintf(fout_test,"gsms2=%f, gwmw2=%f, grmr2=%f (all in fm^2) \n",gsms2f,gwmw2f,grmr2f);
	fprintf(fout_test,"\n");
        fprintf(fout_test,"g_sigma=%f, g_omega=%f, g_rho=%f (all dimesionless)\n",gsigsol,gomegsol,grho);
	fprintf(fout_test,"\n");
        fprintf(fout_test," b=%f, c=%f (both dimensionless)\n",bgle,cgle); //see Glend, PRL 1991)
     //   fprintf(fout_plot,"grmr2=%f, g_rho=%f \n",grmr2,grho);

// INPUT FOR EOS
   //     fprintf(fout_plot,"%f %f %f %f %f  %f %f %f  %f %f %f\n",rho0in,lasat0in,ksat0in,jsym0in,effmin,gsigsol,gomegsol,grho,bgle,cgle,lamomeg);
fprintf(fout_plot,"%f %f %f %f %f %f  %f %f %f  %f %f %f\n",rho0in,lasat0in,ksat0in,jsym0in,lsym0in,effmin,gsigsol,gomegsol,grho,bgle,cgle,lamomeg);


  fclose(fout_plot);
  fclose(fout_test);

return 0;
}


// ========================= FUNCTIONS ===========================
// function to calculate EoS couplings
 struct Couplings Calc_Couplings( struct Empirical satdata ){
 struct Couplings eos_couplings ;
 struct EOS eos_snm;
// float rho0,lasat0,ksat0,effm,jsym0,lsym0,ksym0;
// float rho,del ;

//====== CALCULATION OF COUPLINGS ========================      

// MULTIROOT SOLVE// MULTIROOT SOLVER    
  const gsl_multiroot_fsolver_type *T;
  gsl_multiroot_fsolver *s;

  int status;
  size_t i, iter = 0;

  const size_t n = 4;

  double x_init[4] = {8.,10.,0.005,-0.002};
///    double x_init[4] = {8.782,8.712,0.00865,-0.00242}; // TEST
  gsl_vector *x = gsl_vector_alloc (n);

  gsl_vector_set (x, 0, x_init[0]);
  gsl_vector_set (x, 1, x_init[1]);
  gsl_vector_set (x, 2, x_init[2]);
  gsl_vector_set (x, 3, x_init[3]);

//   printf("params rho0,lasat0,ksat0,effm0=%f %f %f %f\n",satdata.rho0,satdata.lasat0,satdata.ksat0,satdata.effm) ;
  struct rparams p = {satdata.rho0,satdata.lasat0,satdata.ksat0,satdata.jsym0,satdata.lsym0,satdata.ksym0,satdata.effm};
//  struct rparams p = {satdata};
  gsl_multiroot_function f = {&any_f, n, &p};

  T = gsl_multiroot_fsolver_hybrids;
//  T = gsl_multiroot_fsolver_dnewton;
//  T = gsl_multiroot_fsolver_broyden;
  s = gsl_multiroot_fsolver_alloc (T, 4);
  gsl_multiroot_fsolver_set (s, &f, x);

  //print_state (iter, s);

  do
    {
      iter++;
      status = gsl_multiroot_fsolver_iterate (s);

 //     print_state (iter, s);

      if (status)   /* check if solver is stuck */
        break;

      status =
        gsl_multiroot_test_residual (s->f, 1e-7);
    }
  while (status == GSL_CONTINUE && iter < 1000);

  printf ("status = %s\n", gsl_strerror (status));

          eos_couplings.gsigm = gsl_vector_get (s->x, 0),
          eos_couplings.gomeg = gsl_vector_get (s->x, 1),
          eos_couplings.bsig = gsl_vector_get (s->x, 2),
          eos_couplings.csig = gsl_vector_get (s->x, 3),

  gsl_multiroot_fsolver_free (s);
  gsl_vector_free (x);

   return eos_couplings ;
}

//==================================================================
// function to calculate EoS parameters
 struct EOS Calc_EOS( struct Empirical satdata, struct Couplings eos_couplings ){
 struct EOS eos_any;
// float rho0,lasat0,ksat0,effm,jsym0,lsym0,ksym0;
 float rhob,del ;

// double gsig,gomeg,bsig,csig,grho ;
	 //,grwn,gw2n ;
// double endens,dedx,pres,comp,jsym,lsym,ksym ;
 double rmn,pfp;

//  printf("%f %f %f %f\n",pi,pi2,rmn,hc1);
//  printf("rmsig,rmome=%f %f\n",rmsig,rmome);

 
 /*
   printf("%f %f %f \n",satdata.rho0,satdata.lasat0,satdata.ksat0);
   printf("%f %f %f \n",satdata.jsym0,satdata.lsym0,satdata.ksym0);
   printf("%f \n",satdata.effm);
   printf("%f %f \n",eos_couplings.gsigm,eos_couplings.gomeg);
   printf("%f %f \n",eos_couplings.bsig,eos_couplings.csig);
 */


// calculation of EoS parameters of SNM
          double gsms2 = pow(eos_couplings.gsigm,2.)/rmsig/rmsig ; 
	  double gwmw2 = pow(eos_couplings.gomeg,2.)/rmome/rmome ;
	  double b1 = eos_couplings.bsig*rmp*pow(eos_couplings.gsigm,3);
          double c1 = eos_couplings.csig*pow(eos_couplings.gsigm,4);
          double bcomp = eos_couplings.bsig*rmp*eos_couplings.gsigm ;
          double ccomp = eos_couplings.csig*pow(eos_couplings.gsigm,2) ;
        //  printf("gsms2,gwmw2,b1,c1=%f %f %f %f \n",gsms2,gwmw2,b1,c1);
        //  printf("bcomp,ccomp=%f %f \n",bcomp,ccomp);
	//  printf("\n");


          rhob=satdata.rho0*hc3 ;
          rmn = satdata.effm*rmp ;
          pfp = pow(1.5*pi2*rhob,1./3.) ;
          double pfp2 = pfp*pfp ;
          double rmn2 = rmn*rmn ;
        double x1p=sqrt(pfp2+rmn2) ;
        double y1p=pfp*x1p ;
        double y2p=rmn2*log( (pfp+x1p)/rmn ) ;
//	printf("rho0,pfp,rmn=%f %f %f\n",satdata.rho0,pfp,rmn);
//	printf("rhob,x1p,y1p,y2p=%f %f %f %f\n ",rhob,x1p,y1p,y2p);
//	printf("\n");

//  SCALAR DENSITY =========
        double rhop=pow(pfp,3.)/3./pi2 ; // MeV^3
	double omeg0=sqrt(gwmw2)*2.*rhop/rmome ;
        double rnss=rmn*(y1p-y2p)/pi2 ;
        double sigm = sqrt(gsms2)*rnss*rmsig ;
//	printf("rhop,omeg0,rnss,sigm=%f %f %f %f\n",rhop,omeg0,rnss,sigm);
//	printf("\n");

// SOLVE CUBIC EQUATION FOR SIGMA	
         double x0,x1,x2; 
	 int roots;
// aa x^3 + bb * x^2 + cc * x + dd = 0
// Solve Cubic Equation 

	 double aa = c1 ; 
	 double bb = b1/3.;
	 double cc = rmsig*rmsig/3.;
	 double dd = -sigm;

//	printf("aa,bb,cc,dd=%f %f %f %f\n",aa,bb,cc,dd);

    //    roots = gsl_poly_solve_cubic(bb/aa, cc/aa, dd/aa, &x0, &x1, &x2);
    //    printf("roots %d, x0,x1,x2 =%g %g %g\n",roots,x0,x1,x2);
	 double hh=aa*cc-pow(bb,2);
         double gg=dd*pow(aa,2) - 3.*aa*bb*cc + 2.*pow(bb,3);
         double gh= pow(gg,2) + 4.*pow(hh,3);
//	 printf("hh,gg,gh=%g %g %g\n",hh,gg,gh);
         if(gh < 0.0e0) {
          double phi = acos(-gg/2./(pow(-hh,1.5))) ;
          double h1 = sqrt(-hh) ;
          double y1= 2.*h1*cos(phi/3.) ;
          double y2 = -2.*h1*cos((phi+pi)/3.) ;
          double y3 = -2.*h1*cos((phi-pi)/3.) ;
//	  printf("phi,h1,y1,y2,y3=%g %g %g %g %g\n",phi,h1,y1,y2,y3);
          double x0 = (y1-bb)/aa ;
          double x1 = (y2-bb)/aa ;
          double x2 = (y3-bb)/aa ;
//	  printf("x0,x1,x2=%f %f %f\n",x0,x1,x2);
         if ( x0 < 0.e0 ) { x0=1.e10; } 
         if ( x1 < 0.e0 ) { x1=1.e10; }
         if ( x2 < 0.e0 ) { x2=1.e10; } 
//	  printf("x0,x1,x2=%f %f %f\n",x0,x1,x2);
         double xx = MIN(x0,x1) ;
         eos_any.sigmo = MIN(xx,x2) ;
//	printf("sol 1 aa,bb,cc,dd,sigmo=%f %f %f %f %f\n",aa,bb,cc,dd,eos_any.sigmo);
	  }
	 else {
          double p1 = (-gg + sqrt(gh))/2.  ;
          double pp = pow(p1,0.333333) ;
          double xx = pp - hh/pp ;
          eos_any.sigmo = (xx-bb)/aa ;
//	printf("sol 2 bb/aa,cc/aa,dd/aa,sigmo=%f %f %f %f\n",bb/aa, cc/aa, dd/aa,eos_any.sigmo);
	 }

   
// BARYON DENSITY =========
        double rhobsym = 2.*rhop ; // MeV^3
        double tpe = pow(rmome*omeg0,2)/2. + pow(rmsig*eos_any.sigmo,2)/2.
            + b1*pow(eos_any.sigmo,3)/3. + c1*pow(eos_any.sigmo,4)/4.  ; // MeV^4
        double tkep = ( pfp*pow(x1p,3)/4. -rmn2*(y1p + y2p)/8. )/pi2 ;
        double tke = 2.*tkep ; // MeV^4
        double tend0 = tpe + tke ; // MeV^4
        double tend = tend0/rhobsym - rmp ; // MeV
        double tend1 = tend0/hc3 ; // MeV/fm^3 
//	eos_any.endens = tend0/rhobsym ; // E/A in MeV
	eos_any.endens = tend ; // E/A in MeV
//	printf("tpe,tkep,tend0,tend,tend1=%f %f %f %f %f\n",tpe,tkep,tend0,tend,tend1);
//	printf("\n");

//==  PRESSURE CALCULATION ========= 
         double prpe = pow(rmome*omeg0,2)/2. - pow(rmsig*eos_any.sigmo,2)/2.
             - b1*pow(eos_any.sigmo,3)/3. - c1*pow(eos_any.sigmo,4)/4.  ;
         double prkep = ( pfp*pow(x1p,3)/4. -rmn2*(5.*y1p - 3.*y2p)/8. )/3./pi2 ;
         double prke = 2.*prkep ;
         double tpres = prpe + prke ;
         eos_any.pres = tpres/hc3 ; // P in MeV/fm^3
//	printf("prpe,prkep,tpres=%f %f %f \n",prpe,prkep,tpres);
//	printf("\n");

// COMPRESS CALCULATION  ================
///         double  comp = 9.*(tpe + prke)/rhob ;
           double compke = 2.*(y1p/2. + rmn2*pfp/x1p - 1.5*y2p)/pi2 ;
          double comg = gwmw2*6.0*pow(pfp,3)/pi2 + 3.*pfp2/x1p
             - 6.*pow(pfp,3)*rmn2*gsms2/pi2/pow(x1p,2)/( 1. + gsms2*
               (2.*bcomp*eos_any.sigmo + 3.*ccomp*eos_any.sigmo*eos_any.sigmo +compke) );
          eos_any.comp = comg ;
//	  printf("compke,comg=%f %f\n",compke,comg);
//	  printf("\n");

//       printf("Calc_EOS sigmo=%f, endens=%f, pres=%f, comp=%f\n",eos_any.sigmo,eos_any.endens,eos_any.pres,eos_any.comp);
//	printf("\n");
//   printf("%f %f %f \n",satdata.rho0,satdata.rhob,eos_any.endens);

   return eos_any ;
}
// ==================================================================

