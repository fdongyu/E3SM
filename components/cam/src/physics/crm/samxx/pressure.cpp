
#include "pressure.h"

using yakl::SArray;

void pressure() {
  auto &p             = :: p;
  auto &rhow          = :: rhow;
  auto &adz           = :: adz;
  auto &adzw          = :: adzw;
  auto &dz            = :: dz;
  auto &dx            = :: dx;
  auto &dy            = :: dy;
  auto &rho           = :: rho;
  auto &ncrms         = :: ncrms;

  int npressureslabs = nsubdomains;
  int nzslab = max(1,nzm/npressureslabs); 
  int nx2 = nx+2;
  int ny2 = ny+2*YES3D;
  int constexpr n3i=3*nx_gl/2+1;
  int constexpr n3j=3*ny_gl/2+1;
  int constexpr fftySize = ny > 4 ? ny : 4;

  real4d f ("f" , nzslab, ny2, nx2, ncrms);
  real4d ff("ff", nzm,ny2,nx+1,ncrms);
  real2d a ("a" , nzm, ncrms);
  real2d c ("c" , nzm, ncrms);

  int iwall = 0;
  int nypp, jwall;

  if (RUN2D) {
    nypp = 1;
    jwall = 0;
  } else {
    nypp = ny+2;
  }

  real2d eign("eign",nypp,nx+1);

  press_rhs();

  // for (int k=0; k<nzslab; k++) {
	//   for (int j=0; j<ny; j++) {
	//     for (int i=0; i<nx; i++) {
	//       for (int icrm=0; icrm<ncrms; icrm++) {
  parallel_for( Bounds<4>(nzslab,ny,nx,ncrms) , YAKL_LAMBDA (int k, int j, int i, int icrm) {
    f(k,j,i,icrm) = p(k,j+offy_p,i+offx_p,icrm);
	});

  yakl::FFT<nx> fftx;
  yakl::FFT<fftySize> ffty;

  // for (int k=0; k<nzslab; k++) {
  //  for (int j=0; j<ny; j++) {
	//      for (int icrm=0; icrm<ncrms; icrm++) {
  parallel_for( Bounds<3>(nzslab,ny,ncrms) , YAKL_LAMBDA (int k, int j, int icrm) {
    real ftmp[nx+2];
    real tmp [nx];

    for (int i=0; i<nx ; i++) { ftmp[i] = f(k,j,i,icrm); }

    fftx.forward(ftmp, tmp);

    for (int i=0; i<nx2; i++) { f(k,j,i,icrm) = ftmp[i]; }
  });

  if (RUN3D) {
    // for (int k=0; k<nzslab; k++) {
    //  for (int i=0; j<nx+1; i++) {
    //    for(int l=0; l<ny2; l++) {
    //      for (int icrm=0; icrm<ncrms; icrm++) {
    parallel_for( Bounds<3>(nzslab,nx+1,ncrms) , YAKL_LAMBDA (int k, int i, int icrm) {
      real ftmp[ny+2];
      real tmp [ny];

      for (int j=0; j<ny ; j++) { ftmp[j] = f(k,j,i,icrm); }

      ffty.forward(ftmp, tmp);

      for (int j=0; j<ny2; j++) { f(k,j,i,icrm) = ftmp[j]; }
	  });
  }

  // for(int k=0; k<nzslab; k++) {
  //  for(int j=0; j<nypp; j++) {
  //    for(int i=0; i<nx+1; i++) {
  //      for(int icrm=0; icrm<ncrms; icrm++) {
  parallel_for( Bounds<4>(nzslab,nypp,nx+1,ncrms) , YAKL_LAMBDA (int k, int j, int i, int icrm) {
    ff(k,j,i,icrm) = f(k,j,i,icrm);
  });

  // for (int k=0; k<nzm; k++) {
  //  for (int icrm=0; icrm<ncrms; icrm++) {
  parallel_for( Bounds<2>(nzm,ncrms) , YAKL_LAMBDA (int k, int icrm) {
    a(k,icrm)=rhow(k,icrm)/(adz(k,icrm)*adzw(k,icrm)*dz(icrm)*dz(icrm));
    c(k,icrm)=rhow(k+1,icrm)/(adz(k,icrm)*adzw(k+1,icrm)*dz(icrm)*dz(icrm));
	});

	//   for (int j=0; j<nypp; j++) {
	//     for (int i=0; i<nx+1; i++) {
  parallel_for( Bounds<2>(nypp,nx+1) , YAKL_LAMBDA (int j, int i) {
    int jt = 0;
    int it = 0;

    real ddx2=1.0/(dx*dx);
    real ddy2=1.0/(dy*dy);
    real pii = 3.14159265358979323846;
    real xnx=pii/nx;
    real xny=pii/ny;
    int jd=((j+1)+jt-0.1)/2.0;
    real facty = 2.0;
    real xj=jd;
    int id=((i+1)+it-0.1)/2.0;
    real factx = 2.0;
    real xi=id;
    eign(j,i)=(2.0*cos(factx*xnx*xi)-2.0)*ddx2+(2.0*cos(facty*xny*xj)-2.0)*ddy2;
	});

	// for (int j=0; j<nypp; j++) {
	//  for (int i=0; i<nx+1; i++) {
	//    for (int icrm=0; icrm<ncrms; icrm++) {
  parallel_for( Bounds<3>(nypp,nx+1,ncrms) , YAKL_LAMBDA (int j, int i, int icrm) {
    SArray<real,nzm-1> alfa;
    SArray<real,nzm-1> beta;

    int jt = 0;
    int it = 0;
    int jd=((j+1)+jt-0.1)/2.0;
    int id=((i+1)+it-0.1)/2.0;
    real b;
    if(id+jd == 0) {
      b=1.0/(eign(j,i)*rho(0,icrm)-a(0,icrm)-c(0,icrm));
      alfa(0)=-c(0,icrm)*b;
      beta(0)=ff(0,j,i,icrm)*b;
    }
    else {
      b=1.0/(eign(j,i)*rho(0,icrm)-c(0,icrm));
      alfa(0)=-c(0,icrm)*b;
      beta(0)=ff(0,j,i,icrm)*b;
    }

    real e;
    for(int k=1; k<nzm-1; k++) {
      e=1.0/(eign(j,i)*rho(k,icrm)-a(k,icrm)-c(k,icrm)+a(k,icrm)*alfa(k-1));
      alfa(k)=-c(k,icrm)*e;
      beta(k)=(ff(k,j,i,icrm)-a(k,icrm)*beta(k-1))*e;
    }
    ff(nzm-1,j,i,icrm)=(ff(nzm-1,j,i,icrm)-a(nzm-1,icrm)*beta(nzm-2))/(eign(j,i)*rho(nzm-1,icrm)-a(nzm-1,icrm)+a(nzm-1,icrm)*alfa(nzm-2));
    for(int k=nzm-2; k>=0; k--) {
      ff(k,j,i,icrm)=alfa(k)*ff(k+1,j,i,icrm)+beta(k);
    }
  });

  // for (int k=0; k<nzslab; k++) {
	//   for (int j=0; j<nypp; j++) {
	//     for (int i=0; i<nx+1; i++) {
	//       for (int icrm=0; icrm<ncrms; icrm++) {
  parallel_for( Bounds<4>(nzslab,nypp,nx+1,ncrms) , YAKL_LAMBDA (int k, int j, int i, int icrm) {
    f(k,j,i,icrm) = ff(k,j,i,icrm);
	});

  if (RUN3D) {
    // for (int k=0; k<nzslab; k++) {
		//  for (int i=0; i<nx+1; i++) {
    //    for(int l=0; l<ny+2; l++) { 
		//      for (int icrm=0; icrm<ncrms; icrm++) {
    parallel_for( Bounds<3>(nzslab,nx+1,ncrms) , YAKL_LAMBDA (int k, int i, int icrm) {
      real ftmp[ny+2];
      real tmp [ny];
      
      for(int j=0; j<ny+2; j++) { ftmp[j] = f(k,j,i,icrm); }

      ffty.inverse(ftmp,tmp);

      for(int j=0; j<ny  ; j++) { f(k,j,i,icrm) = ftmp[j]; } 
	  });
  }

  // for (int k=0; k<nzslab; k++) {
	//  for (int j=0; i<ny; i++) {
  //    for(in l=0; l<nx; l++) { 
	//      for (int icrm=0; icrm<ncrms; icrm++) {
  parallel_for( Bounds<3>(nzslab,ny,ncrms) , YAKL_LAMBDA (int k, int j, int icrm) {
    real ftmp[nx+2];
    real tmp [nx];

    for(int i=0; i<nx+2; i++) { ftmp[i] = f(k,j,i,icrm); }

    fftx.inverse(ftmp,tmp);

    for(int i=0; i<nx  ; i++) { f(k,j,i,icrm) = ftmp[i]; }
  });

  parallel_for( Bounds<4>(nzslab,dimy_p,nx+1,ncrms) , YAKL_LAMBDA (int k, int j, int i, int icrm) {
    int jj, ii;

    if (YES3D) {
      if (j == 0) {
        jj = ny-1; 
      } else {
        jj = j-1;
      }
    } else {
      jj = j;
    }

    if (i == 0) {
      ii = nx-1;
    } else {
      ii = i-1;
    }

    p(k,j,i,icrm) = f(k,jj,ii,icrm);
  });

  press_grad();
}


