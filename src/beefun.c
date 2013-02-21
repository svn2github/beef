#define _XOPEN_SOURCE 500

#include <math.h>
#include <stdlib.h>
#include <string.h>

#include "beefleg.h"

#include "pbecor.h"

// evaluate bee exchange energy and its derivatives de/drho and ( de/d|grad rho| ) / |grad rho|
void beefx_(double *r, double *g, double *e, double *dr, double *dg, int *addlda)
{
    double s2,t,r43,r83,s,sx,dx,fx,dl,dfx;
    const int n=nmax;
    const int i1=1;
    const int i2=1;

    r43 = pow(*r, 4./3.);
    r83 = r43*r43;
    sx = r2e * r43;
    dx = 4./3. * sx / (*r);

    s2 = *g*pix / r83;
    s = sqrt(s2);
    t = 2.*s2/(4.+s2)-1.;

    if(beeforder==-1)
    {
	calclegdleg(t);
	
	if(!(*addlda))
	    fx = ddot_(&n, mi, &i1, L, &i2) - 1.;
	else
	    fx = ddot_(&n, mi, &i1, L, &i2);
	dl = ddot_(&n, mi, &i1, dL, &i2);
	
	dfx = dl*( 4.*s / (4.+s2) - 4.*s2*s/pow((4.+s2), 2) );
	*dr = dx*fx - 4./3.*s2/(s*(*r))*sx*dfx;
	*dg = sx*dfx*pix/(s*r83);
	*e = sx*fx;
	return;
    }
    
    if(beeforder>=0)
    {
	(*LdLn[beeforder])(t, &fx, &dl);

	dfx = dl*( 4.*s / (4.+s2) - 4.*s2*s/pow((4.+s2), 2) );
	*dr = dx*fx - 4./3.*s2/(s*(*r))*sx*dfx;
	*dg = sx*dfx*pix/(s*r83);
	*e = sx*fx;
    }
    else
    {
	*dr = 0.;
	*dg = 0.;
	*e = 0.;
    }
}


// evaluate local part of bee correlation and its derivatives de/drho and ( de/d|grad rho| ) / |grad rho|
void beeflocalcorr_(double *r, double *g, double *e, double *dr, double *dg, int *addlda)
{
    double rs, ldac, ldadr, pbec, pbedr, pbed2rho;

    if(beeforder>=0)
    {
	*e = 0.;
	*dr = 0.;
	*dg = 0.;
	return;
    }
    
    rs = invpi075tothird / pow(*r,1./3.);
    corpbe(rs, 0.5/r2k * sqrt(*g*rs) / (*r),
	(beeforder!=-2), 1, &ldac, &ldadr, &pbec, &pbedr, &pbed2rho);

    if(beeforder==-1)
    {
	if(!(*addlda))
	{
	    *e = beefpbecfrac*pbec*(*r);
	    *dr = beefpbecfrac*pbedr;
	}
	else
	{
	    *e = (beefpbecfrac*pbec+ldac)*(*r);
	    *dr = beefpbecfrac*pbedr+ldadr;
	}
	*dg = beefpbecfrac*pbed2rho / (*r);
	return;
    }
    
    if(beeforder==-2)
    {
	*e = 0.;
	*dr = 0.;
	*dg = 0.;
    }
    else
    {
	if(!(*addlda))
	{
	    *e = (pbec-ldac)*(*r);
	    *dr = pbedr-ldadr;
	    *dg = pbed2rho / (*r);
	}
	else
	{
	    *e = pbec*(*r);
	    *dr = pbedr;
	    *dg = pbed2rho / (*r);
	}
    }
}

// evaluate bee exchange energy only
void beefxpot_(double *r, double *g, double *e, int *addlda)
{
    double s2,t,s,r43;
    const int n=nmax;
    const int i1=1;
    const int i2=1;

    r43 = pow(*r, 4./3.);

    s2 = *g*pix / (r43*r43);
    t = 2.*s2/(4.+s2)-1.;

    if(beeforder==-1)
    {
	calcleg(t);
	
	if(!(*addlda))
	    *e = (ddot_(&n, mi, &i1, L, &i2) - 1.) * r2e * r43;
	else
	    *e = ddot_(&n, mi, &i1, L, &i2) * r2e * r43;
	return;
    }
    
    if(beeforder>=0)
	*e = (*Ln[beeforder])(t) * r2e * r43;
    else
	*e = 0.;
}

// evaluate local part of bee correlation - energy only
void beeflocalcorrpot_(double *r, double *g, double *e, int *addlda)
{
    double rs, ldac, ldadr, pbec, pbedr, pbed2rho;
    
    if(beeforder>=0)
    {
	*e = 0.;
	return;
    }
    
    rs = invpi075tothird / pow(*r,1./3.);
    corpbe(rs, 0.5/r2k * sqrt(*g*rs) / (*r),
	(beeforder!=-2), 0, &ldac, &ldadr, &pbec, &pbedr, &pbed2rho);

    if(beeforder==-1)
    {
	if(!(*addlda))
	    *e = beefpbecfrac*pbec*(*r);
	else
	    *e = (beefpbecfrac*pbec+ldac)*(*r);
	return;
    }
    
    if(beeforder==-2)
	*e = 0.;
    else
    {
	if(!(*addlda))
	    *e = (pbec-ldac)*(*r);
	else
	    *e = pbec*(*r);
    }
}



// evaluate local part of bee correlation for spin polarized system
void beeflocalcorrspin_(double *r, double *z, double *g, double *e,
    double *drup, double *drdown, double *dg, int *addlda) {
    double rs, ldac, ldadrup, ldadrdown, pbec, pbedrup, pbedrdown, pbed2rho;
    
    if(beeforder>=0)
    {
	*e = 0.;
	*drup = 0.;
	*drdown = 0.;
	*dg = 0.;
	return;
    }
    
    rs = invpi075tothird / pow(*r,1./3.);
    corpbespin(rs, 0.5/r2k * sqrt(*g*rs) / (*r), *z,
	(beeforder!=-2), 1, &ldac, &ldadrup, &ldadrdown, &pbec,
	&pbedrup, &pbedrdown, &pbed2rho);

    if(beeforder==-1)
    {
	if(!(*addlda))
	{
	    *e =  beefpbecfrac*pbec*(*r);
	    *drup = beefpbecfrac*pbedrup;
	    *drdown = beefpbecfrac*pbedrdown;
	}
	else
	{
	    *e =  (beefpbecfrac*pbec+ldac)*(*r);
	    *drup = beefpbecfrac*pbedrup+ldadrup;
	    *drdown = beefpbecfrac*pbedrdown+ldadrdown;
	}
	*dg = beefpbecfrac*pbed2rho / (*r);
	return;
    }
    
    if(beeforder==-2)
    {
	*e = 0.;
	*drup = 0.;
	*drdown = 0.;
	*dg = 0.;
    }
    else
    {
	if(!(*addlda))
	{
	    *e = (pbec-ldac)*(*r);
	    *drup = pbedrup-ldadrup;
	    *drdown = pbedrdown-ldadrdown;
	}
	else
	{
	    *e = pbec*(*r);
	    *drup = pbedrup;
	    *drdown = pbedrdown;
	}
	*dg = pbed2rho / (*r);
    }
}

// evaluate local part of bee correlation for spin polarized system - energy only
void beeflocalcorrpotspin_(double *r, double *z, double *g, double *e, int *addlda)
{
    double rs, ldac, ldadrup, ldadrdown, pbec, pbedrup, pbedrdown, pbed2rho;
    
    if(beeforder>=0)
    {
	*e = 0.;
	return;
    }
    
    rs = invpi075tothird / pow(*r,1./3.);
    corpbespin(rs, 0.5/r2k * sqrt(*g*rs) / (*r), *z,
	(beeforder!=-2), 0, &ldac, &ldadrup, &ldadrdown, &pbec,
	&pbedrup, &pbedrdown, &pbed2rho);

    if(beeforder==-1)
    {
	if(!(*addlda))
	    *e = beefpbecfrac*pbec*(*r);
	else
	    *e = (beefpbecfrac*pbec+ldac)*(*r);
	return;
    }
    
    if(beeforder==-2)
	*e = 0.;
    else
    {
	if(!(*addlda))
	    *e = (pbec-ldac)*(*r);
	else
	    *e = pbec*(*r);
    }
}



// mode >= 0: for perturbed parameters --- calc Legendre order mode only
// -1: standard beefxc expansion coefficients
// -2: no exchange, LDA correlation only
// else: no exchange, PBE correlation (without LDA contribution)
void beefsetmode_(int *mode)
{
    beeforder = *mode;
}

// initialize pseudo random number generator
void beefrandinit_(unsigned int *seed)
{
    srandom(*seed);
}

// initialize pseudo random number generator with default seed
void beefrandinitdef_()
{
    srandom(defaultseed);
}

// calculate ensemble energies
void beefensemble_(double *beefxc, double *ensemble)
{
    double vec[nmax+2],randvec[nmax+1];
    const double alpha=1.;
    const double beta=0.;
    const int m=nmax+1;
    const int n=nmax+1;
    const int la=nmax+1;
    const int ix=1;
    const int iy=1;
    const int n2=nmax+2;

    for(int i=0;i<nsamples;i++)
    {
	for(int j=0;j<nmax+1;j++) randvec[j] = normrand();
	dgemv_("T", &m, &n, &alpha, beefmat, &la, randvec, &ix,
	    &beta, vec, &iy);
	vec[nmax+1] = -vec[nmax];
	ensemble[i] = ddot_(&n2, vec, &ix, beefxc, &iy);
    }
}
