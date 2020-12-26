//Function to calculate the values of shape function at gauss points in each cell and face.

#include "functions.h"
#include "twod_header.h"

void gauss_point_values()
{
	
	double *x,*y,xc,yc,xmax,ymax,xmin,ymin,xgp,ygp;
	int *ip;
	double **L;
	
	x=new double[4];
	y=new double[4];
	ip=new int[4];
	L=new double* [4];
	
	for(int i=0;i<=3;i++) L[i]=new double[4];
	
	
//-----------------------------------------------------------------------------
//	shape function for each cell at gauss points
//	Bxquad[ielem][1:3]= shape function(Bx) at three gauss points for ielem
//	Byquad[ielem][1:3]= shape function(By) at three gauss points for ielem
//-----------------------------------------------------------------------------
	deltax=new double[nelem+nface+1];
	deltay=new double[nelem+nface+1];
	
	Bxquad=new double* [nelem+nface+1];
	Byquad=new double* [nelem+nface+1];
	for(int ielem=1;ielem<=nelem+nface;ielem++)
	{
	Bxquad[ielem]=new double[4];
	Byquad[ielem]=new double[4];
	}
	Wequad=new double[4];

	L[1][1]=0.5;		L[1][2]=0.0;		L[1][3]=0.5;
	L[2][1]=0.5;		L[2][2]=0.5;		L[2][3]=0.0;
	L[3][1]=0.0;		L[3][2]=0.5;		L[3][3]=0.5;
	Wequad[1]=1.0/3.0;	Wequad[2]=1.0/3.0;  Wequad[3]=1.0/3.0;	
	
for(int ielem=1; ielem<=nelem+nface;ielem++)
{
		
	for(int ipoint=1;ipoint<=3;ipoint++)
	ip[ipoint]=inpoel[ielem][ipoint];

	for(int ipoint=1;ipoint<=3;ipoint++)
    {
		x[ipoint]=coord[ip[ipoint]][1];
		y[ipoint]=coord[ip[ipoint]][2];
	}

	xc=(x[1]+x[2]+x[3])/3.0;
	yc=(y[1]+y[2]+y[3])/3.0;
	
	xmax=max(x[1],x[2]);
	xmax=max(xmax,x[3]);
	ymax=max(y[1],y[2]);
	ymax=max(ymax,y[3]);

	xmin=min(x[1],x[2]);
	xmin=min(xmin,x[3]);
	ymin=min(y[1],y[2]);
	ymin=min(ymin,y[3]);
	
	deltax[ielem]=(xmax-xmin)/2.0;
	deltay[ielem]=(ymax-ymin)/2.0;
	
	for(int igp=1;igp<=3;igp++)
	{
		xgp=L[igp][1]*x[1]+L[igp][2]*x[2]+L[igp][3]*x[3];
		ygp=L[igp][1]*y[1]+L[igp][2]*y[2]+L[igp][3]*y[3];
		
		Bxquad[ielem][igp]=(xgp-xc)/deltax[ielem];
		Byquad[ielem][igp]=(ygp-yc)/deltay[ielem];
	}
}

//-----------------------------------------------------------------------------
//					End of shape function for each cell at gauss points
//-----------------------------------------------------------------------------

cout<<"End of calculaitng quad shape funcions for domain"<<endl;
//-----------------------------------------------------------------------------
//	shape function for each cell at gauss points
//	Bxcubic[ielem][0:3]= shape function(Bx) at four gauss points for ielem
//	Bycubic[ielem][0:3]= shape function(By) at four gauss points for ielem
//-----------------------------------------------------------------------------
	Bxcubic=new double* [nelem+nface+1];
	Bycubic=new double* [nelem+nface+1];
	for(int ielem=1;ielem<=nelem+nface;ielem++)
	{
	Bxcubic[ielem]=new double[4];
	Bycubic[ielem]=new double[4];
	}
	Wecubic=new double[4];

	L[0][1]=(1.0/3.0);		L[0][2]=(1/3.0);		L[0][3]=(1.0/3.0);
	L[1][1]=0.6;		L[1][2]=0.2;		L[1][3]=0.2;
	L[2][1]=0.2;		L[2][2]=0.6;		L[2][3]=0.2;
	L[3][1]=0.2;		L[3][2]=0.2;		L[3][3]=0.6;
	
	Wecubic[0]=-(27.0/48.0);
	Wecubic[1]=+(25.0/48.0);
	Wecubic[2]=+(25.0/48.0);
	Wecubic[3]=+(25.0/48.0);
	
for(int ielem=1; ielem<=nelem+nface;ielem++)
{
		
	for(int ipoint=1;ipoint<=3;ipoint++)
	ip[ipoint]=inpoel[ielem][ipoint];

	for(int ipoint=1;ipoint<=3;ipoint++)
    {
		x[ipoint]=coord[ip[ipoint]][1];
		y[ipoint]=coord[ip[ipoint]][2];
	}

	xc=(x[1]+x[2]+x[3])/3.0;
	yc=(y[1]+y[2]+y[3])/3.0;
	
	for(int igp=0;igp<=3;igp++)
	{
		xgp=L[igp][1]*x[1]+L[igp][2]*x[2]+L[igp][3]*x[3];
		ygp=L[igp][1]*y[1]+L[igp][2]*y[2]+L[igp][3]*y[3];
		
		Bxcubic[ielem][igp]=(xgp-xc)/deltax[ielem];
		Bycubic[ielem][igp]=(ygp-yc)/deltay[ielem];
	}
}

//-----------------------------------------------------------------------------
//					End of cubic shape function for each cell at gauss points
//-----------------------------------------------------------------------------
cout<<"End of calculaitng cubic shape funcions for domain"<<endl;

//-----------------------------------------------------------------------------
//	shape function for each face at gauss points
//	Bxline_lin[iface][1]= shape function(Bx) for linear  order for left cell
//	Bxline_lin[iface][2]= shape function(Bx) for linear  order for right cell
//	Byline_lin[iface][1]= shape function(By) for linear  order for left cell
//	Byline_lin[iface][2]= shape function(By) for linear  order for right cell
//-------------------------------------------------------------------------------
//	Bxline_quad[iface][1][1:2]= shape function(Bx) for quadratic order for left cell at two gauss points
//	Bxline_quad[iface][2][1:2]= shape function(Bx) for quadratic order for right cell at two gauss points
//	Byline_quad[iface][1][1:2]= shape function(By) for quadratic order for left cell at two gauss points
//	Byline_quad[iface][2][1:2]= shape function(By) for quadratic order for right cell at two gauss points
//-----------------------------------------------------------------------------

	Bxline_lin=new double* [naface+1];
	Byline_lin=new double* [naface+1];
	Weline_lin=2.0;
	
	Bxline_quad=new double** [naface+1];
	Byline_quad=new double** [naface+1];
	Weline_quad=new double[3];
	
	Weline_quad[1]=1.0;
	Weline_quad[2]=1.0;

	Bxline_cubic=new double** [naface+1];
	Byline_cubic=new double** [naface+1];
	Weline_cubic=new double[3];
	
	Weline_cubic[0]=(8.0/9.0);
	Weline_cubic[1]=(5.0/9.0);
	Weline_cubic[2]=(5.0/9.0);
	
	cout<<"check0"<<endl;

	for(int iface=1;iface<=naface;iface++)
	{
	Bxline_lin[iface]=new double[3];
	Byline_lin[iface]=new double[3];

	Bxline_quad[iface]=new double* [3];
	Byline_quad[iface]=new double* [3];

	Bxline_cubic[iface]=new double* [3];
	Byline_cubic[iface]=new double* [3];
	}
	cout<<"check1"<<endl;
	
	for(int iface=1;iface<=naface;iface++) 
	{
	Bxline_quad[iface][1]=new double[3];
	Bxline_quad[iface][2]=new double[3];
	
	Byline_quad[iface][1]=new double[3];
	Byline_quad[iface][2]=new double[3];

	Bxline_cubic[iface][1]=new double[3];
	Bxline_cubic[iface][2]=new double[3];
	
	Byline_cubic[iface][1]=new double[3];
	Byline_cubic[iface][2]=new double[3];
	}
	cout<<"check2"<<endl;

	for(int iface=1;iface<=naface;iface++)
	{
		Bxline_cubic[iface][1][0]=0.0;
		Bxline_cubic[iface][1][1]=0.0;
		Bxline_cubic[iface][1][2]=0.0;
		Bxline_cubic[iface][2][0]=0.0;
		Bxline_cubic[iface][2][1]=0.0;
		Bxline_cubic[iface][2][2]=0.0;

		Byline_cubic[iface][1][0]=0.0;
		Byline_cubic[iface][1][1]=0.0;
		Byline_cubic[iface][1][2]=0.0;
		Byline_cubic[iface][2][0]=0.0;
		Byline_cubic[iface][2][1]=0.0;
		Byline_cubic[iface][2][2]=0.0;

		Bxline_quad[iface][1][1]=0.0;
		Bxline_quad[iface][1][2]=0.0;
		Bxline_quad[iface][2][1]=0.0;
		Bxline_quad[iface][2][2]=0.0;

		Byline_quad[iface][1][1]=0.0;
		Byline_quad[iface][1][2]=0.0;
		Byline_quad[iface][2][1]=0.0;
		Byline_quad[iface][2][2]=0.0;
		
		Bxline_lin[iface][1]=0.0;
		Bxline_lin[iface][2]=0.0;
	}
	cout<<"check3"<<endl;

	
int ip1,ip2,ileft,iright;
double xcl,ycl,xp1,xp2,yp1,yp2,tx,ty;
for(int iface=1;iface<=naface;iface++)
{
	ileft=intface[iface][1];
	iright=intface[iface][2];
	ip1=intface[iface][3];
	ip2=intface[iface][4];

	xp1=coord[ip1][1];
	yp1=coord[ip1][2];
	xp2=coord[ip2][1];
	yp2=coord[ip2][2];
	
	xcl=0.5*(xp1+xp2);
	ycl=0.5*(yp1+yp2);
	
	Bxline_lin[iface][1]=(xcl-centerx[ileft])/deltax[ileft];
	Byline_lin[iface][1]=(ycl-centery[ileft])/deltay[ileft];

	Bxline_lin[iface][2]=(xcl-centerx[iright])/deltax[iright];
	Byline_lin[iface][2]=(ycl-centery[iright])/deltay[iright];


	tx=-ny[iface];
	ty=nx[iface];
		
	xgp=xcl+tx*length[iface]/(2.0*sqrt(3.0));
	ygp=ycl+ty*length[iface]/(2.0*sqrt(3.0));

	Bxline_quad[iface][1][1]=(xgp-centerx[ileft])/deltax[ileft];
	Byline_quad[iface][1][1]=(ygp-centery[ileft])/deltay[ileft];

	Bxline_quad[iface][2][1]=(xgp-centerx[iright])/deltax[iright];
	Byline_quad[iface][2][1]=(ygp-centery[iright])/deltay[iright];

	xgp=xcl-tx*length[iface]/(2.0*sqrt(3.0));
	ygp=ycl-ty*length[iface]/(2.0*sqrt(3.0));

	Bxline_quad[iface][1][2]=(xgp-centerx[ileft])/deltax[ileft];
	Byline_quad[iface][1][2]=(ygp-centery[ileft])/deltay[ileft];

	Bxline_quad[iface][2][2]=(xgp-centerx[iright])/deltax[iright];
	Byline_quad[iface][2][2]=(ygp-centery[iright])/deltay[iright];
	
	
	xgp=xcl;
	ygp=ycl;

	Bxline_cubic[iface][1][0]=(xgp-centerx[ileft])/deltax[ileft];
	Byline_cubic[iface][1][0]=(ygp-centery[ileft])/deltay[ileft];

	Bxline_cubic[iface][2][0]=(xgp-centerx[iright])/deltax[iright];
	Byline_cubic[iface][2][0]=(ygp-centery[iright])/deltay[iright];
	
	xgp=xcl+(sqrt(15.0)/5)*tx*length[iface]/2.0;
	ygp=ycl+(sqrt(15.0)/5)*ty*length[iface]/2.0;

	Bxline_cubic[iface][1][1]=(xgp-centerx[ileft])/deltax[ileft];
	Byline_cubic[iface][1][1]=(ygp-centery[ileft])/deltay[ileft];

	Bxline_cubic[iface][2][1]=(xgp-centerx[iright])/deltax[iright];
	Byline_cubic[iface][2][1]=(ygp-centery[iright])/deltay[iright];

	xgp=xcl-(sqrt(15.0)/5)*tx*length[iface]/2.0;
	ygp=ycl-(sqrt(15.0)/5)*ty*length[iface]/2.0;

	Bxline_cubic[iface][1][2]=(xgp-centerx[ileft])/deltax[ileft];
	Byline_cubic[iface][1][2]=(ygp-centery[ileft])/deltay[ileft];

	Bxline_cubic[iface][2][2]=(xgp-centerx[iright])/deltax[iright];
	Byline_cubic[iface][2][2]=(ygp-centery[iright])/deltay[iright];

}

//-----------------------------------------------------------------------------
//	End of shape function for each face at gauss points
//-----------------------------------------------------------------------------
cout<<"End of calculaitng shape funcions for faces"<<endl;

	delete [] x;
	delete [] y;
	delete [] ip;
	for(int i=0;i<4;i++) delete [] L[i];
	delete [] L;

}
