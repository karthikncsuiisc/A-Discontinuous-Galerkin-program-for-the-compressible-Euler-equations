//Subroutine for calculating the right hand side of euler equations
#include "functions.h"
#include "twod_header.h"

void flux2D_eul_calc_p1()
{
	
int ileft,iright,ip1,ip2;
double u1left,u2xleft,u2yleft,u3left;
double u1right,u2xright,u2yright,u3right;
double rho0,velx0,vely0,pres0,c0,vdotn;
double fiplus0,fiplus1x,fiplus1y,fiplus2;
double fiminus0,fiminus1x,fiminus1y,fiminus2;
double Mach;

double rhoin,velxin,velyin,u3in,presin,cin;
double velxg,velyg,rhog,u3g,presg;
double velxref,velyref,Machn;
double velnref,veltref,velnb,veltb,cref,r1,r2,r3,r4,velng,veltg,cg;
double normx,normy;
double xcl,ycl,xgp,ygp,tx,ty;

//-----------------------------------------------------------
// Start of flux calculation
//-----------------------------------------------------------

for(int iface=1;iface<=naface;iface++)
{
	
for(int igp=1;igp<3;igp++)
{ 
//------------------------------Left cell calculation--------------------------	
	ileft=intface[iface][1];
	iright=intface[iface][2];
	ip1=intface[iface][3];
	ip2=intface[iface][4];
	

	u1left=unkel[ileft][0][0]+unkel[ileft][0][1]*Bxline_quad[iface][1][igp]
							 +unkel[ileft][0][2]*Byline_quad[iface][1][igp];
							 
   u2xleft=unkel[ileft][1][0]+unkel[ileft][1][1]*Bxline_quad[iface][1][igp]
							 +unkel[ileft][1][2]*Byline_quad[iface][1][igp];
							 
   u2yleft=unkel[ileft][2][0]+unkel[ileft][2][1]*Bxline_quad[iface][1][igp]
							 +unkel[ileft][2][2]*Byline_quad[iface][1][igp];

	u3left=unkel[ileft][3][0]+unkel[ileft][3][1]*Bxline_quad[iface][1][igp]
							 +unkel[ileft][3][2]*Byline_quad[iface][1][igp];
	

	if(iright>4*nelem)
	{
	xcl=0.5*(coord[ip1][1]+coord[ip2][1]);
	ycl=0.5*(coord[ip1][2]+coord[ip2][2]);
	
	tx=-ny[iface];
	ty=nx[iface];
	
	xgp=xcl+tx*(-igp)*pow(-1,igp%2)*length[iface]*sqrt(3.0/5.0)/2.0;
	ygp=ycl+ty*(-igp)*pow(-1,igp%2)*length[iface]*sqrt(3.0/5.0)/2.0;
	
	if(sqrt(xgp*xgp+ygp*ygp)>5.0)
	{
	normx=xgp/sqrt(pow(xgp,2)+pow(ygp,2));
	normy=ygp/sqrt(pow(xgp,2)+pow(ygp,2));
	}
	else
	{
	normx=-xgp/sqrt(pow(xgp,2)+pow(ygp,2));
	normy=-ygp/sqrt(pow(xgp,2)+pow(ygp,2));
//	cout<<iface<<" "<<igp<<" "<<nx[iface]*normy-ny[iface]*normx<<" "<<sqrt(xgp*xgp+ygp*ygp)<<endl;
	}
	}
	else
	{
	normx=nx[iface];
	normy=ny[iface];
	}
	
	rho0=u1left;
	velx0=u2xleft/u1left;
	vely0=u2yleft/u1left;
	pres0=(gama-1)*(u3left-0.5*(velx0*velx0+vely0*vely0)*rho0);
	c0=sqrt(gama*max(pres0,hrd_lmtr_cin)/rho0);
	
	vdotn=velx0*normx+vely0*normy;

	Mach=vdotn/c0;
	
	if(Mach<-1.0)
	{
	fiplus0=0;
	fiplus1x=0;
	fiplus1y=0;
	fiplus2=0;
    }
	else if(Mach>=-1.0&&Mach<=1.0)
	{
	double factor;
	
	factor=rho0*c0*(Mach+1)*(Mach+1)*0.25;
	
	fiplus0=factor;
	fiplus1x=factor*(velx0+normx*(-vdotn+2*c0)/gama);
	fiplus1y=factor*(vely0+normy*(-vdotn+2*c0)/gama);
	
	fiplus2=0.5*(velx0*velx0+vely0*vely0)-0.5*vdotn*vdotn;
	fiplus2=fiplus2+((gama-1)*vdotn+2.0*c0)*((gama-1)*vdotn+2.0*c0)/(2.0*gama*gama-2.0);
	fiplus2=factor*fiplus2;
    }
	else if(Mach>1.0)
	{
	fiplus0=rho0*vdotn;
	fiplus1x=rho0*vdotn*velx0+pres0*normx;
	fiplus1y=rho0*vdotn*vely0+pres0*normy;
	fiplus2=vdotn*(u3left+pres0);
    }
//------------------------End of Left cell calculation--------------------------	
    
//-------------Fj-(right cell) calcualtion---------------
if(iright<=nelem)
{
	u1right=unkel[iright][0][0]+unkel[iright][0][1]*Bxline_quad[iface][2][igp]
       						   +unkel[iright][0][2]*Byline_quad[iface][2][igp];

   u2xright=unkel[iright][1][0]+unkel[iright][1][1]*Bxline_quad[iface][2][igp]
       						   +unkel[iright][1][2]*Byline_quad[iface][2][igp];

   u2yright=unkel[iright][2][0]+unkel[iright][2][1]*Bxline_quad[iface][2][igp]
       						   +unkel[iright][2][2]*Byline_quad[iface][2][igp];

    u3right=unkel[iright][3][0]+unkel[iright][3][1]*Bxline_quad[iface][2][igp]
       						   +unkel[iright][3][2]*Byline_quad[iface][2][igp];
	
	rho0=u1right;
	velx0=u2xright/u1right;
	vely0=u2yright/u1right;
	pres0=(gama-1)*(u3right-0.5*(velx0*velx0+vely0*vely0)*rho0);
	c0=sqrt(gama*max(pres0,hrd_lmtr_cin)/rho0);
	
	vdotn=velx0*normx+vely0*normy;

	Mach=vdotn/c0;
		
		if(Mach>1.0)
		{
		fiminus0=0;
		fiminus1x=0;
		fiminus1y=0;
		fiminus2=0;
		}
		else if(Mach>=-1.0&&Mach<=1.0)
		{
			double factor;
		factor=-rho0*c0*(Mach-1)*(Mach-1)*0.25;
	
		fiminus0=factor;
		fiminus1x=factor*(velx0+normx*(-vdotn-2*c0)/gama);
		fiminus1y=factor*(vely0+normy*(-vdotn-2*c0)/gama);
	
		fiminus2=0.5*(velx0*velx0+vely0*vely0)-0.5*vdotn*vdotn;
		fiminus2=fiminus2+((gama-1)*vdotn-2.0*c0)*((gama-1)*vdotn-2.0*c0)/(2*gama*gama-2.0);
		fiminus2=factor*fiminus2;
		}
		else if(Mach<-1.0)
		{
		fiminus0=rho0*vdotn;
		fiminus1x=rho0*vdotn*velx0+pres0*normx;
		fiminus1y=rho0*vdotn*vely0+pres0*normy;
		fiminus2=vdotn*(u3right+pres0);
		}

}
else
{

//		cout<<iface<<" "<<igp<<" "<<nx[iface]<<" "<<ny[iface]<<" "<<sqrt(xgp*xgp+ygp*ygp)<<" "<<normx<<" "<<normy<<" "<<sqrt(normx*normx+normy*normy)<<endl;

		rhoin=u1left;
		velxin=u2xleft/u1left;
		velyin=u2yleft/u1left;
		u3in=u3left;
		
		presin=(gama-1)*(u3in-0.5*(velxin*velxin+velyin*velyin)*rhoin);
		cin=sqrt(gama*max(presin,hrd_lmtr_cin)/rhoin);
	
		vdotn=velxin*normx+velyin*normy;
		
		if(bound_cond[iface][1]==2.0) //Wall boundary condition
		{
				velxg=velxin-2*vdotn*normx;
				velyg=velyin-2*vdotn*normy;
		
				rhog=rhoin;
				presg=presin;
		}
		else if(bound_cond[iface][1]==4.0) //Far field boundary condition
		{
				Machn=vdotn/cin;
			
			    velxref=velref*cos(alpha*3.14159265/180.0);
				velyref=velref*sin(alpha*3.14159265/180.0);
				if(Machn<=-1.0)
				{
						velxg=velxref;
						velyg=velyref;
		
						rhog=rhoref;
						presg=presref;
				}
				else if(Machn>=1.0)
				{
						velxg=velxin;
						velyg=velyin;
		
						rhog=rhoin;
						presg=presref;
				}
				
				else if(Machn<=0&&Machn>-1.0)
				{
					velnref=velxref*(normx)+velyref*(normy);
					veltref=velxref*(-normy)+velyref*(normx);
			
					velnb=velxin*normx+velyin*normy;
					veltb=velxin*(-normy)+velyin*(normx);
					cref=sqrt(gama*presref/rhoref);
					
					r1=velnref-2.0*cref/(gama-1);
					r2=veltref;
					r3=presref/pow(rhoref,gama);
					r4=velnb+2.0*cin/(gama-1);
			
					velng=0.5*(r4+r1);
					cg=(gama-1)*0.25*(r4-r1);
					veltg=r2;
			
					rhog=pow(cg*cg/(r3*gama),1.0/(gama-1.0));
					presg=r3*pow(rhog,gama);
			
					velxg=velng*normx-veltg*normy;
					velyg=velng*normy+veltg*normx;
			
				}
				else if(Machn<1&&Machn>0.0)
				{
					
					velnref=velxref*(normx)+velyref*(normy);
					veltref=velxref*(-normy)+velyref*(normx);
			
					velnb=velxin*normx+velyin*normy;
					veltb=velxin*(-normy)+velyin*(normx);
			
					cref=sqrt(gama*presref/rhoref);
					
					r1=velnref-2.0*cref/(gama-1);
					r2=veltb;
					r3=presin/pow(rhoin,gama);
					r4=velnb+2*cin/(gama-1);
			
					velng=0.5*(r4+r1);
					cg=(gama-1)*0.25*(r4-r1);
					veltg=r2;
			
					rhog=pow(cg*cg/(r3*gama),1.0/(gama-1.0));
					presg=r3*pow(rhog,gama);
			
					velxg=velng*normx-veltg*normy;
					velyg=velng*normy+veltg*normx;
			
//					u3g=presg/(gama-1)+0.5*rhog*(velxg*velxg+velyg*velyg);
				 }
		}

	rho0=rhog;
	velx0=velxg;
	vely0=velyg;
	
	pres0=presg;
	c0=sqrt(gama*max(pres0,hrd_lmtr_cin)/rho0);
	
	vdotn=velx0*normx+vely0*normy;

	Mach=vdotn/c0;
		
		if(Mach>1.0)
		{
		fiminus0=0;
		fiminus1x=0;
		fiminus1y=0;
		fiminus2=0;
		}
		else if(Mach>=-1.0&&Mach<=1.0)
		{
			double factor;
		factor=-rho0*c0*(Mach-1)*(Mach-1)*0.25;
	
		fiminus0=factor;
		fiminus1x=factor*(velx0+normx*(-vdotn-2*c0)/gama);
		fiminus1y=factor*(vely0+normy*(-vdotn-2*c0)/gama);
	
		fiminus2=0.5*(velx0*velx0+vely0*vely0)-0.5*vdotn*vdotn;
		fiminus2=fiminus2+((gama-1)*vdotn-2.0*c0)*((gama-1)*vdotn-2.0*c0)/(2*gama*gama-2.0);
		fiminus2=factor*fiminus2;
		}
		else if(Mach<-1.0)
		{
		fiminus0=rho0*vdotn;
		fiminus1x=rho0*vdotn*velx0+pres0*normx;
		fiminus1y=rho0*vdotn*vely0+pres0*normy;
		fiminus2=vdotn*(u3g+pres0);
		}
}
		
    
 //---------------------Calculating rhspo------------------------------
    
   	if(ileft<=nelem)
   	{
    rhspo[ileft][0][0]=rhspo[ileft][0][0]-Weline_quad[igp]*(fiplus0+fiminus0)*length[iface]/2.0;
	rhspo[ileft][1][0]=rhspo[ileft][1][0]-Weline_quad[igp]*(fiplus1x+fiminus1x)*length[iface]/2.0;
    rhspo[ileft][2][0]=rhspo[ileft][2][0]-Weline_quad[igp]*(fiplus1y+fiminus1y)*length[iface]/2.0;
    rhspo[ileft][3][0]=rhspo[ileft][3][0]-Weline_quad[igp]*(fiplus2+fiminus2)*length[iface]/2.0;

    rhspo[ileft][0][1]=rhspo[ileft][0][1]-Weline_quad[igp]*(fiplus0+fiminus0)*Bxline_quad[iface][1][igp]*length[iface]/2.0;
	rhspo[ileft][1][1]=rhspo[ileft][1][1]-Weline_quad[igp]*(fiplus1x+fiminus1x)*Bxline_quad[iface][1][igp]*length[iface]/2.0;
    rhspo[ileft][2][1]=rhspo[ileft][2][1]-Weline_quad[igp]*(fiplus1y+fiminus1y)*Bxline_quad[iface][1][igp]*length[iface]/2.0;
    rhspo[ileft][3][1]=rhspo[ileft][3][1]-Weline_quad[igp]*(fiplus2+fiminus2)*Bxline_quad[iface][1][igp]*length[iface]/2.0;

    rhspo[ileft][0][2]=rhspo[ileft][0][2]-Weline_quad[igp]*(fiplus0+fiminus0)*Byline_quad[iface][1][igp]*length[iface]/2.0;
	rhspo[ileft][1][2]=rhspo[ileft][1][2]-Weline_quad[igp]*(fiplus1x+fiminus1x)*Byline_quad[iface][1][igp]*length[iface]/2.0;
    rhspo[ileft][2][2]=rhspo[ileft][2][2]-Weline_quad[igp]*(fiplus1y+fiminus1y)*Byline_quad[iface][1][igp]*length[iface]/2.0;
    rhspo[ileft][3][2]=rhspo[ileft][3][2]-Weline_quad[igp]*(fiplus2+fiminus2)*Byline_quad[iface][1][igp]*length[iface]/2.0;

	}
    
    if(iright<=nelem)
    {
    rhspo[iright][0][0]=rhspo[iright][0][0]+Weline_quad[igp]*(fiplus0+fiminus0)*length[iface]/2.0;
    rhspo[iright][1][0]=rhspo[iright][1][0]+Weline_quad[igp]*(fiplus1x+fiminus1x)*length[iface]/2.0;
    rhspo[iright][2][0]=rhspo[iright][2][0]+Weline_quad[igp]*(fiplus1y+fiminus1y)*length[iface]/2.0;
    rhspo[iright][3][0]=rhspo[iright][3][0]+Weline_quad[igp]*(fiplus2+fiminus2)*length[iface]/2.0;

    rhspo[iright][0][1]=rhspo[iright][0][1]+Weline_quad[igp]*(fiplus0+fiminus0)*Bxline_quad[iface][2][igp]*length[iface]/2.0;
    rhspo[iright][1][1]=rhspo[iright][1][1]+Weline_quad[igp]*(fiplus1x+fiminus1x)*Bxline_quad[iface][2][igp]*length[iface]/2.0;
    rhspo[iright][2][1]=rhspo[iright][2][1]+Weline_quad[igp]*(fiplus1y+fiminus1y)*Bxline_quad[iface][2][igp]*length[iface]/2.0;
    rhspo[iright][3][1]=rhspo[iright][3][1]+Weline_quad[igp]*(fiplus2+fiminus2)*Bxline_quad[iface][2][igp]*length[iface]/2.0;

    rhspo[iright][0][2]=rhspo[iright][0][2]+Weline_quad[igp]*(fiplus0+fiminus0)*Byline_quad[iface][2][igp]*length[iface]/2.0;
    rhspo[iright][1][2]=rhspo[iright][1][2]+Weline_quad[igp]*(fiplus1x+fiminus1x)*Byline_quad[iface][2][igp]*length[iface]/2.0;
    rhspo[iright][2][2]=rhspo[iright][2][2]+Weline_quad[igp]*(fiplus1y+fiminus1y)*Byline_quad[iface][2][igp]*length[iface]/2.0;
    rhspo[iright][3][2]=rhspo[iright][3][2]+Weline_quad[igp]*(fiplus2+fiminus2)*Byline_quad[iface][2][igp]*length[iface]/2.0;
	}
}
	
}

//-----------------------------------------------------------
// End of flux calculation
//-----------------------------------------------------------

//-----------------------------------------------------------
// Start of domain source contribution 
//-----------------------------------------------------------

double u1,u2x,u2y,u3;

for(int ielem=1;ielem<=nelem;ielem++)
{
	for(int igp=1; igp<4;igp++)
	{
		
		u1=unkel[ielem][0][0]+unkel[ielem][0][1]*Bxquad[ielem][igp]
       						 +unkel[ielem][0][2]*Byquad[ielem][igp];

       u2x=unkel[ielem][1][0]+unkel[ielem][1][1]*Bxquad[ielem][igp]
       						 +unkel[ielem][1][2]*Byquad[ielem][igp];

	   u2y=unkel[ielem][2][0]+unkel[ielem][2][1]*Bxquad[ielem][igp]
       						 +unkel[ielem][2][2]*Byquad[ielem][igp];
       					
  		u3=unkel[ielem][3][0]+unkel[ielem][3][1]*Bxquad[ielem][igp]
       						 +unkel[ielem][3][2]*Byquad[ielem][igp];
       						   
       	rho0=u1;
		velx0=u2x/u1;
		vely0=u2y/u1;
		pres0=(gama-1)*(u3-0.5*(velx0*velx0+vely0*vely0)*rho0);
       						   
       	rhspo[ielem][0][1]=rhspo[ielem][0][1]+Wequad[igp]*rho0*velx0*area[ielem]/deltax[ielem];
       	rhspo[ielem][0][2]=rhspo[ielem][0][2]+Wequad[igp]*rho0*vely0*area[ielem]/deltay[ielem];

       	rhspo[ielem][1][1]=rhspo[ielem][1][1]+Wequad[igp]*(rho0*velx0*velx0+pres0)*area[ielem]/deltax[ielem];
       	rhspo[ielem][1][2]=rhspo[ielem][1][2]+Wequad[igp]*rho0*velx0*vely0*area[ielem]/deltay[ielem];

       	rhspo[ielem][2][1]=rhspo[ielem][2][1]+Wequad[igp]*rho0*velx0*vely0*area[ielem]/deltax[ielem];
       	rhspo[ielem][2][2]=rhspo[ielem][2][2]+Wequad[igp]*(rho0*vely0*vely0+pres0)*area[ielem]/deltay[ielem];

       	rhspo[ielem][3][1]=rhspo[ielem][3][1]+Wequad[igp]*(u3+pres0)*velx0*area[ielem]/deltax[ielem];
       	rhspo[ielem][3][2]=rhspo[ielem][3][2]+Wequad[igp]*(u3+pres0)*vely0*area[ielem]/deltay[ielem];

	}
}


}


