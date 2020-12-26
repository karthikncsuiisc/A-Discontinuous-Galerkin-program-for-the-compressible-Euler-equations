 //Function to implement boundary condition.

#include "functions.h"
#include "twod_header.h"

void boundary_condition_euler_2d()
{

int ileft,iright,jleft;
double u1,u2x,u2y,u3;
double u1p,u2px,u2py,u3p;
double rhoin,velxin,velyin;
double rhog,velxg,velyg,u3g,presg,cg;
double vdotn;
double presin,cin,Machn,cref;
double velxref,velyref;
double velnref,veltref;
double velnb,veltb;
double velng,veltg;
double r1,r2,r3,r4;


for(int iface=1;iface<=nbface;iface++)
{
	ileft=intface[iface][1];
	iright=intface[iface][2];
	
	u1=unkel[ileft][0][0];
	u2x=unkel[ileft][1][0];
	u2y=unkel[ileft][2][0];
	u3=unkel[ileft][3][0];
	
	rhoin=u1;
	velxin=u2x/u1;
	velyin=u2y/u1;
//	presin=max((gama-1)*(u3-0.5*(velxin*velxin+velyin*velyin)*rhoin),1e-16);
//	cin=sqrt(gama*presin/rhoin);

	presin=(gama-1)*(u3-0.5*(velxin*velxin+velyin*velyin)*rhoin);
	cin=sqrt(gama*max(presin,hrd_lmtr_cin)/rhoin);

	vdotn=velxin*nx[iface]+velyin*ny[iface];

	if(bound_cond[iface][1]==2.0) //Wall boundary condition
	{
		velxg=velxin-2*vdotn*nx[iface];
		velyg=velyin-2*vdotn*ny[iface];
		
		rhog=rhoin;
		u3g=u3;
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
		u3g=presref/(gama-1)+0.5*rhog*(velxg*velxg+velyg*velyg);
		}
		else if(Machn>=1.0)
		{
		velxg=velxin;
		velyg=velyin;
		
		rhog=rhoin;
		u3g=u3;
		}
		else if(Machn<=0&&Machn>-1.0)
		{
			velnref=velxref*(nx[iface])+velyref*(ny[iface]);
			veltref=velxref*(-ny[iface])+velyref*(nx[iface]);
			
			velnb=vdotn;
			veltb=velxin*(-ny[iface])+velyin*(nx[iface]);
			
			cref=sqrt(gama*presref/rhoref);
			r1=velnref-2*cref/(gama-1);
			r2=veltref;
			r3=presref/pow(rhoref,gama);
			r4=velnb+2*cin/(gama-1);
			
			velng=0.5*(r4+r1);
			cg=(gama-1)*0.25*(r4-r1);
			veltg=r2;
			
			rhog=pow(cg*cg/(r3*gama),1.0/(gama-1.0));
			presg=r3*pow(rhog,gama);
			
			velxg=velng*nx[iface]-veltg*ny[iface];
			velyg=velng*ny[iface]+veltg*nx[iface];
			
			u3g=presg/(gama-1)+0.5*rhog*(velxg*velxg+velyg*velyg);
	
		}
		else if(Machn<1&&Machn>0.0)
		{
			velnref=velxref*(nx[iface])+velyref*(ny[iface]);
			veltref=velxref*(-ny[iface])+velyref*(nx[iface]);
			
			velnb=vdotn;
			veltb=velxin*(-ny[iface])+velyin*(nx[iface]);
			
			cref=sqrt(gama*presref/rhoref);
			r1=velnref-2*cref/(gama-1);
			r2=veltb;
			r3=presin/pow(rhoin,gama);
			r4=velnb+2*cin/(gama-1);
			
			velng=0.5*(r4+r1);
			cg=(gama-1)*0.25*(r4-r1);
			veltg=r2;
			
			rhog=pow(cg*cg/(r3*gama),1.0/(gama-1.0));
			presg=r3*pow(rhog,gama);
			
			velxg=velng*nx[iface]-veltg*ny[iface];
			velyg=velng*ny[iface]+veltg*nx[iface];
			
			u3g=presg/(gama-1)+0.5*rhog*(velxg*velxg+velyg*velyg);
		}
		
	}
	else if(bound_cond[iface][1]==6.0) //Periodic boundary condition
	{
		for(int jface=1;jface<=nbface;jface++)
		{
			if(jface!=iface)
			{
			if(bound_cond[jface][3]==bound_cond[iface][3])
			{
				jleft=intface[jface][1];
				
				u1p=unkel[jleft][0][0];
				u2px=unkel[jleft][1][0];
				u2py=unkel[jleft][2][0];
				u3p=unkel[jleft][3][0];
				
				rhog=u1p;
				velxg=u2px/u1p;
				velyg=u2py/u1p;
				u3g=u3p;
			}
			}
		}
	}

unkel[iright][0][0]=rhog;
unkel[iright][1][0]=rhog*velxg;
unkel[iright][2][0]=rhog*velyg;
unkel[iright][3][0]=u3g;
}

}
