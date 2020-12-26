//Subroutine for finding elements surrounding point.

#include "functions.h"
#include "twod_header.h"

void findfaces()
{
cout<<"-----------------------------------------"<<endl;
cout<<"Finding interface connectvity matrix"<<endl;

int maface,mbface,nnode,ip1,ip2,je;

nbface=0;
maface=1000000;
mbface=maface;
nnode=3;

intface=new int* [maface];


for(int ielem=1;ielem<=nelem;ielem++)
{
	for(int inode=1;inode<=nnode;inode++)
	{
	ip1=inode%nnode+1;
	ip2=ip1%nnode+1;
	je=elsuel[ielem][inode];
	
	if(je==0)
	{
		nbface++;
		if(nbface>=mbface)
		{
		cout<<"Please increase the value of mbface= "<<nbface;
		exit(1);
		}
		intface[nbface]=new int[5];
		intface[nbface][1]=ielem;
		intface[nbface][2]=nelem+nbface;
		intface[nbface][3]=inpoel[ielem][ip1];
		intface[nbface][4]=inpoel[ielem][ip2];
		elsuel[ielem][inode]=nelem+nbface;
	}
    }
}

naface=nbface;

for(int ielem=1;ielem<=nelem;ielem++)
{
	for(int inode=1;inode<=3;inode++)
	{
		ip1=inode%3+1;
		ip2=ip1%nnode+1;
		je=elsuel[ielem][inode];
		
		if(je>ielem&&je<=nelem)
		{
			naface++;
			if(naface>maface)
			{
			cout<<"Please increase size of maface in subroutine findface.cpp"<<naface<<endl;
			exit(1);
			}
    		intface[naface]=new int[5];
			intface[naface][1]=ielem;
			intface[naface][2]=je;
			intface[naface][3]=inpoel[ielem][ip1];
			intface[naface][4]=inpoel[ielem][ip2];
				
		}
	}
}

/*cout<<"Boundary Faces"<<endl;
for(int iface=1;iface<=nbface;iface++)
cout<<iface<<" "<<intface[iface][1]<<" "<<intface[iface][2]<<" "<<intface[iface][3]<<" "<<intface[iface][4]<<endl;

cout<<"Rest of the Faces"<<endl;
for(int iface=nbface+1;iface<=naface;iface++)
cout<<iface<<" "<<intface[iface][1]<<" "<<intface[iface][2]<<" "<<intface[iface][3]<<" "<<intface[iface][4]<<endl;*/

cout<<"number of boundary face="<<nbface<<endl;
cout<<"Number of all faces="<<naface<<endl;

//-----------------------------------------------------------------------------
//				Calculating length of each face and area of each element
//-----------------------------------------------------------------------------

int ileft,iright;
int *ip;
ip =new int[nnode+1];
double *x,*y;
x=new double[nnode+1];
y=new double[nnode+1];

length=new double[naface+1];
nx=new double[naface+1];
ny=new double[naface+1];

for(int iface=1;iface<=naface;iface++)
{
	ip[1]=intface[iface][3];
	ip[2]=intface[iface][4];
    
    x[1]=coord[ip[1]][1];
    y[1]=coord[ip[1]][2];
    x[2]=coord[ip[2]][1];
    y[2]=coord[ip[2]][2];

	length[iface]=sqrt((x[2]-x[1])*(x[2]-x[1])+(y[2]-y[1])*(y[2]-y[1]));

	nx[iface]=(y[2]-y[1])/length[iface];
	ny[iface]=(-x[2]+x[1])/length[iface];
}

area=new double[nelem+nbface+1];
centerx=new double[nelem+nbface+1];
centery=new double[nelem+nbface+1];

for(int ielem=1;ielem<=nelem;ielem++)
{
	ip[1]=inpoel[ielem][1];
	ip[2]=inpoel[ielem][2];
	ip[3]=inpoel[ielem][3];
	
	x[1]=coord[ip[1]][1];
    y[1]=coord[ip[1]][2];
    
    x[2]=coord[ip[2]][1];
    y[2]=coord[ip[2]][2];
    
    x[3]=coord[ip[3]][1];
    y[3]=coord[ip[3]][2];
    
    centerx[ielem]=(x[1]+x[2]+x[3])/3.0;
    centery[ielem]=(y[1]+y[2]+y[3])/3.0;
   
    area[ielem]=triangle_area(x[1],y[1],x[2],y[2],x[3],y[3]);	
}
//--------------------------------------------------------------------------
//		End of Calculating length of each face and area of each element
//--------------------------------------------------------------------------

//---------------------------------------------------------------------
//							Adding ghost cells
//---------------------------------------------------------------------

double height;

for(int iface=1;iface<=nbface;iface++)
{
	ileft=intface[iface][1];
	iright=intface[iface][2];
	ip[1]=intface[iface][3];
	ip[2]=intface[iface][4];
	

	if(iright>nelem)
	{
		for(int inode=1;inode<=nnode;inode++)
		if(inpoel[ileft][inode]!=ip[1]&&inpoel[ileft][inode]!=ip[2])	
			{
				ip[3]=inpoel[ileft][inode];
				break;
			}
	}
	else
	{
		cout<<"right in find faces"<<endl; exit(1);		
		for(int inode=1;inode<=nnode;inode++)
		if(inpoel[iright][inode]!=ip[1]&&inpoel[iright][inode]!=ip[2])	
			{
				ip[3]=inpoel[iright][inode]; 
				break;
			}
	}
 
	inpoel[nelem+iface][1]=ip[2];
	inpoel[nelem+iface][2]=ip[1];
	inpoel[nelem+iface][3]=npoint+iface;
	
    height=area[ileft]/(0.5*length[iface]);

	coord[npoint+iface][1]=coord[ip[3]][1]+nx[iface]*2*height;
	coord[npoint+iface][2]=coord[ip[3]][2]+ny[iface]*2*height;
	
	x[1]=coord[inpoel[nelem+iface][1]][1];
	y[1]=coord[inpoel[nelem+iface][1]][2];
	
	x[2]=coord[inpoel[nelem+iface][2]][1];
	y[2]=coord[inpoel[nelem+iface][2]][2];

	x[3]=coord[npoint+iface][1];
	y[3]=coord[npoint+iface][2];

    centerx[nelem+iface]=(x[1]+x[2]+x[3])/3.0;
    centery[nelem+iface]=(y[1]+y[2]+y[3])/3.0;

	area[nelem+iface]=triangle_area(x[1],y[1],x[2],y[2],x[3],y[3]);
    
}
nelemtot=nelem+nbface;

//Finding minimum and maximum values of x and y
double xmin,ymin,xmax,ymax;
xmin=1e10;
ymin=1e10;
xmax=-1e10;
ymax=-1e10;

for(int ipoint=1; ipoint<npoint;ipoint++)
{
	xmin=min(xmin,coord[ipoint][1]);
	ymin=min(ymin,coord[ipoint][2]);
	xmax=max(xmax,coord[ipoint][1]);
	ymax=max(ymax,coord[ipoint][2]);
}

cout<<"xmin "<<xmin<<" ymin "<<ymin<<" xmax "<<xmax<<" ymax "<<ymax<<endl;

//--------------------------------------------------------------------------
//							End of Adding ghost cells
//--------------------------------------------------------------------------

//--------------------------------------------------------------------------
//		Assigning boundary conditions
//--------------------------------------------------------------------------

double xi1,yi1,xi2,yi2;
double xj1,yj1,xj2,yj2;
int *jp;
jp =new int[nnode];
int *mark;
mark=new int[nelem+nbface+1];
for(int ielem=1;ielem<=nelem+nbface;ielem++) mark[ielem]=0;

bound_cond=new double* [nface+1]; 
for (int j=1; j <=nface; j++) bound_cond[j] = new double[4];

cout<<"check1"<<endl;
for(int iface=1;iface<=nbface;iface++)
{
	ip[1]=intface[iface][3];
	ip[2]=intface[iface][4];
	xi1=coord[ip[1]][1];
	yi1=coord[ip[1]][2];
	xi2=coord[ip[2]][1];
	yi2=coord[ip[2]][2];

	if(mark[iface]==0)
	{
	for(int jface=1;jface<=nface;jface++)
	{
		jp[1]=(int)bface[jface][1];
		jp[2]=(int)bface[jface][2];
		
		xj1=coord[jp[1]][1];
		yj1=coord[jp[1]][2];
		xj2=coord[jp[2]][1];
		yj2=coord[jp[2]][2];
		
		if((xi1==xj1&&xi2==xj2&&yi1==yj1&&yi2==yj2)||(xi1==xj2&&xi2==xj1&&yi1==yj2&&yi2==yj1))
		{
			bound_cond[iface][1]=bface[jface][3];
			bound_cond[iface][2]=bface[jface][4];
			bound_cond[iface][3]=bface[jface][5];
			mark[iface]=1;
		}
	}
	}
}

/*for(int ielem=1;ielem<=nelem+nbface;ielem++)
cout<<ielem<<" "<<centerx[ielem]<<" "<<centery[ielem]<<endl;*/
/*for(int iface=1;iface<=nbface;iface++)
cout<<iface<<" point1 "<<intface[iface][3]<<" point2: "<<intface[iface][4]<<" "<<bound_cond[iface][1]<<" "<<bound_cond[iface][2]<<" "<<bound_cond[iface][3]<<endl;
*/

//Deallocating memory
delete [] ip;
delete [] jp;
delete [] x;
delete [] y;
delete [] mark;

for(int ielem=1;ielem<=nface;ielem++) delete [] bface[ielem];

delete [] bface;

}
