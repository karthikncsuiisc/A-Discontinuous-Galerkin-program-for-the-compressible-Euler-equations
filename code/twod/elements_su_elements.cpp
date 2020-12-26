//Subroutine for finding elements surrounding point.

#include "functions.h"
#include "twod_header.h"

void elements_su_elements()
{
cout<<"-----------------------------------------"<<endl;
cout<<"Finding Elements surrounding Elements"<<endl;

int **help1;
int *poin1;

help1=new int* [4];
for(int i=0;i<=3;i++)
help1[i]=new int[3];
help1[1][1]=2;
help1[1][2]=3;
help1[2][1]=3;
help1[2][2]=1;
help1[3][1]=1;
help1[3][2]=2;
 
poin1=new int[npoint+1];
for(int ipoint=0; ipoint<=npoint;ipoint++) poin1[ipoint]=0;

//Allocating memory and Initializing
elsuel=new int* [nelem+1];
for(int ielem=0;ielem<=nelem; ielem++)
elsuel[ielem]=new int[4];

for(int ielem=0;ielem<=nelem;ielem++)
for(int ifael=0;ifael<=3;ifael++)
elsuel[ielem][ifael]=0;

int ip1,ip2,jelem,icoun,ieadj;

for(int ielem=1;ielem<=nelem;ielem++)
{
	for(int ifael=1;ifael<=3;ifael++)
	{
		ip1=inpoel[ielem][help1[ifael][1]];
		ip2=inpoel[ielem][help1[ifael][2]];
		poin1[ip1]=1;
		poin1[ip2]=1;
		
		for(int istor=esup2[ip1]+1;istor<=esup2[ip1+1];istor++)
		{
			jelem=esup1[istor];
			if(jelem!=ielem)
			{
				icoun=poin1[inpoel[jelem][1]]+poin1[inpoel[jelem][2]]+poin1[inpoel[jelem][3]];
				if(icoun==3-1) ieadj=jelem;
			}
		}
		elsuel[ielem][ifael]=ieadj;
		ieadj=0;
		poin1[ip1]=0;
		poin1[ip2]=0;
	}
}

/*for(int ielem=1;ielem<=nelem;ielem++)
cout<<ielem<<" "<<elsuel[ielem][1]<<" "<<elsuel[ielem][2]<<" "<<elsuel[ielem][3]<<endl;*/

//Deleting the element surround point arrays as it is no longer used.
delete [] esup1;
delete [] esup2;

for(int i=0;i<=3;i++)delete [] help1[i];
delete [] help1;

//cout<<"-----------------------------------------"<<endl;
}
