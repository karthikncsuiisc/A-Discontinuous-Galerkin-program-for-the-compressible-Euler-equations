//Subroutine for finding elements surrounding point.

#include "functions.h"
#include "twod_header.h"

void elements_su_points()
{
cout<<"-----------------------------------------"<<endl;
cout<<"Finding Elements surrounding points"<<endl;

//Allocating memory and Initializing
esup2=new int[npoint+2];
for(int ipoint=0; ipoint<=npoint+1;ipoint++) esup2[ipoint]=0;

//---------------------------------------------------------
// Counting the number of elements surrounding each point
//---------------------------------------------------------
for(int ielem=1; ielem<=nelem; ielem++)
for(int inode=1; inode<=3; inode++)
esup2[inpoel[ielem][inode]+1]=esup2[inpoel[ielem][inode]+1]+1;

for( int ipoint =2; ipoint<=npoint+1;ipoint++)
esup2[ipoint]=esup2[ipoint]+esup2[ipoint-1];

int mesup=esup2[npoint+1];
esup1=new int[mesup+1];
for(int ielem=1;ielem<=nelem;ielem++)
for(int inode=1;inode<=3;inode++)
{
	int istore,ipoint;
	ipoint=inpoel[ielem][inode];
	istore=esup2[ipoint]+1;
	esup2[ipoint]=istore;
	esup1[istore]=ielem;
}

for(int ipoint=npoint+1; ipoint>=2;ipoint--)
esup2[ipoint]=esup2[ipoint-1];

esup2[1]=0;

/*cout<<"esup1 ";
for(int i=1;i<=mesup;i++) cout<<" "<<esup1[i];
cout<<endl;

cout<<"esup2 ";
for(int inode=1;inode<=npoint+1;inode++) cout<<" "<<esup2[inode];*/

}
