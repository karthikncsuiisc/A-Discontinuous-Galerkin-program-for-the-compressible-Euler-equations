//Subroutine for finding points surrounding points

#include "functions.h"
#include "twod_header.h"

void points_su_points()
{
cout<<"-----------------------------------------"<<endl;
cout<<"Finding points surrounding points"<<endl;

int *lpoin;
lpoin=new int[npoint+1];
for(int ipoint=0; ipoint<=npoint;ipoint++) lpoin[ipoint]=0;

//Allocating memory and Initializing
int istor =0;
int ielem,jpoint;
int totcount,swap;
int maxall=100000;
psup1=new int[maxall];
psup2=new int[npoint+2];

for(int ipoint=0; ipoint<=npoint+1;ipoint++)
psup2[ipoint]=0;

for(int ipoint=1; ipoint<=npoint; ipoint++)
{
	totcount=0;

//-------------------------------------------------------------------------
//   finding the point surrounding points including repitative points
//-------------------------------------------------------------------------
	for(int iesup=esup2[ipoint]+1; iesup<=esup2[ipoint+1]; iesup++)
	{
		ielem=esup1[iesup];
		for(int inode=1; inode<=3;inode++)
		{
			jpoint=inpoel[ielem][inode];
			if(jpoint!=ipoint)
			{
				lpoin[totcount]=jpoint;
				totcount++;
			}
		}
	        	     if (istor>=maxall)
        	      {
					  cout<<"hi"<<endl;
                   cout<<"Maxall is more than "<<maxall<<".Increase the size in function points_su_points.cpp"<<endl;
                   exit(1);
                  }
	}

//-------------------------------------------------------------------------
//   Sorting lpoin in increasing order matrix
//-------------------------------------------------------------------------
	 for (int i=0;i<totcount-1;i++)
     for (int j=0;j<totcount-i-1;j++)
 	 if (lpoin[j] > lpoin[j+1])
      {
        swap       = lpoin[j];
        lpoin[j]   = lpoin[j+1];
        lpoin[j+1] = swap;
      }
//-------------------------------------------------------------------------
//   Deleting the repitative points
//-------------------------------------------------------------------------
	istor=istor+1;
    psup1[istor]=lpoin[0];
    int x1=lpoin[0];
    for(int i=1 ;i<totcount;i++)
	if(x1 != lpoin[i])  
	{ 
   	  istor++;
      psup1[istor]=lpoin[i];
	  x1=lpoin[i];
	}
	psup2[ipoint+1]=istor;
}

/*cout<<"psup1 ";
for(int i=1;i<=istor;i++) cout<<" "<<psup1[i];
cout<<endl;

cout<<"psup2 ";
for(int ipoint=1; ipoint<=npoint+1; ipoint++) 
cout<<" "<<psup2[ipoint];
cout<<endl;*/

//delete [] lpoin;

//cout<<"-----------------------------------------"<<endl;
}
