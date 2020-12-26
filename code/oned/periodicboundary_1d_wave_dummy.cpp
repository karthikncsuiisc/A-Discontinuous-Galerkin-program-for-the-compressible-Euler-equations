//Function to decide which function to choose.

#include "functions.h"
#include "oned_header.h"

void periodicboundary_1d_wave()
{	

//Periodic boundary condition

unkel[0][0][0]=unkel[nelem][0][0];
unkel[1][0][0]=unkel[nelem+1][0][0];

unkel[nelem+2][0][0]=unkel[2][0][0];
unkel[nelem+3][0][0]=unkel[3][0][0];



if(ndegr>=2)
{
unkel[0][0][1]=unkel[nelem][0][1];
unkel[1][0][1]=unkel[nelem+1][0][1];

unkel[nelem+2][0][1]=unkel[2][0][1];
unkel[nelem+3][0][1]=unkel[3][0][1];
}

if(ndegr==3)
{
unkel[0][0][2]=unkel[nelem][0][2];
unkel[1][0][2]=unkel[nelem+1][0][2];

unkel[nelem+2][0][2]=unkel[2][0][2];
unkel[nelem+3][0][2]=unkel[3][0][2];
}

}
