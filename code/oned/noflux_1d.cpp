//Function to decide which function to choose.

#include "functions.h"
#include "oned_header.h"

void noflux_1d()
{

for(int ivar=0;ivar<nvar;ivar++)
for(int ideg=0;ideg<ndegr;ideg++)
{
unkel[0][ivar][ideg]=unkel[2][ivar][ideg];
unkel[1][ivar][ideg]=unkel[2][ivar][ideg];

unkel[nelem+2][ivar][ideg]=unkel[nelem+1][ivar][ideg];
unkel[nelem+3][ivar][ideg]=unkel[nelem+1][ivar][ideg];
}
	
}
