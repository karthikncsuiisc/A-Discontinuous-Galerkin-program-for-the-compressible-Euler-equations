//Function to calculate minmode of two variable a,b

#include "functions.h"
#include "oned_header.h"
#include "twod_header.h"

double minmod(double a, double b)
{
if(abs(a)<abs(b)&&(a*b)>0.0)	return a;
else if(abs(b)<abs(a)&&(a*b)>0.0) return b;
else return 0;
}
