//Function to calculate minmode of two variable a,b

#include "functions.h"
#include "oned_header.h"
#include "twod_header.h"

double MClimiter(double a, double b,double c)
{

double val=1e10;
if(a*b>0.0&&a*c>0.0)
{
	val=min(2*abs(a),abs(b));
	val=min(val,2*abs(c));
	return val*(b/abs(b+1e-16));
}
else return 0.0;

}
