//Program to calculate the area of a triangle

#include "functions.h"

double triangle_area(double x1,double y1,double x2,double y2,double x3, double y3)
{
double area;
area=x1*(y2-y3)+y1*(x3-x2)+(x2*y3-x3*y2);
area=abs(0.5*area);
return(area);
}
