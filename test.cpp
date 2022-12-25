//test.cpp
//This file is created in order to exhibit the funcionalities of the rk4s class to solve a system of differential equations using the runge-kutta of order 4 method.

#include <iostream>
#include "rk4s.h"

using namespace std;

// // Definition of the function to call the system of differential equations // //
// // You must always write it in this way, as is it not implemented inside the rk4s.cc file // //
double rk4s::function( const int & num , const double * values )
{
  // -                                     -                                    - //
  // *(values+0) = independent variable       //
  // *(values+i) = i-th variable                           //
  // -                                     -                                   - //
  double res;
  double x = *(values + 0);
  double y1 = *(values + 1);
  double y2 = *(values + 2);
  double y3 = *(values + 3);
  double y4 = *(values + 4);
  switch (num)
    {
    case 1:  /* res = y1' = f1(x,yj) */
      res = y2;
      break;
    case 2:  /* res = y2' = f2(x,yj) */
      res = -(y1-cos(x))/(pow(pow(y1-cos(x),2.)+pow(y3-sin(x),2.),3/2));
      break;
    case 3:  /* res = y3' = f3(x,yj) */
      res = y4;
      break;
    case 4:  /* res = y4' = f4(x,yj) */
      res = -(y3-sin(x))/(pow(pow(y1-cos(x),2.)+pow(y3-sin(x),2.),3/2));
      break;
    }
  return res;
}


// // Main function // //
int main()
{
  double initial_conditions1[] = { 1.05 , -0.25 , 0 , 0.5 }; //Writing the initial conditions
  
  double * ptr_initial = initial_conditions1; //Pointer to initial conditions

  rk4s system(4,ptr_initial,0,20,1000);
  
  system.solve();
  system.print_to_file("output1");

  system.clear(); // Set the system as unsolved to make changes to it

  double initial_conditions2[] = { 8.34709 ,0.19942 ,1.64137, -0.0534173 }; //Writing the new initial conditions
  ptr_initial = initial_conditions2;
  
  system.set_initials(ptr_initial);
  system.set_boundaries(10,20);
  system.set_steps(1000/2);
  
  system.solve();
  system.print_to_file("output2");
  
  return 0;
}
