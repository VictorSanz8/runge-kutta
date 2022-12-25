// rk4s.h
// Header file to declare class Runge-Kutta of order 4 for solving sistems of diferential equiations

#ifndef rk4s_h_
#define rk4s_h_

#include <iostream>
#include <string>
#include <cmath>
#include <fstream>

using namespace std;


class rk4s
{
public:
  // // Constructors // //
  rk4s( const long & dim , const  double * const ptrinitial , const double & a , const double & b , const long & N ); /*Usual constructor*/
  
  // // // // // // // // // //    I N F O R M A T I O N   // // // // // // // // // // // // // // // // // // // //
  /* (dim) is the number of diff equations that are in your sistem                                                 */
  /* (ptrinitial) is a pointer to the array of initial conditions of your dependent variables     */
  /* (a) and (b)  the lower and upper boundaries of the independent variabe, respectively   */
  /* (N) is the number of steps you want the methods to do                                                           */
  // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // 
  
  // // Destructor // //
  ~rk4s();
  

  // // Function to call the system of differential equations // //
  double function( const int & num , const double * values );



  // // Setters // //
  void set_initials( const double * const ptrinitial ); //Sets the pointer to the initial conditions if the system is unsolved
  void set_steps( const long & N ); //Sets the number of steps if the system is unsolved
  void set_boundaries( const double & min , const double & max ); //Sets the boundaries if the system is unsolved



  
  // // Method to solve the system // //
  void solve();

  // // Method to unsolve the system // // (Clears everything except the initial conditions)
  void clear();
  
  // // Method to print in to a file named as name.txt in columns // //
  void print_to_file( const string & name );

  
protected:
  // // K-functions of the runge-kutta 4 method // //
  double k1( const int & num , const double * const values )
  {
    double* aux_values;
    aux_values= new double[dim_+1];
    
    for( int i = 0 ; i < dim_+1 ; i++ )
      {
	*(aux_values+i) = *(values+i );
      }
    
    double res =  function( num , aux_values );
    delete aux_values;
    return res;
  };

  double k2( const int & num , const double * const values )
  {
    double* aux_values;
    aux_values = new double[dim_+1];
      
    for( int i = 0 ; i < dim_+1 ; i++ )
      {
	if( i == 0 )
	  {
	    *(aux_values+i) = *(values+i) + 0.5*h_;
	  }
	else
	  {
	    *(aux_values+i) = *(values+i) + 0.5*h_*k1( i , values  );
	  }
      }
    
    double res =  function( num , aux_values );
    delete aux_values;
    return res;
  };
    
  double k3( const int & num , const double * const values )
  {
    double* aux_values;
    aux_values = new double[dim_+1];
      
    for( int i = 0 ; i < dim_+1 ; i++ )
      {
	if( i == 0 )
	  {
	    *(aux_values+i) = *(values+i) + 0.5*h_;
	  }
	else
	  {
	    *(aux_values+i) = *(values+i)+ 0.5*h_*k2( i , values  );
	  }
      }
    
    double res =  function( num , aux_values );
    delete aux_values;
    return res;
  };
    
  double k4( const int & num , const double * const values )
  {
    double* aux_values;
    aux_values = new double[dim_+1];
      
    for( int i = 0 ; i < dim_+1 ; i++ )
      {
	if( i == 0 )
	  {
	    *(aux_values+i) = *(values+i) + h_;
	  }
	else
	  {
	    *(aux_values+i) = *(values+i) + h_ * k3( i , values  );
	  }
      }
    
    double res =  function( num , aux_values );
    delete aux_values;
    return res;
  };


  
  // // Protected methods of control // //
  void show_message( const char & mode ,  const int & code )
  {
    string msg,header;

    switch (mode)
      {
      case 'i':
	header = "Info: ";
	switch (code)
	  {
	  case 0:
	    msg = "The system is not valid.";
	    break;
	  case 1:
	    msg ="Cannot change the conditions if the system is already solved.";
	    break;
	  case 2:
	    msg ="Final value of the independent variable must be exclusively greater than the initial value.";
	    break;
	  case 3:
	    msg = "System of differential equations not solved yet.";
	    break;
	  case 4:
	    msg = "System of differential equations is already solved.";
	    break;
	  default:
	    msg = "Something went wrong.";
	  }
	cout << header+msg << endl;
	break;
	
      case 'e':
	header = "Error: ";
	switch (code)
	  {
	  case 0:
	    msg = "The system is not valid.";
	    break;
	  case 1:
	    msg = "Not valid number of differential equations, it must be positive.";
	    break;
	  case 2:
	    msg ="Final value of the independent variable must be exclusively greater than the initial value.";
	    break;
	  case 3:
	    msg = "The maximum numbers of steps must be at least greater than 1.";
	    break;
	  default:
	    msg = "Something went wrong.";
	  }
	cerr << header+msg << endl;
	break;
      }
  };


  
private:
  // // Private members // //
  int dim_;
  double h_;
  long steps_;
  double* values_;
  bool solved_;
  
  // // Private members to control if everything went ok // //
  bool ok_=true;
  int errorcode_= 0;
};

#endif
