// rk4s.cc
// Implementation of the rk4s class

#include "rk4s.h"
#include <iostream>
#include <string>
#include <cmath>
#include <fstream>

using namespace std;

// // constructors // //
rk4s::rk4s( const long & dim , const  double * const ptrinitial , const double & a , const double & b , const long & N)/*Usual constructor*/
{
  if(dim <= 0) /*Checking validity of dim*/
    {
      ok_ = false;
      errorcode_ = 1;
      show_message( 'e' , errorcode_ );
    }
  else
    {
      dim_ = dim;
    }
  
  if( b <= a ) /*Checking the validity of a and b*/
    {
      ok_ = false;
      errorcode_ = 2;
      show_message( 'e' , errorcode_  );
    }
  else if( N < 1) /*Checking the validity of N*/
    {
      ok_ = false;
      errorcode_ = 3;
      show_message( 'e' , errorcode_  );
    }
  else
    {
      h_ = (b-a)/(N);
      steps_ = N;
    }

  if( ok_ == true ) /*Set the remaining members if at first everything went ok*/
    {
      
      values_ = new double[ (steps_ + 1) * (dim_ + 1) ]; /*Create the values pointer that will add the values of the variables of each step*/

      for( long i = 0 ; i < dim+1 ; i++ ) /*Write the initial conditions in the row 0*/
	{
	  if( i == 0 )
	    {
	      *( values_ + (0)*(dim_+1) + i ) = a;
	    }
	  else
	    {
	      *( values_+(0)*(dim_+1)+i) = *(ptrinitial+i-1);
	    }
	}
      solved_ = false;
    }
}










// // Setters // //
void rk4s::set_initials( const double * const ptrinitial )//Sets the pointer to the initial conditions if the system is unsolved
{
  if( ok_ == false )
    {
      show_message('i',0);
	}
  else
    {
      if( solved_ == false )
	{
	  double aux = *values_;
	  delete values_;
	  values_ = new double[ (steps_ + 1) * (dim_ + 1) ]; /*Create the values pointer that will add the values of the variables of each step*/

	  for( long i = 0 ; i < dim_+1 ; i++ ) /*Write the initial conditions in the row 0*/
	    {
	      if( i == 0 )
		{
		  *( values_ + (0)*(dim_+1) + i ) = aux;
		}
	      else
		{
		  *( values_+(0)*(dim_+1)+i) = *(ptrinitial+i-1);
		}
	    }
	}
      else
	{
	  show_message('i',1);
	}
    }
}

void rk4s::set_steps( const long & N ) //Sets the number of steps if the system is unsolved
{
  if( ok_ == false )
    {
      show_message('i',0);
    }
  else
    {
      if( solved_ == false )
	{
	  if( N < 1) /*Checking the validity of N*/
	    {
	      show_message('e',3);
	    }
	  else
	    {
	      h_ = h_*steps_/N;
	      steps_ = N;
	    }
	}
      else
	{
	  show_message('i',1);
	}
    }
}
void rk4s::set_boundaries( const double & min , const double & max ) //Sets the minimum boundary if the system is unsolved
{
  if( ok_ == false )
    {
      show_message('i',0);
    }
  else
    {
      if( solved_ == false )
	{
	  if( min >= max ) /*Checking the validity of min and max*/
	    {
	      show_message('e',2);
	    }
	  else
	    {
	      h_ = (max-min)/steps_;
	      *values_ = min;
	    }
	}
      else
	{
	  show_message('i',1);
	}
    }
}








// // Destructors // //
rk4s::~rk4s()
{
  delete values_;
}





// // Method to solve the system // //
void rk4s::solve()
{
  if( ok_ == true)
    {
      if( solved_ == false )
	{
	  for( long i = 1 ; i < steps_+1 ; i++ )
	    {
	      double* previous_values = new double[dim_+1]; // store the line previous to the writing one to send to the functions
	      for( long j = 0 ; j < dim_+1 ; j++ )
		{
		  *(previous_values + j) = *(values_+(i-1)*(dim_+1)+j);
		}
	      for( long j = 0 ; j < dim_+1 ; j++ )
		{
		  if( j == 0 )
		    {
		      *(values_ + (i)*(dim_+1) + j) = *(previous_values + j) + h_;
		    }
		  else
		    {
		      *(values_ + (i)*(dim_+1) + j) = *(previous_values + j) + (h_/6.)*( k1( j , previous_values )+2*k2( j , previous_values)+2*k3(j,previous_values)+k4(j,previous_values) );
		    }
		}
	      delete previous_values;
	    }
	  cout << "Info: System solved successfully." << endl;
	  solved_ = true;
	}
      else
	{
	  show_message('i',4);
	}
    }
  else
    {
      show_message('i',0);
    }
}





// // Method to unsolve the system // //
void rk4s::clear()
{
  solved_ = false;
}





// // Method to print in to a file named as name.txt in columns // //
void rk4s::print_to_file( const string & name )
{
  if( ok_ == true )
    {
      if( solved_ == false )
	{
	  show_message('i',3);
	}
      else
	{
	  //Create object for output file
	  ofstream outfile;

	  //Create the nameof the file
	  string extension = ".txt";

	  //output file
	  outfile.open(name+extension);

	  for( long i = 0 ; i < steps_+1 ; i++ )
	    {
	      for( long j = 0 ; j < dim_+1 ; j++ )
		{
		  outfile << *( values_ + ( i )*(dim_+1) + j ) << " ";
		}
	      outfile << endl;
	    }
	  outfile.close();
	  cout <<"Info: Solution saved in "<< name+extension << " successfully."<<endl;
	}
    }
  else
    {
      show_message('i',0);
    }
}
