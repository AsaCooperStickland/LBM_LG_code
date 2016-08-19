//=======================================================
//
// = LIBRARY
//     Timer
//
// = FILENAME
//     RT_Timer.cc
//
// = AUTHOR(S)
//     Frederic Guidec (EPFL)
//     Alexandre Dupuis
//
// = VERSION
//     $Revision: 1.1 $ 
//
// = DATE RELEASED
//     $Date: 2000/05/29 14:33:44 $
//
// = COPYRIGHT
//     University of Geneva, Switzerland
//
//=======================================================

#include "RT_Timer.hh"
#include "conditions.hh"

#include <sys/time.h>
#include <iostream>


//=======================================================
// = CONSTRUCTOR
//=======================================================

RT_Timer::RT_Timer()
{

}

//=======================================================
RT_Timer::RT_Timer(const RT_Timer &other)
{
  *this = other;
}

//=======================================================
// = DESTRUCTOR
//=======================================================

RT_Timer::~RT_Timer()
{

}

//=======================================================
// = OPERATORS
//=======================================================
   
RT_Timer& RT_Timer::operator=(const RT_Timer& other)
{
  if (this == &other)
    return *this;

  Timer::operator=(other);

  return *this;
  
}

//=======================================================
// = ACCESSOR
//=======================================================

double RT_Timer::currentTime()
{
  double result;

  timeval timeVal;
  gettimeofday(&timeVal, 0);

  result = (double)timeVal.tv_sec
    + ((double)timeVal.tv_usec / 1e6);

  return result;
}
