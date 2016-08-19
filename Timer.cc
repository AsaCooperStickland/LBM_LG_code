//=======================================================
//
// = LIBRARY
//     Timer
//
// = FILENAME
//     Timer.cc
//
// = AUTHOR(S)
//     Frederic Guidec (EPFL)
//     Alexandre Dupuis
//
// = VERSION
//     $Revision: 1.1 $
//
// = DATE RELEASED
//     $Date: 2000/05/29 14:33:49 $
//
// = COPYRIGHT
//     University of Geneva, Switzerland
//
//=======================================================

#include "Timer.hh"
#include "conditions.hh"

#include <iostream>

//=======================================================
// = CONSTRUCTORS
//=======================================================

Timer::Timer()
{
  isRunning_ = false;
  elapsed_   = 0.0;
  startTime_ = 0.0;
  stopTime_  = 0.0;
}

//=======================================================
Timer::Timer(const Timer &other)
{
  *this = other;
}

//=======================================================
// = DESTRUCTOR
//=======================================================

Timer::~Timer()
{

}

//=======================================================
// = OPERATORS
//=======================================================
   
Timer& Timer::operator=(const Timer& other)
{

  if (this == &other)
    return *this;

  isRunning_ = other.isRunning_;
  startTime_ = other.startTime_;
  stopTime_  = other.stopTime_;
  elapsed_   = other.elapsed_;

  return *this;
}

//=======================================================
// = ACCESSORS
//=======================================================
   
void Timer::reset()
{
  elapsed_ = 0.0;
  isRunning_ = false;
}

//=======================================================
double Timer::elapsed()
{
  double result;

  if (isRunning_)
    result = elapsed_ + currentTime() - startTime_;
  else
    result = elapsed_;

  return result;
}

//=======================================================
bool Timer::isRunning()
{
  return isRunning_;
}

//=======================================================
// = ACTORS
//=======================================================

void Timer::start()
{
  require("Not running", !isRunning_);

  startTime_ = currentTime();
  isRunning_ = true;
}

//=======================================================
void Timer::stop()
{
  require("Is running", isRunning_);

  stopTime_ = currentTime();
  isRunning_ = false;
  elapsed_ += stopTime_ - startTime_;
}


