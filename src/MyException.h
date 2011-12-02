#ifndef _MYEXCEPTION_H_
#define _MYEXCEPTION_H_

#include <iostream>
#include <cstring>

class MyException
{
public:
  std::string message;
  int  type;

  MyException()
  { 
    message.clear(); 
    type = 0; 
  }
  
  MyException(const char* s, int e) 
  { 
    message = s; 
    type = e;
  }

};

#endif
