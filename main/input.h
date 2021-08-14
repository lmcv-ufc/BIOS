// -------------------------------------------------------------------------
// input.h - input data functions.
// -------------------------------------------------------------------------
// Copyright (c) 2021 LMCV/UFC
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without 
// modification, are permitted provided that the following conditions
// are met:
//
// 1. Redistributions of source code must retain the above copyright 
//    notice, this list of conditions and the following disclaimer.
// 2. Redistributions in binary form must reproduce the above copyright 
//    notice, this list of conditions and the following disclaimer in the
//    documentation and/or other materials provided with the distribution.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS 
// "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT 
// LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
// PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT 
// HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, 
// SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED 
// TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR 
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF 
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING 
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS 
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
// -------------------------------------------------------------------------
//
// The input operations are made using Neutral File (text file with tags). The
// class cInpMap is responsible for store a map with string as key and pointer
// to read functions (pointer to void func(void)) as value. The global method
// ReadNeutralFile read an input file using a given cInpMap object. Global
// cInpMap objects are defined using Singletorn pattern, these are CtrlMap ->
// store inputs for .opt file, ProbMap -> store inputs for problem file.
// -------------------------------------------------------------------------
//
// bool Insert(string label, ReadFunc func ) 
//
//   label - read function label                                   (in)
//   func  - pointer to read function                              (in)
//
//
// This method insert a pair (label,func) into input map.
// -------------------------------------------------------------------------
//
// ReadFunc GetFunc(string label)
//
//   label - read function label                                   (in)
//
//
// This method return a pointer to function based on a given label.
// -------------------------------------------------------------------------

#ifndef _INPUT_H
#define _INPUT_H

#include <map>
#include <string>
#include <istream>


class cAbsReadEntry
{
  public:
  virtual void Read(std::istream&) = 0;
};

class cReadEntry : public cAbsReadEntry
{
  void (*rfunc)(std::istream&);

 public:
  cReadEntry(void (*rf)(std::istream&)) : rfunc(rf) { }
  void Read(std::istream &in) { rfunc(in); }
};

template <typename T>
class cReadObjEntry : public cAbsReadEntry
{
  T *obj;
  void (T::*rfunc)(std::istream&);

  public:
  cReadObjEntry(T *o, void (T::*rf)(std::istream&))
  {
    obj   = o;
    rfunc = rf;
  }

  void Read(std::istream &in) 
  { 
   (obj->*(rfunc))(in);
  }
};

#define makeRead(a)      new cReadEntry(a)
#define makeReadObj(a,b) new cReadObjEntry<a>(this,&a::b)

class cInpMap
{
 private:
  bool Feedback;
  std::map<std::string,cAbsReadEntry*> Map;

 public:
                 cInpMap(void);
                ~cInpMap(void);
  bool           Insert(std::string,cAbsReadEntry*); 
  cAbsReadEntry* GetFunc(std::string);

  void           SetFeedback(bool fb) { Feedback = fb; }

  bool           Insert(std::string,void (*) (void)) {return true;}
};

void     ReadNeutralFile(cInpMap&,std::ifstream&);    // Read from a neutral file
void     CopyFileData(std::ifstream&,std::ofstream&); // Copy the input data to the output file

// Singleton global variables based on [Meyer1996] 

cInpMap&  CtrlMap(void); 
cInpMap&  ProbMap(void); 

// Reference

//   [Meyer1996] - Scott Meyers, More Effective C++: 35 New Ways to Improve Your
//                 Programs and Designs. Addison-Wesley Pub Co, 1995

#endif
