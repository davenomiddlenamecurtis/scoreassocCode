/* DCERROR.CXX */ 
/* Copyright Dave Curtis 1994 */ 
/* dcurtis@hgmp.mrc.ac.uk */ 
/* no warranty or liability of any kind is accepted, expressed or implied */ 
 
/* DCERROR.CPP */ 
/* Copyright Dave Curtis 1991 */ 
#include "dcerror.hpp"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#ifdef __ZTC__
#include <bios.h>
#endif
#if 0
#include <fg.h>
#endif
#if defined wx_msw || defined wx_xview || defined wx_motif
#include "allwx.hpp"
int default_wx_error_func(char *format,va_list arg_ptr)
{
char s[2000];
vsprintf(s,format,arg_ptr);
return wxMessageBox(s,"Error",wxOK,mainframeptr);
}
#endif

extern "C" {
int same_address(const void *a, const void *b)
{
return a==b;
}

}

//void error_object::show_message(char *format)
//   { vprintf(format,arg_ptr); }

int error_object::operator()(int e,char *format,...)
  {
  va_list arg_ptr;
  if (!format)
#ifndef SUN_CC
    return operator()(e,"%s\n",strerror(e));
#else
    return operator()(e,"unknown type of error\n");
#endif
  status=e;
#if 0
  if (fg_inited) fg_term();
#endif
  va_start(arg_ptr,format);
  if (showing)
    {
    if (display_func) 
      display_func(format,arg_ptr);
    else
      {
#ifndef wx_x
      vfprintf(stderr,format,arg_ptr);
//      vprintf(format,arg_ptr);
#else
      default_wx_error_func(format,arg_ptr);
#endif
      }
    }
  va_end(arg_ptr);
#ifndef wx_x
  if (fatal) exit(e);
#else
  if (fatal) wxExit();
#endif
#if 0
  if (fg_inited) fg_init();
#endif
  return 0;
  }

error_object::error_object() { showing=fatal=1; status=0; display_func=0;}

error_object::~error_object() {;}
void error_object::set_display(void (*df)(char *,va_list)) { display_func=df; }
int error_object::stat() { return status; }
void error_object::clear() { status=0; }
void error_object::hide() { showing=0; }
void error_object::show() { showing=1; }
void error_object::warn() { fatal=0; }
void error_object::kill() { fatal=1; }


error_object dcerror;

//int _allocerr = 0;

void default_error(int n,char *s)
  { dcerror(n,s); }
// for old CPP toolkit

