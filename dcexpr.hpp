/* DCEXPR.HPP */
/* Copyright Dave Curtis 1994 */
/* dcurtis@hgmp.mrc.ac.uk */
/* no warranty or liability of any kind is accepted, expressed or implied */

#ifndef DCEXPRESSHPP
#define DCEXPRESSHPP 1

#include "dcerror.hpp"
#include <math.h>
#include <stdio.h>

#define TOKENMAXLEN 100

class dcexpr_val {
public:
virtual int is_string_really()=0;
virtual operator char*()=0;
virtual operator double()=0;
// virtual ~dcexpr_val()=0;
// above doesn't seem to work
virtual ~dcexpr_val();
};

class dcexpr_double:public dcexpr_val {
double val;
char buff[50];
public:
virtual int is_string_really();
virtual operator char*();
virtual operator double();
~dcexpr_double();
dcexpr_double(double v=0.0);
};

class dcexpr_string:public dcexpr_val {
char *buff;
public:
virtual int is_string_really();
virtual operator char*();
virtual operator double();
~dcexpr_string();
dcexpr_string(char *t,int len=0);
};

class vnode {
char *str;
public:
int nbranches;
vnode *branch[4];
vnode(int nb);
virtual ~vnode(); // destroy correct object
vnode& operator=(vnode &old);
int add(vnode *vn,int b);
void wipe();
virtual dcexpr_val *eval()=0;
virtual vnode *copy()=0;
int matches(const char *s);
};

extern char *vprimitive(char*,vnode **);
extern char *vun_op(char*,vnode **);
extern char *vbin_op(char*,vnode **,int);
extern char *vbracket(char*,vnode **);
extern int vcompute(char,vnode **,vnode **);
extern int un_vcompute(char,vnode **);

class express {
protected:
vnode *head;
char token[TOKENMAXLEN ];
char *vbin_op(char *s,vnode **br,int level);
char *vun_op(char*s,vnode **br);
char *vbracket(char*s,vnode **br);
virtual char *vprimitive(char *s,vnode **br);
int vcompute(char o,vnode **i,vnode **j);
int un_vcompute(char op,vnode **i);
virtual char *get_next(char *s);
public:
void wipe();
express();
virtual ~express(); // only virtual so no warning
dcexpr_val *eval();
int parse(char *s);
};

class variable : public vnode {
public:
variable();
virtual ~variable();
dcexpr_val *eval()=0;
};

class vconstant : public vnode {
double value;
public:
vconstant(double v):vnode(0){ value=v;}
vnode *copy();
~vconstant();
dcexpr_val *eval();
};

class vstrconstant : public vnode {
char *value;
public:
vstrconstant(char *v);
vnode *copy();
~vstrconstant();
dcexpr_val *eval();
};

struct new_op_t { char *str; vnode *inst; };

extern char dc_expr_buff[];

#define MAX_N_OPS 30
#define MAX_OP_LEVEL 15
void add_un_op(char *lab,dcexpr_val *(*f)(vnode *)) ;
void add_bin_op_same(char *lab,dcexpr_val *(*f)(vnode *,vnode *));
void add_bin_op_next(char *lab,dcexpr_val *(*f)(vnode *,vnode *));
void add_bin_op(char *lab,dcexpr_val *(*f)(vnode *,vnode *),int level);

#endif

