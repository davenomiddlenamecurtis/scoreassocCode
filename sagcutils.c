/* gcutils.c */

#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <limits.h>
#include <math.h>
#include <ctype.h>

#include "sagcutils.h"
#include "cdflib.h"

#ifndef MSDOS
#include <unistd.h>
#define NUL "NUL"
#else 
#define NUL "/dev/null"
#define _CRTIMP
#endif

int ISTEMP;

#ifdef USEDCASSERT
// catch assertion failures so Windows does not do something clever with them
#ifdef  __cplusplus
extern "C" {
#endif

_CRTIMP void __cdecl _assert(void *str, void *fn, unsigned line)
{
fprintf(stderr,"Assertion failed: %s, file %s, line %d\n",(char *)str,(char *)fn,line);
exit(1);
}

#ifdef  __cplusplus
}
#endif

#endif

int error_func(unsigned l,char *f,char *s1,char *s2)
{
fprintf(stderr,"Error on line %d of source file %s:\n %s%s\n\n",l,f,s1,s2);
exit(1);
return 0;
}

double chistat(double x,double df)
{
double p,q,d1=1.0,bound;
int status,which=1;
if (x==0.0) return 1.0;
// if (x<0.0 && x>-1.0) return 1.0; 
if (x<0.0) return 1.0; 
/* do not worry about negative lrt values */
cdfchi(&which,&p,&q,&x,&df,&status,&bound);
if (status!=0)
	error("cdfchi failed","");
return q;
}

double tstat(double t,double df)
{
	int which,status;
	double p,q,f,dfn,dfd,bound;
	f=t*t;
	which=1;
	dfn=1;
	dfd=df;
	cdff(&which,&p,&q,&f,&dfn,&dfd,&status,&bound);
	return q;
}

float get_perm_lrt(subject **s,int nsub,par_info *pi)
{
float lrt;
FILE *fo;
gc_res *case_res,*cont_res,*all_res;
assert ((case_res=(gc_res *)malloc(sizeof(gc_res)))!=0);
assert ((cont_res=(gc_res *)malloc(sizeof(gc_res)))!=0);
assert ((all_res=(gc_res *)malloc(sizeof(gc_res)))!=0);
run_gc("permall",s,nsub,2,0,pi);
assert((fo=fopen("permall.gco","r"))!=0);
get_likes(all_res,fo);
fclose(fo);
run_gc("permcase",s,nsub,1,0,pi);
assert((fo=fopen("permcase.gco","r"))!=0);
get_likes(case_res,fo);
fclose(fo);
run_gc("permcont",s,nsub,0,0,pi);
assert((fo=fopen("permcont.gco","r"))!=0);
get_likes(cont_res,fo);
fclose(fo);
lrt=-2*(all_res->ll_h1-case_res->ll_h1-cont_res->ll_h1);
free(case_res);
free(cont_res);
free(all_res);
return lrt;
}

int get_haps(gc_res *res,FILE *fp,int nloci,int maxhap)
{
char line[241],rest[241];
int i,j;
hap_res *rhap;
assert((rhap=(hap_res*)malloc(sizeof(hap_res)))!=0);
fseek(fp,0L,SEEK_SET);
do
  {
  if (fgets(line,240,fp)==0)
    { error("Could not find number of haplotypes line in gene-counting output file",""); return 0; }
  } while (strncmp("Number of nonzero haplotypes",line,strlen("Number of nonzero haplotypes")));
if (sscanf(line,"Number of nonzero haplotypes = %d",&res->nhap)!=1)
    { error("Could not read number of haplotypes in gene-counting output file",""); return 0; }
fseek(fp,0L,SEEK_SET);
for (i=0;i<maxhap;++i)
  {
  res->hap[i].h0_freq=res->hap[i].h1_freq=0;
  res->hap[i].id=i+1;
  }
do
  {
  if (fgets(line,240,fp)==0)
    { error("Could not find haplotype frequency lines in gene-counting output file",""); return 0; }
  } while (strncmp("       h1",line,strlen("       h1")));
fgets(line,240,fp);
for (i=0;i<res->nhap;++i)
  {
  if (fgets(line,240,fp)==0)
    { error("Could not find all haplotype frequency lines in gene-counting output file",""); return 0; }
  if (sscanf(line,"%*f %f %f %*f %[^\n]",
        &rhap->h1_freq,&rhap->h0_freq,rest)!=3)
    { error("Could not read from haplotype frequency line in gene-counting output file: ",line); return 0; }
  for (j=0;j<nloci;++j)
    {
    strcpy(line,rest);
    *rest='\0';
    if (sscanf(line,"%d %[^\n]",
        &rhap->all[j],rest)!=2)
      { error("Could not read alleles in haplotype in gene-counting output file",""); return 0; }
    }
  if (sscanf(rest,"%d",&rhap->id)!=1)
    { error("Could not read id of haplotype in gene-counting output file",""); return 0; }
  res->hap[rhap->id-1]=*rhap;
  }
free(rhap);
return 0;
}

int get_likes(gc_res *res,FILE *fp)
{
char line[241];
int full,partial;
do
  {
  if (fgets(line,240,fp)==0)
    { error("Could not find line with number of subjects in gene-counting output file",""); return 0; }
  } while (sscanf(line,"%d/%d",&full,&partial)!=2);
res->nsub=full+partial;
do
  {
  if (fgets(line,240,fp)==0)
    { error("Could not find null log-likelihood line in gene-counting output file",""); return 0; }
  } while (strncmp("log-likelihood",line,strlen("log-likelihood")));
if (sscanf(line,"log-likelihood assuming linkage equilibrium %f",&res->ll_h0)!=1)
    { error("Could not read H0 log-likelihood in gene-counting output file",""); return 0; }
do
  {
  if (fgets(line,240,fp)==0)
    { error("Could not find H1 log-likelihood line in gene-counting output file",""); return 0; }
  } while (strncmp("log-likelihood",line,strlen("log-likelihood")));
if (sscanf(line,"log-likelihood assuming linkage disequilibrium =%f",&res->ll_h1)!=1)
    { error("Could not read H1 log-likelihood in gene-counting output file",""); return 0; }
do
  {
  if (fgets(line,240,fp)==0)
    { error("Could not find number of haplotypes line in gene-counting output file",""); return 0; }
  } while (strncmp("Number",line,strlen("Number")));
if (sscanf(line,"Number of nonzero haplotypes = %d",&res->nhap)!=1)
    { error("Could not read number of haplotypes in gene-counting output file",""); return 0; }
return 1;
}

int read_subject(char *line,subject *s,par_info *pi)
{
char *ptr;
int i,found_error,n_to_skip;
found_error=0;
switch (pi->use_cc)
  {
  case 3:
  if (sscanf(line,"%s %d %d",s->id,&s->cc,&s->group)!=3)
    { error("Syntax error in subject line: ",line); return 0; }
  n_to_skip=3;
  break;
  case 2:
  if (sscanf(line,"%s %d",s->id,&s->group)!=2)
    { error("Syntax error in subject line: ",line); return 0; }
  n_to_skip=2;
  break;
  case 1:
  if (sscanf(line,"%s %d",s->id,&s->cc)!=2)
    { error("Syntax error in subject line: ",line); return 0; }
  n_to_skip=2;
  break;
  case 0:
  if (sscanf(line,"%s",s->id)!=1)
    { error("Syntax error in subject line: ",line); return 0; }
  n_to_skip=1;
  break;
  }
for (ptr=line,i=0;i<n_to_skip;++i)
{
	while (isspace(*ptr++)) ;
	while (!isspace(*ptr++)) ;
}
for (i=0;i<pi->nloci;++i)
  {
  if (sscanf(ptr,"%d %d",&s->all[i][0],&s->all[i][1])<2)
    { 
    printf("Not enough alleles in subject line:\n%s\n",line);
    found_error=1;
    break;
    }
  n_to_skip=2;
  if ((s->all[i][0]==0)!=(s->all[i][1]==0))
    {
    printf("In subject %s there is only one zero allele for locus %d\n",
      s->id,i+1);
    found_error=1;
    }
  if (s->all[i][0]<0 || s->all[i][1]<0 || s->all[i][0]>pi->n_alls[i] || s->all[i][1]>pi->n_alls[i])
    {
    printf("In subject %s bad allele number for locus %d\n",
      s->id,i+1);
    found_error=1;
    }
	while (isspace(*ptr++)) ;
	while (!isspace(*ptr++)) ;
	while (isspace(*ptr++)) ;
	while (!isspace(*ptr++)) ;
  }
return (!found_error);
}

int read_all_subjects(FILE *fi,subject **s,int *nsub,par_info *pi)
{
char id[MAX_ID_LENGTH+1];
int found_error;
found_error=0;
*nsub=0;
while (fgets(long_line,LONG_LINE_LENGTH,fi) && sscanf(long_line,"%s",id)==1)
  {
  if (*nsub==MAX_SUB)
    { error("Number of subjects exceeds MAX_SUB",""); return 0; }
  if (!read_subject(long_line,s[(*nsub)++],pi)) found_error=1;
  }
if (found_error)
  error("Error[s] found in data file","");
return !found_error;
}

int alls_to_genotype(int *all,int n_alls)
{
int geno;
if (all[0]>all[1]) INT_SWAP(all[0],all[1]);
/* I know this is inefficient but I do not mind */
if (all[0]<=0 && all[1]<=0)
  geno=n_alls*n_alls+1;
else if (all[1]>0 && all[0]<=0)
  { error("Locus with only one allele missing",""); return 0; }
else if (all[1]>n_alls)
  { error("Locus with allele number greater than number of alleles in ",""); return 0; }
else 
  geno=(all[0]-1)*n_alls+all[1];
return geno-1;
}

int get_max_geno(int n_loci_to_use,int *loci_to_use,int *n_alls)
{
long geno;
int i;
for (geno=1L,i=0;i<n_loci_to_use;++i)
  {
  geno*=n_alls[loci_to_use[i]]*n_alls[loci_to_use[i]]+1;
  }
if (geno>INT_MAX)
  { error("Maximum possible genotypes exceeds integer limit on this machine",""); return 0;}
return (int)geno;
}

long get_max_long_geno(int n_loci_to_use,int *loci_to_use,int *n_alls)
{
long geno;
int i;
float old_geno;
for (geno=1L,i=0;i<n_loci_to_use;++i)
  {
  old_geno=geno;
  geno*=n_alls[loci_to_use[i]]*n_alls[loci_to_use[i]]+1;
  if (geno<old_geno)
    { error("Maximum possible genotypes exceeds long integer limit on this machine",""); return 0L;}
/* hope that if overflow occurs geno will go negative or at least get smaller */
  /* seems this may not work with 14 biallelic loci */
  /*
  		geno	1808548329	long
		old_geno	1.2207031e+009	float
		but geno should be 5 * old_geno


  */
  }
return geno;
}

long get_genotype(subject *s,int n_loci_to_use,int *loci_to_use,int *n_alls)
{
int i;
long g;
g=0;
for (i=0;i<n_loci_to_use;++i)
  {
  g+=alls_to_genotype(s->all[loci_to_use[i]],n_alls[loci_to_use[i]]);
  if (i!=n_loci_to_use-1)
    g*=n_alls[loci_to_use[i+1]]*n_alls[loci_to_use[i+1]]+1;
  }
return g;
}

int comp_geno(const void *s1,const void *s2)
{
subject **sub1,**sub2;
sub1=(subject**)s1;
sub2=(subject**)s2;
return (*sub1)->geno-(*sub2)->geno;
}

int fill_genotypes(subject **s,int nsub,par_info *pi)
{
int i;
long g;
for (i=0;i<nsub;++i)
  {
  g=get_genotype(s[i],pi->n_loci_to_use,pi->loci_to_use,pi->n_alls);
  s[i]->geno=g;
  }
qsort(s,nsub,sizeof(subject*),comp_geno);
return 0;
}

void output_alls(FILE *fo,long g,int n_loci_to_use,int *loci_to_use,int *n_alls)
{
int all[MAX_LOCI][2],this_geno,i;
for (i=n_loci_to_use-1;i>=0;--i)
  {
  this_geno=g%(n_alls[loci_to_use[i]]*n_alls[loci_to_use[i]]+1)+1;
  if (this_geno==(n_alls[loci_to_use[i]]*n_alls[loci_to_use[i]]+1))
    {
    all[i][0]=all[i][1]=0;
    g-=n_alls[loci_to_use[i]]*n_alls[loci_to_use[i]];
    }
  else
    {
    all[i][1]=(this_geno-1)%n_alls[loci_to_use[i]]+1;
    all[i][0]=(this_geno-1)/n_alls[loci_to_use[i]]+1;
    }
  g/=n_alls[loci_to_use[i]]*n_alls[loci_to_use[i]]+1;
  }
for (i=0;i<n_loci_to_use;++i)
  {
  fprintf(fo,"  %2d %2d",all[i][0],all[i][1]);
  if ((all[i][0]==0)!= (all[i][1]==0))
     error("Failure in output_alls()","");
  }

}

int output_genotypes(FILE *fo,subject **s,int nsub,int cc,int group,par_info *pi)
{
int i,j,gcount,which_g,found_one;
long last_geno,this_geno;
this_geno=last_geno=-1L;
found_one=0;
gcount=1;
which_g=1;
for (i=0;i<nsub;++i)
  {
  s[i]->skip=0;
  if (pi->skip_missing)
    {
    for (j=0;j<pi->n_loci_to_use;++j)
      if (s[i]->all[pi->loci_to_use[j]][0]==0)
        {
        s[i]->skip=1;
        break;
        }
    }
  if (s[i]->skip==0  && (cc==2 || s[i]->cc==cc || (cc==-1 && group>0 && s[i]->group==group)  || (cc==-1 && group<0 && s[i]->group!=-group) ))
  /* cc==-1 is signal to look at group status rather than cc status */
    {
    this_geno=s[i]->geno;
    if (this_geno==last_geno) ++gcount;
    found_one=1;
    }
  if (this_geno!=last_geno && last_geno!=-1)
    {
    fprintf(fo,"%6d  %6d",which_g++,gcount);
    output_alls(fo,last_geno,pi->n_loci_to_use,pi->loci_to_use,pi->n_alls);
    fprintf(fo,"\n");
    gcount=1;
    }
  if (i==nsub-1 && this_geno!=-1)
    {
    fprintf(fo,"%6d  %6d",which_g,gcount);
    output_alls(fo,this_geno,pi->n_loci_to_use,pi->loci_to_use,pi->n_alls);
    fprintf(fo,"\n");
    gcount=1;
    }
  if (s[i]->skip==0 && (cc==2 || s[i]->cc==cc || (cc==-1 && group>0 && s[i]->group==group)  || (cc==-1 && group<0 && s[i]->group!=-group) ))
    {
	last_geno=this_geno;
	s[i]->gc_geno=which_g;
	}
  }
return found_one;
}

char long_line[LONG_LINE_LENGTH+1];
int read_par(FILE *fp,par_info *pi)
{
char *ptr;
int i,use;
pi->skip_missing=0;
if (!fgets(long_line,LONG_LINE_LENGTH,fp)
       || sscanf(long_line,"%d %d %d",&pi->nloci,&pi->use_cc,&pi->skip_missing)<2)
  { error("Could not read first long_line of parameter file",""); return 0; }
if (pi->use_cc<0||pi->use_cc>3)
  { error("Bad value for whether to use case-control and/or group label in first long_line of parameter file:\n",long_line); return 0; }
if (pi->nloci<1)
  { error("Less than one marker in parameter file: ",long_line); return 0; }
if (pi->nloci>MAX_LOCI)
  { error("Number of loci exceeds MAX_LOCI - need to change MAX_LOCI and recompile:\n",long_line); return 0; }
if (!fgets(long_line,LONG_LINE_LENGTH,fp))
  { error("Could not read second line of parameter file",""); return 0; }
for (i=0,ptr=long_line;i<pi->nloci;++i)
  {
  if (sscanf(ptr,"%d",&pi->n_alls[i])<1)
    { error("Not enough allele numbers in parameter file",""); return 0; }
  while(isspace(*ptr++)) ;
  while(!isspace(*ptr++)) ;
  }
if (!fgets(long_line,LONG_LINE_LENGTH,fp))
  { error("Could not read third line of parameter file",""); return 0; }
for (pi->n_loci_to_use=0,i=0,ptr=long_line;i<pi->nloci;++i)
  {
  if (sscanf(ptr,"%d",&use)<1)
    { error("Not enough alleles to use in parameter file",""); return 0; }
  while(isspace(*ptr++)) ;
  while(!isspace(*ptr++)) ;
  if (use==1) pi->loci_to_use[(pi->n_loci_to_use)++]=i;
  }
return 1;
}

int read_perm_par(FILE *fp,par_info *pi)
{
char line[2000];
pi->do_perms=0;
pi->target_for_perms=-1;
if (fgets(line,1999,fp) && 
     sscanf(line,"%d %ld %ld",&pi->do_perms,&pi->n_perms,&pi->target_for_perms)>=1 &&
     pi->do_perms!=0)
  {
  if (pi->n_perms<0)
    { error("Bad number of permutations in parameter file",""); return 0; }
  if (pi->target_for_perms==-1)
    pi->target_for_perms=pi->n_perms; /* use this as default, behaviour should be OK */
  if (pi->target_for_perms<0 || pi->target_for_perms>pi->n_perms)
    { error("Bad number for target for number of permutations in parameter file",""); return 0; }
  pi->permseed=43;
  if (fgets(line,1999,fp)!=0)
    sscanf(line,"%d",&pi->permseed);
  srand(pi->permseed);
  }
return 1;
}

int run_gc(char *root,subject **s,int nsub,int cc,int group,par_info *pi)
{
char line[200];
FILE *fi;
int i;
sprintf(line,"%s.gci",root);
fi=fopen(line,"w");
if (fi==0)
  { error("Could not open .gci file for writing ",root); return 0; }
for (i=0;i<pi->n_loci_to_use;++i)
  fprintf(fi,"%3d ",pi->n_alls[pi->loci_to_use[i]]);
fprintf(fi,"\n");
output_genotypes(fi,s,nsub,cc,group,pi);
fclose(fi);
sprintf(line,"%s.gco",root);
unlink(line);
/*
I could insert here:
sprintf(line,"gc %s.gci %s.gco > gc.log",root,root);
but it might slow permutations to write to a real file
*/
#if 0
ftest=fopen(line,"r");
if (ftest!=0)
  {
  fclose(ftest);
  ftest=fopen(line,"w");
  fclose(ftest);
  }
#endif
sprintf(line,"gc %s.gci %s.gco > %s",root,root,NUL);
system(line);
return 1;
}

int get_geno_probs(geno_probs *gp,FILE *fi,par_info *pi)
{
char line[2][241],rest[241];
int thisg,lastg,whichpair,i,j,wrote_message;
do {
   if (fgets(line[0],240,fi)==0)
     { error("Could not find haplotype assignments in GC output file",""); return 0; }
   } while (strncmp(line[0],"Assignment",strlen("Assignment")));
assert(fgets(line[0],240,fi));
assert(fgets(line[0],240,fi));
lastg=-1;
while (fgets(line[0],240,fi) && sscanf(line[0],"%d [%*d] %[^\n]",&thisg,rest)==2)
  {
  if (thisg!=lastg)
    {
	wrote_message=0;
	whichpair=0;
	lastg=thisg;
	}
  else
    {
	++whichpair;
	if (whichpair>=MAXHAPSPERGENO)
	  {
	  if (!wrote_message)
	    {
		printf("\n\nNote: genotype %d has more than MAXHAPSPERGENO=%d combinations\n - ignoring excess\n\n",
	      thisg,MAXHAPSPERGENO);
		wrote_message=1;
		}
	  fgets(line[1],240,fi);
	  continue;
	  }
	}
  gp[thisg-1].npairs=whichpair+1;
  assert(fgets(line[1],240,fi));
  for (i=0;i<2;++i)
    {
	sscanf(line[i],"%*d [%*d] %[^\n]",rest);
	for (j=0;j<pi->n_loci_to_use;++j)
	  {
	  strcpy(line[i],rest);
	  assert(sscanf(line[i],"%*d %[^\n]",rest)==1);
	  }
    sscanf(rest,"%f %d",&gp[thisg-1].hap_pair_prob[whichpair],&gp[thisg-1].hap[whichpair][i]);
	}
  }
return 1;
}

void perm_subs(subject **s,int nsub)
{
int i,k;
for (i=nsub-1;i>=1;--i)
  {
  k=rand()%(i+1);
  INT_SWAP(s[i]->cc,s[k]->cc);
  }
}

void get_means_ll(float *means_ll,int *found_hap,gc_res *cont_res,gc_res *case_res,gc_res *all_res,int maxhap)
{
gc_res *mean_sample;
int i;
float sigma_p_logp_case,sigma_p_logp_cont,sigma_p_logp_all;
assert ((mean_sample=(gc_res *)malloc(sizeof(gc_res)))!=0);
assert(mean_sample->hap=(hap_res*)calloc(maxhap,sizeof(hap_res)));
sigma_p_logp_case=sigma_p_logp_cont=sigma_p_logp_all=0;
for (i=0;i<maxhap;++i)
  {
  if (cont_res->hap[i].h1_freq!=0)
    sigma_p_logp_cont+=cont_res->hap[i].h1_freq*log(cont_res->hap[i].h1_freq);
  if (case_res->hap[i].h1_freq!=0)
    sigma_p_logp_case+=case_res->hap[i].h1_freq*log(case_res->hap[i].h1_freq);
  if (all_res->hap[i].h1_freq!=0)
    sigma_p_logp_all+=all_res->hap[i].h1_freq*log(all_res->hap[i].h1_freq);
  }
cont_res->nprime=cont_res->ll_h1/sigma_p_logp_cont;
case_res->nprime=case_res->ll_h1/sigma_p_logp_case;
all_res->nprime=all_res->ll_h1/sigma_p_logp_all;
*means_ll=0;
*found_hap=0;
for (i=0;i<maxhap;++i)
  {
  mean_sample->hap[i].h1_freq=
    (cont_res->hap[i].h1_freq*cont_res->nprime+
     case_res->hap[i].h1_freq*case_res->nprime)/
     (cont_res->nprime+case_res->nprime);
  if (mean_sample->hap[i].h1_freq!=0)
    {
	  *means_ll+=mean_sample->hap[i].h1_freq*log(mean_sample->hap[i].h1_freq)*(cont_res->nprime+case_res->nprime);
    ++*found_hap;
    }
  else if (all_res->hap[i].h1_freq!=0)
    ++*found_hap;
    /* this is when a haplotype is found in whole sample but not in cases or controls */
  }
free(mean_sample->hap);
free(mean_sample);
}

int read_window_par(FILE *fp,par_info *pi)
{
char line[2000];
pi->wt=SLIDE;
if (fgets(line,1999,fp)==0 || sscanf(line,"%d %d",&pi->n_window,&pi->wt)<1)
    { error("Bad window definition in parameter file",""); return 0; }
return 1;
}

void fill_res(char *root,gc_res *res,subject **s,int nsub,int cc,int gr,par_info *pi,int maxhap)
{
char fname[100];
FILE *fo;
run_gc(root,s,nsub,cc,gr,pi);
sprintf(fname,"%s.gco",root);
assert((fo=fopen(fname,"r"))!=0);
get_likes(res,fo);
get_haps(res,fo,pi->n_loci_to_use,maxhap);
fclose(fo);
}



