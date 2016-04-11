/* sagcutils.h */
#ifndef SAGCUTILSH
#define SAGCUTILSH

#include <stdio.h>

#ifndef VERSION
#define VERSION "2.5"
#endif

#ifndef MAX_LOCI
#define MAX_LOCI 12000
#endif
#ifndef MAX_ALL
#define MAX_ALL 2
#endif
#ifndef MAX_SUB
#define MAX_SUB 10000
#endif
#ifndef MAX_LINE
#define MAX_LINE 5000
#endif
#ifndef MAX_ID_LENGTH 
#define MAX_ID_LENGTH 100
#endif


extern int ISTEMP;
#define INT_SWAP(x,y) (ISTEMP=x,x=y,y=ISTEMP)


#define error(s,t) error_func(__LINE__,__FILE__,s,t)

enum WINDOW_TYPE {SLIDE=1,COMB};
enum TEST_TYPE {TEST_CC=1,TEST_GROUP};


struct par_info_t {
int loci_to_use[MAX_LOCI],n_alls[MAX_LOCI],n_loci_to_use,nloci,use_cc,skip_missing;
int do_perms,permseed;
long n_perms,target_for_perms;
int n_window;
int which_gr,last_gr,*ngr;
enum WINDOW_TYPE wt;
enum TEST_TYPE tt;
};

typedef struct par_info_t par_info;

struct hap_res_t {
float h1_freq,h0_freq;
int all[MAX_LOCI];
int id;
};

typedef struct hap_res_t hap_res;

struct gc_res_t {
float a_freq[MAX_LOCI][MAX_ALL];
float ll_h0,ll_h1;
hap_res *hap;
int nhap,nsub;
float nprime;
} ;

typedef struct gc_res_t gc_res;

#define MAXHAPSPERGENO 200 /* 32 */
#define MAXLINELENGTH 50000

struct geno_probs_t {
int hap[MAXHAPSPERGENO][2],npairs;
float hap_pair_prob[MAXHAPSPERGENO];
};

typedef struct geno_probs_t geno_probs;

struct subject_t { 
	char id[35]; int cc, group, skip; union { int all[MAX_LOCI][2]; float prob[MAX_LOCI][3]; }; long geno; int gc_geno;
};

typedef struct subject_t subject;

int error_func(unsigned l,char *f,char *s1,char *s2);
double chistat(double x,double df);
double tstat(double t,double df);
int read_subject(char *line,subject *s,par_info *pi);
int get_haps(gc_res *res,FILE *fp,int nloci,int maxhap);
int get_likes(gc_res *res,FILE *fp);
int read_all_subjects(FILE *fi,subject **s,int *nsub,par_info *pi);
int alls_to_genotype(int *all,int n_alls);
long get_genotype(subject *s,int n_loci_to_use,int *loci_to_use,int *n_alls);
int fill_genotypes(subject **s,int nsub,par_info *pi);
void output_alls(FILE *fo,long g,int n_loci_to_use,int *loci_to_use,int *n_alls);
int output_genotypes(FILE *fo,subject **s,int nsub,int cc,int group,par_info *pi);
int read_par(FILE *fp,par_info *pi);
int read_perm_par(FILE *fp,par_info *pi);
int run_gc(char *root,subject **s,int nsub,int cc,int group,par_info *pi);
float get_perm_lrt(subject **s,int nsub,par_info *pi);
int get_max_geno(int n_loci_to_use,int *loci_to_use,int *n_alls);
long get_max_long_geno(int n_loci_to_use,int *loci_to_use,int *n_alls);
int comp_geno(const void *s1,const void *s2);
int get_geno_probs(geno_probs *gp,FILE *fi,par_info *pi);
void perm_subs(subject **s,int nsub);
void get_means_ll(float *means_ll,int *found_hap,gc_res *cont_res,gc_res *case_res,gc_res *all_res,int maxhap);
void fill_res(char *root,gc_res *res,subject **s,int nsub,int cc,int gr,par_info *pi,int maxhap);

#define LONG_LINE_LENGTH MAX_LOCI*100
extern char long_line[]; /* for thousands of markers */
#endif





