#include <ctype.h>
#include "scoreassoc.hpp"

void write_scores(FILE *fs,subject **sub,int nsub,float *score)
{
	int s;
	for (s=0;s<nsub;++s)
		fprintf(fs,"%20s %d %8.4f\n",sub[s]->id,sub[s]->cc,score[s]);
}

float do_score_onetailed_ttest(FILE *fo,float *score,subject **sub,int nsub,par_info *pi,sa_par_info *spi,float cc_freq[2][MAX_LOCI],float cc_count[2][MAX_LOCI],int max_cc[2],float *weight,float *missing,int *rarer)
{
	int s,i,n[2],cc,l;
	float sigma_x[2],sigma_x2[2],mean[2],var[2],tval,SE,s2,rfreq,fscore;
	double p;
	for (i=0;i<2;++i)
		sigma_x[i]=sigma_x2[i]=n[i]=0;
	for (s=0;s<nsub;++s)
    {
		 cc=sub[s]->cc;
		 ++n[cc];
		 sigma_x[cc]+=score[s];
		 sigma_x2[cc]+=score[s]*score[s];
    }
    for (i=0;i<2;++i)
    {
		if (spi->use_cc_freqs[i])
		{
			float total_score=0;
			n[i]=max_cc[i];
			var[i]=0;
			for (l=0;l<pi->n_loci_to_use;++l)
			{
				float tempvar=0;
				float ngen[3]; // number of typed subjects who are AA,AB,BB
				sigma_x2[i]=0;
				sigma_x[i]=0; 
				if (rarer[l]==2)
					rfreq=cc_freq[i][pi->loci_to_use[l]]; 
				// assume supplied frequency is of alt allele, i.e. allele 2
				// but score will be added to by rarer allele, even if this is allele 1
				else
					rfreq=1-cc_freq[i][pi->loci_to_use[l]];
				// as actual counts missing, assume HWE
				ngen[0]=(1-rfreq)*(1-rfreq)*cc_count[i][pi->loci_to_use[l]];
				ngen[1]=2*rfreq*(1-rfreq)*cc_count[i][pi->loci_to_use[l]];
				ngen[2]=rfreq*rfreq*cc_count[i][pi->loci_to_use[l]];
				fscore=weight[l]; // average score for each subject for this locus
				sigma_x[i]+=ngen[1]*fscore; // AB
				sigma_x2[i]+=ngen[1]*fscore*fscore;
				sigma_x[i]+=ngen[2]*2*fscore; // BB
				sigma_x2[i]+=ngen[2]*4*fscore*fscore;
				fscore=missing[l]; // as if missing scores were also independent
				sigma_x[i]+=(max_cc[i]-cc_count[i][pi->loci_to_use[l]])*fscore;
				sigma_x2[i]+=(max_cc[i]-cc_count[i][pi->loci_to_use[l]])*fscore*fscore;
				tempvar=(sigma_x2[i]-sigma_x[i]*sigma_x[i]/n[i])/(n[i]-1);
				var[i]+=tempvar; // add variances due to each marker
				total_score+=sigma_x[i];
			}
			mean[i]=total_score/n[i];
		}
		else
		{
			var[i]=(sigma_x2[i]-sigma_x[i]*sigma_x[i]/n[i])/(n[i]-1);
			mean[i]=sigma_x[i]/n[i];
		}
    }
	s2=((n[0]-1)*var[0]+(n[1]-1)*var[1])/(n[0]+n[1]-2);
    SE=sqrt(s2*(1/(float)n[0]+1/(float)n[1]));
	if (SE==0)
		tval=0;
	else
		tval=(mean[1]-mean[0])/SE;
    p=tstat(tval,n[0]+n[1]-2.0)/2;
	if (mean[0]>mean[1])
		p=1.0-p;
	if (fo!=NULL)
		fprintf(fo,"             Controls  Cases     \n"
		           "N            %9d %9d\n"
				   "Mean score   %9.3f %9.3f\n"
				   "t (%d df) = %6.3f\n"
				   "p = %10.8f\n"
				   "-log(p) = %8.2f\n",
				   n[0],n[1],mean[0],mean[1],n[0]+n[1]-2,tval,p,-log10(p));
	return p;
}

/* treats male subjects as homozygous females for X loci*/
void get_scores(float *score,float *weight,float *missing,int *rarer,subject **sub,int nsub,par_info *pi,sa_par_info *spi)
{
	int l,ll,s,a,p;
	for (s=0;s<nsub;++s)
	{
		score[s]=0;
        for (l=0;l<pi->n_loci_to_use;++l)
			{
				ll=pi->loci_to_use[l];
				if (spi->use_probs)
				{
					if (sub[s]->prob[ll][0]+sub[s]->prob[ll][1]+sub[s]->prob[ll][2]==0)
						score[s] += missing[l]*2; 
					else if (rarer[l]==2)
						score[s]+=sub[s]->prob[ll][1]+sub[s]->prob[ll][2]*2;
					else
						score[s]+=sub[s]->prob[ll][1]+sub[s]->prob[ll][0]*2;
				}
				else
				{
					for (a = 0; a < 2; ++a)
						if (sub[s]->all[ll][a] == rarer[l])
							score[s] += weight[l];
						else if (sub[s]->all[ll][a] == 0)
							score[s] += missing[l];
				}
			}
	}
}

// treats male subjects as homozygous females for X loci
// this function is here to allow us to exclude loci with wildly different allele frequencies, etc.
void get_freqs(subject **sub,int nsub,par_info *pi,sa_par_info *spi,float cc_freq[2][MAX_LOCI],float cc_count[2][MAX_LOCI],float cc_genocount[2][3][MAX_LOCI])
{
	int l,s,nh[2],cc,i,g;
	float ccfreq[2],gencount[2][3],vcount[2];
	for (l=0;l<pi->n_loci_to_use;++l)
	{
		nh[0]=nh[1]=vcount[0]=vcount[1]=0;
		for (i=0;i<3;++i)
			for (cc=0;cc<2;++cc)
				gencount[cc][i]=0;
		if (spi->use_probs)
			for (s = 0; s < nsub; ++s)
			{
			for (g=0;g<3;++g)
				gencount[sub[s]->cc][g]+=sub[s]->prob[pi->loci_to_use[l]][g];
			vcount[sub[s]->cc]+=sub[s]->prob[pi->loci_to_use[l]][1]+sub[s]->prob[pi->loci_to_use[l]][2]*2;
			// allele count
			nh[sub[s]->cc]+=(sub[s]->prob[pi->loci_to_use[l]][0]+sub[s]->prob[pi->loci_to_use[l]][1]+sub[s]->prob[pi->loci_to_use[l]][2])*2;
			// allow for possibility that unknowns could be coded as 0 0 0
			}
		else
		for (s=0;s<nsub;++s)
			if (sub[s]->all[pi->loci_to_use[l]][0]!=0)
			{
			i=(sub[s]->all[pi->loci_to_use[l]][0]==2)+(sub[s]->all[pi->loci_to_use[l]][1]==2);
			++gencount[sub[s]->cc][i];
			vcount[sub[s]->cc]+=(sub[s]->all[pi->loci_to_use[l]][0]==2)+(sub[s]->all[pi->loci_to_use[l]][1]==2);
			nh[sub[s]->cc]+=2;
			}
		for (cc=0;cc<2;++cc)
			if (!spi->use_cc_freqs[cc])
			{
				if (nh[cc]==0)
					ccfreq[cc]=0;
				else
					ccfreq[cc]=vcount[cc]/nh[cc];
				cc_freq[cc][pi->loci_to_use[l]]=ccfreq[cc];
				// this line is here because I will use this for filtering bad loci
				cc_count[cc][pi->loci_to_use[l]]=nh[cc]/2;
				for (g=0;g<3;++g)
					cc_genocount[cc][g][pi->loci_to_use[l]]=gencount[cc][g];
			}
	}
}

float get_quadratic_weight(float freq,float wfactor)
{
	float wt;
	wt=(4*wfactor-4)*freq*freq-(4*wfactor-4)*freq+wfactor;
	return wt;
}

float get_quartic_weight(float freq,float wfactor)
{
	static float a,b,kept_wfactor=-1;
	float wt;
	if (wfactor!=kept_wfactor)
	{
		a=(1/pow(0.4,2)-wfactor/pow(0.5,2))/(pow(0.4,2)-pow(0.5,2));
		b=(1/pow(0.4,4)-wfactor/pow(0.5,4))/(1/pow(0.4,2)-1/pow(0.5,2));
		kept_wfactor=wfactor;
	}
	wt=a*pow(freq-0.5,4)+b*pow(freq-0.5,2);
	return wt;
}

float get_zero_based_quadratic_weight(float freq,float wfactor)
/* wfactor is gradient when freq==0 or freq==1 - a=wfactor */
{
	float wt;
	wt=wfactor*(freq-0.5)*(freq-0.5);
	return wt;
}

void set_weights(FILE *f,float *weight,float *missing_score,int *rarer,subject **sub,int nsub,par_info *pi,sa_par_info *spi,float *func_weight,float cc_freq[2][MAX_LOCI],float cc_count[2][MAX_LOCI],int max_cc[2],char names[MAX_LOCI][20],char comments[MAX_LOCI][MAX_COMMENT_LENGTH])
{
	int l,ll,s,nh[2],cc,i,g;
	float freq,ccfreq[2],vcount[2],gencount[2][3];

	for (l=0;l<pi->n_loci_to_use;++l)
	{		
		nh[0]=nh[1]=vcount[0]=vcount[1]=0;
		for (i=0;i<3;++i)
			for (cc=0;cc<2;++cc)
				gencount[cc][i]=0;
		ll=pi->loci_to_use[l];
		for (s = 0; s < nsub; ++s)
		{
			if (spi->use_probs)
			{
				for (g=0;g<3;++g)
					gencount[sub[s]->cc][g]+=sub[s]->prob[ll][g];
				vcount[sub[s]->cc] += sub[s]->prob[ll][1]+sub[s]->prob[ll][2]*2;
				nh[sub[s]->cc]+=(sub[s]->prob[ll][0]+sub[s]->prob[ll][1]+sub[s]->prob[ll][2])*2;
			}
			else
			{
				if (sub[s]->all[ll][0] != 0)
				{
					i = (sub[s]->all[ll][0] == 2) + (sub[s]->all[ll][1] == 2);
					++gencount[sub[s]->cc][i];
					vcount[sub[s]->cc] += (sub[s]->all[ll][0] == 2) + (sub[s]->all[ll][1] == 2);
					nh[sub[s]->cc] += 2;
				}
			}
		}
		for (cc=0;cc<2;++cc)
			if (spi->use_cc_freqs[cc])
			{
				ccfreq[cc]=cc_freq[cc][ll];
				nh[cc]=2*cc_count[cc][ll];
			}
			else
			{
				if (nh[cc]==0)
					ccfreq[cc]=0;
				else
					ccfreq[cc]=vcount[cc]/nh[cc];
				cc_freq[cc][ll]=ccfreq[cc];
				// this line is here because I will use this for filtering bad loci
			}
		if (nh[0]+nh[1]==0)
			freq=0.0;
		else
			freq=(ccfreq[0]*nh[0]+ccfreq[1]*nh[1])/(nh[0]+nh[1]);
		if (freq>0.5)
		{
			freq=1.0-freq;
			rarer[l]=1;
		}
		else
			rarer[l]=2;
		if (freq==0.0) /* monomorphic */
			weight[l]=0;
		else
			weight[l]=get_quadratic_weight(freq,spi->wfactor);
		if (spi->use_func_weights)
			weight[l]*=func_weight[pi->loci_to_use[l]];
		missing_score[l]=weight[l]*freq; // I think
	if (f!=0)
	{
		if (!spi->use_locus_names)
			sprintf(names[ll],"LOC%05d",ll+1);
		fprintf(f,"%-20s",names[ll]);
		fprintf(f,
			spi->use_probs ?
			"%6.2f : %6.2f : %6.2f  %8.6f  %6.2f : %6.2f : %6.2f  %8.6f  %8.6f  %d      %5.2f  %s\n" :
			"%6.0f : %6.0f : %6.0f  %8.6f  %6.0f : %6.0f : %6.0f  %8.6f  %8.6f  %d      %5.2f  %s\n",
			gencount[0][0],gencount[0][1],gencount[0][2],
			ccfreq[0],
			gencount[1][0],gencount[1][1],gencount[1][2],
			ccfreq[1],
			freq,rarer[l],weight[l],
			spi->use_comments?comments[ll]:"");

	}
	}
}

int read_sa_par(FILE *fp, par_info *pi, sa_par_info *spi, float *func_weight, float cc_freq[2][MAX_LOCI], float cc_count[2][MAX_LOCI], int max_cc[2], char names[MAX_LOCI][20], char comments[MAX_LOCI][MAX_COMMENT_LENGTH])
{
char *ptr;
int l,cc,c;
if (!fgets(long_line,LONG_LINE_LENGTH,fp) || sscanf(long_line,"%f",&spi->wfactor)<1)
  {  
  error("Could not read weight factor\n",long_line);
  return 0;
  }
spi->use_func_weights=spi->use_cc_freqs[0]=spi->use_cc_freqs[1]=spi->use_locus_names=spi->use_comments=spi->do_recessive_test=0;
spi->weight_threshold=0;
spi->LD_threshold=1;
spi->use_haplotypes=0;
if (fgets(long_line,LONG_LINE_LENGTH,fp)) // extra info is optional
{
	sscanf(long_line,"%d %d %d %d %d %d %f %f %d",
		&spi->use_func_weights,&spi->use_cc_freqs[0],&spi->use_cc_freqs[1],&spi->use_locus_names,&spi->use_comments,&spi->do_recessive_test,&spi->weight_threshold,&spi->LD_threshold,&spi->use_haplotypes);
}
if (spi->use_func_weights)
{
	if (!fgets(long_line,LONG_LINE_LENGTH,fp))
	{  
		error("Could not read functional weights\n","");
		return 0;
	}
	for (l=0,ptr=long_line;l<pi->nloci;++l)
	{
		if (sscanf(ptr,"%f",&func_weight[l])<1)
		{
			error("Could not read all functional weights\n",long_line);
			return 0;
		}
		while(isspace(*ptr++)) ;
		while(!isspace(*ptr++)) ;

	}
}
else
	for (l=0;l<pi->nloci;++l)
		func_weight[l]=1.0;
for (cc=0;cc<2;++cc) if (spi->use_cc_freqs[cc])
{
	if (!fgets(long_line,LONG_LINE_LENGTH,fp))
	{  
		error("Could not read supplied frequencies\n","");
		return 0;
	}
	for (l=0,ptr=long_line;l<pi->nloci;++l)
	{
		if (sscanf(ptr,"%f",&cc_freq[cc][l])<1)
		{
			error("Could not read all supplied frequencies\n",long_line);
			return 0;
		}
		while(isspace(*ptr++)) ;
		while(!isspace(*ptr++)) ;
	}
	if (!fgets(long_line,LONG_LINE_LENGTH,fp))
	{  
		error("Could not read subject counts\n","");
		return 0;
	}
	max_cc[cc]=0;
	for (l=0,ptr=long_line;l<pi->nloci;++l)
	{
		if (sscanf(ptr,"%f",&cc_count[cc][l])<1)
		{
			error("Could not read all subject counts\n",long_line);
			return 0;
		}
		if (max_cc[cc]<cc_count[cc][l])
			max_cc[cc]=cc_count[cc][l];
		while(isspace(*ptr++)) ;
		while(!isspace(*ptr++)) ;
	}
}
if (spi->use_locus_names)
{
	if (!fgets(long_line,LONG_LINE_LENGTH,fp))
	{  
		error("Could not read locus names\n","");
		return 0;
	}
	for (l=0,ptr=long_line;l<pi->nloci;++l)
	{
		if (sscanf(ptr,"%s",names[l])<1)
		{
			error("Could not read all locus names\n",long_line);
			return 0;
		}
		while(isspace(*ptr++)) ;
		while(!isspace(*ptr++)) ;

	}
}
if (spi->use_comments)
{
	if (!fgets(long_line,LONG_LINE_LENGTH,fp))
	{  
		error("Could not read locus comments\n","");
		return 0;
	}
	for (l=0,ptr=long_line;l<pi->nloci;++l)
	{
		while(isspace(*ptr)) 
			++ptr;
		if (*ptr=='\0')
		{
			error("Could not read all locus comments\n",long_line);
			return 0;
		}
		c=0;
		while(!isspace(*ptr))
		{
			if (c==MAX_COMMENT_LENGTH-1)
				comments[l][c]='\0';
			else
				comments[l][c++]=*ptr;
			++ptr;
		}
		comments[l][c]='\0';
	}
}
return 1;
}

int read_score_assoc_par(FILE *fp,par_info *pi,float *wfactor,int *use_func_weights,int use_cc_freqs[2],float *func_weight,float cc_freq[2][MAX_LOCI],float cc_count[2][MAX_LOCI],int max_cc[2],int *use_locus_names,char names[MAX_LOCI][20],int *use_comments,char comments[MAX_LOCI][MAX_COMMENT_LENGTH],int *do_recessive_test,float *weight_threshold,float *LD_threshold,int *use_haplotypes)
{
char *ptr;
int l,cc,c;
if (!fgets(long_line,LONG_LINE_LENGTH,fp) || sscanf(long_line,"%f",wfactor)<1)
  {  
  error("Could not read weight factor\n",long_line);
  return 0;
  }
*use_func_weights=use_cc_freqs[0]=use_cc_freqs[1]=*use_locus_names=*use_comments=*do_recessive_test=0;
*weight_threshold=0;
*LD_threshold=1;
*use_haplotypes=0;
if (fgets(long_line,LONG_LINE_LENGTH,fp)) // extra info is optional
{
	sscanf(long_line,"%d %d %d %d %d %d %f %f %d",use_func_weights,&use_cc_freqs[0],&use_cc_freqs[1],use_locus_names,use_comments,do_recessive_test,weight_threshold,LD_threshold,use_haplotypes);
}
if (*use_func_weights)
{
	if (!fgets(long_line,LONG_LINE_LENGTH,fp))
	{  
		error("Could not read starting weights\n","");
		return 0;
	}
	for (l=0,ptr=long_line;l<pi->nloci;++l)
	{
		if (sscanf(ptr,"%f",&func_weight[l])<1)
		{
			error("Could not read all starting weights\n",long_line);
			return 0;
		}
		while(isspace(*ptr++)) ;
		while(!isspace(*ptr++)) ;

	}
}
else
	for (l=0;l<pi->nloci;++l)
		func_weight[l]=1.0;
for (cc=0;cc<2;++cc) if (use_cc_freqs[cc])
{
	if (!fgets(long_line,LONG_LINE_LENGTH,fp))
	{  
		error("Could not read supplied frequencies\n","");
		return 0;
	}
	for (l=0,ptr=long_line;l<pi->nloci;++l)
	{
		if (sscanf(ptr,"%f",&cc_freq[cc][l])<1)
		{
			error("Could not read all supplied frequencies\n",long_line);
			return 0;
		}
		while(isspace(*ptr++)) ;
		while(!isspace(*ptr++)) ;
	}
	if (!fgets(long_line,LONG_LINE_LENGTH,fp))
	{  
		error("Could not read subject counts\n","");
		return 0;
	}
	max_cc[cc]=0;
	for (l=0,ptr=long_line;l<pi->nloci;++l)
	{
		if (sscanf(ptr,"%f",&cc_count[cc][l])<1)
		{
			error("Could not read all subject counts\n",long_line);
			return 0;
		}
		if (max_cc[cc]<cc_count[cc][l])
			max_cc[cc]=cc_count[cc][l];
		while(isspace(*ptr++)) ;
		while(!isspace(*ptr++)) ;
	}
}
if (*use_locus_names)
{
	if (!fgets(long_line,LONG_LINE_LENGTH,fp))
	{  
		error("Could not read locus names\n","");
		return 0;
	}
	for (l=0,ptr=long_line;l<pi->nloci;++l)
	{
		if (sscanf(ptr,"%s",names[l])<1)
		{
			error("Could not read all locus names\n",long_line);
			return 0;
		}
		while(isspace(*ptr++)) ;
		while(!isspace(*ptr++)) ;

	}
}
if (*use_comments)
{
	if (!fgets(long_line,LONG_LINE_LENGTH,fp))
	{  
		error("Could not read locus comments\n","");
		return 0;
	}
	for (l=0,ptr=long_line;l<pi->nloci;++l)
	{
		while(isspace(*ptr)) 
			++ptr;
		if (*ptr=='\0')
		{
			error("Could not read all locus comments\n",long_line);
			return 0;
		}
		c=0;
		while(!isspace(*ptr))
		{
			if (c==MAX_COMMENT_LENGTH-1)
				comments[l][c]='\0';
			else
				comments[l][c++]=*ptr;
			++ptr;
		}
		comments[l][c]='\0';
	}
}
return 1;
}
