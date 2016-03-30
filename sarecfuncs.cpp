/* sarecfuncs.cpp */

#include "scoreassoc.hpp"

#define TESTFRACTION 0.1

#define MIN(X,Y) (((X) < (Y)) ? (X) : (Y))

float occurs_with[MAX_LOCI][MAX_LOCI],occurs[MAX_LOCI],occ[2],MAF[MAX_LOCI],contMAF[MAX_LOCI];

double cumulBinom(int N,int k,double p)
// return cumulative probability of getting k successes or fewer with N attempts at probability p of successes
// return cumulative probability of getting k successes or fewer with N attempts at probability p of successes
// there is a real problem with rounding here, even using doubles
// to avoid adding big numbers to little numbers, begin by calculating the probability for the expected number of successes
// then go up to k and down to 0
{
	double logExactP,cumulP,logP,log1mP,binom,x,startLogExactP;
	int i,starti;
	if (p==0 || k==N)
		return 1.0;
	if (p==1)
		return 0.0;
	logP=log(p);
	x=exp(logP);
	log1mP=log(1.0-p);
	x=exp(log1mP);
	binom=0;
	starti=(int)(p*N+0.5); // most likely number or successes
	for (i=1;i<=starti;++i)
		binom+=log((double)N-i+1)-log((double)i); // log of binomial coeffcient for most likely number of successes
	logExactP=startLogExactP=starti*logP+(N-starti)*log1mP+binom;
	cumulP=exp(logExactP);
	for (i=starti;i<k;++i)
	{
		logExactP+=logP-log1mP+log((double)N-i)-log((double)i+1);
		cumulP+=exp(logExactP);
	}
	logExactP=startLogExactP; // now go back and work down t0 0
	for (i=starti;i>=1;--i)
	{
		logExactP+=log1mP-logP-log((double)N-i+1)+log((double)i);
		cumulP+=exp(logExactP);
	}
	return cumulP;
}

void do_recessive_HWE_test(FILE *fo,float *score,subject **sub,int nsub,par_info *pi,sa_par_info *spi,float cc_freq[2][MAX_LOCI],float cc_count[2][MAX_LOCI],int max_cc[2],float *weight,float *missing,int *old_rarer,char names[MAX_LOCI][20])
{
	float tab[2][2],ex[2][2],col_tot[3],row_tot[2],N,counts[2][3],hom_counts[2],all_count[MAX_LOCI],exp_hom_freq,homoz;
	int r,c,s,l,ll,lll,a,genotype,first,all_genotype,rec_rarer[MAX_LOCI],g;
	int n_loci_to_remove,loci_to_remove[MAX_LOCI];
	double p,df,dchi;
	float freq,tot,geno_prob[3],geno_count[MAX_LOCI][3];
	int status,which=1;
	par_info rec_pi; // use this to keep track of which loci we are using
	// first thing we will do is remove loci which fail the weight threshold
	// note that weights and rarer are indexed up to pi->n_loci_to_use
	for (l=0,ll=0;l<pi->n_loci_to_use;++l)
	{
		if (weight[l]>=spi->weight_threshold)
		{
			rec_pi.loci_to_use[ll]=pi->loci_to_use[l];
			rec_rarer[rec_pi.loci_to_use[ll]]=rarer[l]; // here is where we go back to indexing rec_rarer in the same way as other locus attributes
			if (spi->use_cc_freqs[0])
				contMAF[ll]=rarer[l]==2?cc_freq[0][pi->loci_to_use[l]]:(1-cc_freq[0][pi->loci_to_use[l]]);
			// we think we are looking at the original frequencies as they were read in 
			// but now we are only keeping those for the qualifying loci and they are stored in contMAF
			++ll;
		}
	}
	rec_pi.n_loci_to_use=ll;

	// from now on it impossible to marry up the weights with the loci listed in rec_pi
	// hope that is OK
	N=0;
	for (r=0;r<2;++r)
	{
		for (c=0;c<3;++c)
			counts[r][c]=0;
		row_tot[r]=0;
		hom_counts[r]=0;
	}
	for (l=0;l<rec_pi.n_loci_to_use;++l)
		for (g=0;g<3;++g)
			geno_count[l][g]=0;
	// begin by just listing all subjects with two or more variants at these loci before excluding on LD
	fprintf(fo,"\nRecessive analyses based on loci with weight reaching %f\n",spi->weight_threshold);
	fprintf(fo,"Using the following loci:\n");
	for (l=0;l<rec_pi.n_loci_to_use;++l)
		fprintf(fo,"%3d %s\n",l+1,names[rec_pi.loci_to_use[l]]);
	fprintf(fo,"\nSubjects with two or more minor alleles:\n");
	for (s=0;s<nsub;++s)
	{
		genotype=0;
        for (l=0;l<rec_pi.n_loci_to_use;++l)
			{
				for (a=0;a<2;++a)
					if (sub[s]->all[rec_pi.loci_to_use[l]][a]==rec_rarer[rec_pi.loci_to_use[l]])
					{
						if (genotype<2)
							++genotype;
					}
			}

		if (genotype==2)
		{
			first=1;
			fprintf(fo,"%-10s %d   ",sub[s]->id,sub[s]->cc);
			for (l=0;l<rec_pi.n_loci_to_use;++l)
				for (a=0;a<2;++a)
					if (sub[s]->all[rec_pi.loci_to_use[l]][a]==rarer[l])
					{
						fprintf(fo,"%s%d",first?"":"-",l+1);
						first=0;
					}
			fprintf(fo,"\n");
		}
	}
	// now do the exclusions based on LD
	for (l=0;l<rec_pi.n_loci_to_use;++l)
	{
		for (ll=0;ll<rec_pi.n_loci_to_use;++ll)
			if (l!=ll)
				occurs_with[l][ll]=0;
		occurs[l]=0;
	}
	for (s=0;s<nsub;++s)
	{
		for (l=0;l<rec_pi.n_loci_to_use;++l)
		{
			occ[0]=(sub[s]->all[rec_pi.loci_to_use[l]][0]==rec_rarer[rec_pi.loci_to_use[l]])+(sub[s]->all[rec_pi.loci_to_use[l]][1]==rec_rarer[rec_pi.loci_to_use[l]]);
			if (occ[0]!=0)
			{
				occurs[l]+=occ[0];
				for (ll=0;ll<rec_pi.n_loci_to_use;++ll)
					if (l!=ll)
					{
					occ[1]=(sub[s]->all[rec_pi.loci_to_use[ll]][0]==rec_rarer[rec_pi.loci_to_use[ll]])+(sub[s]->all[rec_pi.loci_to_use[ll]][1]==rec_rarer[rec_pi.loci_to_use[ll]]);
					occurs_with[l][ll]+=MIN(occ[0],occ[1]); // I think this is right - number of variants occurring at both positions
					}
			}
		}
	}
	n_loci_to_remove=0;
	for (l=0;l<rec_pi.n_loci_to_use;++l)
		for (ll=0;ll<rec_pi.n_loci_to_use;++ll)
			if (l!=ll)
			{
				if (occurs_with[ll][l]>1 && occurs_with[ll][l]/occurs[ll]>spi->LD_threshold && occurs_with[ll][l]/occurs[ll]>=occurs_with[l][ll]/occurs[l])
				{
					loci_to_remove[n_loci_to_remove++]=ll;
					if (occurs_with[ll][l]==occurs[l] && occurs_with[ll][l]==occurs[ll])
						fprintf(fo,"Removing locus %d because it always occurs together with locus %d (%.0f times)\n",ll+1,l+1,occurs_with[l][ll]);
					else
						fprintf(fo,"Removing locus %d because it usually occurs with locus %d\n",ll+1,l+1);
					for (lll=0;lll<rec_pi.n_loci_to_use;++lll)
						occurs_with[ll][lll]=occurs_with[lll][ll]=0; // this stops ll being added to the list again and stops it being used to remove l
				}
			}
	if (n_loci_to_remove==0)
		fprintf(fo,"No loci need to be removed for being in LD\n");
	for (l=0;l<n_loci_to_remove;++l)
	{
		for (ll=loci_to_remove[l];ll<rec_pi.n_loci_to_use-1;++ll)
			rec_pi.loci_to_use[ll]=rec_pi.loci_to_use[ll+1];
		--rec_pi.n_loci_to_use;
		for (ll=l;ll<n_loci_to_remove;++ll)
			if (loci_to_remove[ll]>loci_to_remove[l])
				--loci_to_remove[ll]; // loci with larger numbers will have to be renumbered down
	}
	fprintf(fo,"Will retain the following loci:\n");
		for (l=0;l<rec_pi.n_loci_to_use;++l)
			fprintf(fo,"%3d %s\n",l+1,names[rec_pi.loci_to_use[l]]);
	for (c=0;c<3;++c)
		col_tot[c]=0;
    for (l=0;l<rec_pi.n_loci_to_use;++l)
		all_count[l]=0;
	for (s=0;s<nsub;++s)
	{
		genotype=0;
		homoz=0;
        for (l=0;l<rec_pi.n_loci_to_use;++l)
			{
				g=0;
				for (a=0;a<2;++a)
					if (sub[s]->all[rec_pi.loci_to_use[l]][a]==rec_rarer[rec_pi.loci_to_use[l]])
					{
						++g;
						if (genotype<2)
							++genotype;
						++all_count[l];
					}
				if (g==2)
					homoz=1;
				if (sub[s]->all[rec_pi.loci_to_use[l]][0]!=0)
					++geno_count[l][g];
			}
		++counts[sub[s]->cc][genotype];
		if (homoz)
			++hom_counts[sub[s]->cc];
		if (genotype==2)
		{
			first=1;
			fprintf(fo,"%-10s %d   ",sub[s]->id,sub[s]->cc);
			for (l=0;l<rec_pi.n_loci_to_use;++l)
				for (a=0;a<2;++a)
					if (sub[s]->all[rec_pi.loci_to_use[l]][a]==rec_rarer[rec_pi.loci_to_use[l]])
					{
						fprintf(fo,"%s%d",first?"":"-",l+1);
						first=0;
					}
			fprintf(fo,"\n");
		}
	}

for (r=0;r<2;++r)
  for (c=0;c<3;++c)
    {
    row_tot[r]+=counts[r][c];
    col_tot[c]+=counts[r][c];
    N+=counts[r][c];
    }
exp_hom_freq=0;
for (l=0;l<rec_pi.n_loci_to_use;++l)
{
	tot=0;
	for (g=0;g<3;++g)
		tot+=geno_count[l][g];
	MAF[l]=all_count[l]/(2*tot);
	if (spi->use_cc_freqs[0])
	{
		MAF[l]=(MAF[l]*N+contMAF[l]*max_cc[0])/(N+max_cc[0]); 
		// weighted average of control and case frequencies
		// just use max number of controls because I can no longer keep track of cc_count
	}
	exp_hom_freq=1-((1-exp_hom_freq)*(1-MAF[l]*MAF[l]));
}
// trying to get expected total proportion of cases who will be true homozygotes
// at at least one locus

// now work out probs for a sub to have 0, 1 or more B alleles, i.e. include compound heterozygotes
	geno_prob[0]=1;
	geno_prob[1]=0;
	for (l=0;l<rec_pi.n_loci_to_use;++l)
	{
		geno_prob[1]*=(1-MAF[l])*(1-MAF[l]); // probability that this locus is AA so does not change overall allele count
		geno_prob[1]+=geno_prob[0]*2*(1-MAF[l])*MAF[l]; // probability that this locus is AB
		geno_prob[0]*=(1-MAF[l])*(1-MAF[l]); 
	}
	geno_prob[2]=1-geno_prob[0]-geno_prob[1];
	fprintf(fo,"\nGenotype counts for cases and controls\n");
	fprintf(fo,"\n(AB have one minor allele, BB have two or more minor alleles)\n");
	fprintf(fo,"            AA      AB      BB   \n");
	for (r=0;r<2;++r)
	{
		fprintf(fo,"%8s ",r==0?"Controls":"Cases");
		for (c=0;c<3;++c)
			fprintf(fo,"%7.2f ",counts[r][c]);
		fprintf(fo,"\n");
	}
	for (r=0;r<2;++r)
	{
		tab[r][0]=counts[r][0]+counts[r][1];
		tab[r][1]=counts[r][2];
    }
	fprintf(fo,"%10s","\nObserved:\n");
	for (r=0;r<2;++r)
	{
		for (c=0;c<2;++c)
			fprintf(fo,"%7.2f ",tab[r][c]);
		fprintf(fo,"\n");
	}
	col_tot[0]=col_tot[0]+col_tot[1];
	col_tot[1]=col_tot[2];
	for (r=0;r<2;++r)
		for (c=0;c<2;++c)
			ex[r][c]=row_tot[r]*col_tot[c]/N;
	fprintf(fo,"%10s","\nExpected:\n");
	for (r=0;r<2;++r)
	{
		for (c=0;c<2;++c)
			fprintf(fo,"%7.2f ",ex[r][c]);
		fprintf(fo,"\n");
	}

	dchi=0;
	if (spi->use_cc_freqs[0]==0)
	{
	for (r=0;r<2;++r)
		for (c=0;c<2;++c)
			if (ex[r][c]!=0)
//				dchi+=(fabs(tab[r][c]-ex[r][c])-0.5)*(fabs(tab[r][c]-ex[r][c])-0.5)/(ex[r][c]<1?1:ex[r][c]);
// this was with Yates correction - leave out for now
				dchi+=(tab[r][c]-ex[r][c])*(tab[r][c]-ex[r][c])/(ex[r][c]<1?1:ex[r][c]);
// but do not allow expected values less than one
	}
	if (col_tot[0]==0 || col_tot[1]==0)
		dchi=0; // will not be otherwise because of Yates correction
	p=chistat(dchi,1.0)/2;
	if (ex[1][1]>tab[1][1])
		p=1-p;
	fprintf(fo,"\n\nRecessive chi-squared = %f, 1 df, p = %f\n",dchi,p);	
	fprintf(fo,"-log(p) = %.2f\n",-log10(p));


 	fprintf(fo,"\nExpected probabilities for subjects to have total of 0, 1 or more minor alleles:\n%8.6f %8.6f %8.6f \n",
		geno_prob[0],geno_prob[1],geno_prob[2]);
	// here we do the recessive HWE test for all loci combined as compound heterozygotes
	tab[0][1]=counts[1][2];
	ex[0][1]=geno_prob[2]*row_tot[1];
	tab[0][0]=row_tot[1]-tab[0][1];
	ex[0][0]=row_tot[1]-ex[0][1];
	fprintf(fo,"%10s","\nObserved: ");
	for (c=0;c<2;++c)
		fprintf(fo,"%7.2f ",tab[0][c]);
	fprintf(fo,"%10s","\nExpected: ");
	for (c=0;c<2;++c)
		fprintf(fo,"%7.2f ",ex[0][c]);
	if (ex[0][1]>5)
	{
		dchi=0;
		for (c=0;c<2;++c)
			if (ex[0][c]!=0 && tab[0][c]!=0)
				dchi+=(fabs(tab[0][c]-ex[0][c])-0.5)*(fabs(tab[0][c]-ex[0][c])-0.5)/(ex[0][c]<1?1:ex[0][c]);
		p=chistat(dchi,1.0)/2;
		if (ex[0][1]>tab[0][1])
			p=1-p;
		fprintf(fo,"\nRecessive HWE for cases chi-squared = %f, 1 df, p = %f\n",dchi,p);	
	}
	else
	{
		if (tab[0][1]==0)
			p=1.0;
		else
			p=1-cumulBinom(row_tot[1],tab[0][1]-1,geno_prob[2]);
		if (p<=0)
			p=pow((double)10,(double)-20);
		fprintf(fo,"\n\nRecessive HWE for cases exact test, p = %f\n",p);	
	}
	fprintf(fo,"-log(p) = %.2f\n",-log10(p));
	tab[0][1]=counts[0][2];
	ex[0][1]=geno_prob[2]*row_tot[0];
	tab[0][0]=row_tot[0]-tab[0][1];
	ex[0][0]=row_tot[0]-ex[0][1];
	fprintf(fo,"%10s","\nObserved: ");
	for (c=0;c<2;++c)
		fprintf(fo,"%7.2f ",tab[0][c]);
	fprintf(fo,"%10s","\nExpected: ");
	for (c=0;c<2;++c)
		fprintf(fo,"%7.2f ",ex[0][c]);
	if (ex[0][1]>5)
	{
		dchi=0;
		for (c=0;c<2;++c)
			if (ex[0][c]!=0 && tab[0][c]!=0)
				dchi+=(fabs(tab[0][c]-ex[0][c])-0.5)*(fabs(tab[0][c]-ex[0][c])-0.5)/(ex[0][c]<1?1:ex[0][c]);
		p=chistat(dchi,1.0)/2;
		if (ex[0][1]>tab[0][1])
			p=1-p;
		fprintf(fo,"\n\nRecessive HWE for controls chi-squared = %f, 1 df, p = %f\n",dchi,p);	
	}
	else
	{
		if (tab[0][1]==0)
			p=1.0;
		else
			p=1-cumulBinom(row_tot[0],tab[0][1]-1,geno_prob[2]);
		if (p<=0)
			p=pow((double)10,(double)-20);
		fprintf(fo,"\n\nRecessive HWE for controls exact test, p = %f\n",p);	
	}
	fprintf(fo,"-log(p) = %.2f\n",-log10(p));

	// Now do staandard chi-squared test for true homozygotes
	col_tot[0]=col_tot[1]=0;
	for (r=0;r<2;++r)
	{
		col_tot[1]+=tab[r][1]=hom_counts[r];
		col_tot[0]+=tab[r][0]=row_tot[r]-hom_counts[r];
    }
		fprintf(fo,"\nCounts for cases and controls homozygous for at least one variant\n");
	fprintf(fo,"          NotHOM    HOM   \n");
	for (r=0;r<2;++r)
	{
		fprintf(fo,"%8s ",r==0?"Controls":"Cases");
		for (c=0;c<2;++c)
			fprintf(fo,"%7.2f ",tab[r][c]);
		fprintf(fo,"\n");
	}
	for (r=0;r<2;++r)
		for (c=0;c<2;++c)
			ex[r][c]=row_tot[r]*col_tot[c]/N;
	fprintf(fo,"\nExpected:\n");
	for (r=0;r<2;++r)
	{
		fprintf(fo,"%8s ",r==0?"Controls":"Cases");
		for (c=0;c<2;++c)
			fprintf(fo,"%7.2f ",ex[r][c]);
		fprintf(fo,"\n");
	}

	dchi=0;
	for (r=0;r<2;++r)
		for (c=0;c<2;++c)
			if (ex[r][c]!=0)
//				dchi+=(fabs(tab[r][c]-ex[r][c])-0.5)*(fabs(tab[r][c]-ex[r][c])-0.5)/(ex[r][c]<1?1:ex[r][c]);
// this was with Yates correction - leave out for now
				dchi+=(tab[r][c]-ex[r][c])*(tab[r][c]-ex[r][c])/(ex[r][c]<1?1:ex[r][c]);
// but do not allow expected values less than one
	if (col_tot[0]==0 || col_tot[1]==0)
		dchi=0; // will not be otherwise because of Yates correction
	p=chistat(dchi,1.0)/2;
	if (ex[1][1]>tab[1][1])
		p=1-p;
	fprintf(fo,"\n\nRecessive chi-squared for homozygotes = %f, 1 df, p = %f\n",dchi,p);	
	fprintf(fo,"-log(p) = %.2f\n",-log10(p));


	 // now we do the recessive HWE test for all loci considered invidually as true homozygotes then summated 
	tab[0][1]=hom_counts[1];
	ex[0][1]=exp_hom_freq*row_tot[1];
	tab[0][0]=row_tot[1]-tab[0][1];
	ex[0][0]=row_tot[1]-ex[0][1];
	fprintf(fo,"%10s","\nObserved: ");
	for (c=0;c<2;++c)
		fprintf(fo,"%7.2f ",tab[0][c]);
	fprintf(fo,"%10s","\nExpected: ");
	for (c=0;c<2;++c)
		fprintf(fo,"%7.2f ",ex[0][c]);
	if (ex[0][1]>5)
	{
		dchi=0;
		for (c=0;c<2;++c)
			if (ex[0][c]!=0 && tab[0][c]!=0)
				dchi+=(fabs(tab[0][c]-ex[0][c])-0.5)*(fabs(tab[0][c]-ex[0][c])-0.5)/(ex[0][c]<1?1:ex[0][c]);
		p=chistat(dchi,1.0)/2;
		if (ex[0][1]>tab[0][1])
			p=1-p;
		fprintf(fo,"\n\nRecessive HWE for true homozygote cases chi-squared = %f, 1 df, p = %f\n",dchi,p);	
	}
	else
	{
		if (tab[0][1]==0)
			p=1.0;
		else
			p=1-cumulBinom(row_tot[1],tab[0][1]-1,exp_hom_freq);
		if (p<=0)
			p=pow((double)10,(double)-40);
		fprintf(fo,"\n\nRecessive HWE for true homozygote cases exact test, p = %f\n",p);	
	}
	fprintf(fo,"-log(p) = %.2f\n",-log10(p));
	tab[0][1]=hom_counts[0];
	ex[0][1]=exp_hom_freq*row_tot[0];
	tab[0][0]=row_tot[0]-tab[0][1];
	ex[0][0]=row_tot[0]-ex[0][1];
	fprintf(fo,"%10s","\nObserved: ");
	for (c=0;c<2;++c)
		fprintf(fo,"%7.2f ",tab[0][c]);
	fprintf(fo,"%10s","\nExpected: ");
	for (c=0;c<2;++c)
		fprintf(fo,"%7.2f ",ex[0][c]);
	if (ex[0][1]>5)
	{
		dchi=0;
		for (c=0;c<2;++c)
			if (ex[0][c]!=0 && tab[0][c]!=0)
				dchi+=(fabs(tab[0][c]-ex[0][c])-0.5)*(fabs(tab[0][c]-ex[0][c])-0.5)/(ex[0][c]<1?1:ex[0][c]);
		p=chistat(dchi,1.0)/2;
		if (ex[0][1]>tab[0][1])
			p=1-p;
		fprintf(fo,"\n\nRecessive HWE for true homozygote controls chi-squared = %f, 1 df, p = %f\n",dchi,p);	
	}
	else
	{
		if (tab[0][1]==0)
			p=1.0;
		else
			p=1-cumulBinom(row_tot[0],tab[0][1]-1,exp_hom_freq);
		if (p<=0)
			p=pow((double)10,(double)-40);
		fprintf(fo,"\n\nRecessive HWE for true homozygote controls exact test, p = %f\n",p);	
	}
	fprintf(fo,"-log(p) = %.2f\n",-log10(p));


	return;
}

