/* scoreassoc.cpp */

#include "scoreassoc.h"

// Throughout, I am going to calculate a one-tailed p but then report SLP for the two-tailed p

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

double one_tailed_binomial_p(int N, int k, double p)
{
	double pval;
	if (k>N*p) // use right tail
		pval=1-cumulBinom(N,k-1,p);
	else if (k+1<N*p)
		pval=cumulBinom(N,k,p);
	else
		pval=0.5;
	return pval;
}

void do_recessive_HWE_test(FILE *fo,float *score,subject **sub,int nsub,par_info *pi,sa_par_info *spi,float cc_freq[2][MAX_LOCI],float cc_count[2][MAX_LOCI],int max_cc[2],float *weight,float *missing,int *rarer,char names[MAX_LOCI][20])
{
	float tab[2][2],ex[2][2],col_tot[3],row_tot[2],N,counts[2][3],hom_counts[2],all_count[MAX_LOCI],exp_hom_freq,homoz;
	int r,c,s,l,ll,lll,pl,pll,plll,a,genotype,first,all_genotype,g;
	int n_loci_to_remove,loci_to_remove[MAX_LOCI];
	double p,df,dchi;
	float freq,tot,geno_prob[3],geno_count[MAX_LOCI][3];
	int status,which=1;
	par_info rec_pi; // use this to keep track of which loci we are using
	// first thing we will do is remove loci which fail the weight threshold
	// note that weights and rarer are indexed up to pi->n_loci_to_use
	for (pl=0,ll=0;pl<pi->n_loci_to_use;++pl)
	{
		l=pi->loci_to_use[pl];
		if (weight[l]>=spi->weight_threshold)
		{
			rec_pi.loci_to_use[ll]=l;
			if (spi->use_cc_freqs[0])
				contMAF[l]=rarer[l]==2?cc_freq[0][l]:(1-cc_freq[0][l]);
			++ll;
		}
	}
	rec_pi.n_loci_to_use=ll;

	N=0;
	for (r=0;r<2;++r)
	{
		for (c=0;c<3;++c)
			counts[r][c]=0;
		row_tot[r]=0;
		hom_counts[r]=0;
	}
	for (pl=0;pl<rec_pi.n_loci_to_use;++pl)
		for (g=0;g<3;++g)
			geno_count[pl][g]=0;
	// begin by just listing all subjects with two or more variants at these loci before excluding on LD
	fprintf(fo,"\nRecessive analyses based on loci with weight reaching %f\n",spi->weight_threshold);
	fprintf(fo,"Using the following loci:\n");
	for (pl = 0; pl < rec_pi.n_loci_to_use; ++pl)
	{
		l = rec_pi.loci_to_use[pl];
		fprintf(fo, "%3d %s\n", pl + 1, names[l]);
	}
	fprintf(fo,"\nSubjects with two or more minor alleles:\n");
	for (s=0;s<nsub;++s)
	{
		genotype=0;
        for (pl=0;pl<rec_pi.n_loci_to_use;++pl)
			{
				l = rec_pi.loci_to_use[pl];
				for (a=0;a<2;++a)
					if (sub[s]->all[l][a]==rarer[l])
					{
						if (genotype<2)
							++genotype;
					}
			}

		if (genotype==2)
		{
			first=1;
			fprintf(fo,"%-10s %d   ",sub[s]->id,sub[s]->cc);
        for (pl=0;pl<rec_pi.n_loci_to_use;++pl)
			{
				l = rec_pi.loci_to_use[pl];
				for (a=0;a<2;++a)
					if (sub[s]->all[l][a] == rarer[l])
					{
					fprintf(fo, "%s%d", first ? "" : "-", pl + 1);
					first = 0;
					}
			}
			fprintf(fo,"\n");
		}
	}
	// now do the exclusions based on LD
	for (pl=0;pl<rec_pi.n_loci_to_use;++pl)
	{
		l=rec_pi.loci_to_use[pl];
		for (pll = 0; pll < rec_pi.n_loci_to_use; ++pll)
		{
		ll=rec_pi.loci_to_use[pll];
			if (l != ll)
				occurs_with[l][ll] = 0;
		}
		occurs[l]=0;
	}
	for (s=0;s<nsub;++s)
	{
		for (pl=0;pl<rec_pi.n_loci_to_use;++pl)
		{
			l=rec_pi.loci_to_use[pl];
			occ[0]=(sub[s]->all[l][0]==rarer[l])+(sub[s]->all[l][1]==rarer[l]); // count number of rare alleles at first locus
			if (occ[0]!=0)
			{
				occurs[l]+=occ[0];
				for (pll=0;pll<rec_pi.n_loci_to_use;++pll)
					if (pl!=pll)
					{
					ll=rec_pi.loci_to_use[pll];
					occ[1]=(sub[s]->all[ll][0]==rarer[ll])+(sub[s]->all[ll][1]==rarer[ll]);
					occurs_with[l][ll]+=MIN(occ[0],occ[1]); // number of variants occurring together at both positions, 2 if both homozygous
					}
			}
		}
	}
	n_loci_to_remove=0;
	for (pl=0;pl<rec_pi.n_loci_to_use;++pl)
		for (pll=0;pll<rec_pi.n_loci_to_use;++pll)
			if (pl!=pll)
			{
			l=rec_pi.loci_to_use[pl];
			ll=rec_pi.loci_to_use[pll];
				if (occurs_with[ll][l]>1 && occurs_with[ll][l]/occurs[ll]>spi->LD_threshold && occurs_with[ll][l]/occurs[ll]>=occurs_with[l][ll]/occurs[l])
				{
					loci_to_remove[n_loci_to_remove++]=pll;
					if (occurs_with[ll][l]==occurs[l] && occurs_with[ll][l]==occurs[ll])
						fprintf(fo,"Removing locus %d because it always occurs together with locus %d (%.0f times)\n",pll+1,pl+1,occurs_with[l][ll]);
					else
						fprintf(fo,"Removing locus %d because it usually occurs with locus %d\n",pll+1,pl+1);
					for (plll = 0; plll < rec_pi.n_loci_to_use; ++plll)
					{
						lll = rec_pi.loci_to_use[plll];
						occurs_with[ll][lll] = occurs_with[lll][ll] = 0; // this stops pll being added to the list again and stops it being used to remove pl
					}
				}
			}
	if (n_loci_to_remove==0)
		fprintf(fo,"No loci need to be removed for exceeding the LD threshold of %f\n",spi->LD_threshold);
	for (pl=0;pl<n_loci_to_remove;++pl)
	{
		for (pll=loci_to_remove[pl];pll<rec_pi.n_loci_to_use-1;++pll)
			rec_pi.loci_to_use[pll]=rec_pi.loci_to_use[pll+1];
		--rec_pi.n_loci_to_use;
		for (pll=l;pll<n_loci_to_remove;++pll)
			if (loci_to_remove[pll]>loci_to_remove[pl])
				--loci_to_remove[pll]; // loci with larger numbers will have to be renumbered down
	}
	fprintf(fo,"Will retain the following loci:\n");
	for (pl = 0; pl < rec_pi.n_loci_to_use; ++pl)
	{
		l=rec_pi.loci_to_use[pl];
		fprintf(fo, "%3d %s\n", pl + 1, names[l]);
	}
	for (c=0;c<3;++c)
		col_tot[c]=0;
	for (pl = 0; pl < rec_pi.n_loci_to_use; ++pl)
	{
		l=rec_pi.loci_to_use[pl];
		all_count[l] = 0;
	}
	fprintf(fo,"\nSubjects with two or more minor alleles:\n");
	for (s=0;s<nsub;++s)
	{
		genotype=0;
		homoz=0;
        for (pl=0;pl<rec_pi.n_loci_to_use;++pl)
			{
			l=rec_pi.loci_to_use[pl];
				g=0;
				for (a=0;a<2;++a)
					if (sub[s]->all[l][a]==rarer[l])
					{
						++g;
						if (genotype<2)
							++genotype;
						++all_count[l];
					}
				if (g==2)
					homoz=1;
				if (sub[s]->all[l][0]!=0)
					++geno_count[pl][g];
			}
		++counts[sub[s]->cc][genotype];
		if (homoz)
			++hom_counts[sub[s]->cc];
		if (genotype==2)
		{
			first=1;
			fprintf(fo,"%-10s %d   ",sub[s]->id,sub[s]->cc);
			for (pl = 0; pl < rec_pi.n_loci_to_use; ++pl)
			{
				l=rec_pi.loci_to_use[pl];
				for (a = 0; a < 2; ++a)
					if (sub[s]->all[l][a] == rarer[l])
					{
					fprintf(fo, "%s%d", first ? "" : "-", pl + 1);
					first = 0;
					}
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
for (pl=0;pl<rec_pi.n_loci_to_use;++pl)
{
	l=rec_pi.loci_to_use[pl];
	tot=0;
	for (g=0;g<3;++g)
		tot+=geno_count[pl][g];
	MAF[l]=all_count[l]/(2*tot);
	if (spi->use_cc_freqs[0])
	{
		MAF[l]=(MAF[l]*N+contMAF[l]*max_cc[0])/(N+max_cc[0]); 
		// weighted average of control and case frequencies
		// just use max number of controls 
	}
	exp_hom_freq=1-((1-exp_hom_freq)*(1-MAF[l]*MAF[l]));
}
// trying to get expected total proportion of cases who will be true homozygotes
// at at least one locus

// now work out probs for a sub to have 0, 1 or more B alleles, i.e. include compound heterozygotes
	geno_prob[0]=1;
	geno_prob[1]=0;
	for (pl=0;pl<rec_pi.n_loci_to_use;++pl)
	{
		l=rec_pi.loci_to_use[pl];
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
	p=chistat(dchi,1.0)/2; // one-tailed
	fprintf(fo,"\n\nRecessive chi-squared = %f, 1 df, p = %f\n",dchi,2*p);	
	fprintf(fo,"SLP = %.2f(signed log10(p), positive if cases more frequently have two minor alleles than controls)\n",log10(2*p)*(ex[1][1]>tab[1][1]?1:-1));
	if (ex[1][1]>tab[1][1])
		p=1-p;

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
		fprintf(fo,"\nRecessive HWE for cases chi-squared = %f, 1 df, p = %f\n",dchi,2*p);	
		if (ex[0][1]>tab[0][1])
			p=1-p;
	}
	else
	{
		if (tab[0][1]==0)
			p=0.5;
		else
			p=one_tailed_binomial_p(row_tot[1],tab[0][1],geno_prob[2]);
		if (tab[0][1]/row_tot[1]<geno_prob[2])
			p=1-p;
		if (p<=0)
			p=pow((double)10,(double)-20);
		fprintf(fo,"\n\nRecessive HWE for cases exact test, p = %f\n",2*p);	
	}
	fprintf(fo,"SLP = %.2f (signed log10(p), positive if cases more frequently have two minor alleles than expected)\n",p<0.5?-log10(2*p):log10(2*(1-p)));
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
		fprintf(fo,"\n\nRecessive HWE for controls chi-squared = %f, 1 df, p = %f\n",dchi,2*p);	
		if (ex[0][1]>tab[0][1])
			p=1-p;
	}
	else
	{
		if (tab[0][1]==0)
			p=0.5;
		else
			p=one_tailed_binomial_p(row_tot[1],tab[0][1],geno_prob[2]);
		if (tab[0][1]/row_tot[1]<geno_prob[2])
			p=1-p;
		if (p<=0)
			p=pow((double)10,(double)-20);
		fprintf(fo,"\n\nRecessive HWE for controls exact test, p = %f\n",2*p);	
	}
	fprintf(fo,"SLP = %.2f (signed log10(p), positive if controls more frequently have two minor alleles than expected)\n",p<0.5?-log10(2*p):log10(2*(1-p)));

	// Now do standard chi-squared test for true homozygotes
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
	p=chistat(dchi,1.0)/2; // one-tailed
	fprintf(fo,"\n\nRecessive chi-squared for homozygotes = %f, 1 df, p = %f\n",dchi,2*p);	
	fprintf(fo,"SLP = %.2f (signed log10(p), positive if homozygotes more frequent in cases than controls)\n",log10(2*p)*(ex[1][1]>tab[1][1]?1:-1));
	if (ex[1][1]>tab[1][1])
		p=1-p;


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
		fprintf(fo,"\n\nRecessive HWE for true homozygote cases chi-squared = %f, 1 df, p = %f\n",dchi,2*p);	
		if (ex[0][1]>tab[0][1])
			p=1-p;
	}
	else
	{
		if (tab[0][1]==0)
			p=0.5;
		else
			p=one_tailed_binomial_p(row_tot[1],tab[0][1]-1,exp_hom_freq);
		if (tab[0][1]/row_tot[1]<exp_hom_freq)
			p=1-p;
		if (p<=0)
			p=pow((double)10,(double)-40);
		fprintf(fo,"\n\nRecessive HWE for true homozygote cases exact test, p = %f\n",2*p);	
	}
	fprintf(fo,"SLP = %.2f (signed log10(p), positive if homozygote cases more frequent than expected under HWE)\n",p<0.5?-log10(2*p):log10(2*(1-p)));
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
		fprintf(fo,"\n\nRecessive HWE for true homozygote controls chi-squared = %f, 1 df, p = %f\n",dchi,2*p);	
		if (ex[0][1]>tab[0][1])
			p=1-p;
	}
	else
	{
		if (tab[0][1]==0)
			p=0.5;
		else
			p=one_tailed_binomial_p(row_tot[1],tab[0][1]-1,exp_hom_freq);
		if (tab[0][1]/row_tot[1]<exp_hom_freq)
			p=1-p;
		if (p<=0)
			p=pow((double)10,(double)-40);
		fprintf(fo,"\n\nRecessive HWE for true homozygote controls exact test, p = %f\n",p);	
	}
	fprintf(fo,"SLP = %.2f (signed log10(p), positive if homozygote controls more frequent than expected under HWE)\n",p<0.5?-log10(2*p):log10(2*(1-p)));


	return;
}

