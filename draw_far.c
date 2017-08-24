#include <math.h>
#include <stdio.h>
#define TORADS 57.29577951
#define SPACE 2
#define STEP 4.
#define SIZE 8.
draw_far(t,d,tr,pl,k,fpout)
float t,d,tr[],pl[];
int k;
FILE *fpout;

{
	float b[3];
	double z,z2,z3;
	float dist[2000],x[3][2000],print[2000];
	int flag;
	float hold;
	int i,j,l,m;
	float min;  /* minimum dot product in a swath */
	float angle,angle2,angle1;
	float firsttr,firstpl;
	int first;
	float cosp;
	int nout;

	cosp=SPACE;
	cosp/=TORADS;
	cosp=cos(cosp);
	/* compute best vector */
	z=t/TORADS;
	z2=d/TORADS;
	b[0]=sin(z)*cos(z2);
	b[1]=cos(z)*cos(z2);
	b[2]=sin(z2);

	/* loop over other vectors to compute their parameters */
	for(i=0;i<k;++i){
		/* compute vector */
		z=tr[i]/TORADS;
		z2=pl[i]/TORADS;
		x[0][i]=sin(z)*cos(z2);
		x[1][i]=cos(z)*cos(z2);
		x[2][i]=sin(z2);
		dist[i]=x[0][i]*b[0]+x[1][i]*b[1]*x[2][i]*b[2];
		if(dist[i]<0)dist[i]= -dist[i]; /* two ended axis */
		print[i]=1; /* do print this one */
	}
	/* delete vectors that are closer to b and within SPACE of another */
	for(i=0;i<k;++i){
		for(l=0;l<k;++l){
			if(l==i)continue;
			if(print[l]==0)continue;
			z=x[0][i]*x[0][l]+x[1][i]*x[1][l]+x[2][i]*x[2][l];
			if(z<0)z= -z;
			if(z>=cosp && dist[i]>=dist[l]){
				print[i]=0;
				break;
			}
		}
	}
	for(i=0;i<k;++i)if(print[i])fprintf(fpout,"l %6.2f %6.2f\n",tr[i],pl[i]);
	return;
}
