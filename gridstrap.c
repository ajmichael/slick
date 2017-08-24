/* version with 1-cos(beta) misfit criterion */
#include <math.h>
#include <stdio.h>
#include "controls"
#define TORADS 57.29577951
#define MAXDATA 100
#define MAX6 600
/* COORDINATES ARE EAST,NORTH,UP */

main(argc,argv)  /* slickenside inversion program */
int argc;  /* argument count */
char **argv; /* argument string */
{
	float ddir[MAXDATA];  /* dip direction for data */
	float dip[MAXDATA];   /* dip of data */
	float rake[MAXDATA];  /* rake of data */
	int   set[MAXDATA];  /* 0 if unknown, 1 if correct, 2 if backwards */
	short nobs;  /* number of observations */
	float amat[MAX6][5];  /* coefficient matrix for normal equation */
	float stress[6];  /* stress tensor in vector form, element order is: */
	/* xx,xy,xz,yy,yz,zz */
	float strten[3][3];  /* stress tensor in tensor form */
	float slick[MAX6];    /* slickenside vector elements vector */
	float n1,n2,n3;    /* normal vector elements */
	float norm[MAX6];  /* storage of n1,n2,n3 */
	short i,j,k,l;        /* dummy variables */
	short first=1;
	float z,z2,z3,z4;     /* more dummy variables */
	char name[20];      /* output file name */
	FILE *fpin;   /* input file pointer */
	FILE *fpout;  /* output file pointer */
	FILE *fplot;  /* plot file pointer */
	short whfile=0; /* decides which backup file to use */
	float a2[5][5],cc[5],sigma;  /* for use with leasq subr */
	float a2i[5][5];  /* to get covariance mtrix */
	float lam[3];  /* eigenvalues */
	float vecs[3][3];  /* eigenvectors */
	char line[80];  /* character line */
	float t[3];  /* shear stress vector */
	float iso;  /* isotropic stress mag */
	float angavg,angstd;  /* average and standard deviation of fit angle */
	float isoavg,isostd;  /* same for isotropic stress size */
	float magavg,magstd;  /* same for tangential stress size */
	float tf[3],tnorm;  /* full traction vector  */
	/* and normal traction */
	/* variables from mkten */
	float theta,beta,gamma,phi;
	float x1[3],x2[3],x3[3];
	float p1[3],p2[3];
	float smat[3][3];
	float eigen[3];
	char which[MAXDATA],bestwh[MAXDATA];
	float best[5];
	float error;
	float *ps,*ps3,*pst,*pa,*pa3,tn,tr;

	/* get file pointers */
	-- argc;  
	++argv;
	if(argc == 0){
		printf("usage: gridstrap data_file\n");
		return;
	}
	fpin=fopen(*argv,"r");
	if(fpin==NULL){
		printf("unable to open %s.\n",*argv);
		return;
	}
	sprintf(name,"%s.gboot",*argv);
	fpout=fopen(name,"a");
	if(fpout==NULL){
		printf("unable to open %s.\n",name);
		return;
	}

	/* read and write comment line from data file to output file */
	fgets(line,80,fpin);

	/* loop to get data and make up equation */
	nobs=0;
	while(fscanf(fpin,"%f%f%f%d",&ddir[nobs],&dip[nobs],&rake[nobs],&set[nobs])
	    != EOF )
	{
		i=nobs;
		j=6*nobs;
		++nobs;
		z=ddir[i]/TORADS;
		z2=dip[i]/TORADS;
		z3=rake[i]/TORADS;

		n1=sin(z)*sin(z2);  /* normal vector to fault plane */
		n2=cos(z)*sin(z2);
		n3=cos(z2);

		norm[j+0]=n1;
		norm[j+1]=n2;
		norm[j+2]=n3;

		/* slickenside vector calculation */
		slick[j]= -cos(z3)*cos(z)-sin(z3)*sin(z)*cos(z2);
		slick[j+1]= cos(z3)*sin(z)-sin(z3)*cos(z)*cos(z2);
		slick[j+2]= sin(z3)*sin(z2);

		/* find the matrix elements */
		amat[j][0]= n1-n1*n1*n1+n1*n3*n3;
		amat[j][1]= n2-2.*n1*n1*n2;
		amat[j][2]= n3-2.*n1*n1*n3;
		amat[j][3]= -n1*n2*n2+n1*n3*n3;
		amat[j][4]= -2.*n1*n2*n3;

		amat[j+1][0]= -n2*n1*n1+n2*n3*n3;
		amat[j+1][1]= n1-2.*n1*n2*n2;
		amat[j+1][2]= -2.*n1*n2*n3;
		amat[j+1][3]= n2-n2*n2*n2+n2*n3*n3;
		amat[j+1][4]= n3-2.*n2*n2*n3;

		amat[j+2][0]= -n3*n1*n1-n3+n3*n3*n3;
		amat[j+2][1]= -2.*n1*n2*n3;
		amat[j+2][2]= n1-2.*n1*n3*n3;
		amat[j+2][3]= -n3*n2*n2-n3+n3*n3*n3;
		amat[j+2][4]= n2-2.*n2*n3*n3;

		/* find reversed equations */
		norm[j+3]=slick[j];
		norm[j+4]=slick[j+1];
		norm[j+5]=slick[j+2];

		slick[j+3]= norm[j];
		slick[j+4]= norm[j+1];
		slick[j+5]= norm[j+2];

		n1=norm[j+3];
		n2=norm[j+4];
		n3=norm[j+5];

		/* find the matrix elements */
		amat[j+3][0]= n1-n1*n1*n1+n1*n3*n3;
		amat[j+3][1]= n2-2.*n1*n1*n2;
		amat[j+3][2]= n3-2.*n1*n1*n3;
		amat[j+3][3]= -n1*n2*n2+n1*n3*n3;
		amat[j+3][4]= -2.*n1*n2*n3;

		amat[j+4][0]= -n2*n1*n1+n2*n3*n3;
		amat[j+4][1]= n1-2.*n1*n2*n2;
		amat[j+4][2]= -2.*n1*n2*n3;
		amat[j+4][3]= n2-n2*n2*n2+n2*n3*n3;
		amat[j+4][4]= n3-2.*n2*n2*n3;

		amat[j+5][0]= -n3*n1*n1-n3+n3*n3*n3;
		amat[j+5][1]= -2.*n1*n2*n3;
		amat[j+5][2]= n1-2.*n1*n3*n3;
		amat[j+5][3]= -n3*n2*n2-n3+n3*n3*n3;
		amat[j+5][4]= n2-2.*n2*n3*n3;

		/* check to see if all possible data has been read */
		if(nobs==MAXDATA){
			fprintf(fpout,"NOT ALL DATA COULD BE READ.\n");
			break;
		}
	}  /* end of data read loop */
	/*printf("theta  beta   gamma  phi\n");*/

	for(theta=THETASTART;theta<THETASTOP;theta+= THETASTEP){
		for(beta= BETASTART;beta<= BETASTOP;beta+= BETASTEP){
			z=theta/TORADS;
			z2=beta/TORADS;
			x1[0]=sin(z)*cos(z2);
			x1[1]=cos(z)*cos(z2);
			x1[2]= -sin(z2);

			if(z2==90.){
				p1[0]=1.;
				p1[1]=0.;
				p1[2]=0.;
				p2[0]=0.;
				p2[1]=1.;
				p2[2]=0.;
			}
			else {
				p1[0]= -cos(z);
				p1[1]= sin(z);
				p1[2]= 0.;
				p2[0]= sin(z)*sin(z2);
				p2[1]= cos(z)*sin(z2);
				p2[2]= cos(z2);
			}

			for(gamma= GAMMASTART;gamma<= GAMMASTOP;gamma+= GAMMASTEP){
				/* construct x2,x3 */
				z=gamma/TORADS;
				z2=sin(z);
				z=cos(z);
				for(i=0;i<3;++i){
					x2[i]= z*p1[i] + z2*p2[i];
					x3[i]= -z2*p1[i] + z*p2[i];
				}

				/* fill up s matrix */
				for(i=0;i<3;++i){
					smat[i][0]=x1[i];
					smat[i][1]=x2[i];
					smat[i][2]=x3[i];
				}

				for(phi=PHISTART; phi<=PHISTOP; phi+= PHISTEP){
					/* correct phi equations */
					/* find eigenvalues */
					eigen[0]= 1;
					eigen[1]= 1-2*phi;
					eigen[2]= -1;
					/* remove isotropic component */
					z=(eigen[0]+eigen[1]+eigen[2])/3.;
					eigen[0]-= z;
					eigen[1]-= z;
					eigen[2]-= z;
					/* set maximum tangential traction to 1 */
					z= (eigen[0] - eigen[2])/2.;
					for(i=0;i<3;++i)eigen[i]/= z;

					/* now find stress tensor */
					for(i=0;i<3;++i){
						for(j=i;j<3;++j){
							strten[i][j]=0;
							for(k=0;k<3;++k)
								strten[i][j]+= smat[i][k]*eigen[k]*smat[j][k];
						}
					}
					stress[0]=strten[0][0];
					stress[1]=strten[0][1];
					stress[2]=strten[0][2];
					stress[3]=strten[1][1];
					stress[4]=strten[1][2];
					/* end of finding stress tensor */
					/* time to test it */
					/*printf("%6f %6f %6f %6f\n",theta,beta,gamma,phi);*/
					error=0;

					for(i=0;i<nobs;++i){
						/* find error measures */
						j=i*6;
						z=0;
						z2=0;
						ps = &(slick[j]);
						ps3 = &(slick[j+3]);
						z3=0;
						z4=0;
						ps= &slick[j];
						ps3= &slick[j+3];
						for(k=0;k<3;++k){
							tn=0;
							tr=0;
							pst=stress;
							pa= &amat[j+k][0];
							pa3= &amat[j+k+3][0];
							for(l=0;l<5;++l){
								tn+= (*(pa++)) * (*pst);
								tr+= (*(pa3++)) * (*(pst++));
							}
							z+= tn*(*(ps++));
							z3+= tn*tn;
							z2+= tr*(*(ps3++));
							z4+= tr*tr;
						}
						z3=sqrt(z3);
						z4=sqrt(z4);
						if(z3==0.)z=2.;
						else z=1.-(z/z3);
						if(z4==0.)z2=2.;
						else z2=1.-(z2/z4);
						if(set[i]==0){
							if(z>z2){ /* backwards case */
								error+= z2;
							}
							else { /* fowards case */
								error+= z;
							}
						}
						else if(set[i]==1){
							error+= z;
						}
						else if(set[i]==2){
							error+= z2;
						}
						if(first)continue;
						if(error>best[0])break;
					} /* end of observation loop */
					if(error>best[0] && !first)continue;
					best[0]=error;
					best[1]=theta;
					best[2]=beta;
					best[3]=gamma;
					best[4]=phi;
					first=0;
				} /* end of phi loop */
			} /* end of gamma loop */
		} /* end of beta loop */
	} /* end of theta loop */
	/* reconstruct best stress tensor axes */
	z=best[1]/TORADS;
	z2=best[2]/TORADS;
	x1[0]=sin(z)*cos(z2);
	x1[1]=cos(z)*cos(z2);
	x1[2]= -sin(z2);

	if(z2==90.){
		p1[0]=1.;
		p1[1]=0.;
		p1[2]=0.;
		p2[0]=0.;
		p2[1]=1.;
		p2[2]=0.;
	}
	else {
		p1[0]= -cos(z);
		p1[1]= sin(z);
		p1[2]= 0.;
		p2[0]= sin(z)*sin(z2);
		p2[1]= cos(z)*sin(z2);
		p2[2]= cos(z2);
	}

	/* construct x2,x3 */
	z=best[3]/TORADS;
	z2=sin(z);
	z=cos(z);
	for(i=0;i<3;++i){
		x2[i]= z*p1[i] + z2*p2[i];
		x3[i]= -z2*p1[i] + z*p2[i];
	}

	/* fill up s matrix */
	for(i=0;i<3;++i){
		smat[i][0]=x1[i];
		smat[i][1]=x2[i];
		smat[i][2]=x3[i];
	}

	eigen[0]= 1;
	eigen[1]= 1-2*best[4];
	eigen[2]= -1;
	/* remove isotropic component */
	z=(eigen[0]+eigen[1]+eigen[2])/3.;
	eigen[0]-= z;
	eigen[1]-= z;
	eigen[2]-= z;
	/* set maximum tangential traction to 1 */
	z= (eigen[0] - eigen[2])/2.;
	for(i=0;i<3;++i)eigen[i]/= z;

	/* now find stress tensor */
	for(i=0;i<3;++i){
		for(j=0;j<3;++j){
			strten[i][j]=0;
			for(k=0;k<3;++k)
				strten[i][j]+= smat[i][k]*eigen[k]*smat[j][k];
		}
	}


	fprintf(fpout,"best,theta,gamma,beta,phi=  %6f %6f %6f %6f %6f\n",
	best[0],best[1],best[2],best[3],best[4]);
	fprintf(fpout,"s3= %6f %6f %6f\n",x1[0],x1[1],x1[2]);
	fprintf(fpout,"s2= %6f %6f %6f\n",x2[0],x2[1],x2[2]);
	fprintf(fpout,"s1= %6f %6f %6f\n",x3[0],x3[1],x3[2]);
	for(i=0;i<3;++i){
		for(j=0;j<3;++j)fprintf(fpout,"%g ",strten[i][j]);
		fprintf(fpout,"\n");
	}
}
