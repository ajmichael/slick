#include <stdio.h>
#include <math.h>
main(argc,argv)
int argc; 
char **argv;
{
	char line[100];
	float tr1[2000],tr2[2000],tr3[2000];
	float pl1[2000],pl2[2000],pl3[2000];
	double t1,p1,t2,p2,t3,p3,tr,pl;
	double e,n,u;
	int flag;
	FILE *fpin,*fpout;
	char namein[20],nameout[20];
	float best[3][3],bestmag;
	float stress[3][3],mag;
	float dot[2000];
	int i,j,k;
	float level,level80;
	float phi,phimin,phimax;
	float conf;
	float z;

	phimin=1;
	phimax=0;

	if(argc!=4){
		printf("usage: plotbootg *.gboot output_file confidence_level\n");
		return;
	}

	++argv;
	sscanf(*argv,"%s",namein);
	fpin=fopen(namein,"r");
	if(fpin==NULL){
		printf("unable to open %s\n",namein);
		return;
	}
	++argv;
	sscanf(*argv,"%s",nameout);
	fpout=fopen(nameout,"w");
	if(fpout==NULL){
		printf("unable to open %s\n",nameout);
		return;
	}
	++argv;
	sscanf(*argv,"%f",&conf);

	/* first go through file to find confidence level */
	fprintf(fpout,"title %s %g %%\nr 6.\nsize line 6\n",namein,conf);
	fgets(line,100,fpin);
	fgets(line,100,fpin);
	sscanf(&line[3],"%lf %lf %lf",&e,&n,&u);
	dirplg(e,n,u,&t3,&p3);
	fprintf(fpout,"symbol line 3 3\nl %5.1f %5.1f\n",t3,p3);
	fgets(line,100,fpin);
	sscanf(&line[3],"%lf %lf %lf",&e,&n,&u);
	dirplg(e,n,u,&t2,&p2);
	fprintf(fpout,"symbol line 2 2\nl %5.1f %5.1f\n",t2,p2);
	fgets(line,100,fpin);
	sscanf(&line[3],"%lf %lf %lf",&e,&n,&u);
	dirplg(e,n,u,&t1,&p1);
	fprintf(fpout,"symbol line 1 1\nl %5.1f %5.1f\n",t1,p1);
	fgets(line,100,fpin);
	sscanf(line,"%f %f %f",&best[0][0],&best[0][1],&best[0][2]);
	fgets(line,100,fpin);
	sscanf(line,"%f %f %f",&best[1][0],&best[1][1],&best[1][2]);
	fgets(line,100,fpin);
	sscanf(line,"%f %f %f",&best[2][0],&best[2][1],&best[2][2]);
	tenmag(best,&bestmag);
	fprintf(fpout,"size line 3\n");
	i=0;
	while(fgets(line,100,fpin)!=NULL){
		if(line[0]!='b')fprintf(stderr,"OOPS i=%d\n",i);
		fgets(line,100,fpin);
		fgets(line,100,fpin);
		fgets(line,100,fpin);
		fgets(line,100,fpin);
		sscanf(line,"%f %f %f",&stress[0][0],&stress[0][1],&stress[0][2]);
		fgets(line,100,fpin);
		sscanf(line,"%f %f %f",&stress[1][0],&stress[1][1],&stress[1][2]);
		fgets(line,100,fpin);
		sscanf(line,"%f %f %f",&stress[2][0],&stress[2][1],&stress[2][2]);
		tenmag(stress,&mag);
		tendot(best,stress,bestmag,mag,&dot[i]);
		i++;
	}
	sort(dot,i);
	j=i*((100.-conf)/100.);
	level80=dot[j];
	fclose(fpin);
	fpin=fopen(namein,"r");
	fgets(line,100,fpin);
	fgets(line,100,fpin);
	fgets(line,100,fpin);
	fgets(line,100,fpin);
	fgets(line,100,fpin);
	fgets(line,100,fpin);
	fgets(line,100,fpin);

	k=0;
	for(j=0;j<i;++j){
		fgets(line,100,fpin);
		if(line[0]!='b')fprintf(stderr,"OOPS j=%d\n",j);
		sscanf(&line[27],"%f %f %f %f %f",&z,&z,&z,&z,&phi);
		fgets(line,100,fpin);
		sscanf(&line[3],"%lf %lf %lf",&e,&n,&u);
		dirplg(e,n,u,&tr,&pl);
		tr3[k]=tr;
		pl3[k]=pl;
		fgets(line,100,fpin);
		sscanf(&line[3],"%lf %lf %lf",&e,&n,&u);
		dirplg(e,n,u,&tr,&pl);
		tr2[k]=tr;
		pl2[k]=pl;
		fgets(line,100,fpin);
		sscanf(&line[3],"%lf %lf %lf",&e,&n,&u);
		dirplg(e,n,u,&tr,&pl);
		tr1[k]=tr;
		pl1[k]=pl;
		fgets(line,100,fpin);
		sscanf(line,"%f %f %f",&stress[0][0],&stress[0][1],&stress[0][2]);
		fgets(line,100,fpin);
		sscanf(line,"%f %f %f",&stress[1][0],&stress[1][1],&stress[1][2]);
		fgets(line,100,fpin);
		sscanf(line,"%f %f %f",&stress[2][0],&stress[2][1],&stress[2][2]);
		tenmag(stress,&mag);
		tendot(best,stress,bestmag,mag,&level);
		if(level>=level80){
			k++;
			if(phi<phimin)phimin=phi;
			if(phi>phimax)phimax=phi;
		}
	}
	fprintf(fpout,"# phirange = %g %g\n",phimin,phimax);
	fprintf(fpout,"symbol line 1 1\n");
	draw_far(t1,p1,tr1,pl1,k,fpout);
	fprintf(fpout,"symbol line 2 2\n");
	draw_far(t2,p2,tr2,pl2,k,fpout);
	fprintf(fpout,"symbol line 3 3\n");
	draw_far(t3,p3,tr3,pl3,k,fpout);
}

tendot(ten1,ten2,mag1,mag2,pdot)
float ten1[3][3];
float ten2[3][3];
float mag1,mag2;
float *pdot;
{
	int i,j;

	*pdot=0;
	for(i=0;i<3;++i)for(j=0;j<3;++j)*pdot+=ten1[i][j]*ten2[i][j];
	*pdot/=mag1*mag2;
	return;
}

tenmag(ten,pmag)
float ten[3][3];
float *pmag;
{
	int i,j;
	double z;

	z=0;
	for(i=0;i<3;++i)for(j=0;j<3;++j)z+=ten[i][j]*ten[i][j];
	z=sqrt(z);
	*pmag=z;
	return;
}
