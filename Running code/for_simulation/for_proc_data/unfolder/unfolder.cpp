// DLCA_new.cpp : main project file.

//#include "stdafx.h"



#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include <float.h>   

//definition should be checked before each run 

// "n" is number of particles 
#define n 1000000

// "L" is box length
#define L 374

//"Lc" is cell length. Make sure "L" is divisible by "Lc"!!
#define Lc 2

// "LL" is equal to L/Lc or the lenght of bos in cells
#define LL L/Lc

void Unfold(struct p p[], int clustnext[], int clustnum[], int clustlead[], int Nc, int *next, int ***first);
void clust_unfold(int l,struct p p[],int clustnext[], int clustlead[],int f[], int *next, int ***first);


double rcalc(double X,double Y, double Z,double X1, double Y1,double Z1,int i,int j);
void PRINT( struct p p[], int clustnext[], int clustlead[] );
void READ( struct p p[], int clustnext[], int clustlead[] );
int initialize(struct p p[],int *next,int *prev,int ***first,int ***last,int clustlead[], int clustnext[],int clustnum[]);

struct p
{
	double x,y,z;
	double ux,uy,uz;
	int cx,cy,cz;
};

int main()
{
	
	struct p *p;
	int *clustnum, *clustnext, *clustlead;
	int i,Nc; 
	int *next,*prev,j;
	int ***first;
	int ***last;
	
	
	/////////////////////////////////////////////////////////////////
	p=(struct p *)malloc(n*sizeof(struct p));                      //
	if(p== NULL)                                                   //
    {                                                              //
      printf("p");                                                 //
    }                                                              //
	                                                               //
	clustnum=(int*)malloc(n*sizeof(int));                          //
	if(clustnum==NULL)                                             //
	{                                                              //
		printf("clustnum");                                        //
	}                                                              //
	clustnext=(int*)malloc(n*sizeof(int));                         //
	if(clustnext==NULL)                                            //
	{                                                              //
		printf("clustnext");                                       //
	}                                                              //
	clustlead=(int*)malloc(n*sizeof(int));                         //
	if(clustlead==NULL)                                            //
	{                                                              //
		printf("clustlead");									   //
	}															   //
																   //
																   //
																   //
	next=(int *)malloc(n*sizeof(int));                             //   
	if(next== NULL)                                                //      
	{                                                              //   
		printf("next");                                            //    
	}                                                              //      
                                                                   //     
	prev=(int *)malloc(n*sizeof(int));                             //            
	if(prev== NULL)                                                //          
	{                                                              //     
		printf("prev");                                            //              
	}                                                              //            
	                                                               //        
	first=(int***)malloc(LL*sizeof (int**));                       //          
	if(first==NULL)                                                //        
	{                                                              //              
		printf("first");                                           //            
	}                                                              //                
	for(i=0; i<LL; i++)                                            //           
	{                                                              //            
		first[i]=(int**)malloc(LL*sizeof(int*));                   //                
		if (first[i]==NULL)                                        //           
		{                                                          //       
			printf("first[%d]",i);                                 //           
		}                                                          //             
		for(j=0; j<LL; j++)                                        //              
		{                                                          //         
			first[i][j]=(int*)malloc(LL*sizeof(int));              //             
			if(first[i][j]==NULL)                                  //                 
			{                                                      //            
				printf("first[%d][%d]",i,j);                       //         
			}                                                      //        
		}                                                          //       
	}                                                              //                      
                                                                   //      
	last=(int***)malloc(LL*sizeof (int**));                        //               
	if(last==NULL)                                                 //         
	{                                                              //       
		printf("last");                                            //            
	}                                                              //     
	for(i=0; i<LL; i++)                                            //            
	{                                                              //         
		last[i]=(int**)malloc(LL*sizeof(int*));                    //                         
		if (last[i]==NULL)                                         //         
		{                                                          //             
			printf("last[%d]",i);                                  //         
		}                                                          //          
		for(j=0; j<LL; j++)                                        //                 
		{                                                          //         
			last[i][j]=(int*)malloc(LL*sizeof(int));               //              
			if(last[i][j]==NULL)                                   //      
			{                                                      //
				printf("last[%d][%d]",i,j);                        //          
			}                                                      //
		}                                                          //
	}                                                              // 
																   //
																   //
///////////////////////////////////////////////////////////////////// 
	READ(p, clustnext, clustlead );

	printf(" 1 \n");

	Nc=initialize( p, next, prev, first, last, clustlead, clustnext, clustnum);

	
	
	Unfold(  p ,  clustnext,  clustnum, clustlead,  Nc,  next, first);

	PRINT(  p, clustnext, clustlead );



    return 0;
}


void Unfold(struct p p[], int clustnext[], int clustnum[], int clustlead[], int Nc, int *next, int ***first)
{

	int j,a,A,part,nx,mass, *unfold;


	unfold=(int *)malloc(n*sizeof(int));                           //                                                         
	if(unfold== NULL)                                              //                 
	{                                                              //        
		printf("unfold");                                          //                 
	}                                                              // 
	/////////////////////////////////////////////////////////////////


	/// here is whrer unfloding happens/////////////////////////////////////////////////////////////////////////
	for(j=0; j<n; j++)
	{
		unfold[j]=-1;
	

	
		p[j].ux=p[j].x;
		p[j].uy=p[j].y;
		p[j].uz=p[j].z;

	

	}
	
				
	for(j=0; j<Nc; j++)
	{
		printf(" %d/%d \n", (j+1),Nc);
		part=clustnum[j];

		nx=part;
		mass=0;
		while(nx!=-1)
		{
			mass++;
			nx=clustnext[nx];
		}
		nx=part;
		printf("		%d \n",mass);



		if(part==-1){break;}
		a=0;
		while(1)
		{	
			A=0;
			//unfolding(p,clustlead,clustnext,clustnum,unfold);
			clust_unfold( part, p, clustnext, clustlead, unfold,  next, first);
			nx=part;
			while(nx!=-1)
			{
				if(unfold[nx]==-1)
				{
					A++;
				}
				nx=clustnext[nx];
			}
			printf("		%d \n",A);
			if(A==0)
			{break;}
			
			if (a==A)
			{
				printf("unfold error(  PRATICLE: %d) \n\n\n\n",j);
				break;
			}
			a=A;
						
			
		}
	//printf("unfold_%d \n",j);
	}

	

}
void clust_unfold(int l, struct p p[], int clustnext[], int clustlead[], int f[], int *next, int ***first)
{
	int nx,nb, i,j,k;
	int cxl,cyl,czl,cxn,cyn,czn;
	int cxp,cxm,cyp,cym,czp,czm;

	int cx_home, cy_home, cz_home;
	int cx_next, cy_next, cz_next, cand_part; 

	double X,Y,Z,X1,Y1,Z1;
	double r;
	p[l].ux=p[l].x; p[l].uy=p[l].y; p[l].uz=p[l].z; f[l]=1;
	nb=l;
	//loop cycles through paritcles in cluster

	i=0;
	while(nb!=-1)
	{
		
		
		
		//X=p[nb].x; Y=p[nb].y; Z=p[nb].z; cxl=p[nb].cx; cyl=p[nb].cy; czl=p[nb].cz; 
		//cxp=cxl+1; cxm=cxl-1; 
		//cyp=cyl+1; cym=cyl-1; 
		//czp=czl+1; czm=czl-1;


		cx_home=p[nb].cx;
		cy_home=p[nb].cy;
		cz_home=p[nb].cz;


		
		for(i=-1; i<=1; i++)
		{
			for(j=-1; j<=1; j++)
			{
				for(k=-1; k<=1; k++)
				{

					cx_next= cx_home + i;
					cy_next= cy_home + j;
					cz_next= cz_home + k;

					if(cx_next>=L/Lc)
					{cx_next=0;}
					if(cx_next<0)
					{cx_next=(L/Lc)-1;}
		
					if(cy_next>= L/Lc)
					{cy_next=0;}
					if(cy_next<0)
					{cy_next=(L/Lc)-1;}
		
					if(cz_next>=L/Lc)
					{cz_next=0;}
					if(cz_next<0)
					{cz_next=(L/Lc)-1;}



					cand_part=first[cx_next][cy_next][cz_next];


					while(cand_part !=-1)
					{
						if( clustlead[cand_part] == clustlead[nb] )
						{
							// here is where unfolding occurs
									
						

							if(f[nb]==1 && f[cand_part]!=1) 
							{
								// campare and adjust to unfolded corrd:mark cand_part as complete 
								if( (p[nb].ux - p[cand_part].ux)>= L/2 )
								{
									p[cand_part].ux=p[cand_part].ux+L;
								}
								if( (p[cand_part].ux - p[nb].ux) >= L/2 )
								{
									p[cand_part].ux=p[cand_part].ux-L;
								}

								if( (p[nb].uy - p[cand_part].uy)>= L/2 )
								{
									p[cand_part].uy=p[cand_part].uy+L;
								}
								if( (p[cand_part].uy - p[nb].uy) >= L/2 )
								{
									p[cand_part].uy=p[cand_part].uy-L;
								}

								if( (p[nb].uz - p[cand_part].uz)>= L/2 )
								{
									p[cand_part].uz=p[cand_part].uz+L;
								}
								if( (p[cand_part].uz- p[nb].uz) >= L/2 )
								{
									p[cand_part].uz=p[cand_part].uz-L;
								}
								f[cand_part]=1;
							}
							if(f[cand_part]==1 && f[nb]!=1)
							{
								// campare and adjust to unfolded corrd:mark l as complete
								if( (p[cand_part].ux - p[nb].ux)>= L/2 )
								{
									p[nb].ux=p[nb].ux+L;
								}
								if( (p[nb].ux - p[cand_part].ux) >= L/2 )
								{
									p[nb].ux=p[nb].ux-L;
								}

								if( (p[cand_part].uy-p[nb].uy)>= L/2 )
								{
									p[nb].uy=p[nb].uy+L;
								}
								if( (p[nb].uy-p[cand_part].uy) >= L/2 )
								{
									p[nb].uy=p[nb].uy-L;
								}

								if( (p[cand_part].uz-p[nb].uz)>= L/2 )
								{
									p[nb].uz=p[nb].uz+L;
								}
								if( (p[nb].uz-p[cand_part].uz) >= L/2 )
								{
									p[nb].uz=p[nb].uz-L;
								}
								f[nb]=1;

							}

						

						}

						cand_part=next[cand_part];
					}




				}
			}
		}
















		/*

		if(cxp>=L/Lc)
		{cxp=0;}
		if(cxm<0)
		{cxm=(L/Lc)-1;}
		
		if(cyp>= L/Lc)
		{cyp=0;}
		if(cym<0)
		{cym=(L/Lc)-1;}
		
		if(czp>=L/Lc)
		{czp=0;}
		if(czm<0)
		{czm=(L/Lc)-1;}

		nx=l;
		// loop campares particle 'nb' to all other particles in cluster
		while(nx!=-1)
		{
			//printf("%d \n",nx);
				
			X1=p[nx].x; Y1=p[nx].y; Z1=p[nx].z; cxn=p[nx].cx; cyn=p[nx].cy; czn=p[nx].cz;
			// checks if 'l' and 'nx' are in same or neigoring cell 
				

			if(	f[nb]==-1 || f[nx]==-1 )
			{

				if( (cxn==cxl||cxn==cxp||cxn==cxm) && (cyn==cyl||cyn==cyp||cyn==cym) && (czn==czl||czn==czp||czn==czm))
				{
				 
									
					// here is where unfolding occurs
									
						

					if(f[nb]==1 && f[nx]!=1) 
					{
						// campare and adjust to unfolded corrd:mark nx as complete 
						if( (p[nb].ux - p[nx].ux)>= L/2 )
						{
							p[nx].ux=p[nx].ux+L;
						}
						if( (p[nx].ux - p[nb].ux) >= L/2 )
						{
							p[nx].ux=p[nx].ux-L;
						}

						if( (p[nb].uy - p[nx].uy)>= L/2 )
						{
							p[nx].uy=p[nx].uy+L;
						}
						if( (p[nx].uy - p[nb].uy) >= L/2 )
						{
							p[nx].uy=p[nx].uy-L;
						}

						if( (p[nb].uz - p[nx].uz)>= L/2 )
						{
							p[nx].uz=p[nx].uz+L;
						}
						if( (p[nx].uz- p[nb].uz) >= L/2 )
						{
							p[nx].uz=p[nx].uz-L;
						}
						f[nx]=1;
					}
					if(f[nx]==1 && f[nb]!=1)
					{
						// campare and adjust to unfolded corrd:mark l as complete
						if( (p[nx].ux - p[nb].ux)>= L/2 )
						{
							p[nb].ux=p[nb].ux+L;
						}
						if( (p[nb].ux - p[nx].ux) >= L/2 )
						{
							p[nb].ux=p[nb].ux-L;
						}

						if( (p[nx].uy-p[nb].uy)>= L/2 )
						{
							p[nb].uy=p[nb].uy+L;
						}
						if( (p[nb].uy-p[nx].uy) >= L/2 )
						{
							p[nb].uy=p[nb].uy-L;
						}

						if( (p[nx].uz-p[nb].uz)>= L/2 )
						{
							p[nb].uz=p[nb].uz+L;
						}
						if( (p[nb].uz-p[nx].uz) >= L/2 )
						{
							p[nb].uz=p[nb].uz-L;
						}
						f[nb]=1;

					}
									
				
								
				}
			}
			nx=clustnext[nx];

			
		}
		*/
		nb=clustnext[nb];
	}
	
	return;
}



double rcalc(double X,double Y, double Z,double X1, double Y1,double Z1,int i,int j)
 {
	double r,dx,dy,dz;
	
	if(fabs(X-X1) <= L/2 )
	{
		dx=fabs(X-X1);
	}
	if( (X-X1) >= L/2)
	{
		dx=fabs(X-X1-L);
	}
	if( (X1-X) >= L/2)
	{
		dx=fabs(X1-X-L);
	}

	if(fabs(Y-Y1) <= L/2 )
	{
		dy=fabs(Y-Y1);
	}
	if( (Y-Y1) >= L/2)
	{
		dy=fabs(Y-Y1-L);
	}
	if( (Y1-Y) >= L/2)
	{
		dy=fabs(Y1-Y-L);
	}

	if(fabs(Z-Z1) <= L/2 )
	{
		dz=fabs(Z-Z1);
	}
	if( (Z-Z1) >= L/2)
	{
		dz=fabs(Z-Z1-L);
	}
	if( (Z1-Z) >= L/2)
	{
		dz=fabs(Z1-Z-L);
	}

	r=(dx*dx)+(dy*dy)+(dz*dz);
	r=sqrt(r);

	return r;
 }
void PRINT( struct p p[], int clustnext[], int clustlead[] )
{
	int i;
	
	FILE*output;
	output=fopen("output","w");

	for(i=0; i<n; i++)
	{
		fprintf(output,"%lf %lf %lf %d %d \n ", p[i].ux, p[i].uy, p[i].uz, clustlead[i], clustnext[i]);
	}


	return; 
}

void READ( struct p p[], int clustnext[], int clustlead[] )
{
	int i,x;
	double y;

	FILE*input;
	input=fopen("input","r");


	for(i=0; i<n ; i++)
	{
		fscanf(input,"%d %lf %lf %lf %lf %d %d %lf %lf %lf %lf %lf ",&x, &p[i].x, &p[i].y, &p[i].z, &y, &clustlead[i], &clustnext[i], &y, &y, &y, &y, &y);
		
		//printf("%d \n",i);	
	}

	return; 
}

int initialize(struct p p[],int *next,int *prev,int ***first,int ***last,int clustlead[] ,int clustnext[],int clustnum[])
{
	int lastp, Nc;
	int i,j,k;
	
	printf(" 2 \n");
	
	//set all of clustnum to -1
	for(i=0; i<n; i++)
	{ 
		clustnum[i]=-1;

		/////// set cells ////////////////////////////////////////
		p[i].cx= int( (p[i].x)/Lc );
		p[i].cy= int( (p[i].y)/Lc );
		p[i].cz= int( (p[i].z)/Lc );

	} 
	Nc=0;		
	printf(" 3 \n");
	//assings each clust a number
	for(i=0; i<n; i++)
	{
		if(i==clustlead[i])
		{
			clustnum[Nc]=i;
			Nc++;		
		}
	}
	///////////////////////

	printf(" 4 \n");
	
	printf("# of clust:%d\n",Nc);



	/* sets all first vaules to -1 */
	for (i=0; i<L/Lc; i++)
	{
		for(k=0; k<L/Lc; k++)
		{
			for(j=0; j<L/Lc; j++)
			{	
				first[i][k][j]=-1;
			}
		}
	}
	printf("First\n");

	/* sets all last vaules to -1 */
	for (i=0; i<L/Lc; i++)
	{
		for(k=0; k<L/Lc; k++)
		{
			for(j=0; j<L/Lc; j++)
			{	
				last[i][k][j]=-1;
			}
		}
	}
	printf("LAST\n");

	
	for(i=0; i<n; i++)
	{

		/////// set cells ////////////////////////////////////////
		p[i].cx= int( (p[i].x)/Lc );
		p[i].cy= int( (p[i].y)/Lc );
		p[i].cz= int( (p[i].z)/Lc );

		////////// initialize link list /////////////////////////////
		if(first[p[i].cx][p[i].cy][p[i].cz]==-1)
		{
			first[p[i].cx][p[i].cy][p[i].cz]=i;
			last[p[i].cx][p[i].cy][p[i].cz]=i;
			prev[i]=-1;
			next[i]=-1;
		}
		else
		{
			lastp=last[p[i].cx][p[i].cy][p[i].cz];
			next[lastp]=i;
			next[i]=-1;
			prev[i]=lastp;
			last[p[i].cx][p[i].cy][p[i].cz]=i;
		}

		
	}
	printf("Link list set!!\n");

	return Nc;
}