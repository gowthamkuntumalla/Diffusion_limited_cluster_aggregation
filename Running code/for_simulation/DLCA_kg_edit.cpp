// DLCA_new.cpp : main project file.

/**************************NOTES****************************************/
/*
Stop the iterations when Nc =1 happens.
Lesser L is computationally preferred
change names of DLCA files(if possible) for each run.  This command gives working directory $ pwdx <PID> in shell
*/


//#include "stdafx.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include <float.h>   

//definition should be checked before each run 

/***********SIMULATION DEFINITIONS**************/

// "n" is number of particles 
#define n 1000000

// "L" is box length
#define L 374

// "seed" is the seed for random numbers . It can be anything. Change seed (large odd number)for different output 
#define seed (int) 32103255

//"Lc" is cell length. Make sure "L" is divisible by "Lc"!!
#define Lc 2

// record each cluster mass at with Kinetic data : 0 is on and 1 is off
#define mode2 1 

// "dt" how often to take kinetics data in mode 0 
#define dt 1000

// "dt_start"  
#define dt_start 10

// unfold time/ 
/**starts after large time i.e virtually unfolding is not done**/
#define dt2 100000000000

// "Stop_time" is the time the simulation stops
#define Stop_time 10000000

// "Start_time" is the the simulation starts used only for restarts. Must include the unfold file
#define Start_time 0





//min overlap allowed between monomers
#define r_min 0.99999

//max gap allowed between monomers
#define r_max 1.00001

//"mode" is data taking mode: if 0 at regular intervel of dt
#define mode 0 

// "LL" is equal to L/Lc or the length of bos in cells
#define LL L/Lc

//"PI" is PI
#define PI 4.0*atan(1.0)


//fv, initial volume fraction of monomers
#define fv n*(4/3)*PI*0.5*0.5*0.5*(1/(L*L*L));

//struct for particles 
struct p
{
	double  x, y, z; // actual position
	double ux,uy,uz;// unfolded coordinates
	double vx,vy,vz; // random velocity (assigned later in initialize)
	double P;
	int    cx,cy,cz; // cell coordinates
};


// FUNCTION LIST //

//random number generator 
double psdrand(int );

//finds memory for array's
int start_array(struct p *p);

//initialize the arrays
void initialize(struct p p[],int *next,int *prev,int ***first,int ***last,int clustlead[],int clustlast[],int clustnext[],int clustnum[]);

// check to see if initalize placed particles in overlap 
void initialcheck(struct p p[],int *next,int *prev,int ***first,int ***last);

// measures starting particles radius
int init_part_check(int i,struct p p[],int *next,int *prev,int ***first,int ***last);

//moves particle that are in overlap
void init_part_move(int i,struct p p[],int *next,int *prev,int ***first,int ***last);

// checks PBC for particle "i" ; PBC is periodic boundary condition
void PBC(struct p p[], int i);

//checks and sets link list of particle "i" when it moves into a new cell
void link_list(struct p p[], int i,int cxo,int cyo,int czo,int *next,int *prev,int ***first,int ***last);

// finds r between two particles using PBC
double rcalc(double X,double Y, double Z,double X1, double Y1,double Z1,int i,int j);

// here is the engine that runs the aggregation 
void iteration(struct p p[],double *time,int *Nc,int clustlead[],int clustnext[],int clustnum[],int clustlast[],int *next,int *prev,int ***first,int ***last,int unfold[]);

// checks to see if new cluster is formed 
int cluster_check(int c,struct p p[], int *next, int *prev, int ***first, int ***last, int clustlead[], int clustlast[], int clustnext[], int clustnum[],int new_part_list[],int *new_part_count);

// moves monomers that are in overlap 
void overlap(int part,int j,double r,struct p p[],int clustlead[],int clustlast[],int clustnext[],int *next,int *prev,int ***first,int ***last);

// takes when clust and sets up the clust list 
void make_clust(struct p p[], int i, int j, int clustlead[], int clustlast[], int clustnext[]);

//
//measures Nc Rg with time; NC = number of clusters and Rg = radius of gyration 
void Kinetics(struct p p[],int clustnum[],int clustlead[],int clustnext[],double *avRg, int Rg_check, int Nc, double *avRnn,  double time);
//

//finds PROB for DLCA 
void PROB(int c,struct p p[],int clustnext[], int clustlead[],int unfold[]);
//

//unfolds clust 
void clust_unfold(int l,struct p p[],int clustnext[],int f[]);

// finds a cluster center of mass and returns cluster mass
int R_CM(int clust_lead, double r_cm[],int clustnext[],struct p p[]);

//Radius of gyration
double Rg(int clust_lead,struct p p[], int clustnext[]);

double PROB(int part,struct p p[],int clustnext[]);

void print_files(int Rg_check, struct p p[],int clustlead[], int clustlast[], int clustnext[], double time, int Nc, double avRg);

void Unfold(int unfold[], struct p p[], int clustnext[], int clustnum[], int Nc);

/***********MAIN FUNCTION**************/
int main()
{
	//varibles
	struct p *p;
	int *next,*prev,i,j;
	int ***first;
	int ***last;
	int *clustlead,*clustlast,*clustnext;
	int *clustnum;
	int Nc, Nc_count[500], nc, update;
	
	double time,t,t2;

	int *unfold, Rg_check;
	double avRg, avRnn;
	int index=1;

	FILE *output;
	FILE *update_time;
	
	char filename[30],filename_time[30];

	nc=0;

	/*Nc_count[0]=50000;
	Nc_count[1]=10000;
	Nc_count[2]=8000;
	Nc_count[3]=2000;
	Nc_count[4]=1000;
	Nc_count[5]=500;
	Nc_count[6]=300;
	Nc_count[7]=200;
	Nc_count[8]=150;
	Nc_count[9]=100;
	Nc_count[10]=60;
	Nc_count[11]=35;
	Nc_count[12]=20;
	Nc_count[13]=10;
	Nc_count[14]=2;*/
	Nc_count[0]=50000;
	Nc_count[1]=10002;
	for(int foo=2;foo<18;foo++)
	{
		Nc_count[foo]=Nc_count[index]-500;
		index=foo;
	}
	for(int foo=18;foo<500;foo++)
	{
		Nc_count[foo]=Nc_count[index]-5;
		index=foo;
	}	
	/*In paper_runs data_i goes from 1 and 2 to 18 (500 steps) and 19 to 418 (5 steps) */
	
	// Edit it with taking it every time steps 5 or less
	// script something to simulate IMPORT function in python or C++
	// Structure factor calc for Single cluster fractal dimension (2.5). Reference: Sorensen paper. small angle , Xray scattering, Dissertations
	
//////////////////////// memory allocation///////////////////////////
	//array of particle data
	p=(struct p *)malloc(n*sizeof(struct p));                      //  
	if(p== NULL)                                                   //  
	{                                                              // 
		printf("p");                                               //  
	}                                                              // 
	                                                               //
	//next[i]=next particle to 'i'th particle
	next=(int *)malloc(n*sizeof(int));                             //   
	if(next== NULL)                                                //      
	{                                                              //   
		printf("next");                                            //    
	}                                                              //      
                                                                   //     
	//prev[i]=previous particle to 'i'th particle
	prev=(int *)malloc(n*sizeof(int));                             //            
	if(prev== NULL)                                                //          
	{                                                              //     
		printf("prev");                                            //              
	}                                                              //            
	                                                               //        
	//Array of clusters first to last
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
	//Array of clusters last to first
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
	//Lead particle of each cluster
	clustlead=(int*)malloc(n*sizeof(int));                         //               
	if(clustlead==NULL)                                            //       
	{                                                              //
		printf("clustlead");                                       //         
	}                                                              //  
                                                                   //              
	//Last particle of each cluster
	clustlast=(int *)malloc(n*sizeof(int));                        //                      
	if(clustlast==NULL)                                            //          
	{                                                              //             
		printf("clustlast");                                       //                
	}                                                              //                 
                                                                   //                
	clustnext=(int *)malloc(n*sizeof(int));                        //                         
	if(clustnext==NULL)                                            //                
	{                                                              //
		printf("clustnext");                                       //                  
	}                                                              //   
                                                                   // 
	clustnum=(int*)malloc(n*sizeof(int));                          //                                
	if(clustlead==NULL)                                            //                      
	{                                                              //                
		printf("clustnum");                                        //                  
	}                                                              //          
	                                                               //                             
	unfold=(int *)malloc(n*sizeof(int));                           //                                                         
	if(unfold== NULL)                                              //                 
	{                                                              //        
		printf("unfold");                                          //                 
	}                                                              //           
///////////////////////////////////////////////////////////////////// 

	initialize(p,next,prev,first,last,clustlead,clustlast,clustnext,clustnum);
	if(Start_time== 0)
	{
		FILE *kinetics; 
		kinetics=fopen("time","w");
		
		fprintf(kinetics,"time Nc  avRg  \n");
		
		fclose(kinetics);

		time=0.0;
		printf("System pre-initialize ");
		initialcheck(p,next,prev,first,last);

		 Nc=n;  t=dt_start;  t2=dt2;
	}
	else
	{
		
		
		time=Start_time;
		t=Start_time+dt;  t2=Start_time+dt2;
		Nc=0; 
		for(i=0; i<n; i++)
		{
			clustnum[i]=-1;
		}
		for(i=0; i<n; i++)
		{
			if( clustlead[i]==i )
			{
				clustnum[Nc]=i;
				Nc++;
				
			}	
		}
	}

    /** can possibly change this to .1 or something low to save time **/
	
	update=(int)(0.1 * Nc);// initializing the 'update' variable
	/*EDIT DONE by KG in the while condition*/
	while((time<=Stop_time)&&(Nc>1))
	{
		
		iteration(p,&time,&Nc,clustlead,clustnext,clustnum,clustlast,next,prev,first,last,unfold);

		printf("Nc:%d Time:%lf \n",Nc,time);
		

		
		
		
		/// here is where data is written 
		if(mode==0)
		{
			if(time>=t)
			{
				t=t+dt;
				
				Rg_check=0;
				if(time>=t2)
				{
					t2=t2+dt2;
					Rg_check=1;
					
					
					/// here is where unfloding happens/////////////////////////////////////////////////////////////////////////
					Unfold(unfold, p , clustnext, clustnum, Nc);
				}
				
				///////////////// here is where <Rg> is found ///////////////
				
				avRg=0;  avRnn=0; 
				Kinetics(p, clustnum, clustlead, clustnext, &avRg, Rg_check, Nc, &avRnn, time);
				/////////////////////////////////////////////////////////////

				
				print_files( Rg_check,  p ,clustlead, clustlast, clustnext, time, Nc, avRg);

				
				///////////////////////////////////////////////////////////////////////////////////////////////////////////////

			}
		}

		if(Nc<update)
		{
			output=fopen("unfold","w");
			for(i=0; i<n; i++)
			{
				fprintf(output, "%d %lf %lf %lf %lf %d %d %d %lf %lf %lf %lf \n",i, p[i].x, p[i].y, p[i].z, 0.5, clustlead[i], clustnext[i], clustlast[i], p[i].P, p[i].vx, p[i].vy, p[i].vz ); 
			}
			fclose(output);
            /** can possibly change this to .1 or something low to save time **/
			update=(int)(0.1*Nc); // updating the 'update' variable
			update_time=fopen("update_time","w");
			fprintf(update_time,"%lf \n",time);
			fclose(update_time); 
			
		}


		/** GETTING THE SNAPSHOT OF THE SYSTEM**/	
		if(Nc< Nc_count[nc])
		{
		
			nc++;		
			/// here is where unfloding happens/////////////////////////////////////////////////////////////////////////
			Unfold(unfold, p , clustnext, clustnum, Nc);
			
			sprintf(filename, "data%d",nc); ///Add time to file name
			//creating a new file with name data%nc
			output=fopen(filename,"w");
			for(i=0; i<n; i++)
			{
				fprintf(output, "%d %lf %lf %lf %lf %d %d %d %lf %lf %lf %lf \n",i, p[i].ux, p[i].uy, p[i].uz, 0.5, clustlead[i], clustnext[i], clustlast[i], p[i].P, p[i].vx, p[i].vy, p[i].vz ); 
			}
			fclose(output);
			// ADDED THIS ON 15 June to get time, avg Nc, Avg Rg 
			sprintf(filename_time, "data_time%d",nc); ///Add time to file name
			output=fopen(filename_time,"w");
			fprintf(output,"%lf \n",time);
			fclose(output);
		}




	}

    return 0;
}

//
//
//
//

/***********AUXILIARY FUNCTIONS**************/

double psdrand(int iseed)
{
	int i, j, k, inx;
	double ran_num;
	static const int ndim = 55, m10 = 1000000000, is = 21, ir = 30;
	static const double base = 1.0E9;
	static int jrand, istack[58];
	static bool init = false;

	if((!init) || (iseed < 0))
	{
		iseed = abs(iseed);
 		istack[ndim] = iseed;
		j = iseed;
		k = 1;

		for(i = 1; i <= (ndim - 1); ++i)
		{
			inx = i*is - int((double)(i*is)/(double)(ndim))*ndim;
			istack[inx] = k;
			k = j - k;
			if(k < 0) {k += m10;}
			j = istack[inx];
		}

		for(j = 1; j <= 3; ++j)
		{
			for(i = 1; i <= ndim; ++i)
			{
				inx = i + ir - int((double)(i+ir)/(double)(ndim))*ndim;
				istack[i] -= istack[inx+1];
				if(istack[i] < 0) {istack[i] += m10;}
			}
		}
		jrand = 0;
		init = true;
	}

	jrand += 1;

	if(jrand > ndim)
	{
		for(i = 1; i <= ndim; ++i)
		{
			inx = i + ir - ((int)((double)(i+ir)/(double)(ndim)))*ndim;
		    istack[i] -= istack[inx+1];
			if(istack[i] < 0) {istack[i] += m10;}
		}
		jrand = 1;
	}

	ran_num = ((double)istack[jrand]) / base;

	return ran_num;
}
void initialize(struct p p[],int *next,int *prev,int ***first,int ***last,int clustlead[],int clustlast[],int clustnext[],int clustnum[])
{
	int lastp;
	int i,j,k;
	double a,O,I,Vx,Vy,Vz,x;
	FILE* input;
	



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

	//randomly place particles in box and assings cells 
	if(Start_time !=0)
	{
		input=fopen("unfold","r");
	}
	for(i=0; i<n; i++)
	{

		if(Start_time==0)
		{
			p[i].x =L*psdrand(seed); 
			p[i].y =L*psdrand(seed);
			p[i].z =L*psdrand(seed);

			if(p[i].x==L)
			{
				p[i].x=0;
			}
			if(p[i].y==L)
			{
				p[i].y=0;
			}
			if(p[i].z==L)
			{
				p[i].z=0;
			}
			//set unfold coord//
			p[i].ux=p[i].x;
			p[i].uy=p[i].y;
			p[i].uz=p[i].z;
			p[i].P=1.0000;
			
			//initalizes clustarrays
			clustlead[i]=i; clustlast[i]=i; clustnext[i]=-1; clustnum[i]=i;
		}
		else //i.e. restarting a program from unfold file
		{
			fscanf(input,"%lf %lf %lf %lf %lf %d %d %d %lf %lf %lf %lf ",&x, &p[i].ux, &p[i].uy, &p[i].uz, &x, &clustlead[i], &clustnext[i], &clustlast[i], &p[i].P, &p[i].vx, &p[i].vy, &p[i].vz); 
			// put unfolded coord back into the box   
			p[i].x=p[i].ux; 
			p[i].y=p[i].uy; 
			p[i].z=p[i].uz; 
			
			PBC(p,i);

			
		}
		
		
		
		///// sets velocity////////
		if(Start_time==0)
		{   /**Random**/
			a=2*( psdrand(seed))-1;
			O=acos(a);
			I=2*PI* (psdrand(seed));
			
			Vx=sin(O)*cos(I);
			Vy=sin(O)*sin(I);
			Vz=cos(O);

		
		
			p[i].vx=Vx;
			p[i].vy=Vy;
			p[i].vz=Vz; 
		}
		
		

		
	
		/////// set cells ////////////////////////////////////////
		p[i].cx= int( (p[i].x)/Lc );
		p[i].cy= int( (p[i].y)/Lc );
		p[i].cz= int( (p[i].z)/Lc );

		////////// initialize link list for each cell /////////////////////////////
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
	if(Start_time!=0)
	{
		fclose(input);
	}



	return;
}
void init_part_move(int i,struct p p[],int *next,int *prev,int ***first,int ***last)
{
  	int cxo,cyo,czo;
	

   
	  cxo=p[i].cx;
	  cyo=p[i].cy;
	  czo=p[i].cz;
	  p[i].x =L*psdrand(seed); 
	  p[i].y =L*psdrand(seed);
	  p[i].z =L*psdrand(seed);

	  //set unfold coord//
	  p[i].ux=p[i].x;
	  p[i].uy=p[i].y;
	  p[i].uz=p[i].z;

	  p[i].cx= int( (p[i].x)/Lc );
	  p[i].cy= int( (p[i].y)/Lc );
	  p[i].cz= int( (p[i].z)/Lc );
	 
	  link_list(p,i,cxo,cyo,czo,next,prev,first,last);

	  return;
  
}
int init_part_check(int i,struct p p[],int *next,int *prev,int ***first,int ***last)
{
    int j,k,w;
  	int nx,pr,nb;
	int cx1,cy1,cz1;
	
	int marker;
	double X,Y,Z,X1,Y1,Z1,r;

	X=p[i].x; Y=p[i].y; Z=p[i].z;

     marker=0;
    //checks the next's
    nx=next[i];
    while(1)
    {
      if(nx==-1)
      {break;}
      X1=p[nx].x; Y1=p[nx].y; Z1=p[nx].z;
      r=rcalc(X,Y,Z,X1,Y1,Z1,i,nx);
      if(r<r_min )
      {  
        marker=1;
         
		init_part_move(i,p,next,prev,first,last);

     				
      }
      if( marker==1)
       {break;}
       nx=next[nx];
    }
  
	// checks the prev's
   pr=prev[i];
   while(1)
   {
       if( marker==1)
	  {break;}
       if(pr==-1)
		{break;}
       X1=p[pr].x; Y1=p[pr].y; Z1=p[pr].z;
       r=rcalc(X,Y,Z,X1,Y1,Z1,i,pr);
       if( r<r_min)
       {
		 marker=1;
	 
		 init_part_move(i,p,next,prev,first,last);
	   }
       if(marker==1)
		{break;}
       pr=prev[pr];
    }
	

    //neiboring cells check
    for(w=-1; w<=1; w++)
    {
       if(marker==1)
	  {break;}
       for(j=-1; j<=1; j++)
       {
           if(marker==1)
	   {break;}
           for(k=-1; k<=1; k++)
           {
	      if(marker==1)
	       {break;}
	      cx1=p[i].cx+w; cy1=p[i].cy+j; cz1=p[i].cz+k;
					
	      // if statements for PBC
	      if(cx1<0)
	      {
		  cx1=int((L/Lc)-1);
	      }
	      if(cx1>( (L/Lc)-1 ))
	      {
		  cx1=0;
	      }
		
	      if(cy1<0)
	      {
		  cy1=int((L/Lc)-1);
	      }
	      if(cy1>( (L/Lc)-1 ))
	      {
		 cy1=0;
	      }

	      if(cz1<0)
	      {
		 cz1=int((L/Lc)-1);
	      }
	      if(cz1>( (L/Lc)-1 ))
	      {
		 cz1=0;
	      }
	      // so it does not count the home cell agian
	      if(cx1==p[i].cx && cy1==p[i].cy && cz1==p[i].cz)
	      { break;}
						
	      nb=first[cx1][cy1][cz1];
	      while(1)
	      {
		 if(nb==-1)
		 {break;}
		 X1=p[nb].x; Y1=p[nb].y; Z1=p[nb].z;
		 r=rcalc(X,Y,Z,X1,Y1,Z1,i,nb);
		 if( r<r_min)
		 {
		      marker=1;
		     
		      init_part_move(i, p , next, prev,first,last);

		
		 }
		  if(marker==1)
		     {break;}
		  nb=next[nb];
	      }
	   }
       }
    }

       return marker;
}
void initialcheck(struct p p[],int *next,int *prev,int ***first,int ***last)
{
	int i,w;
	
	

	int marker;
	
	
	
	w=0;
       
	//cycles through the particles
	for(i=0; i<n; i++)
	{
	        marker=0;

	       
		
		do
		  {
		    //printf("%d\n",i);
		    marker=  init_part_check( i, p , next, prev,first,last);
			if(marker==1)
			{
				w++;
			}
		  }
		while(marker==1);
	}
	printf("DONE %d \n",w);
	return;
}
void PBC(struct p p[], int i)
{
	/* PBC condition */	
	if ( (p[i].x)>=L )
	{
		p[i].x=(p[i].x)-L;
	}
	if ( (p[i].y)>=L )
	{
		p[i].y=(p[i].y)-L;
	}
	if ( (p[i].z)>=L )
	{
		p[i].z=(p[i].z)-L;
	}
		
	if( (p[i].x)<0 )
	{
		p[i].x= (p[i].x)+L;
	}
	if( (p[i].y)<0 )
	{
		p[i].y= (p[i].y)+L;
	}
	if( (p[i].z)<0 )
	{
		p[i].z= (p[i].z)+L;
	}
	//for computer rounding error
	if( p[i].x==L)
	{
		p[i].x=p[i].x-L;
	}
	if (p[i].y==L)
	{
		p[i].y=p[i].y-L;
	}
	if (p[i].z==L)
	{
		p[i].z=p[i].z-L;
	}
	//Updated cell coordinates
	p[i].cx= int( (p[i].x)/Lc );
	p[i].cy= int( (p[i].y)/Lc );
	p[i].cz= int( (p[i].z)/Lc );

	return;
}
void link_list(struct p p[], int part,int cxo,int cyo,int czo,int *next,int *prev,int ***first,int ***last)
{
	int ofirst,olast,oprev,onext,nfirst,nlast;
	
	nfirst=0;
	nlast=0;
	ofirst=0;
	olast=0;
	oprev=0;
	onext=0;


	
	if(cxo != p[part].cx || cyo != p[part].cy || czo != p[part].cz )
	{
		ofirst=first[cxo][cyo][czo];
		olast=last[cxo][cyo][czo];
		oprev=prev[part];
		onext=next[part];
		nfirst=first[p[part].cx][p[part].cy][p[part].cz];
		nlast=last[p[part].cx][p[part].cy][p[part].cz];
			
		if(nfirst==-1)
		{
			first[p[part].cx][p[part].cy][p[part].cz]=part;
			last[p[part].cx][p[part].cy][p[part].cz]=part;
			next[part]=-1;
			prev[part]=-1;
		}
		else
		{
			last[p[part].cx][p[part].cy][p[part].cz]=part;
			next[part]=-1;
			prev[part]=nlast;
			next[nlast]=part;
		}
		
		if( ofirst==part)
		{
			first[cxo][cyo][czo]=onext;
			if(onext!=-1)
			{
				prev[onext]=-1;
			}
		}
		if(olast==part)
		{
			last[cxo][cyo][czo]=oprev;
			if(oprev!=-1)
			{
				next[oprev]=-1;
			}
		}
		
		if( ofirst != part && olast !=part)
		{
			next[oprev]=onext;
			prev[onext]=oprev;
		}	
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
void iteration(struct p p[],double *time,int *Nc,int clustlead[],int clustnext[],int clustnum[],int clustlast[],int *next,int *prev,int ***first,int ***last,int unfold[])
{
	int i,j,count,NC,N,c;
	int nx,check,part;
	int cxo,cyo,czo;
	int new_part_list[1000],new_part_count;
	double Prob,prob;
	double O,a,I;
	double Vx,Vy,Vz;
	double X,Y,Z,X1,Y1,Z1,r;
	

	NC=*Nc;

	count=0;
	for(i=0; i<NC; i++)
	{

		//counts the number of clusters that have been picked this iteration 
		if(i==count)
		{
			count=count+10000;
			//printf("%d \n",count); 
		}

		// picks a random cluster
		c=(int)(psdrand(seed)*NC);
		//printf("c:%d \n",c);

		
		


		// finds cluster mass N //////////
		nx=clustnum[c];
		
		N=0;
		while(nx!=-1)
		{
			nx=clustnext[nx];
			
			N++;
		}
		//////////////////////////////////
		
		nx=clustnum[c];
		
		Prob=PROB(nx, p, clustnext);
		

		prob=psdrand(seed);
		
		
		
		*time=*time+(1.0/((double) NC));
		
		
	
		
		
		


		// here is where cluster movement happens ///////////
		if(prob<Prob)
		{

			// set velocity(a random unit vector)  of clust for Diffusional Case////////////
			
			
		  		a=2*( psdrand(seed) )-1;
				O=acos(a);
				I=2*PI* (psdrand(seed) );

				Vx=sin(O)*cos(I);
				Vy=sin(O)*sin(I);
				Vz=cos(O); 
				

			
			///////////////////////////////////////////////////////


			/// move the cluster here  ////
			nx=clustnum[c];
			while(1)
			{
				
				cxo=p[nx].cx;
				cyo=p[nx].cy;
				czo=p[nx].cz;

				
				
				p[nx].vx=Vx;
				p[nx].vy=Vy;
				p[nx].vz=Vz; 
				

				p[nx].x = p[nx].x + p[nx].vx;
				p[nx].y = p[nx].y + p[nx].vy;
				p[nx].z = p[nx].z + p[nx].vz;

				


				PBC(p,nx);
				link_list( p, nx, cxo, cyo, czo, next, prev, first, last);
				
				

				nx=clustnext[nx];
				if(nx==-1){break;}

			}

		

			for(j=0; j<1000; j++)
			{
				new_part_list[j]=-1;
			}
			new_part_count=0;
		
			check=cluster_check(c,p, next, prev, first, last, clustlead, clustlast, clustnext, clustnum, new_part_list, &new_part_count);
			if(check==2 || check==1)
			{
				check=1;
				while(check==1)
				{
					//printf("%d %d \n",new_part_count,c);
					check=cluster_check(c,p, next, prev, first, last, clustlead, clustlast, clustnext, clustnum, new_part_list, &new_part_count);
				}
			
		


				for(j=0; j<1000; j=j+2)
				{
					if(j==998)
					{
						printf("probabily going to have a seg fault \n"); 
					}
					if( new_part_list[j]!=-1 || new_part_list[j+1]!=-1 )
					{
						X=p[  new_part_list[j] ].x; Y=p[  new_part_list[j] ].y; Z=p[  new_part_list[j] ].z;
						X1=p[new_part_list[j+1]].x; Y1=p[new_part_list[j+1]].y; Z1=p[new_part_list[j+1]].z;
				
						r=rcalc(X,Y,Z,X1,Y1,Z1,j,j+1);
						if(r<=1.001)
						{
							make_clust(p,new_part_list[j],new_part_list[j+1],clustlead,clustlast,clustnext);
							if(r<0.99 || r>1.01)
							{
								printf("still have overlap with %d, %d \n",new_part_list[j],new_part_list[j+1]);
							}
						}


					}
					else{break;}
				}


				
				part=clustnum[c];
				//set all of clustnum to -1
				for(j=0; j<n; j++)
				{
					clustnum[j]=-1;
				}
				NC=0;

				//assings each clust a number and counts Nc
				for(j=0; j<n; j++)
				{
					if(j==clustlead[j])
					{
						clustnum[NC]=j;
						NC++;
					}
				}
				*Nc=NC;

				
				////////////////////////////////////////
				
				////////////////////////////////////////

			}
		}
		
		
	
	
	
	
	
	
	}

    //free all the memory allocations
	return;


}
int cluster_check(int c,struct p p[], int *next, int *prev, int ***first, int ***last, int clustlead[], int clustlast[], int clustnext[], int clustnum[],int new_part_list[], int *new_part_count)
{
	int i,j,k;
	int nx,pr,nb;
	int cx1,cy1,cz1;
	int part; 
	int ret_val;
	int count;
	double X,Y,Z,X1,Y1,Z1,r;


	
	count=*new_part_count;

	part=clustnum[c];
	ret_val=0;
	
	while(part!=-1)
	{
		
		X=p[part].x; Y=p[part].y; Z=p[part].z;

		
		
		//// checks part's next  //////////////////////////////////////////////////////////////////
		nx=next[part];
		while(nx!=-1)
		{

			// this "if" is to make sure to only check monmomers not in the cluster
			if(clustlead[part]!=clustlead[nx])
			{
				X1=p[nx].x; Y1=p[nx].y; Z1=p[nx].z;
				r=rcalc(X,Y,Z,X1,Y1,Z1,part,nx);
				if(r<=r_max)
				{
					ret_val=2;

					new_part_list[count]=part; new_part_list[count+1]=nx;
					
					count=count+2;
					
					//deal with particle overlap here
					if(r<r_min)
					{
						overlap(part,nx,r,p,clustlead,clustlast,clustnext,next,prev,first,last);
						
						X1=p[nx].x; Y1=p[nx].y; Z1=p[nx].z;
						X=p[part].x; Y=p[part].y; Z=p[part].z;
						r=rcalc(X,Y,Z,X1,Y1,Z1,part,nx);
						ret_val=1;
					}
					
					if(r<r_min)
					{
						printf("r-move error: r=%0.3lf i=%d j=%d \n",r,part,nx);
						printf("");
					}
					////////////////////////////////

					

				}
			}
			
			
			nx=next[nx];
		}
		//////////////////////////////////////////////////////////////////////////////////////////////////////////////////

		//// checks part's prev  //////////////////////////////////////////////////////////////////
		pr=prev[part];
		while(pr!=-1)
		{

			// this "if" is to make sure to only check monmomers not in the cluster
			if(clustlead[part]!=clustlead[pr])
			{
				X1=p[pr].x; Y1=p[pr].y; Z1=p[pr].z;
				r=rcalc(X,Y,Z,X1,Y1,Z1,part,pr);
				if(r<=r_max)
				{
					if(ret_val!=1)
					{
						ret_val=2;
					}

					new_part_list[count]=part; new_part_list[count+1]=pr;
					
					count=count+2;
					
					//deal with particle overlap here
					if(r<r_min)
					{
						overlap(part,pr,r,p,clustlead,clustlast,clustnext,next,prev,first,last);
						
						X1=p[pr].x; Y1=p[pr].y; Z1=p[pr].z;
						X=p[part].x; Y=p[part].y; Z=p[part].z;
						r=rcalc(X,Y,Z,X1,Y1,Z1,part,pr);
						ret_val=1;
					}
					
					if(r<r_min)
					{
						printf("r-move error: r=%0.3lf i=%d j=%d \n",r,part,pr);
						printf("");
					}
					////////////////////////////////////
				}
			}
			
			
			pr=prev[pr];
		}
		//////////////////////////////////////////////////////////////////////////////////////////////////////////////////

		/// checks the neiboring cells /////////////////////////////////////////////////
		for(i=-1; i<=1; i++)
		{
			for(j=-1; j<=1; j++)
			{
				for(k=-1; k<=1; k++)
				{

					cx1=p[part].cx +i; cy1=p[part].cy +j; cz1=p[part].cz +k;
					// if statements for PBC
					if(cx1<0)
					{
						cx1=int((LL)-1);
					}
					if(cx1>( (LL)-1 ))
					{
						cx1=0;
					}
		
					if(cy1<0)
					{
						cy1=int((LL)-1);
					}
					if(cy1>( (LL)-1 ))
					{
						cy1=0;
					}

					if(cz1<0)
					{
						cz1=int((LL)-1);
					}
					if(cz1>( (LL)-1 ))
					{
						cz1=0;
					}

					// so it does not count the home cell agian
					if(cx1==p[part].cx && cy1==p[part].cy && cz1==p[part].cz)
					{ 
					}
					else
					{
						
						
						//// here is the actual check ////////////////////////////////
						nb=first[cx1][cy1][cz1];
						while(nb!=-1)
						{
							if(clustlead[part]!= clustlead[nb])
							{
								X1=p[nb].x; Y1=p[nb].y; Z1=p[nb].z;
								r=rcalc(X,Y,Z,X1,Y1,Z1,part,nb);
								if(r<=r_max)
								{
									if(ret_val!=1)
									{
										ret_val=2;
									}

									new_part_list[count]=part; new_part_list[count+1]=nb;
					
									count=count+2;
					
									//deal with particle overlap here
									if(r<r_min)
									{
										overlap(part,nb,r,p,clustlead,clustlast,clustnext,next,prev,first,last);
										
										X1=p[nb].x; Y1=p[nb].y; Z1=p[nb].z;
										X=p[part].x; Y=p[part].y; Z=p[part].z;
										r=rcalc(X,Y,Z,X1,Y1,Z1,part,nb);
										ret_val=1;
									}
					
									if(r<r_min)
									{
										printf("r-move error: r=%0.3lf i=%d j=%d \n",r,part,nb);
										printf("");
									}
									////////////////////////////////////
								}


							}
							nb=next[nb];

						}

					   /////////////////////////////////////////////////////////////////////////


					}



				}
			}
		}
		//////////////////////////////////////////////////////////////////////////////////////////////

		part=clustnext[part];
	}

	
	*new_part_count=count;
	return ret_val;
}
void overlap(int part,int j,double r,struct p p[],int clustlead[],int clustlast[],int clustnext[],int *next,int *prev,int ***first,int ***last)
{
	double dx,dy,dz;
	double X,Y,Z,X1,Y1,Z1;
	double m,k,l;
	int nx;
	int cxo,cyo,czo;
	
	X=p[part].x; Y=p[part].y; Z=p[part].z;
	X1=p[j].x; Y1=p[j].y; Z1=p[j].z;

	dx=X1- X;
	dy=Y1- Y;
	dz=Z1- Z;

	// PBC check/////////////////
	if(abs((int)dx)>L/2)
	{
		if(X1>X)
		{
			X1=X1-L;
			dx=X1-X;
		}
		else
		{
			X1=X1+L;
			dx=X1- X;
		}
	}
	
	
	if(abs((int)dy)>L/2)
	{
		if(Y1>Y)
		{
			Y1=Y1-L;
			dy=Y1-Y;
		}
		else
		{
			Y1=Y1+L;
			dy=Y1- Y;
		}
	}
	
	if(abs((int)dz)>L/2)
	{
		if(Z1>Z)
		{
			Z1=Z1-L;
			dz=Z1-Z;
		}
		else
		{
			Z1=Z1+L;
			dz=Z1- Z;
		}
	}
///////////////////////////////////

	k=2.0*(dx*p[part].vx + dy*p[part].vy + dz*p[part].vz);
	
	l=-4.0*(-1.0 + dx*dx + dy*dy +dz*dz); l=l+k*k;
	
	m=0.5*(-1.0*k + sqrt(l));

	dx=m*p[part].vx;
	dy=m*p[part].vy;
	dz=m*p[part].vz;

	
	// move the cluster 
	nx=clustlead[part];

	while(nx!=-1)
	{
		
		cxo=p[nx].cx;
		cyo=p[nx].cy;
		czo=p[nx].cz;

		p[nx].x=p[nx].x -dx;
		p[nx].y=p[nx].y -dy;
		p[nx].z=p[nx].z -dz;

		

		PBC(p,nx);
		link_list(p,nx,cxo,cyo,czo,next,prev,first,last);

		nx=clustnext[nx];

	}
	//
	
	//

	return;
}
void make_clust(struct p p[], int part, int j, int clustlead[], int clustlast[], int clustnext[])
{
	int V,B,l;
	double A,O,I; 
	double Vx,Vy,Vz;


	// set velocity of new clust for ballistic only  
	
	A=2*( psdrand(seed) )-1;
	O=acos(A);
	I=2*PI* (psdrand(seed) );

	Vx=sin(O)*cos(I);
	Vy=sin(O)*sin(I);
	Vz=cos(O);
	
	////////////////////////////////////////////////
	
	// if statement makes sure they are not clustered already
	if(clustlead[part]!=clustlead[j])
	{
		
		
		//sets the clustlast of new clust
		B=clustlast[part];
			
		V=clustlead[part];
		while(V!=-1)
		{
			clustlast[V]=clustlast[j];
			V=clustnext[V];		
		}
		//bridges the clustlast of clust1 to the clustlead of clust2
		clustnext[B]=clustlead[j];

		// sets the clustlead of new clust 
			
		V=clustlead[j];
		while(V!=-1)
		{
			clustlead[V]=clustlead[part];
			V=clustnext[V];
		}
		

		l=clustlead[part];
		while(1)
		{
			if(l==-1)
			{break;}
			p[l].vx=Vx;
			p[l].vy=Vy;
			p[l].vz=Vz;
			l=clustnext[l];
		} 
	}

	return;
	
}
void Kinetics(struct p p[],int clustnum[],int clustlead[],int clustnext[],double *avRg, int Rg_check, int Nc, double *avRnn, double time)
{
	int i,part,nx,count,mass;
	double xc,yc,zc,rg;
	double *Xc, *Yc, *Zc, *RG, *Mass;
	//double rnn, r, avRn;
	//int j,k;
	double r_cm[3]={};

	char filename[30];
	FILE* Mass_of_clust;

	Xc=(double *)malloc(Nc*sizeof(double));                      //  
	if(Xc== NULL)                                                   //  
	{                                                              // 
		printf("Xc");                                               //  
	}   
	Yc=(double *)malloc(Nc*sizeof(double));                      //  
	if(Yc== NULL)                                                   //  
	{                                                              // 
		printf("Yc");                                               //  
	}  
	
	Zc=(double *)malloc(Nc*sizeof(double));                      //  
	if(Zc== NULL)                                                   //  
	{                                                              // 
		printf("Zc");                                               //  
	} 

	RG=(double *)malloc(Nc*sizeof(double));                      //  
	if(RG== NULL)                                                   //  
	{                                                              // 
		printf("RG");                                               //  
	}  
	
	Mass=(double *)malloc(Nc*sizeof(double));                      //  
	if(Mass== NULL)                                                   //  
	{                                                              // 
		printf("mass");                                               //  
	}  
	

	for(i=0;i<Nc; i++)
	{
		RG[i]=0.75;
	}


	count=0;
	
	if (mode2==0)
	{
		sprintf(filename, "time_%lf",time);
		Mass_of_clust=fopen(filename,"w");
	}

	for(i=0;i<Nc; i++)
	{
		
		
		part=clustnum[i];

		if(part==-1){break;}
		
		xc=0; yc=0; zc=0; mass=1; rg=0; 
		
		
		nx=part;
		mass=R_CM(nx, r_cm, clustnext, p );


		xc=r_cm[0];
		yc=r_cm[1];
		zc=r_cm[2];
		
		
		

		Mass[i]=mass;
		
		
		Xc[i]=xc;
		Yc[i]=yc;
		Zc[i]=zc;
		
		// <Rg> ..
		nx=part;

		if(mass>1)
		{
			rg=Rg( nx, p , clustnext);

			*avRg=*avRg+rg;
				
			RG[i]=rg;
				
			count++;
		}
		else
		{
			RG[i]=0.75; //for monomer = root(3/5)/a
		}
		
	if (mode2==0)
		{
			fprintf(Mass_of_clust,"%lf %d %lf \n",time,mass,RG[i]);// note rg for monomer is .75
		}

	
	}
	

	
	*avRg=*avRg/count;
	
	
	/*
	//<Rnn>...................................................................
	avRn=0; avla_Rn=0;
	for(i=0; i<Nc; i++)
	{
		
		
		rnn=10*L;
			
		for(j=0; j<Nc; j++)
		{
			if( i != j)
			{
				r=rcalc(Xc[i], Yc[i], Zc[i], Xc[j], Yc[j], Zc[j], i, j);
					
				if(r<rnn)
				{
					rnn=r;
					k=j;
						
				}
					
			}

		}
			
		avRn=avRn+rnn;
		
		
		
		/////////////////////////////////////////
	}
	avRn=avRn/Nc;
	*avRnn=avRn;

	*/


	free(Xc); free(Yc); free(Zc); free(RG); free(Mass);
	//........................................................................


	if (mode2==0)
	{
		fclose(Mass_of_clust);
	}

	return;
}
void clust_unfold(int l,struct p p[],int clustnext[],int f[])
{
	int nx,nb;
	int cxl,cyl,czl,cxn,cyn,czn;
	int cxp,cxm,cyp,cym,czp,czm;
	double X,Y,Z,X1,Y1,Z1;
	double r;
	p[l].ux=p[l].x; p[l].uy=p[l].y; p[l].uz=p[l].z; f[l]=1;
	nb=l;
	//loop cycles through paritcles in cluster
	while(nb!=-1)
	{
		

		X=p[nb].x; Y=p[nb].y; Z=p[nb].z; cxl=p[nb].cx; cyl=p[nb].cy; czl=p[nb].cz; 
		cxp=cxl+1; cxm=cxl-1; 
		cyp=cyl+1; cym=cyl-1; 
		czp=czl+1; czm=czl-1;

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
			
				
			X1=p[nx].x; Y1=p[nx].y; Z1=p[nx].z; cxn=p[nx].cx; cyn=p[nx].cy; czn=p[nx].cz;
			// checks if 'l' and 'nx' are in same or neigoring cell 
				

				
			if( (cxn==cxl||cxn==cxp||cxn==cxm) && (cyn==cyl||cyn==cyp||cyn==cym) && (czn==czl||czn==czp||czn==czm))
			{
				r=rcalc(X,Y,Z,X1,Y1,Z1,nb,nx);
				if(r<=(4*Lc ) )
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
		nb=clustnext[nb];
	}
	
	return;
}
int R_CM(int clust_lead, double r_cm[],int clustnext[],struct p p[])
{
		double xc,yc,zc,x,y,z;
		int nx,mass;


		xc=0; yc=0; zc=0; mass=1; 
		
		
		nx=clust_lead;
		xc=p[nx].x;
		yc=p[nx].y;
		zc=p[nx].z;
		// center of mass..........................................
		while(nx!=-1)
		{
			
			nx=clustnext[nx];
			if(nx==-1)
			{break;}

			x=p[nx].x;
			y=p[nx].y;
			z=p[nx].z;

			if(fabs(x-xc)>L/2 && x>xc)
			{
				x=x-L;
			}
			if(fabs(x-xc)>L/2 && x<xc)
			{
				x=x+L;
			}

			if(fabs(y-yc)>L/2 && y>yc)
			{
				y=y-L;
			}
			if(fabs(y-yc)>L/2 && y<yc)
			{
				y=y+L;
			}

			if(fabs(z-zc)>L/2 && z>zc)
			{
				z=z-L;
			}
			if(fabs(z-zc)>L/2 && z<zc)
			{
				z=z+L;
			}


			xc=mass*xc+x;
			yc=mass*yc+y;
			zc=mass*zc+z;

			mass++;

			xc=xc/(1.0*mass);
			yc=yc/(1.0*mass);
			zc=zc/(1.0*mass);
			
		}

		r_cm[0]=xc;
		r_cm[1]=yc;
		r_cm[2]=zc;

	return mass;
}
double Rg(int clust_lead,struct p p[], int clustnext[])
{
	int m,nx; 
	double rg, dx, dy, dz,r_cm[3]={};

	R_CM( clust_lead, r_cm, clustnext, p);

	m=0; rg=0;

	nx=clust_lead;
	while( nx !=-1)
	{
		dx=fabs(p[nx].x - r_cm[0]);
		if(fabs(dx-L) < dx)
		{
			dx=fabs(dx-L);
		}


		dy=fabs(p[nx].y - r_cm[1]);
		if(fabs(dy-L) < dy)
		{
			dy=fabs(dy-L);
		}

		dz=fabs(p[nx].z - r_cm[2]);
		if(fabs(dz-L) < dz)
		{
			dz=fabs(dz-L);
		}


		rg=rg + dx*dx  + dy*dy + dz*dz ;
		m++;
		
		nx=clustnext[nx];
	}

	rg=rg/(1.0*m); 
	rg=sqrt(rg); 


	return rg; 
}
double PROB(int part,struct p p[],int clustnext[])
{
	int nx,N; 
	double xc,yc,zc,r_cm[3]={};
	double h,x,Prob;

	
	
	// finds cluster mass N //////////
	nx=part;
		
	N=0;
	while(nx!=-1)
	{
		nx=clustnext[nx];	
		N++;
	}
	//////////////////////////////////
	

	
	

	
	//finds center mass
	xc=0; yc=0; zc=0;
	

	nx=part;
	R_CM(nx, r_cm, clustnext, p);

	xc=r_cm[0]; 
	yc=r_cm[1]; 
	zc=r_cm[2];
		
	Prob=0;
	while(nx!=-1)
	{
		h=rcalc(xc,yc,zc,p[nx].x,p[nx].y,p[nx].z,part,nx);
		h=h*h;
		x=(20.0/3.0)*h;
		Prob=Prob+( x )+1;
			
		nx=clustnext[nx];
			
	}

	Prob=Prob/N;
	Prob=sqrt(Prob);
	Prob=1/Prob;

	//sets prob
	nx=part;
	while(nx!=-1)
	{
		p[nx].P=Prob;
		nx=clustnext[nx];
	}


	return Prob;
}
void print_files(int Rg_check, struct p p[],int clustlead[], int clustlast[], int clustnext[], double time, int Nc, double avRg)
{

	int i; 

	
	FILE* kinetics;
	FILE* output;

	if(Rg_check==1)
	{
		output=fopen("unfold","w");
		for(i=0; i<n; i++)
		{
			fprintf(output, "%d %lf %lf %lf %lf %d %d %d %lf %lf %lf %lf \n",i, p[i].ux, p[i].uy, p[i].uz, 0.5, clustlead[i], clustnext[i], clustlast[i], p[i].P, p[i].vx, p[i].vy, p[i].vz ); 
		}
		fclose(output);

		output=fopen("fold","w");
		for(i=0; i<n; i++)
		{
			fprintf(output, "%d   %7.3lf    %7.3lf    %7.3lf    %7.3lf    %d    %d   %d  %lf %lf %lf %lf   \n",i, p[i].x, p[i].y, p[i].z, 0.5, clustlead[i], clustnext[i], clustlast[i], p[i].P, p[i].vx, p[i].vy, p[i].vz ); 
		}
		fclose(output);


	}

	kinetics=fopen("time","a");
				
				
	fprintf(kinetics,"%lf %d %lf \n",time,Nc,avRg);
				
				
	fclose(kinetics);


				
				///////////////////////////////////////////////////////////////////////////////////////////////////////////////
}
void Unfold(int unfold[], struct p p[], int clustnext[], int clustnum[], int Nc)
{

	int j,a,A,part,nx;
	/// here is whrer unfloding happens/////////////////////////////////////////////////////////////////////////
	for(j=0; j<n; j++)
	{
		unfold[j]=-1;
	}

	for(j=0; j<n; j++)
	{
		p[j].ux=p[j].x;
		p[j].uy=p[j].y;
		p[j].uz=p[j].z;

	}
	/*
				
	for(j=0; j<Nc; j++)
	{
		part=clustnum[j];
		if(part==-1){break;}
		a=0;
		while(1)
		{	
			A=0;
			//unfolding(p,clustlead,clustnext,clustnum,unfold);
			clust_unfold(part, p ,clustnext, unfold);
			nx=part;
			while(nx!=-1)
			{
				if(unfold[nx]==-1)
				{
					A++;
				}
				nx=clustnext[nx];
			}
			if(A==0)
			{break;}
			
			if (a==A)
			{
				printf("unfold error(  PRATICLE: %d) \n\n\n\n,",j);
				break;
			}
			a=A;
						
			
		}
	//printf("unfold_%d \n",j);
	}

	*/

}

