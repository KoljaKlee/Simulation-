 /*
  * Simple code for the Langevin dynamics simulation of dimers in 2D based on code provided by Jaeoh Shin
  * 
 */


#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <sys/time.h>
#include <math.h>
#include <string.h>
#include "/Volumes/Macintosh HD/Users/KoljaSchule/Studium/Bachelor/Progamme/Zufallszahl/mt19937-64.c" //path of the random number generator


//system setup
#define tple  2   //number of particles (monomers), only even numbers
#define config_capture 1 // 1=store VMD data



// Control the simulation
int pbc = 1;			//1=periodic boundary condition
int steps= 1e8;			//simulation integration steps
int capture1=2e2;		//interval of data saving



/**
 **  parameters for LJ and harmonic
**/

static int tple2 = tple/2; //number of dimers
double rcut,rcut2; //cutoff distance of LJ potential and its square
static double sigma=1; //particle size [sigma]=1nm
double box; //we use a periodic boundary condition. this is box size
static double T=1.; //temperature
double eps; //strength of LJ potential
double disteq; // harmonic equilibrium distance
double ks; //spring constant 



/**
* neighbour list and related arrays: Reference, M.P. Allen and D. J. Tildesley (1987)
**/


short int nlist[tple][tple];
double displacement[tple]; //displacement after last nlist update
double rskin, rskin2; //neighbor distance
int nlist_updates=0;

int roundToint(double x);

/**
 * Update list
 **/
int update_nlist(double *x, double *y);


/**
 * Check whether updating neighbor list is needed or not
 * */
int check_nlist();


/**
 * Force: Lenneard-Jones
 **/

int lj_force(double *x, double *y, double *fljx, double *fljy);

// Harmonic force
int h_force(double *x, double *y, double *fhx, double *fhy);


// changes the angle based on changes of the coordinates
int dangle(double *xt, double *xto, double *yt, double *yto, double *theta);

//Set coordinate of the center of mass
int CoM(double *xt, double *yt, double *CoMx, double *CoMy);

/** Makes two Gaussian random numbers with zero mean and unit variance.
 * Argument is double array with two elements. **/
void gaussrand(double* grnd);


/** 
 * Generation of random number seed, get a current time in micri-second unit
 * */
int gus() {
	struct timeval tp;
	gettimeofday(&tp,0);
	return tp.tv_usec;
}


int main(){
	
	double dt=0.005 ;  //integration time step
	unsigned long long idum;
	idum=gus();
	init_genrand64(idum);
	time_t start, end ;
	time(&start);
	int simulation_time ;
	int is_overlapped ;	
	int Nint ;
    
	
	int i, j, k ;
    

//below open a few files to save data duing simulation.

	FILE *infoFile ;
	FILE *eqFile;     
		eqFile=fopen("./eqFile.txt","w");
	FILE *angleFile ;
		angleFile=fopen("./angle.txt","w");//angles
    FILE *dFile ;
    dFile=fopen("./dis.txt","w");//distances
	FILE *xFile ;
		xFile=fopen("./coorx.txt", "w"); //coordinates of tracer particle
    FILE *yFile ;
        yFile=fopen("./coory.txt", "w");
	FILE *config;
	  config=fopen("./ini_config.txt","w");	
	


/****
 system parameters
*****/


	//rcut=1.12246205*sigma; //=2^(1/6)*sigma, repulsive LJ
    	rcut = 2*sigma; // attractive + repulsive LJ
	rcut2=rcut*rcut;
	
	rskin=4*sigma; //Verlet-list skin radius. ~4.0 should be fine.
	rskin2=rskin*rskin;
    
    
    eps=1;//set to unity
    box=20*sigma; // size of the box

    
	double grnd[2]; // two gaussian random numbers
	grnd[0]=0.0;
	grnd[1]=0.0;
	
	

	
    	double m=1;
	double gamma=1;  //friction coefficient set to unity
	int need_update_nlist=0;

	double x[tple], y[tple] ;// position
    double xt[tple], yt[tple]; // unbound positions
    double xto[tple], yto[tple]; //previous coordinates
    
    double vx[tple], vy[tple] ; //particle velocities
    
	double xpcm, ypcm; //coordinates of the origin
    
	double dx, dy; // changes to the coordinates
    
	double fx[tple], fy[tple] ; //forces
	double fx1[tple], fy1[tple] ; //previous forces
	
	double fljx[tple],fljy[tple]; //LJ force
    double fhx[tple], fhy[tple];//harmonic force
	
	double drx[tple], dry[tple]; //position noise
	double drvx[tple], drvy[tple];//velocity noise
    
    double theta[tple2];// angle between harmonic rod and x axis
  
    double CoMx[tple2];//center of mass coordinates
    double CoMy[tple2];
    
    double d[tple2]; //relative coordinate
    double reld[tple2];// d/disteq
    
	double vmx, vmy, temp1 ;
	double fs;	
	double xcm;
	
    
	
/** The random displacements drx and drvx are correlated Gaussian random numbers.
     They can be obtained from zero-mean unit-variance distribution by transformation
     drx[i]=cdr1*grnd[0]+cdr2*grnd[1];
     drvx[i]=cdrv*grnd[0];
     given in, e.g., Allen-Tildesley Sec 9.3 and Appendix G.3.**/

	double c0, c1, c2, cx1, cx2, cv1, cv2, cv3 ;
	double varr, varv, covrv;
	double cdr1, cdr2, cdrv;


/**
***  setting up the coefficients 
**/

		c0=exp(-1.0*gamma*dt);
		c1=(1.0-c0)/(gamma*dt);
		c2=(1.0-c1)/(gamma*dt);
		
		cx1=c1*dt; //= (dt-dt*c0)/(gamma*dt) = (1-c0)/gamma
		cx2=c2*dt*dt/m; // = dt/m *(dt-dt*c1)/(gamma*dt)= dt/m*(1-c1)/gamma= 1/(m*gamma)*(dt-dt*c1)= (dt-cx1)/(m*gamma)
		
		cv1=c0;\mathrm
		cv2=(c1-c2)*dt/m; //= (c1-(1-c1)/gamma*dt)*dt/m = (c1*gamma*dt-1+c1)/(gamma*m) = (c1*dt*(gamma+1/dt)-1)/(gamma*m) = (cx1*(gamma+1/dt)-1)/(gamma*m)
		cv3=c2*dt/m; // = (1-cx1/dt)/(m*gamma)
		
		varr=T*dt/(gamma*m)*(2.0-(3.0-4.0*exp(-1.0*gamma*dt)+exp(-2.0*gamma*dt))/(gamma*dt));
		varv=T/m*(1.0-exp(-2.0*gamma*dt));
		covrv=T/(gamma*m)*(1.0-exp(-1.0*gamma*dt))*(1.0-exp(-1.0*gamma*dt))/(sqrt(varr)*sqrt(varv));
		
		cdr1=sqrt(varr)*covrv;
		cdr2=sqrt(varr)*sqrt(1.0-covrv*covrv);
		cdrv=sqrt(varv);
		
	
	   

/** start by distribute particles inside the box **/
    
    int kmax = box/(2.01*sigma+0.01); //number of dimers in a line
    int imax = box/(1.01*sigma+0.01); //number of lines
    j=0;
    
    for(i=0;i<imax;i++){
        if(j<tple2){
            for(k=0;k<kmax;k++){
                if(j<tple2){
                    x[j]=(-box/2)+2.01*sigma*k+0.01;
                    y[j]=(-box/2)+1.01*sigma*i+0.01;
                    x[j+tple2]=x[j]+sigma;
                    y[j+tple2]=y[j];
                    j++;
                }
            }
        }
    }
    
    
    
    for(i=0;i<tple;i++){
   fprintf(config, "%4.3f\t%4.3f\n", x[i],y[i]);}
   fclose(config); 
  
    
    for(j=0;j<tple;j++){ //save the coordinates before they are changed
        xto[j]=x[j];
        yto[j]=y[j];
    }
    
    for(i=0;i<tple2;i++){//set the angle to zero
         theta[i]=0.;
    }
    
 
  
/** initialize the velocity **/

	vmx=0., vmy=0. ;
	for(j=0;j<tple;j++){
		vx[j]=(genrand64_real2()-0.5); //random velocities
		vy[j]=(genrand64_real2()-0.5);
		vmx=vmx+vx[j];
		vmy=vmy+vy[j];
	}
	
	vmx=vmx/(double)tple; //average velocity
	vmy=vmy/(double)tple;
	temp1=0.;
	
	for(j=0;j<tple;j++){ //calculate the variance of the velocities
		vx[j]=vx[j]-vmx; 
		vy[j]=vy[j]-vmy;
		temp1=temp1 +vx[j]*vx[j] +vy[j]*vy[j];	} 

	temp1=temp1/(2.0*(double)tple); 
	fs=sqrt(T/temp1); 
	temp1=0.;

	for(j=0;j<tple;j++){ 
		vx[j]=fs*vx[j]; // new velocities = sqrt(T/variance)* old velocities 
		vy[j]=fs*vy[j];
		temp1=temp1+vx[j]*vx[j]+vy[j]*vy[j]; 
	}
				
	update_nlist(x,y );
    
    for(j= 0; j<tple; j++){
        xt[j]=x[j];// initialize unbounded coordinates
        yt[j]=y[j];
    }
    
	lj_force(x, y, fljx, fljy) ; // calculate forces
        h_force(x,y, fhx, fhy);

	/**** time integration of the Langevin equation *******/			
	for (i=0 ; i< steps ; i++)	{ 

        xpcm=0.0; // origin of the box
	ypcm=0.0;
        
		for(j=0;j<tple;j++){ 
		  if(pbc){
			  	x[j]=x[j]-roundToint((x[j]-xpcm)/box)*(double)box; //applies perodic boundary conditions
				y[j]=y[j]-roundToint((y[j]-ypcm)/box)*(double)box;  }
				      }

		for(j=0 ; j< tple; j++) {
			gaussrand(grnd);
			drx[j]=cdr1*grnd[0]+cdr2*grnd[1]; //calculate noise
			drvx[j]=cdrv*grnd[0];
			
			gaussrand(grnd);
			dry[j]=cdr1*grnd[0]+cdr2*grnd[1];
			drvy[j]=cdrv*grnd[0];}
			
			
		// New value for x and y
		for(j=0;j<tple;j++){
            fx[j]=fljx[j]+fhx[j]; //add forces
            fy[j]=fljy[j]+fhy[j];

			dx=cx1*vx[j]+cx2*fx[j]+drx[j]; 
			x[j]=x[j]+dx;
            		xt[j]= xt[j]+dx; //change unbounded coordinates
			
			dy=cx1*vy[j]+cx2*fy[j]+dry[j];
			y[j]=y[j]+dy;
			yt[j]= yt[j]+dy;//change unbounded coordinates
            
			displacement[j]=displacement[j]+sqrt(dx*dx+dy*dy);
            
						}
        for(j=0;j<tple2;j++){ //determine the relative distance
            d[j]=sqrt((xt[j]-xt[j+tple2])*(xt[j]-xt[j+tple2])+(yt[j]-yt[j+tple2])*(yt[j]-yt[j+tple2]));
            reld[j]=d[j]/(1.5*sigma);
        }
        
        CoM(xt,yt,CoMx,CoMy);//change the center of mass
        
        dangle(xt,xto,yt,yto,theta);//update the orientation
        
        
	//update neighbor list if the displacement is larger than rskin-rcut
		need_update_nlist=check_nlist();
		if(need_update_nlist) {
			update_nlist(x,y);
			}
			
        lj_force(x, y, fljx, fljy) ;//calculate forces
        h_force(x,y, fhx, fhy);
        
	for(j=0;j<tple;j++){
        fx1[j]=fljx[j]+fhx[j];
        fy1[j]=fljy[j]+fhy[j];

			vx[j]=cv1*vx[j]+cv2*fx[j]+cv3*fx1[j]+drvx[j];// change velocities
			vy[j]=cv1*vy[j]+cv2*fy[j]+cv3*fy1[j]+drvy[j];
		}
        
        //data collection, starts after 2000 relaxation steps
        
	if((i%capture1)==0 && (i>capture1*2000)){//center of mass x
        
        fprintf(xFile, " %.2f\t",dt*i);
        
        for(j=0;j<tple2;j++){
            fprintf(xFile, " %.3f\t", CoMx[j]);
            }
         fprintf(xFile,"\n");
        
    }
    
        if((i%capture1)==0 && (i>capture1*2000)){//center of mass <
            
        fprintf(yFile, " %.2f\t",dt*i);
        
        for(j=0;j<tple2;j++){
            fprintf(yFile, " %.3f\t", CoMy[j]);
        }
         fprintf(yFile,"\n");
    }
        
    if((i%capture1)==0 && (i>capture1*2000)){//angle
        
         fprintf(angleFile, " %.2f\t",dt*i);
            
         for(j=0;j<tple2;j++){
            fprintf(angleFile, " %.3f\t", theta[j]);
         }
        
         fprintf(angleFile,"\n");
            
    }
        
    if((i%capture1)==0 && (i>capture1*2000) ){//distance
            
        fprintf(dFile, " %.2f\t",dt*i);
        
        for(j=0;j<tple2;j++){
            fprintf(dFile, " %.3f\t", d[j]);
        }
        
        fprintf(dFile,"\n");
    }
		
		
//  Storing particle coordinate and velocity for VMD (=visulization program.)

	if ( ((i%(capture1))==0) && (config_capture) && (i<capture1*7000) ){
		
			 //Output format for VMD.
			fprintf(eqFile, "ITEM: TIMESTEP\n"); 
			fprintf(eqFile,"%d\n", (int)i); 
			fprintf(eqFile,"ITEM: NUMBER OF ATOMS\n");
			fprintf(eqFile,"%d\n", (int)tple); 
			fprintf(eqFile,"ITEM: BOX BOUNDS pp pp pp\n");
			fprintf(eqFile,"%d %d\n", 0, (int)box);
			fprintf(eqFile,"%d %d\n", 0, (int)box);
			fprintf(eqFile,"%d %d\n", 0, (int)box);
			fprintf(eqFile, "ITEM: ATOMS id type xu yu zu vx vy vz \n");
										
			for(j=0;j<tple;j++){
				if(j<tple/2)
					fprintf(eqFile,"%d %d %.3f %.3f %.3f %.3f %.3f %.3f\n", (int)(j), 1, x[j], y[j], 0., vx[j], vy[j], 0. ); //Koordinaten und Geschwindigkeiten
				else	
				     fprintf(eqFile,"%d %d %.3f %.3f %.3f %.3f %.3f %.3f\n", (int)(j), 2, x[j], y[j], 0., vx[j], vy[j], 0.);
								}
										}
											
												
	} // integration steps end here
	
	fclose(xFile);
    fclose(yFile);
	fclose(eqFile);
	fclose(angleFile);
    fclose(dFile);
    
	infoFile = fopen("./info.txt", "a"); //simulation time
 	time(&end);
 	simulation_time = (int)difftime(end, start);
	fprintf(infoFile, "Simulation time= %d minutes %d sec\n", (int)(simulation_time/60),(simulation_time%60));
	fclose(infoFile);
	
	return 0;
	
 }
  
                 /************** End of main ****************/







int lj_force(double *x, double *y, double *fljx, double *fljy){  /** lennard jones force**/
	int j, k;
	double xdist, ydist, dist2,dist2i,dist6i, lj_factor, sigma2, sigma6 ;
    sigma2=sigma*sigma;
	for(j=0;j<tple;j++){
		fljx[j]=0.0;        
		fljy[j]=0.0;   }


	for (j=0;j<tple-1;j++) {
		for(k=j+1;k<tple;k++){
			
		if(nlist[j][k]==1) {   //check the neighbor list
			xdist=x[j]-x[k];  // interaction only between near neighbours 
			ydist=y[j]-y[k];
			
			if(pbc){
				  xdist=xdist-roundToint(xdist/box)*(double)box; 
				  ydist=ydist-roundToint(ydist/box)*(double)box; }
			
			dist2=xdist*xdist +ydist*ydist; 

			if(dist2<rcut2){ 
                    sigma6 = sigma2*sigma2*sigma2;
                    dist2i=1./dist2;
					dist6i=dist2i*dist2i*dist2i;
					lj_factor=48.0*eps*dist2i*dist6i*sigma6*(sigma6*dist6i-0.5); // calculate the strength of the interaction force
				
					fljx[j]=fljx[j]+lj_factor*xdist; 
					fljy[j]=fljy[j]+lj_factor*ydist;
					fljx[k]=fljx[k]-lj_factor*xdist;
					fljy[k]=fljy[k]-lj_factor*ydist;
						}
			}
		}
	}
	return 0;
}


//Harmonic force

int h_force(double *x, double *y, double *fhx, double *fhy){
    int j;
    double xdist, ydist, dist, disti, h_factor,sigma2 ;
    sigma2 = sigma*sigma;
    disteq = 1.5*sigma;
    ks = 100;//=100eps
    
    for(j=0;j<tple2;j++){
        fhx[j]=0.0;        
        fhy[j]=0.0;
        fhx[j+tple2]=0.0;        
        fhy[j+tple2]=0.0;
    
    xdist=x[j]-x[j+tple2]; //Force between the particle and its partner tple2 away
    ydist=y[j]-y[j+tple2];
    
    if(pbc){
        xdist=xdist-roundToint(xdist/box)*(double)box; 
        ydist=ydist-roundToint(ydist/box)*(double)box; }
    
    dist=sqrt(xdist*xdist +ydist*ydist);
    disti = 1/dist;
    h_factor =-ks*(1-disteq*disti);//calculate the interaction force
    
    fhx[j]=h_factor*xdist;
    fhy[j]=h_factor*ydist;
    fhx[j+tple2]=-h_factor*xdist;
    fhy[j+tple2]=-h_factor*ydist;
    }
    return 0;
}



//Update theta

int dangle(double *xt, double *xto, double *yt, double *yto, double *theta){
    int j, angle_sign;
    double xdist, ydist, dist, dxo, dyo, dxn, dyn;
    for(j=0;j<tple2;j++){
        
        xdist = xto[j + tple2]-xto[j]; //old distances
        ydist = yto[j + tple2]-yto[j];
        dist=sqrt(xdist*xdist +ydist*ydist);
        
        dxo=xdist/dist;
        dyo=ydist/dist;
        
        xdist = xt[j + tple2]-xt[j]; //new distances
        ydist = yt[j + tple2]-yt[j];
        dist=sqrt(xdist*xdist +ydist*ydist);
        
        dxn=xdist/dist;
        dyn=ydist/dist;
        
        if (dxo*dyn-dxn*dyo>0) { //check rotation
            angle_sign=1;
        }
        else angle_sign=-1;
        
        theta[j] = theta[j] + 2* asin(0.5*(sqrt((dxn-dxo)*(dxn-dxo)+(dyn-dyo)*(dyn-dyo))))*angle_sign; //add the angels with the respective signs
        
    }
    for(j=0;j<tple;j++){ //new old coordinates
        xto[j]=xt[j];
        yto[j]=yt[j];
    }
    return 0;
}

//Update the coordinates of the center of mass

int CoM(double *xt, double *yt, double *CoMx, double *CoMy){
    int j;
    double xdist2, ydist2;
    
    
    for (j=0; j<tple2; j++) {
        CoMx[j]=0.0;
        CoMy[j]=0.0;
        
        xdist2 = (xt[j + tple2] - xt[j])/2;
        ydist2 = (yt[j + tple2] - yt[j])/2;
        
        CoMx[j] = xdist2 + xt[j];
        CoMy[j] = ydist2 + yt[j];
    }
    return 0;
}


/**
 * Neighbor list
 **/
int update_nlist(double *x, double *y){
	int j,k;
	double xdist, ydist, dist2;
	
	for(j=0;j<tple-1;j++){
		for(k=j+1;k<tple; k++){
			xdist=x[j]-x[k];
			ydist=y[j]-y[k];
						// Interaction between closed images for crowders
			if(pbc){
				{
					xdist=xdist-roundToint(xdist/box)*(double)box;
					ydist=ydist-roundToint(ydist/box)*(double)box;}
				}
				
			dist2=xdist*xdist +ydist*ydist;
			if(dist2<rskin2){
				nlist[j][k]=1;  
				nlist[k][j]=1;
			}else{
				nlist[j][k]=0;
				nlist[k][j]=0;
			}
		}
	}
	for(j=0;j<tple;j++){
		displacement[j]=0.0;
	}
	nlist_updates++;
	return 0;
}
	

/**
 * Check neighbor list
 * */
int check_nlist(){
	int j;
	
	double maxdisp=0.0;
	double max2disp=0.0;
	
	for(j=0;j<tple;j++){ //Give two largest displacement
		if(displacement[j]>maxdisp){
			max2disp=maxdisp;
			maxdisp=displacement[j];
		}else if(displacement[j]>max2disp){
			max2disp=displacement[j];
		}else {
			//do nothing
		}
	}
	if(maxdisp+max2disp>rskin-rcut){
		return 1;
	}else {
		return 0;
	}
}


/**
 * Gives two gaussian random numbers with zero mean and unit variance.
 * Argument is double array with two elements.
 * Uses the polar form of Box-Muller transform to convert two random
 * variables of [0,1)-uniform distribution into gaussian random numbers.
 */
void gaussrand(double* grnd){
	double r1, r2,y1,y2, w;
		do{
			r1=2.0*genrand64_real2()-1.0;
			r2=2.0*genrand64_real2()-1.0;
			w=r1*r1+r2*r2;
		}while(w >= 1.0);

	w=sqrt((-2.0*log(w))/w);
	y1=r1*w;
	y2=r2*w;
	grnd[0]=y1;
	grnd[1]=y2;
	}

	// Define roundToint
int roundToint(double x){
	if(x>=0) return (int)(x+0.5);
	return (int)(x-0.5);
}
