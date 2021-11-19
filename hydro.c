#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <linux/kd.h>

double rnd();

int main () {
  
  int i,j,k,n,l,a;
  double r1[3],r2[3],rp1[3],rp2[3],acc,p,x,delta,beta,r1R1,r1R2,r2R1,r2R2,r1r2,lnpsi,lnpsip,v,vp,t,tp,epot,ekin,etot,etotmin,evir,R,e11,e12,e21,e22;

  srand((unsigned)time(NULL));

  R=1.0;
  
  beta=2.2;   // parametro variazionale 

  n=1000;// number of monte carlo steps
  k=1000;
  
  delta=0.2;   // size of the move
  
  r1[0]=R/2+delta*(rnd()-0.5);
  r1[1]=delta*(rnd()-0.5);
  r1[2]=delta*(rnd()-0.5);
  
  r2[0]=-R/2+delta*(rnd()-0.5);//-
  r2[1]=delta*(rnd()-0.5);
  r2[2]=delta*(rnd()-0.5);
  
  r1R1=sqrt(pow(r1[0]+R/2,2)+pow(r1[1],2)+pow(r1[2],2));                 //distanza e1 n1
  r1R2=sqrt(pow(r1[0]-R/2,2)+pow(r1[1],2)+pow(r1[2],2));                 //distanza e1 n2
  r2R1=sqrt(pow(r2[0]+R/2,2)+pow(r2[1],2)+pow(r2[2],2));                 //distanza e2 n1
  r2R2=sqrt(pow(r2[0]-R/2,2)+pow(r2[1],2)+pow(r2[2],2));                 //distanza e2 n2

  r1r2=sqrt(pow(r1[0]-r2[0],2)+pow(r1[1]-r2[1],2)+pow(r1[2]-r2[2],2));   //distanza e1 e2

  e11=exp(-beta*r1R1);
  e12=exp(-beta*r1R2);
  e21=exp(-beta*r2R1);
  e22=exp(-beta*r2R2);
 
  lnpsi=log((e11 + e12)*(e21 + e22));
  v=-1.0/r1R1 -1.0/r1R2 -1.0/r2R2 -1.0/r2R1 +1.0/R +1.0/r1r2 ;
  t=-0.5*((beta*(beta-2/r1R1)*e11+beta*(beta-2/r1R2)*e12)/(e11+e12)+(beta*(beta-2/r2R1)*e21+beta*(beta-2/r2R2)*e22)/(e21+e22));

  FILE * fp;
   for(x=2.0; x<3.; x+=0.1){  //al variare di un parametro
     etotmin=10;
     R=0.1;
     //beta=x;               //parametro
     delta=x;
     for(a=0; a<25; a++){
    
       char nome[10];
        
       snprintf(nome, sizeof(char) * 32, "R_%f_p_%f.out",R,x);
       fp = fopen(nome, "w+");
       for (j=0; j<k; j++) {
	 r1[0]=R/2+delta*(rnd()-0.5);
	 r1[1]=delta*(rnd()-0.5);
	 r1[2]=delta*(rnd()-0.5);
  
	 r2[0]=-R/2+delta*(rnd()-0.5);//-
	 r2[1]=delta*(rnd()-0.5);
	 r2[2]=delta*(rnd()-0.5);
	 etot=0.0;
	 epot=0.0;
	 ekin=0.0;
	 acc=0.0;
	 for (i=0; i<n ; i++) {
	   rp1[0]=r1[0]+delta*(rnd()-0.5);
	   rp1[1]=r1[1]+delta*(rnd()-0.5);
	   rp1[2]=r1[2]+delta*(rnd()-0.5);
    
	   rp2[0]=r2[0]+delta*(rnd()-0.5);
	   rp2[1]=r2[1]+delta*(rnd()-0.5);
	   rp2[2]=r2[2]+delta*(rnd()-0.5);
    
	   r1R1= sqrt(pow(rp1[0]+R/2,2)+pow(rp1[1],2)+pow(rp1[2],2));
	   r1R2= sqrt(pow(rp1[0]-R/2,2)+pow(rp1[1],2)+pow(rp1[2],2));
	   r2R1= sqrt(pow(rp2[0]+R/2,2)+pow(rp2[1],2)+pow(rp2[2],2));
	   r2R2= sqrt(pow(rp2[0]-R/2,2)+pow(rp2[1],2)+pow(rp2[2],2));

	   r1r2=sqrt(pow(rp1[0]-rp2[0],2)+pow(rp1[1]-rp2[1],2)+pow(rp1[2]-rp2[2],2));

	   e11=exp(-beta*r1R1);
	   e12=exp(-beta*r1R2);
	   e21=exp(-beta*r2R1);
	   e22=exp(-beta*r2R2);
    
	   lnpsip=log((e11+e12)*(e21+e22));
	   vp=-1.0/r1R1 -1.0/r1R2 -1.0/r2R2 -1.0/r2R1 +1.0/R +1.0/r1r2 ;
	   tp=-0.5*((beta*(beta-2/r1R1)*e11+beta*(beta-2/r1R2)*e12)/(e11+e12)+(beta*(beta-2/r2R1)*e21+beta*(beta-2/r2R2)*e22)/(e21+e22));
	   p=exp(2*(lnpsip-lnpsi));
    
	   if(p>rnd() ){//metropolis
      
	     r1[0]=rp1[0];
	     r1[1]=rp1[1];
	     r1[2]=rp1[2];
      
	     r2[0]=rp2[0];
	     r2[1]=rp2[1];
	     r2[2]=rp2[2];
      
	     lnpsi=lnpsip;
	     v=vp;
	     t=tp;
	     acc++;
	   }

	   epot=epot+v;
	   ekin=ekin+t;
	   etot=etot+v+t;
	 }
	 epot=epot/n;
	 ekin=ekin/n;
	 etot=etot/n;
	 acc=acc/n;
	 evir=epot/2;
	 fprintf(fp,"%f %f %f %f %f %f %f %f %f %f %f\n",r1[0],r1[1],r1[2],r2[0],r2[1],r2[2],etot,epot,ekin,evir,acc);
	 //fprintf(fp,"%f %f %f %f %f\n",etot,epot,ekin,evir,acc);
       }

       //chiudo il file
       fclose(fp);
       //definisco e costruisco il comando
       char command[128]="";
    
       snprintf(command, sizeof(char)*127, "awk '{print $7}' R_%f_p_%f.out | ./statfor> p_%f_S_%f.out",R,x,x,R);
       //mando al terminale e gli faccio applicare statfor ai dati generati
       printf("%s\n",command);
       system(command);
    
       R+=0.1;
    
        }  
      }
   fprintf(stdout, "\aBeep!\n" );
   fclose(fp); 
   return(0);

  
}

double rnd() {
  double r=((double)rand()+1.)/((double)RAND_MAX+2.);
  return r;
}
