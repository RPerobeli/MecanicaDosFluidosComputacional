/*Rafael Alves Bonfim de Queiroz*/
/*Código Revisado 04/08/2020 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

/* definindo um vetor de numeros reais */
typedef float         *vreal;

#define pi	3.142

/*-----------------------------------------------------------------------------------------------*/
vreal  ALLOCVREAL (int m){
	vreal v;    /* ponteiro para um vetor*/
  	int   i;    /* variavel auxiliar*/	
 	v = (vreal) malloc((m)*sizeof(vreal));
  	if (v == NULL) {
    		printf("Vetor 1 : não alocado\n");
    		exit(1);
  	}
  	for (i = 0; i <= m; i++) 
		v[i]= 0.0;
  	return ( v );
}

/*-----------------------------------------------------------------------------------------------*/
int main(){   
	int i,j,aux,bux,yn,xn,esquema;
  	float deltat,deltax,x,t;
	float ubarra1,ubarra2,conv,visc,fiu,numa,dena,numb,denb,a,b,cfl,phiF;
  	float aux1,aux2;
  	FILE *p;
  	vreal u,unew;
  	float eps=1.0e-14;
       
   	deltax = 0.005;
   	xn = 200; /*xn = 1/deltax*/

    	printf("\n Informe a valor da viscosidade: ");
    	scanf("%f",&visc);

    	printf("\n Informe o valor do cfl: ");
    	scanf("%f",&cfl);
  
    	printf("\n Informe o tempo final: ");
    	scanf("%f",&t);

    	deltat = cfl * deltax;
    	yn = (int)(t/deltat);
   
    	printf("\n Esquema (1) - ADBQUICKEST");
    	printf("\n Esquema (2) - Esquema que implementou");
  
    	printf("\n Informe numero do esquema: ");
    	scanf("%d",&esquema);
 
     	t = yn*deltat;
   
     	printf("----------dados:---------------\n");
     	printf("tempo: %f \n",t);
     	printf("iteracoes: %d \n",yn);
     	printf("viscosidade: %f \n",visc);
     	printf("cfl: %f \n",cfl);
     	printf("deltat: %f \n",deltat);
	
     	p = fopen("Burger.dat","w");/*Arquivo que guarda os resultados*/
		
     	u = (vreal) ALLOCVREAL(xn+1);
     	unew = (vreal) ALLOCVREAL(xn+1);
		
        /*parametros do esquema ADBQUICKEST*/
     	numa = (2.0-(3.0*fabs(cfl))+(cfl*cfl));
     	dena = (7.0-(6.0*cfl)-(3.0*fabs(cfl))+(2.0*cfl*cfl));
     	a      = numa/dena;
     	numb = ((cfl*cfl)-4.0-(3.0*fabs(cfl))+(6*cfl));
     	denb = (-5.0-(3.0*fabs(cfl))+(2.0*cfl*cfl)+(6.0*cfl));
     	b    = numb/denb;

    	/*Condição inicial*/
    	for (i = 1; i < xn; i++)
       		u[i] = sin(2.0 * pi * i * deltax);
    	u[0] = 0.0;
    	u[xn] = 0.0;
	
    	for (j=0; j<=yn; j++){
        	for(i=1; i<xn; i++){
	  		ubarra1 = 0.5*(u[i+1]+u[i]);
	  		ubarra2 = 0.5*(u[i]+u[i-1]);
	  		if (ubarra1 > 0.0){
            			if (fabs(u[i+1]-u[i-1])<=eps)
	       				aux1 = ubarra1*u[i];/*FOU*/
	    			else{
	       				fiu = (u[i]-u[i-1])/(u[i+1]-u[i-1]);
               				/*phiU=u[i]; phiR=u[i-1]; phiD=u[i+1]*/
               				if ((fiu < 0.0) || (fiu > 1.0)){
                  				aux1 = ubarra1*u[i];
               				}else{                      
                  				if (esquema == 1){ //ADBQUICKEST            
                     					if (fiu < a)
                                                                phiF = u[i-1] + (u[i+1]-u[i-1]) * ((2.0 - cfl) * fiu);
                     					if ((a <= fiu) && (fiu <= b))
                                                                phiF = u[i-1] + (u[i+1]-u[i-1]) * (fiu + 0.5 * (1.0 -fabs(cfl)) *(1.0-fiu) -(1.0/6.0) *(1.0-pow(fiu,2.0)) * (1.0-2.0*fiu));
            	     					if (b < fiu)
                                                                phiF = u[i-1] + (u[i+1]-u[i-1]) * (1.0 - cfl + cfl * fiu);

                                                        aux1 = ubarra1* phiF;
		  				}
		  				else if (esquema == 2){
		    					//Esquema que implementou
	          				}
                  
               				}						
      	   			}
	 		}
         		if (ubarra1 <= 0.0){
            			bux = i+2;
            			if (bux == xn+1){
	      				aux1 = ubarra1*0.5*(u[i]+u[i+1]); /*Central*/
	    			}
	    			else{
              				if (fabs(u[i]-u[i+2]) <= eps){
                 				aux1 = ubarra1*u[i+1]; /*FOU*/
              				} 
             				else{
                				fiu = (u[i+1]-u[i+2])/(u[i]-u[i+2]);
               					/*phiU=u[i+1]; phiR=u[i+2]; phiD=u[i]*/
                				if ((fiu < 0.0) || (fiu > 1.0)){
                    					aux1 = ubarra1*u[i+1];
                				}else{
		   					if (esquema==1){ //ADBQUICKEST
               	      						if (fiu < a)
                                                                          phiF = u[i+2] + (u[i]-u[i+2]) * ((2.0 - cfl) * fiu);
                      						if ((a <= fiu) && (fiu <= b))
                                                                          phiF = u[i+2] + (u[i]-u[i+2]) * (fiu + 0.5 * (1.0 -fabs(cfl)) *(1.0-fiu) -(1.0/6.0) *(1.0-pow(fiu,2.0)) * (1.0-2.0*fiu));
                      						if (b < fiu)
                                                                          phiF = u[i+2] + (u[i]-u[i+2]) * (1.0 - cfl + cfl * fiu);
                                                                aux1 = ubarra1* phiF;
		   					}
		   					else if (esquema == 2){
		      						//equema que implementou
		   					}
						}
              				}
	    			}
	  		}
			
          		if (ubarra2 > 0.0){
	     			aux = i-2;
             			if (aux == -1.0){
	       				aux2 = ubarra2*0.5*(u[i]+u[i-1]); /*Central*/
	     			}
	     			else{
	       				if (fabs(u[i]-u[i-2]) <= eps){
                  				aux2 = ubarra2*u[i-1];/*FOU*/
               				}else{
                  				fiu = (u[i-1]-u[i-2])/(u[i]-u[i-2]);
                  				/*phiU=u[i-1], phiR=u[i-2], phiD=u[i]*/
                  				if ((fiu < 0.0) || (fiu > 1.0)) {
                     					aux2 = ubarra2*u[i-1];
                  				}else{
                     					if (esquema == 1){ //ADBQUICKEST
                						if (fiu < a)
                                                                        phiF = u[i-2] + (u[i]-u[i-2]) * ((2.0 - cfl) * fiu);
                           					if ((a <= fiu) && (fiu <=  b))
                                                                        phiF = u[i-2] + (u[i]-u[i-2]) * (fiu + 0.5 * (1.0 -fabs(cfl)) *(1.0-fiu) -(1.0/6.0) *(1.0-pow(fiu,2.0)) * (1.0-2.0*fiu));
                           					if (b <  fiu)
                                                                        phiF = u[i-2] + (u[i]-u[i-2]) *  (1.0 - cfl + cfl * fiu);
                                                                aux2 = ubarra2* phiF; 
							}
		     					else if (esquema == 2){
		        					//Esquema que implementou 
		     					}
             	  				}
         				}
	      			}
            		}
	    		if (ubarra2 <= 0.0){
	       			if (fabs(u[i-1]-u[i+1]) <= eps){
		  			aux2 = ubarra2*u[i]; /*FOU*/
      	      			}
	       			else{
                 			fiu = (u[i]-u[i+1])/(u[i-1]-u[i+1]);
                  			/*phiU=u[i], phiR=u[i+1], phiD=u[i-1]*/
                 			if ((fiu < 0.0) || (fiu > 1.0)) {
              	    				aux2 = ubarra2*u[i]; /*FOU*/
                 			} 
		 			else{
		   				if (esquema == 1){ //ADBQUICKEST
                      					if (fiu < a)
                                                                phiF = u[i+1] + (u[i-1]-u[i+1]) * ((2.0 - cfl) * fiu);
                      					if ((a <= fiu) && (fiu <= b))
                                                                 phiF = u[i+1] + (u[i-1]-u[i+1]) * (fiu + 0.5 * (1.0 -fabs(cfl)) *(1.0-fiu) -(1.0/6.0) *(1.0-pow(fiu,2.0)) * (1.0-2.0*fiu));
                      					if (b < fiu)
                                                                 phiF = u[i+1] + (u[i-1]-u[i+1]) * (1.0 - cfl + cfl * fiu);
                                                        aux2 = ubarra2* phiF;
		   				}
		   				else if (esquema == 2){
		        				//Esquema que implementou
		   				}
          	 			}
      	      			}
	    		}
	    		conv = 0.5*(1.0/deltax)*(aux1-aux2);
	    		/*o vetor unew guarda os novos valores de u*/
	    		unew[i] = u[i]-(deltat)*conv+(deltat*visc/(deltax*deltax))*(u[i+1]-2.0*u[i]+u[i-1]);
	  	}
	  	for (i=0; i<=xn; i++){
	     		u[i] = unew[i];
	  	}
       	}

        /*salva a solucao*/       
       	for (i=0; i<=xn; i++){
        	x = i*deltax;
	   	fprintf(p,"%f %lf \n", x, u[i]);
	}
	return(1);   
}
