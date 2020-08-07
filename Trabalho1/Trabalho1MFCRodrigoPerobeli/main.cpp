#include <QCoreApplication>
#include <../Eigen/Eigen/Dense>
#include <stdlib.h>
#include <math.h>
#include <iostream>

#define pi 3.141589

using namespace Eigen;
using namespace std;

/*----------------------------------------------------------------------------------------------------------------------------------------------*/
void ImprimeMatriz(MatrixXf M)
{
    for(int i=0; i<M.rows();i++)
    {
        for(int j=0;j<M.cols();j++)
        {
            cout << M(i,j) <<" ";
        }
        cout << endl;
    }
}

void SalvaArquivo(MatrixXf M, FILE *arquivo)
{
    arquivo = fopen("MatrizSolucaoTrabalho1.txt","w");
    if(!arquivo)
    {
        printf("\nerro ao abrir arquivo");
    }else
    {
        fprintf(arquivo,"%d\n%d\n ", M.rows(), M.cols());

        for(int i=0;i<M.rows();i++)
        {
            for(int j=0;j<M.cols();j++)
            {
               fprintf(arquivo, "%f\n", M(i,j));
            }
        }
        fclose(arquivo);
    }
}
/*----------------------------------------------------------------------------------------------------------------------------------------------*/
class Upwind
{
private:
public:
    //FUNÇOES
    Vector2f FSLS(MatrixXf M, int i, int j, Vector2f v, float beta) //v[0] é u_g e v[1] é u_f /i=cont no tempo e j=cont no espaço
    {

        for(int k =0; k<v.rows();k++)
        {
            if(j == 1 && k==0)
            {
                v(k)= M(i,j-1);
            }else
            {
                float u_chapeu = (M(i,j-1+k)-M(i,j-2+k))/(M(i,j+k)-M(i,j-2+k));
                if(u_chapeu<0 || 1<u_chapeu)
                {
                    v(k)= M(i,j-1+k);
                }else
                {
                    float A = (-2*beta+4)*pow(u_chapeu,4) + (4*beta-8)*pow(u_chapeu,3)+ ((-5*beta+8)/2)*pow(u_chapeu,2)+ ((beta+2)/2)*u_chapeu;
                    v(k) = M(i,j-2+k) + (M(i,j+k)-M(i,j-2+k))*A;
                }
            }
        }
        return v;
    }
};
/*----------------------------------------------------------------------------------------------------------------------------------------------*/
class Burgers_Viscosa
{
  private:
  public:
    //VARIAVEIS
    float visc,
          theta,
          deltaX,
          delta_t,
          comprimento,
          tempo,
          contornoInicio,
          contornoFim,
          beta;

    int contX;
    int cont_t;

    //FUNCOES
    MatrixXf CalculaEquacao()
    {

        cont_t= floor(this->tempo/this->delta_t)+1;
        contX = floor(this->comprimento/this->deltaX)+1;

        VectorXf X(contX);
        MatrixXf U(cont_t,contX); //matriz com as Velocidades


        for(int c=0;c<X.rows();c++)
        {
            //inicia o vetor X com a posição de cada nó do dominio 1D com a borda esquerda sendo a coordenada 0
            if(c==0)
            {
                X[c]=0;

            }else
            {
                X[c]=X[c-1]+deltaX;

            }
        }

        U = this->InsereCondicoesIniciais(U, X);

        //começando o calculo

        float uBarra_f = 0, uBarra_g = 0;

        Upwind upwind;
        Vector2f u_gf = {0,0};

        for(int n=0; n<U.rows()-1;n++)
        {
            for(int i=1; i<U.cols()-1;i++)
            {
                //loop do espaço
                uBarra_f = (U(n,i+1)+U(n,i))*0.5;
                uBarra_g = (U(n,i)+U(n,i-1))*0.5;

                u_gf = upwind.FSLS(U,n,i, u_gf,beta);

                U(n+1,i)= U(n,i) + (delta_t/deltaX)*(uBarra_f* u_gf(1) - uBarra_g * u_gf(0)) + (delta_t*visc/ (deltaX*deltaX))*(U(n,i+1) - 2*U(n,i) + U(n,i-1));

            }
        }
        return U;
    }

    MatrixXf InsereCondicoesIniciais(MatrixXf M, VectorXf X)
    {

        for(int i=0; i<M.rows();i++)
        {
            for(int j=0;j<M.cols();j++)
            {
                if(i==0)
                {
                    //APLICA CONDIÇÃO INICIAL
                    M(i,j) = sin(2.0*pi*X(j));

                }else if(j==0)
                {
                    //APLICA CONDIÇOES DE CONTORNO
                    M(i,j)= this->contornoInicio;
                }else if(j==M.cols()-1)
                {
                    //APLICA CONDIÇOES DE CONTORNO
                    M(i,j)= this->contornoFim;
                }else
                {
                    M(i,j)=0;
                }
            }
        }
        //ImprimeMatriz(M);
        return M;
    }


};

/*----------------------------------------------------------------------------------------------------------------------------------------------*/


int main(int argc, char *argv[])
//int main()
{
    QCoreApplication a(argc, argv);

    FILE *arquivo;
    Burgers_Viscosa Equacao;

    Equacao.visc = 0;
    Equacao.theta = 0;
    Equacao.deltaX=0.1,
    Equacao.delta_t=0.1,
    Equacao.comprimento =1,
    Equacao.tempo = 0.5,
    Equacao.contornoInicio = 0,
    Equacao.contornoFim = 0,
    Equacao.beta = 1.5;

    int tipoDeUpwind=0;

    cout<< "1 = FSLS " <<endl;
    cout<< "2 = ADBQUICKEST"<<endl;
    cin >> tipoDeUpwind;

    MatrixXf solucao = Equacao.CalculaEquacao();

    ImprimeMatriz(solucao);

    SalvaArquivo(solucao, arquivo);



    return a.exec();
    //return 0;
}
