//#include <QCoreApplication>
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

void SalvaArquivo(MatrixXf M, FILE *arquivo, VectorXf tempo, VectorXf X)
{
    arquivo = fopen("MatrizSolucaoTrabalho1.dat","w");
    if(!arquivo)
    {
        printf("\nerro ao abrir arquivo");
    }else
    {
        for(int i=0;i<M.rows();i++)
        {
            for(int j=0;j<M.cols();j++)
            {
               fprintf(arquivo, "%f  %f  %f  \n", tempo(i), X(j), M(i,j));
            }
        }
        cout << "arquivo atualizado" <<endl;
        fclose(arquivo);
    }
}

void SalvaArquivoSemTempo(MatrixXf M, FILE *arquivo, VectorXf X)
{
    arquivo = fopen("SolucaoTrabalho1SemTempo.dat","w");
    if(!arquivo)
    {
        printf("\nerro ao abrir arquivo");
    }else
    {
        int i = M.rows()-1;
        for(int j=0;j<M.cols();j++)
        {
           fprintf(arquivo, "%f  \n", M(i,j));
        }
        cout << "arquivo atualizado" <<endl;
        fclose(arquivo);
    }
}
/*----------------------------------------------------------------------------------------------------------------------------------------------*/
class Upwind
{
private:
public:
    //VARIAVEIS
    int tipoDeUpwind=0;
    float beta=0;
    float a=0,b=0;
    float threshold=1.0e-14;

    //FUNÇOES
    Upwind(float theta)
    {
        //construtor
        cout<< "1 = FSLS " <<endl;
        cout<< "2 = ADBQUICKEST"<<endl;
        cin >> tipoDeUpwind;

        switch(tipoDeUpwind)
        {
        case(0):
        {
            cout<<"Erro ao escolher o tipo de upwind"<<endl;
            break;
        }
        case(1):
        {
            //cout<<"Valor de Beta para FSLS: "<<endl;
            //cin >> this->beta;
            beta= 1.5;
            break;
        }
        case(2):
        {
            float num_a = 2-3*abs(theta)+theta*theta;
            float den_a = 7 - 6*theta- 3*abs(theta)+ 2*theta*theta;
            this->a = num_a/den_a;

            float num_b = -4+6*theta- 3*abs(theta)+ theta*theta;
            float den_b = -5+ 6*theta - 3*abs(theta) + 2*theta*theta;
            this->b = num_b/den_b;
        }
        }

    }

    Vector2f ADBQUICKEST(MatrixXf M, int i, int j, Vector2f v,float u_f, float u_g, float theta) //v[0] é u_g e v[1] é u_f /i=cont no tempo e j=cont no espaço
    {

        //calcula as propriedade convectada
        if(u_g>=0.0)
        {
            //calcula as propriedade convectada qdo a velocidade de transporte é positiva
            if(j == 1)
            {
                v(0)= 0.5*(M(i,j-1)+M(i,j));
            }else
            {
                float phi_chapeu = (M(i,j-1)-M(i,j-2))/(M(i,j)-M(i,j-2));
                if(phi_chapeu<0 || 1<phi_chapeu)
                {
                    v(0)= M(i,j-1);
                }else if(phi_chapeu>=0 && phi_chapeu<a)
                {
                    float A = (2-theta)*phi_chapeu;
                    v(0) = M(i,j-2) + (M(i,j)-M(i,j-2))*A;
                }else if(phi_chapeu>=a && phi_chapeu<=b)
                {
                    float A = phi_chapeu+0.5*(1-abs(theta))-(1/6)*(1-theta*theta)*(1-2*phi_chapeu);
                    v(0) = M(i,j-2) + (M(i,j)-M(i,j-2))*A;
                }else if(phi_chapeu>b && phi_chapeu<=1)
                {
                    float A = 1-theta+theta*phi_chapeu;
                    v(0) = M(i,j-2) + (M(i,j)-M(i,j-2))*A;
                }
            }

        }
        if(u_g<0.0)
        {
            //calcula as propriedade convectada qdo a velocidade de transporte é positiva
            if(fabs(M(i,j-1)-M(i,j+1)) <= threshold)
            {
                v(0) = M(i,j);
            }else
            {
                float phi_chapeu = (M(i,j)-M(i,j+1))/(M(i,j-1)-M(i,j+1));
                if(phi_chapeu<0 || 1<phi_chapeu)
                {
                    v(0)= M(i,j);
                }else if(phi_chapeu>=0 && phi_chapeu<a)
                {
                    float A = (2-theta)*phi_chapeu;
                    v(0) = M(i,j+1) + (M(i,j-1)-M(i,j+1))*A;
                }else if(phi_chapeu>=a && phi_chapeu<=b)
                {
                    float A = phi_chapeu+0.5*(1-abs(theta))-(1/6)*(1-theta*theta)*(1-2*phi_chapeu);
                    v(0) = M(i,j+1) + (M(i,j-1)-M(i,j+1))*A;
                }else if(phi_chapeu>b && phi_chapeu<=1)
                {
                    float A = 1-theta+theta*phi_chapeu;
                    v(0) = M(i,j+1) + (M(i,j-1)-M(i,j+1))*A;
                }
            }

        }
        if(u_f>=0.0)
        {
            //calcula as propriedade convectada qdo a velocidade de transporte é positiva
            if(fabs(M(i,j+1)-M(i,j-1)) <= threshold)
            {
                v(1)= M(i,j);
            }else
            {
                float phi_chapeu = (M(i,j)-M(i,j-1))/(M(i,j+1)-M(i,j-1));
                if(phi_chapeu<0 || 1<phi_chapeu)
                {
                    v(1)= M(i,j);
                }else if(phi_chapeu>=0 && phi_chapeu<a)
                {
                    float A = (2-theta)*phi_chapeu;
                    v(1) = M(i,j-1) + (M(i,j+1)-M(i,j-1))*A;
                }else if(phi_chapeu>=a && phi_chapeu<=b)
                {
                    float A = phi_chapeu+0.5*(1-abs(theta))-(1/6)*(1-theta*theta)*(1-2*phi_chapeu);
                    v(1) = M(i,j-1) + (M(i,j+1)-M(i,j-1))*A;
                }else if(phi_chapeu>b && phi_chapeu<=1)
                {
                    float A = 1-theta+theta*phi_chapeu;
                    v(1) = M(i,j-1) + (M(i,j+1)-M(i,j-1))*A;
                }
            }

        }
        if(u_f<0.0)
        {
            //calcula as propriedade convectada qdo a velocidade de transporte é positiva
            if(j == M.cols()-2)
            {
                v(1)= 0.5*(M(i,j-1)+M(i,j+1));
            }else
            {
                float phi_chapeu = (M(i,j+1)-M(i,j+2))/(M(i,j)-M(i,j+2));
                if(phi_chapeu<0 || 1<phi_chapeu)
                {
                    v(1)= M(i,j+1);
                }else if(phi_chapeu>=0 && phi_chapeu<a)
                {
                    float A = (2-theta)*phi_chapeu;
                    v(1) = M(i,j+2) + (M(i,j)-M(i,j+2))*A;
                }else if(phi_chapeu>=a && phi_chapeu<=b)
                {
                    float A = phi_chapeu+0.5*(1-abs(theta))-(1/6)*(1-theta*theta)*(1-2*phi_chapeu);
                    v(1) = M(i,j+2) + (M(i,j)-M(i,j+2))*A;
                }else if(phi_chapeu>b && phi_chapeu<=1)
                {
                    float A = 1-theta+theta*phi_chapeu;
                    v(1) = M(i,j+2) + (M(i,j)-M(i,j+2))*A;
                }
            }

        }
        return v;
    }

    Vector2f FSLS(MatrixXf M, int i, int j, Vector2f v,float u_f, float u_g)
    {

        //verifica em qual caso a variavel convectada se encontra
        if(u_g>=0.0)
        {
            //calcula as propriedade convectada qdo a velocidade de transporte é positiva
            if(j == 1.0)
            {
                v(0)= 0.5*(M(i,j-1)+M(i,j));
            }else
            {
                if(fabs(M(i,j)-M(i,j-2))<=threshold)
                {
                    v(0)= M(i,j-1);
                }else
                {
                    float phi_chapeu = (M(i,j-1)-M(i,j-2))/(M(i,j)-M(i,j-2));
                    if(phi_chapeu<0.0 || 1.0<phi_chapeu)
                    {
                        v(0)= M(i,j-1);
                    }else
                    {
                        float A = (-2.0*beta+4.0)*pow(phi_chapeu,4.0) + (4.0*beta-8.0)*pow(phi_chapeu,3.0)+ ((-5.0*beta+8.0)/2.0)*pow(phi_chapeu,2.0)+ ((beta+2.0)/2.0)*phi_chapeu;
                        v(0) = M(i,j-2) + (M(i,j)-M(i,j-2))*A;
                    }
                }

            }

        }else if(u_g<0.0)
        {
            //calcula as propriedade convectada qdo a velocidade de transporte é positiva
            if(fabs(M(i,j-1)-M(i,j+1)) <= threshold)
            {
                v(0)= M(i,j);
            }else
            {
                float phi_chapeu = (M(i,j)-M(i,j+1))/(M(i,j-1)-M(i,j+1));
                if(phi_chapeu<0.0 || 1.0<phi_chapeu)
                {
                    v(0)= M(i,j);
                }else
                {
                    float A = (-2.0*beta+4.0)*pow(phi_chapeu,4.0) + (4.0*beta-8.0)*pow(phi_chapeu,3.0)+ ((-5.0*beta+8.0)/2.0)*pow(phi_chapeu,2.0)+ ((beta+2.0)/2.0)*phi_chapeu;
                    v(0) = M(i,j+1) + (M(i,j-1)-M(i,j+1))*A;
                }
            }
        }
        if(u_f>=0.0)
        {
            //calcula as propriedade convectada qdo a velocidade de transporte é positiva
            if(fabs(M(i,j+1)-M(i,j-1))<=threshold)
            {
                v(1)= M(i,j);
            }else
            {
                float phi_chapeu = (M(i,j)-M(i,j-1))/(M(i,j+1)-M(i,j-1));
                if(phi_chapeu<0 || 1<phi_chapeu)
                {
                    v(1)= M(i,j);
                }else
                {
                    float A = (-2.0*beta+4.0)*pow(phi_chapeu,4.0) + (4.0*beta-8.0)*pow(phi_chapeu,3.0)+ ((-5.0*beta+8.0)/2.0)*pow(phi_chapeu,2.0)+ ((beta+2.0)/2.0)*phi_chapeu;
                    v(1) = M(i,j-1) + (M(i,j+1)-M(i,j-1))*A;
                }
            }
        }else if(u_f<0.0)
        {
            //calcula as propriedade convectada qdo a velocidade de transporte é positiva
            if(j == M.cols()-2)
            {
                v(1)= 0.5*(M(i,j)+M(i,j+1));
            }else
            {
                if(fabs(M(i,j)-M(i,j+2)) <= threshold)
                {
                    v(1)= M(i,j+1);
                }else
                {
                    float phi_chapeu = (M(i,j+1)-M(i,j+2))/(M(i,j)-M(i,j+2));
                    if(phi_chapeu<0.0 || 1.0<phi_chapeu)
                    {
                        v(1)= M(i,j+1);
                    }else
                    {
                        float A = (-2.0*beta+4.0)*pow(phi_chapeu,4.0) + (4.0*beta-8.0)*pow(phi_chapeu,3.0)+ ((-5.0*beta+8.0)/2.0)*pow(phi_chapeu,2.0)+ ((beta+2.0)/2.0)*phi_chapeu;
                        v(1) = M(i,j+2) + (M(i,j)-M(i,j+2))*A;
                    }
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
          deltaX,
          delta_t,
          comprimento,
          tempo,
          contornoInicio,
          contornoFim,
          theta=0;


    int contX;
    int cont_t;

    Burgers_Viscosa()
    {
        cout<<"viscosidade: "<<endl;
        cin >> this->visc;
        cout<<"Theta: "<<endl;
        cin >> this->theta;
        //cout<<"Comprimento X total: "<<endl;
        //cin >> this->comprimento;
        this->comprimento = 1;
        //cout<<"Passo em X: "<<endl;
        //cin >> this->deltaX;
        this->deltaX=0.01;
        cout<<"Tempo total: "<<endl;
        cin >> this->tempo;
        //this->tempo = 0.5
        //cout<<"Passo do tempo: "<<endl;
        //cin >> this->delta_t;
        this->delta_t=theta*deltaX;


    }
    //FUNCOES
    MatrixXf CalculaEquacao(VectorXf X)
    {

        MatrixXf U(cont_t,contX); //matriz com as Velocidades

        U = this->InsereCondicoesIniciais(U, X);

        //começando o calculo

        float uBarra_f = 0, uBarra_g = 0;

        Upwind upwind(theta);
        Vector2f u_gf = {0,0};

        switch(upwind.tipoDeUpwind)
        {
            case(0):
            {
                cout<<"erro"<<endl;
                break;
            }
            case(1):
            {
                cout<<"FSLS"<<endl;
                for(int n=0; n<U.rows()-1;n++)
                {
                    for(int i=1; i<U.cols()-1;i++)
                    {
                        //loop do espaço
                        uBarra_f = (U(n,i+1)+U(n,i))*0.5;
                        uBarra_g = (U(n,i)+U(n,i-1))*0.5;

                        u_gf = upwind.FSLS(U,n,i, u_gf,uBarra_f,uBarra_g);
                        float p1 = 0.5*delta_t*(1.0/deltaX)*(uBarra_f* u_gf(1) - uBarra_g * u_gf(0));
                        float p2 = (delta_t*visc/ (deltaX*deltaX))*(U(n,i+1)-2.0*U(n,i) + U(n,i-1));
                        U(n+1,i)= U(n,i) - p1 +  p2;

                    }
                }
                cout << "deu certo fsls" <<endl;
                break;
            }
            case(2):
            {
                cout << "ADBQUICKEST"<<endl;
                for(int n=0; n<U.rows()-1;n++)
                {
                    for(int i=1; i<U.cols()-1;i++)
                    {
                        //loop do espaço
                        uBarra_f = (U(n,i+1)+U(n,i))*0.5;
                        uBarra_g = (U(n,i)+U(n,i-1))*0.5;

                        u_gf = upwind.ADBQUICKEST(U,n,i, u_gf,uBarra_f, uBarra_g,theta);

                        U(n+1,i)= U(n,i) - 0.5*delta_t*(1.0/deltaX)*(uBarra_f* u_gf(1) - uBarra_g * u_gf(0)) + (delta_t*visc/ (deltaX*deltaX))*(U(n,i+1) - 2*U(n,i) + U(n,i-1));

                    }
                }
                cout << "deu certo adb" <<endl;
                break;
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


//int main(int argc, char *argv[])
int main()
{
    QCoreApplication a(argc, argv);

    FILE *arquivo;

    while(true)
    {
        Burgers_Viscosa Equacao;

        Equacao.contornoInicio = 0;
        Equacao.contornoFim = 0;

        Equacao.cont_t= floor(Equacao.tempo/Equacao.delta_t)+1;
        Equacao.contX = floor(Equacao.comprimento/Equacao.deltaX)+1;



        VectorXf X(Equacao.contX);
        VectorXf tempo(Equacao.cont_t);

        for(int c=0;c<X.rows();c++)
        {
            //inicia o vetor X com a posição de cada nó do dominio 1D com a borda esquerda sendo a coordenada 0
            if(c==0)
            {
                X[c]=0;

            }else
            {
                X[c]=X[c-1]+Equacao.deltaX;

            }
        }

        for(int c=0;c<tempo.rows();c++)
        {
            //inicia o vetor tempo com a posição de cada nó no decorrer do tempo
            if(c==0)
            {
                tempo[c]=0;

            }else
            {
                tempo[c]=tempo[c-1]+Equacao.delta_t;

            }
        }

        MatrixXf solucao = Equacao.CalculaEquacao(X);
        cout <<solucao.rows()<<"x"<<solucao.cols()<<endl;
        //ImprimeMatriz(solucao);

        SalvaArquivo(solucao, arquivo, tempo, X);
        SalvaArquivoSemTempo(solucao, arquivo, X);
    }
    //return a.exec();
    return 0;
}
