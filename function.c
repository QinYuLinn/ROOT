#include "function.h"
Double_t function_gauss(Double_t x)
{
    return(exp(-10*pi*x*x));//定义初始注入函数为高斯波包
}
Double_t ** Create2D(Int_t Dimension1,Int_t Dimension2)//
{
    Int_t i=0,j=0;
    Double_t **A;
    A=(Double_t**) malloc(Dimension1*(sizeof(Double_t *)));
    for(int i=0;i<Dimension1;i++)
    {
        A[i]=(Double_t*)malloc(Dimension2*(sizeof(Double_t)));
    }
    for(i=0;i<Dimension1;i++)
    {
        j=0;
        for(;j<Dimension2;j++)
        {
            A[i][j]=0;
        }
    }
    return A;
}
Double_t * Create1D(Int_t Dimension)
{
    Double_t *A;
    A=(Double_t *)malloc(Dimension*(sizeof(Double_t)));
    for(Int_t i=0;i<Dimension;i++)
    {
        A[i]=0;
    }
    return A;
}
