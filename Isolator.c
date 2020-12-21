#include "TStyle.h"
#include "TCanvas.h"
#include "TGraph.h"
#include "TTimer.h"
#include "function.h"
#include "function.c"
#include <time.h>
Double_t xstep=0.1,tstep=0.001,func_delta=0.000484;//步长
clock_t start,finish;
Double_t xmin=0,xmax=2.0,tmin=0.0,tmax=4.0;//自变量范围
Int_t timestochange=(Int_t)(1.0/tstep);
Int_t loop_index=0,loop_jndex=0;//循环指标
Int_t number_in_x=(Int_t)((xmax-xmin)/xstep)+1;
Int_t number_in_t=(Int_t)((tmax-tmin)/tstep)+1;
Double_t factor_second=1.0/3.0*tstep/xstep;
Double_t factor_third=func_delta*tstep/pow(xstep,3);
TF2 *f2;
TGraph *graphtest;
Double_t *xvalues=Create1D(number_in_x);
Double_t **funcvalues=Create2D(2,number_in_x);
Double_t *yyvalue=Create1D(number_in_x);
void Isolator()
{
    for(loop_index=0;loop_index<number_in_x;loop_index++)
    {
        xvalues[loop_index]=xmin+loop_index*xstep;
        funcvalues[0][loop_index]=cos(pi*xvalues[loop_index]);
    }
    /*根据初始条件计算下一个时刻的函数值*/
    funcvalues[1][0]=funcvalues[0][0]-0.5*factor_second*(funcvalues[0][1]+funcvalues[0][0]+funcvalues[0][number_in_x-2])*(funcvalues[0][1]-funcvalues[0][number_in_x-2])-0.5*factor_third*(funcvalues[0][2]-2.0*funcvalues[0][1]+2.0*funcvalues[0][number_in_x-2]-funcvalues[0][number_in_x-3]);//此处乘以0.5来源于时间导数是用相邻两个点计算，而其他的导数则是间隔了一个点取平均，以下几行代码类似
    funcvalues[1][1]=funcvalues[0][1]-0.5*factor_second*(funcvalues[0][2]+funcvalues[0][1]+funcvalues[0][0])*(funcvalues[0][2]-funcvalues[0][0])-0.5*factor_third*(funcvalues[0][3]-2.0*funcvalues[0][2]+2.0*funcvalues[0][0]-funcvalues[0][number_in_x-2]);
    funcvalues[1][number_in_x-1]=funcvalues[1][0];//周期性边界条件
    funcvalues[1][number_in_x-2]=funcvalues[0][number_in_x-2]-0.5*factor_second*(funcvalues[0][number_in_x-1]+funcvalues[0][number_in_x-2]+funcvalues[0][number_in_x-3])*(funcvalues[0][number_in_x-1]-funcvalues[0][number_in_x-3])-0.5*factor_third*(funcvalues[0][1]-2.0*funcvalues[0][0]+2.0*funcvalues[0][number_in_x-3]-funcvalues[0][number_in_x-4]);
    for(loop_jndex=2;loop_jndex<number_in_x-2;loop_jndex++)
    {
        funcvalues[1][loop_jndex]=funcvalues[0][loop_jndex]-0.5*factor_second*(funcvalues[0][loop_jndex+1]+funcvalues[0][loop_jndex]+funcvalues[0][loop_jndex-1])*(funcvalues[0][loop_jndex+1]-funcvalues[0][loop_jndex-1])-0.5*factor_third*(funcvalues[0][loop_jndex+2]-2.0*funcvalues[0][loop_jndex+1]+2.0*funcvalues[0][loop_jndex-1]-funcvalues[0][loop_jndex-2]);
    }
   gStyle->SetCanvasPreferGL(true);
   gStyle->SetFrameFillColor(10);
   TCanvas *c1 = new TCanvas("c1");
   c1->SetFillColor(17);
   TTimer *timer = new TTimer(10);
   //timer->SetCommand("Animate(xvalues,yvalue,yyvalue,number_in_x,factor_left,factor_right,graphtest)");
   timer->SetCommand("Animate()");
   timer->TurnOn();
}
//void Animate(Double_t *xvalues,Double_t * yvalue,Double_t * yyvalue,Int_t  number_in_x,Double_t factor_left,Double_t factor_right,TGraph *graphtest)
void Animate()
{
    start=clock();
    if (!gROOT->GetListOfCanvases()->FindObject("c1")) return;
    loop_index=0;
    for(loop_index=0;loop_index<1;loop_index++)
    {
        loop_jndex=2;
        yyvalue[0]=funcvalues[0][0]-factor_second*(funcvalues[1][1]+funcvalues[1][0]+funcvalues[1][number_in_x-2])*(funcvalues[1][1]-funcvalues[1][number_in_x-2])-factor_third*(funcvalues[1][2]-2.0*funcvalues[1][1]+2.0*funcvalues[1][number_in_x-2]-funcvalues[1][number_in_x-3]);
        yyvalue[1]=funcvalues[0][1]-factor_second*(funcvalues[1][2]+funcvalues[1][1]+funcvalues[1][0])*(funcvalues[1][2]-funcvalues[1][0])-factor_third*(funcvalues[1][3]-2.0*funcvalues[1][2]+2.0*funcvalues[1][0]-funcvalues[1][number_in_x-2]);
        yyvalue[number_in_x-1]=funcvalues[loop_index][0];//周期性边界
        yyvalue[number_in_x-2]=funcvalues[0][number_in_x-2]-factor_second*(funcvalues[1][number_in_x-1]+funcvalues[1][number_in_x-2]+funcvalues[1][number_in_x-3])*(funcvalues[1][number_in_x-1]-funcvalues[1][number_in_x-3])-factor_third*(funcvalues[1][1]-2.0*funcvalues[1][0]+2.0*funcvalues[1][number_in_x-3]-funcvalues[1][number_in_x-4]);
        for(loop_jndex=2;loop_jndex<number_in_x-2;loop_jndex++)
        {
            yyvalue[loop_jndex]=funcvalues[0][loop_jndex]-factor_second*(funcvalues[1][loop_jndex+1]+funcvalues[1][loop_jndex]+funcvalues[1][loop_jndex-1])*(funcvalues[1][loop_jndex+1]-funcvalues[1][loop_jndex-1])-factor_third*(funcvalues[1][loop_jndex+2]-2.0*funcvalues[1][loop_jndex+1]+2.0*funcvalues[1][loop_jndex-1]-funcvalues[1][loop_jndex-2]);
        }
        for(loop_jndex=0;loop_jndex<number_in_x;loop_jndex++)
        {
            funcvalues[0][loop_jndex]=funcvalues[1][loop_jndex];
            funcvalues[1][loop_jndex]=yyvalue[loop_jndex];
        }
    }
    graphtest = new TGraph(number_in_x,xvalues,yyvalue);
    graphtest->GetXaxis()->SetLimits(0,2);
    //graphtest->GetHistogram()->SetMaximum(1.0);
    //graphtest->GetHistogram()->SetMinimum(0);
    graphtest->Draw("AC");
    gPad->Modified();
    gPad->Update();
    finish=clock();
    //printf("程序运行的时间为(单位:秒)\n%lu秒\n",(finish-start)/CLOCKS_PER_SEC);
}
