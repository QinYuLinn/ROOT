#include "TStyle.h"
#include "TCanvas.h"
#include "TGraph.h"
#include "TTimer.h"
#include "function.h"
#include "function.c"
Double_t xstep=0.05,tstep=0.05;//步长
Double_t xmin=-15.0,xmax=1.0,tmin=0.0,tmax=4.0,set_alpha=-1.0;//自变量范围
Int_t loop_index=0,loop_jndex=0;//循环指标
Int_t number_in_x=(Int_t)((xmax-xmin)/xstep)+1;
Int_t number_in_t=(Int_t)((tmax-tmin)/tstep)+1;
Double_t factor_left=factor_left=1.0+set_alpha*tstep/xstep;
Double_t factor_right=-1.0*set_alpha*tstep/xstep;//迎风法的计算下一个值得两个系数
TF2 *f2;
TGraph *graphtest;
Double_t *xvalues=Create1D(number_in_x);
Double_t *yvalue=Create1D(number_in_x);
Double_t *yyvalue=Create1D(number_in_x);
void UpwindWithGauss()
{
    for(loop_index=0;loop_index<number_in_x;loop_index++)
    {
        xvalues[loop_index]=xmin+loop_index*xstep;
        yvalue[loop_index]=function_gauss(xvalues[loop_index]);
    }//初始状态的设定
   gStyle->SetCanvasPreferGL(true);
   gStyle->SetFrameFillColor(10);
   TCanvas *c1 = new TCanvas("c1");
   c1->SetFillColor(17);
   TTimer *timer = new TTimer(20);
   //timer->SetCommand("Animate(xvalues,yvalue,yyvalue,number_in_x,factor_left,factor_right,graphtest)");
   timer->SetCommand("Animate()");
    timer->TurnOn();
}
//void Animate(Double_t *xvalues,Double_t * yvalue,Double_t * yyvalue,Int_t  number_in_x,Double_t factor_left,Double_t factor_right,TGraph *graphtest)
void Animate()
{
    if (!gROOT->GetListOfCanvases()->FindObject("c1")) return;
    loop_index=0;
    loop_jndex=0;
    yyvalue[0]=yvalue[0];
    yyvalue[number_in_x]=yvalue[number_in_x];
    for(loop_index=1;loop_index<number_in_x-1;loop_index++)
    {
        yyvalue[loop_index]=factor_left*yvalue[loop_index]+factor_right*yvalue[loop_index+1];//迎风法计算下一个时刻
    }
    for(loop_index=0;loop_index<number_in_x;loop_index++)
    {
        yvalue[loop_index]=yyvalue[loop_index];
    }
    graphtest = new TGraph(number_in_x,xvalues,yvalue);
    graphtest->GetXaxis()->SetLimits(-15,1);
    graphtest->GetHistogram()->SetMaximum(1.0);
    graphtest->GetHistogram()->SetMinimum(0);
    graphtest->Draw("AC");
    gPad->Modified();
    gPad->Update();
}
