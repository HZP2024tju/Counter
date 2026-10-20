#include "iostream"
#include <iomanip>
#include "thread"
#include <omp.h>
#include "matplotlibcpp.h"
#include "counter.hpp"
#include "my_fun.hpp"
#define SCOUT_OMP_NUM_THREADS 4 //线程个数
#define plot_den 10 //绘图分辨率

//最通透的一次作业
double quarter_integration(double a,double b,std::function<double(double)>& fun,double thegma)
{
//  复合梯形积分        这里需要求到来着 如果要迭代的话
    int n =1;
    double h = (b-a)/static_cast<double>(n);
    while(((b-a)*h*h)/(12*4*pow(a,1.5)) >thegma)
    {
        n *=2;
        h = (b-a)/static_cast<double>(n);
    }//迭代
    std::cout<<"共迭代"<<n<<"次"<<",步长"<<h<<"\n";
    return counter::quarter_adevance_intergration(a,b,fun,n);
}

int main(int argc,char**argv)
{
    if (argc-1 != 4)
    {
        std::cout<<"收到参数不对"<<std::endl;
        return 0;
    }
    double n = std::stod(argv[1]);
    std::cout<<"收到参数个数"<<argc<<" n:"<<n<< std::endl;
    double a = std::stod(argv[2]);//左边界
    double b = std::stod(argv[3]);
    int choose = std::stoi(argv[4]); 

    if (b<a)
    std::swap(a,b);

    target_func fun;
    std::function<double(double)> current_fun = ([&fun,choose](double x){return fun.fun_use(x,choose);});

    std::cout<<"梯形积分结果"<<counter::quarter_adevance_intergration(a,b,current_fun,n)<<std::endl;
    double res = quarter_integration(a,b,current_fun,0.0001);
    std::cout<<"迭代梯形积分结果"<<res<<std::endl;
    std::cout<<"辛普森积分结果"<<counter::sinpcen_adevance_intergration(a,b,current_fun,n)<<std::endl;
    res = counter::dronbag_adevance_intergration(a,b,current_fun,0.0001);
    std::cout<<"龙贝格积分结果"<<res<<std::endl;

    return 0;
}  