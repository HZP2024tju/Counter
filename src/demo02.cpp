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
int main(int argc,char**argv)
{
    if (argc-1 != 8)
    {
        std::cout<<"收到参数不对"<<std::endl;
        return 0;
    }
   
    int n = std::stoi(argv[1]);//先看要几个插值点
    if (n<=0)
    {
        std::cout<<"采样点个数不对"<<std::endl;
        return 0;
    }

    std::cout<<"收到参数个数"<<argc<<" n:"<<n<< std::endl;
    double a = std::stod(argv[2]);//左边界
    double b = std::stod(argv[3]);
    double c = std::stod(argv[4]);
    double d = std::stod(argv[5]);
    double e = std::stod(argv[6]);
    double f = std::stod(argv[7]); //好烦 以后改成launch好了

    int choose = std::stod(argv[8]); 
    if (b<a)
    std::swap(a,b);

    target_func fun(c,d,e,f);//准备拟合目标函数: sin cos
    std::function<double(double)> current_fun = ([fun,choose](double x){return fun.fun_use(x,choose);});

    //制图
    int bigger = 2;
    const int range = static_cast<int>(b-a) +2*bigger;//看前后一格的情况
    std::vector<double> x(range*plot_den);
    std::vector<double> nx(n);//采样点
    std::vector<double> ny(n);//采样点
    std::vector<double> y_target(range*plot_den);
    std::vector<double> y_qiebi(range*plot_den);
    std::vector<double> y_random(range*plot_den);
    std::vector<double> y_legend(range*plot_den);
    double plot_x=a-bigger;//制图用的x
    double plot_step = 1.0/plot_den;//哦真的牛批这里必须是1.0 每次移动多远
    
    for (int i = 0; i < range*plot_den; i++)
    {
        x[i] = plot_x;
        y_target[i] = current_fun(plot_x);
        plot_x += plot_step;
    }

    counter::qieshell_particular(a,b,nx,n);//切比雪夫0点采样
    fun.get_ny(nx,ny,choose);
    counter::newton_inter(nx,n,x,y_qiebi,ny);

    counter::legend_particular(a,b,nx,n);//勒让德采样
    fun.get_ny(nx,ny,choose);
    counter::newton_inter(nx,n,x,y_legend,ny);

    counter::random_particular(a,b,plot_den,nx,n);//随机采样
    fun.get_ny(nx,ny,choose);
    counter::newton_inter(nx,n,x,y_random,ny);

//在一个进程内不可以开两个python线程
    matplotlibcpp::figure_size(1200, 780);

    matplotlibcpp::named_plot("target", x, y_target,"r");
    matplotlibcpp::named_plot("qiebi",x, y_qiebi,"y--");  
    matplotlibcpp::named_plot("legend",x, y_legend,"g--");  
    matplotlibcpp::named_plot("random",x, y_random,"b--"); 
    // 添加标题
    matplotlibcpp::ylim(-4,4);
    matplotlibcpp::xlim(-2,2);
    matplotlibcpp::title("result");
    // 启用图例
    matplotlibcpp::show(); //比较好笑的一点是一旦数据量过多会直接段错误  你开源给我写好了啊
    
    return 0;
}  