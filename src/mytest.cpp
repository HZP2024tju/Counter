#include "matplotlibcpp.h"
#include <iomanip>
#include <vector>
#include <math.h>
#include "counter.hpp"
#include "my_fun.hpp"
#define SCOUT_OMP_NUM_THREADS 4 //线程个数
#define plot_den 12 //绘图分辨率
#define M 10 //测试点
int main(int argc,char**argv) {
    std::srand(1);//随机吧也许
    if (argc-1 != 11)
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
    int devide = std::stod(argv[9]); 
    int up = std::stod(argv[10]); 
    double rand_set = std::stod(argv[11]); //噪声
    if (b<a)
    std::swap(a,b);
    target_func fun(c,d,e,f);//准备拟合目标函数: sin cos

//常规制图
    int bigger = 2;
    const int range = static_cast<int>(b-a) +2*bigger;//看前后一格的情况
    std::vector<double> x(range*plot_den);
    std::vector<double> nx(n);//采样点
    std::vector<double> ny(n);//采样点
    std::vector<double> y_target(range*plot_den);
    std::vector<double> y_legend(range*plot_den);
    std::vector<double> y_base(range*plot_den);
    std::vector<double> y_ransac(range*plot_den);

    double plot_x=a-bigger;//制图用的x
    double plot_step = 1.0/static_cast<double>(plot_den);//哦真的牛批这里必须是1.0 每次移动多远
    
    for (int i = 0; i < range*plot_den; i++)
    {
        x[i] = plot_x;
        y_target[i] = fun.fun_use(plot_x,choose);
        plot_x += plot_step;
    }//准备原函数制图


      //分段制图
    double devide_step=(b-a)/static_cast<int>(devide);
    int per_seg = std::max(up+1,static_cast<int>(devide_step * plot_den));   // 每段点数，一次算好
    if (per_seg < 2) per_seg = 2;    
    int n_seg = std::max(up+1,1+n/devide);                          // 保护

    int total = devide * per_seg;
    
    std::vector<double> nx_devide(n_seg);
    std::vector<double> ny_devide(n_seg);//拟合采样
    std::vector<double> x_devide_final(total);
    std::vector<double> y_devide_final(total);//制图
    std::vector<double> x_devide(per_seg);
    std::vector<double> y_devide(per_seg);//采样
    std::cout << "每段采样点个数" << nx_devide.size()<<std::endl;
    for (int i = 0; i < devide; i++)
    {
        plot_x = a+i*devide_step;
        for (int j = 0; j < per_seg; j++)
        {
            x_devide[j] = plot_x;
            plot_x += plot_step;//没毛
        }

        counter::qieshell_particular(a+i*devide_step,a+(i+1)*devide_step,nx_devide,n_seg);//切比雪夫0点采样
        fun.get_ny(nx_devide,ny_devide,choose);
        for (int i = 0; i < ny_devide.size(); i++)
        {
            ny_devide[i]+=static_cast<double>(10.0-std::rand()%20)*rand_set;//噪声影响采样点
        }
        counter::base_like(a+i*devide_step,a+(i+1)*devide_step,up,nx_devide,ny_devide,x_devide,y_devide);

        for (int j = 0; j < per_seg; j++)
        {
            x_devide_final[j+i*per_seg] = x_devide[j];
            y_devide_final[j+i*per_seg] = y_devide[j];
        }
    }//
    
    matplotlibcpp::figure_size(1200, 780);
    matplotlibcpp::ylim(*std::min_element(y_target.begin(), y_target.end())-2,*std::max_element(y_target.begin(), y_target.end()) +2);
    matplotlibcpp::xlim(a-bigger/2,b+bigger/2);
    matplotlibcpp::title("result");

    matplotlibcpp::named_plot("devide",x_devide_final, y_devide_final,"b--"); 
    matplotlibcpp::named_plot("target", x, y_target,"r"); 
    matplotlibcpp::show();
    _exit(0); 
}