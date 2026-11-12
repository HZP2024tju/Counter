#include "iostream"
#include <iomanip>
#include "thread"
#include <omp.h>
#include "matplotlibcpp.h"
#include "counter.hpp"
#include "my_fun.hpp" 
#define SCOUT_OMP_NUM_THREADS 4 //线程个数
#define plot_den 20 //绘图分辨率
#define M 10 //测试点
//最通透的一次作业
int main(int argc,char**argv)
{
    std::srand(std::time(nullptr));//随机吧也许
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
    std::function<double(double)> current_fun = ([fun,choose](double x){return fun.fun_use(x,choose);});

    //制图
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
        y_target[i] = current_fun(plot_x);
        plot_x += plot_step;
    }//准备原函数制图

    

    counter::legend_particular(a,b,nx,n);//0点采样
    fun.get_ny(nx,ny,choose);
    for (int i = 0; i < ny.size(); i++)
    {
        ny[i]+=static_cast<double>(10.0-std::rand()%20)*rand_set;//噪声影响采样点
    }
    counter::base_like(a,b,up,nx,ny,x,y_base); //基础拟合
    counter::legend_like(a,b,up,nx,ny,x,y_legend); //勒让德拟合
    counter::RANSAC_like(a,b,up,nx,ny,x,y_ransac,nx.size()*2000,0.9,0.1);//RANSAC

    //==求二乘误差==
    double res_2;
    std::vector<double> x_2(M);
    std::vector<double> y_2(M);
    counter::random_particular(a,b,plot_den,x_2,M);//随机采样测试误差
    counter::base_like(a,b,up,nx,ny,x_2,y_2);
    res_2=0;
    for (int i = 0; i < x_2.size(); i++)
    {
        res_2 += y_2[i]*y_2[i] - pow(fun.fun_use(x_2[i],choose),2);
    }
    res_2 = pow(std::abs(res_2),0.5);
    std::cout<<'\n'<<"标准2范数为"<<res_2<<std::endl;

    counter::legend_like(a,b,up,nx,ny,x_2,y_2);
    res_2=0;
    for (int i = 0; i < x_2.size(); i++)
    {
        res_2 += y_2[i]*y_2[i] - pow(fun.fun_use(x_2[i],choose),2);
    }
    res_2 = pow(std::abs(res_2),0.5);
    std::cout<<'\n'<<"勒让德2范数为"<<res_2<<std::endl;

    counter::RANSAC_like(a,b,up,nx,ny,x_2,y_2,nx.size()*2000,0.9,0.1);
    res_2=0;
    for (int i = 0; i < x_2.size(); i++)
    {
        res_2 += y_2[i]*y_2[i] - pow(fun.fun_use(x_2[i],choose),2);
    }
    res_2 = pow(std::abs(res_2),0.5);
    std::cout<<'\n'<<"RANSAC2范数为"<<res_2<<std::endl;


    matplotlibcpp::figure_size(1200, 780);
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

    //M个测试点内测量范数
    res_2=0;
    int rand_temp = rand()%x_devide_final.size();
    for (int i = 0; i < M; i++)
    {
        res_2 += y_devide_final[rand_temp]*y_devide_final[rand_temp] - pow(fun.fun_use(x_devide_final[rand_temp],choose),2);
        rand_temp = rand()%x_devide_final.size();
    }
    res_2 = pow(std::abs(res_2),0.5);
    std::cout<<'\n'<<"分段采样2范数为"<<res_2<<std::endl;

//在一个进程内不可以开两个python线程
    matplotlibcpp::named_plot("devide",x_devide_final, y_devide_final,"b--"); 
    matplotlibcpp::named_plot("target", x, y_target,"r"); 
    matplotlibcpp::named_plot("legend",x, y_legend,"g--");  
    matplotlibcpp::named_plot("ransac",x, y_ransac,"y--"); 
    matplotlibcpp::named_plot("base",x, y_base,"b"); 

    // 添加标题
    matplotlibcpp::ylim(*std::min_element(y_target.begin(), y_target.end())-2,*std::max_element(y_target.begin(), y_target.end()) +2);
    matplotlibcpp::xlim(a-bigger/2,b+bigger/2);
    matplotlibcpp::title("result");
    // 启用图例
    matplotlibcpp::show(); //比较好笑的一点是一旦数据量过多会直接段错误  你开源给我写好了啊
    
    return 0;
}  