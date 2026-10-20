#ifndef __COUNTER_HPP__
#define __COUNTER_HPP__

#include "eigen3/Eigen/Eigen"
#include <random>
#include <unordered_set>
#include <omp.h>

namespace counter
{

//切比雪夫采样 输入区间 0点结果 采样个数
void qieshell_particular(double a,double b,std::vector<double> &res,int n)
{
    double pi = acos(-1.0);
    res.resize(n);//不可以resize否则push_back会错误
    for (double i = 1; i <= n; i++)
    {
        res[i-1] =((2.0*cos( pi*(2*i-1)/(2*n)) *(b -a))+b+a)/2.0;//获取[-1,1]内0点 并扩增
    }
    return;
}

//勒让德插值 p_n在[-1,1]有 n个0点
void legend_particular(double a,double b,std::vector<double> &res, int n)
{
    n++;//实际求n个插值点就要求到P_n 也就是需要求到n+1
    Eigen::MatrixXd vP = Eigen::MatrixXd::Zero(n, n);;//勒让德多项式的系数矩阵//P_0->P_n
    vP(0,0) = 1;
    vP(1,1) = 1;//第0 1行已完成
    Eigen::RowVectorXd row_temp;
    for (int i = 2; i <= n-1; i++)
    {//递推这一块
        row_temp = vP.row(i-1);
        vP.row(i)<<0.0,row_temp.head(n-1);//将上一行向右移动一格后注入这里
        vP.row(i) *=2*(i-1) +1.0;
        vP.row(i) -=(vP.row(i-2)*(i-1));
        vP.row(i) /=i;
    }

    //下面求0点
    Eigen::MatrixXd with = Eigen::MatrixXd::Zero(n-1,n-1);//伴随矩阵
    for (int i = 0; i < n-1; i++)
    {//没招了 最不优雅的方式
        with(0,i) = -vP(n-1,n-2-i)/vP(n-1,n-1);
    }
        
    for (int i = 1; i < n-1; i++)
    {
        with(i,i-1) = 1;
    }

    Eigen::EigenSolver<Eigen::MatrixXd> solver(with);
    Eigen::VectorXd eig = solver.eigenvalues().real();//求解特征值 也就是点 这里不可能为虚数//eig里面就是[-1,1]内的零点
    res.reserve(n-1);
    for (int i = 0; i < eig.size(); i++)
    {
        res[i] = eig(i)*(b-a)/2 +(b+a)/2;
    }
    

/*完成作业特别强
    std::cout<<'\n'<<"解得零点"<<eig.size()<<"个,分别为"<< std::endl;
    for (int i = 0; i < eig.size(); i++)
    {
        std::cout<<std::left  << std::setw(5)<<eig(i)<<' ';
    }
    std::cout<<std::endl;

    std::cout<<"\n"<<"下面开始输出勒让德矩阵的系数"<<std::endl;
    for (int i = 1; i <= n; i++)
    {
        for (int j = 0; j < i; j++)
        {
            std::cout<<std::left  << std::setw(8)<<vP(i-1,j)<<' ';
        }
        std::cout<<"\n";
    }
*/ 
}

//求P_n的值
double legend_P(double x,int n)
{//递推 避免展开P从而被浮点误差影响
    if (n == 0)
    {
        return 1;
    }
    else if(n == 1)
    {
        return x;
    }
    double P1 = 1;
    double P2 = x;
    double P3;
    for (int i = 2; i <=n; i++)
    {
        P3 = ((2*i-1)*x*P2 - (i-1)*P1)/i;
        P1 = P2;
        P2 = P3;
    }
    return P3;
}
//勒让德拟合
void legend_like(double a,double b,int den,std::vector<double> x,std::vector<double> y,std::vector<double> plot_x,std::vector<double> &res)
{           //  拟合多项式次数 采样点数据                                                     制图           
    if (den<1)
    {
        return;
    }
    for (int i=0 ;i< x.size();i++)
    {
        x[i] = ( (2*x[i]-a-b)/(b-a) );//变换采样点
    }
    den++;//求P_den 也就是需要求到den+1
    if (x.size() < den) {
        std::cerr << "采样点太少，无法拟合 " << den-1 << " 次多项式"<<std::endl;
        return;
    }
    
    Eigen::MatrixXd vP(den,den);
    Eigen::VectorXd vb(den);
    double count_temp;
    for (int i = 0; i < den; i++)
    {
        for (int j = 0; j <= i ; j++)
        {
            //计算函数内积的值
            count_temp = 0;
            for (int k = 0; k < x.size(); k++)
            {
                count_temp += legend_P(x[k],i) * legend_P(x[k],j) ;
            }
            vP(i,j) = count_temp;
            if (j != i)
            {
                vP(j,i) = count_temp;//填入矩阵
            }
            
        }

        //计算函数内积的值
        count_temp = 0;
        for (int k = 0; k < x.size(); k++)
        {
            count_temp += legend_P(x[k],i) * y[k] ;
        }
        vb(i) = count_temp;        
    }    

    res.resize(plot_x.size());  
    Eigen::VectorXd va =  vP.colPivHouseholderQr().solve(vb);
    for (int i = 0; i < plot_x.size(); i++)
    {
        res[i]=0;
        for (int j = 0; j < den; j++)
        {
            res[i] +=  va(j)*legend_P((2*plot_x[i]-a-b)/(b-a),j);
        }
        
    }
}//问心无愧了

//回归基本功拟合(1,x,x^2...)
void base_like(double a,double b,int den,std::vector<double> x,std::vector<double> y,std::vector<double> plot_x,std::vector<double> &res)
{           //  拟合多项式次数 采样点数据                                                     制图           
    if (den<1)
    {
        return;
    }

    den++;
    if (x.size() < den) {
        std::cout << "采样点太少，无法拟合 " << den-1 << " 次多项式"<<std::endl;
        return;
    }

    Eigen::MatrixXd vP(den,den);
    Eigen::VectorXd vb(den);
    double count_temp;
    for (int i = 0; i < den; i++)
    {
        for (int j = 0; j <= i ; j++)
        {
            //计算函数内积的值
            count_temp = 0;
            for (int k = 0; k < x.size(); k++)
            {
                count_temp += pow(x[k],i) * pow(x[k],j) ;
            }
            vP(i,j) = count_temp;
            if (j != i)
            {
                vP(j,i) = count_temp;//填入矩阵
            }
        }

        //计算函数内积的值
        count_temp = 0;
        for (int k = 0; k < x.size(); k++)
        {
            count_temp += pow(x[k],i) * y[k] ;
        }
        vb(i) = count_temp;        
    }    

    Eigen::VectorXd va =  vP.colPivHouseholderQr().solve(vb);
    for (int i = 0; i < plot_x.size(); i++)
    {
        res[i]=0;
        for (int j = 0; j < den; j++)
        {
            res[i] +=  va(j)*pow(plot_x[i],j);
        }
        
    }
    
}

void RANSAC_like(double a,double b,int den,std::vector<double> x_in,std::vector<double> y_in,std::vector<double> plot_x,std::vector<double> &res ,int max_step,double ransac,double judge_in_point)
{//已知区间  采样点 制图点 最大迭代步长和内点几率 内点判定阈值 使用RANSAC方法进行迭代
    if (ransac >1 || ransac <0 ||judge_in_point<0||max_step<0)
    {
        std::cout<<"RANSAC错误"<<std::endl;
        return;
    }
    
    std::random_device rd;
    std::mt19937 gen(rd());
    // 创建均匀整数分布
    std::uniform_int_distribution<> dist(0,x_in.size()-1);//防止最为索引时越界
    int in_count_max=0;
    int in_count;
    std::vector<int> in_node_best;
    for (int z = 0;; z++)
    {//这里多线程处理可以先算出多个 最小 比较 最后统一制图
        in_count = 0;
        std::unordered_set<double> seen;
        std::vector<double> x;
        std::vector<double> y;
        std::vector<double> y_res;
        std::vector<int> in_node;
        int i=0;
        int random;
        while(x.size()!=den+1 && i<x_in.size()*3)
        {
            random=dist(gen);
            if (seen.insert(x_in[random]).second) 
            {  // 插入成功说明没出现过
                x.emplace_back(x_in[random]);
                y.emplace_back(y_in[random]);//去重采集拟合点
            }
            i++;
        }//取点
        if (x.size() != den + 1) 
        {
            continue;   // 跳过本次迭代
        }
        legend_like(a,b,den,x,y,x_in,y_res);//随机不重复采样 拟合

        //算内点
        for(i = 0; i < x_in.size(); i++)
        {
            if (std::abs(y_res[i] - y_in[i]) < judge_in_point)
            {
                in_count++;
                in_node.emplace_back(i);//记录所有内点的标签
            }
        }            
        if (in_count > in_count_max)
        {
            in_count_max = in_count;  
            in_node_best = in_node;
            if (static_cast<double>(in_count)/x_in.size() > ransac)
            {
                std::vector<double> best_nx;
                std::vector<double> best_ny;
                for (int i = 0; i < in_count; i++)
                {
                    best_nx.emplace_back(x_in[in_node_best[i]]);
                    best_ny.emplace_back(y_in[in_node_best[i]]);
                }
                legend_like(a,b,den,best_nx,best_ny,plot_x,res);//拟合    
                break;
            }
        }

        if(z == max_step)
        {
            std::vector<double> best_nx;
            std::vector<double> best_ny;
            for (int i = 0; i < in_count_max; i++)
            {
                best_nx.emplace_back(x_in[in_node_best[i]]);
                best_ny.emplace_back(y_in[in_node_best[i]]);
            }
            legend_like(a,b,den,best_nx,best_ny,plot_x,res);//拟合    
            break;
        }
    }
}
//随机采样
void random_particular(double a,double b,int plot_den,std::vector<double> &res,int n)
{
    std::random_device rd;
    std::mt19937 gen(rd());
    // 创建均匀整数分布
    std::uniform_int_distribution<> dist( a*plot_den , b*plot_den );
    res.reserve(n);//不可以resize否则push_back会错误
    for (int i = 0; i <n; i++)
    {
       res[i] =  (static_cast<double>(dist(gen)))/static_cast<double>(plot_den);//好舒服啊
    }
    std::sort(res.begin(),res.end());
    return;
}

//平均采样
void avg_particular(double a,double b,std::vector<double> &res,int n)
{
    if (n == 1)
    {
        res[0] = a;
        return;
    }
    res.reserve(n);//不可以resize否则push_back会错误
    double step = (b-a)/(n-1);
    for (double i = 0; i <n; i++)
    {
       res[i] = a+step*i;
    }
    return;
}
//知道区间 采样点数 目标函数 就能自动拟合并返回参数向量
void vandmod_inter(std::vector<double> x,int n,std::vector<double>& plot_x,std::vector<double>& res ,std::vector<double> y)
{
    Eigen::MatrixXd vA(n,n);
    Eigen::VectorXd vb(n);
    double p;
    for (int i = 0; i < n; i++)
    {//n个采样点 相应的最高次就是n-1次
        p = 1.0;
        for (int j = 0; j < n; j++)
        {//其实这里引用列的话等效内积   后面打表看看耗时
            vA(i,j) = p;//构造范德蒙德阵
            p *= x[i];
        }
        vb(i) =  y[i];//在这里填入采样点
    }

    Eigen::VectorXd va =  vA.colPivHouseholderQr().solve(vb);//忽略稀疏部分 解出va参数向量  SVD很慢

    Eigen::VectorXd vx(n);
    res.resize(plot_x.size());
    double ven_x;
    for (int i = 0; i < plot_x.size(); i++)
    {//static应该不会有伪共享   画出区间内的情况
        ven_x = 1;
        for (int j = 0; j < n; j++)
        {
            vx(j) = ven_x;
            ven_x *= plot_x[i];
            
        }
        res[i] = (va.dot(vx)); //va和vx作内积
    } 
}

//区间 采样点 拟合方程 结果要求点集个数 返回结果
void lagrange_inter(std::vector<double> x,int n,std::vector<double>& plot_x,std::vector<double>& res ,std::vector<double> y)
{
    Eigen::VectorXd vy(n);
    for (int i = 0; i < n; i++)
    {
        vy(i) = y[i];
    }
    
    //qieshell_particular(a,b,vx,n);    

    //计算制图所需
    res.resize(plot_x.size());
    for (int i = 0; i < plot_x.size(); i++)
    {
        res[i] = 0;
        double temp_lk;
        for (int j = 0; j < n; j++)
        {
            temp_lk =1 ;
            for (int k = 0; k < n; k++)
            {
                if (j == k)
                {
                    continue;
                }
                else
                {
                    temp_lk *= (plot_x[i] -x[k])/(x[j] -  x[k]); //计算基函数
                }
            }
            res[i] += vy(j)*temp_lk;
        }
    }
   // #pragma omp parallel num_threads(SCOUT_OMP_NUM_THREADS)
}

//plot_x:制图用点  res:结果
void newton_inter(std::vector<double> x_in,int n,std::vector<double>& plot_x,std::vector<double>&res, std::vector<double> y_in)
{               //这里不太好引用x
    if (n == 1)
    {
        return;
    }
    std::unordered_set<double> seen;
    std::vector<double> x;
    std::vector<double> y;
    for (int i=0 ;i< x_in.size();i++)
    {
        if (seen.insert(x_in[i]).second) 
        {  // 插入成功说明没出现过
            x.emplace_back(x_in[i]);
            y.emplace_back(y_in[i]);//去重
        }
    }
    n = x.size();

    std::vector<double> avg_loss;//均差
    avg_loss.resize(n* (n + 1) /2);
    for (int i = 0; i < n; i++)
    {
        avg_loss[i] = y[i];//先计算f(x)
    }

    for (int i = n-1; i > 0; i--)
    {//循环 n-1次
        for (int k = 0; k < i; k++)
        {
            avg_loss[ ((n-i) *(n+i+1)/2) + k ] = (-avg_loss[ (n-i-1) *(n + i+2)/2 +k ] + avg_loss[ (n-i-1) *(n + i+2)/2 +k+1 ])
            /(x[n-i + k] - x[k]) ; //计算均差
        }
    }

    //计算制图所需
    res.resize(plot_x.size());
    for (int j = 0; j < plot_x.size(); j++)
    {
        res[j] = avg_loss[0];
        for (int i = 1; i < n; i++)
        {
            double temp_ = avg_loss[i * (n + n-i+1)/2] ;
            for (int k = 0; k < i; k++)
            {
                temp_ *= (plot_x[j] - x[k]);
            }
            
            res[j] += temp_;
        }
    }
    return;
}

double quarter_adevance_intergration(double a,double b,std::function<double(double)>& fun,int n)
{//n:分段数
    double h = (b-a)/static_cast<double>(n);
    double x_k = a+h;
    double res=0;
    for (int i = 1; i < n; i++)
    {//循环n-1次
        res +=fun(x_k);
        x_k += h;
    }
    res *=2;
    res = h*(res+fun(a) + fun(b))/2;
    return res;
}

double sinpcen_adevance_intergration(double a,double b,std::function<double(double)>& fun,int n)
{//n:分段数
    double h = (b-a)/static_cast<double>(n); //n=1是另一种情况
    double x_k = a+h;
    double res = fun(x_k-h/2) * 4;
    for (int i = 1; i < n; i++)
    {//循环n-1次
        res +=fun(x_k)*2;
        x_k +=h/2;
        res +=fun(x_k)*4;
        x_k+=h/2;
    }
    res = h*(res+fun(a)+fun(b))/6;
    return res;
}

double dronbag_adevance_intergration(double a,double b,std::function<double(double)>& fun,double thegema)
{
    double h = b-a;//当前步长
    std::vector<double> T_n;        //这里其实不用保存全部的
    T_n.emplace_back((fun(a)+fun(b))*h/2); //最新的结果
    double T_n_last = (fun(a)+fun(b))*h/2;//上一次结果      
    int row = 0;//求到第几行了
    int col = 0;//求到第几个T_n了
    do
    {
        col++;
        if (col>row)
        {
            row++;//步进一行
            if (row >= 2)
            {
                T_n.erase(T_n.begin(), T_n.begin() + row-1);//神了
            }
            
            T_n_last = T_n.back();
            col = 0;
            //求行首元素
            int count_x = 1<<row-1;//总共有几个中点
            double x_temp;
            double T_2n = 0;
            double x_k = a +h/2.0;
            for (int i = 0; i < count_x; i++)
            {
                T_2n += fun(x_k);
                x_k += h;
            }
            T_2n = T_2n *h/2;
            T_2n += T_n[0] /2.0;//没毛
            T_n.emplace_back(T_2n);
            h /=2.0;
        }
        else
        {
            double temp;
            int count_T = T_n.size();
            temp = (T_n[T_n.size()-1]*pow(4,col) - T_n[T_n.size() - row-1])/(pow(4,col)-1);
            //检查
            if (T_n.size() - row-1<0)
            {
                std::cout<<"错误"<<std::endl;
                return 0;
            }
            
            T_n.emplace_back(temp);
        }       
        if (row > 30)
        {
            std::cerr << "达到最大迭代次数，可能未收敛！" << std::endl;
            break;
        } 
    }
    while( !(col == row &&std::abs(T_n.back() - T_n_last)<thegema) );
    return T_n.back();
}
}
#endif