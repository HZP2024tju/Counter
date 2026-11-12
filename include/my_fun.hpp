#ifndef __MY_FUN_HPP__
#define __MY_FUN_HPP__
class target_func
{
public:
    target_func(double c=0,double d=0,double e=0,double f=0)
    : c_(c), d_(d), e_(e), f_(f)
    {
    }

    inline double fun_use(double x,int choose) const
    {//不用引用是线程安全的表现
        switch (choose)
        {
        case 0:
            return this->fun_sincos(x);
            break;
        case 1:
            return this->fun_xn(x);
            break;
        case 2:
            return this->fun_longe(x);
            break;
        case 3:
            return this->fun_x2sincos(x);
            break;
        case 4:
            return this->fun_x2sincos(x);
            break;
        }
    }

    void get_ny(std::vector<double> nx,std::vector<double> &ny,int choose)
    {
        ny.resize(nx.size());
        for (int i = 0; i < nx.size(); i++)
        {
            ny[i] = this->fun_use(nx[i],choose);
        }
    }

private:
    double c_,d_,e_,f_;
    inline double fun_sincos(double x) const
    {//感觉自己这样写老帅了
        return c_* sin(d_*x) + e_*cos(f_*x);
    } 

    inline double fun_xn(double x) const
    {
        return x*x;
    } 

    inline double fun_longe(double x) const
    {
        return 1.0/(1+25*x*x);
    }
    inline double fun_x2sincos(double x) const
    {//感觉自己这样写老帅了
        return x*x+c_* sin(d_*x) + e_*cos(f_*x);
    } 
    inline double fun_2xinx(double x) const
    {//感觉自己这样写老帅了
        return sqrt(x) * log(x);
    } 
};
#endif