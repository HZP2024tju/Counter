#include "iostream"
#include <iomanip>
#include "thread"
#include <omp.h>
#include "mitrx.hpp"

int main(int argc,char**argv)
{
    if (argc < 3)
    {
        std::cout << "参数太少" <<std::endl;
        return 0;
    }
    
    int max_step = std::stoi(argv[1]);
    double err = std::stod(argv[2]);

    Eigen::Matrix<double, 9, 9> A;
    Eigen::Matrix<double, 9, 1> b; 

    A <<  31, -13,   0,   0,   0, -10,   0,   0,   0,
         -13,  35,  -9,   0, -11,   0,   0,   0,   0,
           0,  -9,  31, -10,   0,   0,   0,   0,   0,
           0,   0, -10,  79, -30,   0,   0,   0,  -9,
           0, -11,   0, -30,  57,  -7,   0,  -5,   0,
           0,   0,   0,   0,  -7,  47, -30,   0,   0,
           0,   0,   0,   0,   0, -30,  41,   0,   0,
           0,   0,   0,   0,  -5,   0,   0,  27,  -2,
           0,   0,   0,  -9,   0,   0,   0,  -2,  29;

    b << -15, 27, -23, 0, -20, 12, -7, 7, 10;

    Eigen::VectorXd x_0 = Eigen::VectorXd::Zero(9);
    
    mitrx::SOR_step_solve(b,x_0,A,max_step,err,1);
    
    if(mitrx::SOR_step_solve(b,x_0,A,max_step,err,1))
    {
        std::cout << "高斯法解得:\n" << x_0 << "\n";
    }
    else
    std::cout <<"残差过大 无效\n"<<std::endl;

    std::cout << "\n================================" << "\n\n";
    double w;    
    for (int i = 1; i <= 99; i++)
    {
        w = i/50.0;
        x_0 = Eigen::VectorXd::Zero(9);
        if(mitrx::SOR_step_solve(b,x_0,A,max_step,err,w))
        {
            std::cout <<"w = "<<w<<' ' <<"SOR解得:\n" << x_0 <<std::endl;
        }
        else
        std::cout <<"w = "<<w<<' ' <<"残差过大 无效\n"<<std::endl;
    }
    
    Eigen::VectorXd x = A.partialPivLu().solve(b);
    std::cout << "标准答案 x:\n" << x << "\n";
}