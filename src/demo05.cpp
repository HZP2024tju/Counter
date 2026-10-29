#include "iostream"
#include <iomanip>
#include "thread"
#include <omp.h>
#include "mitrx.hpp"
#define SCOUT_OMP_NUM_THREA 4

int main(int argc,char**argv)
{
    Eigen::MatrixXd A(9, 9);
    Eigen::VectorXd b(9);

    A <<  31, -13,   0,   0,   0, -10,   0,   0,   0,
         -13,  35,  -9,   0, -11,   0,   0,   0,   0,
           0,  -9,  31, -10,   0,   0,   0,   0,   0,
           0,   0, -10,  79, -30,   0,   0,   0,  -9,
           0,   0,   0, -30,  57,  -7,   0,  -5,   0,
           0,   0,   0,   0,  -7,  47, -30,   0,   0,
           0,   0,   0,   0,   0, -30,  41,   0,   0,
           0,   0,   0,   0,  -5,   0,   0,  27,  -2,
           0,   0,   0,  -9,   0,   0,   0,  -2,  29;

    // 填入向量 b 的数据
    b << -15, 
          27, 
         -23, 
           0, 
         -20, 
          12, 
          -7, 
           7, 
          10;

    Eigen::MatrixXd A2(9, 9);
    Eigen::VectorXd b2(9);

    // 按行填入矩阵 A 的数据
    A2 <<   30,   33,   -43,   -11,   -38,   -29,    37,    28,    23,
         -480, -523,   644,   128,   621,   480,  -618,  -489,  -329,
           60,  266, -1862, -1991,   464,   546,  -968, -1567,  1652,
          540,  624,  -782,   290,  -893,   123,   567,     5,  -122,
         -450, -675,  2245,  2326, -1512,  1230,  -822,   129,  -189,
         -300, -120, -1114, -1295,  1946,   302,  -376, -1540,  -609,
         1080,  998,   508,  2460, -1628, -1358,  2896,  2828, -2002,
        -1080, -1408, 3340,  2267,    21, -1202,   866, -2690, -1351,
         -300, -435,  1594,  1685,   340,  2279,   -27,  2917, -2336;

    // 填入向量 b 的数据
    b2 <<   188, 
         -3145, 
         -4994, 
           680, 
          7845, 
          1876, 
          9712, 
        -11599, 
         10127;

    Eigen::VectorXd res,res1,res2,x,x1,x2;
    Eigen::MatrixXd L,U;
    mitrx::mitrx_main_solve(A,res,b);
    x = A.partialPivLu().solve(b);

    std::cout << "要求一:\n";
    std::cout << "高斯消元解得 x:\n" << res << "\n";
    std::cout << "标准答案 x:\n" << x << "\n";
    
    std::cout << "\n================================" << "\n\n";
    A=Eigen::MatrixXd::Random(20,20);
    b=Eigen::VectorXd::Random(20);
    while (mitrx::mitrx_main_solve(A,res,b)==0)
    {
        A=Eigen::MatrixXd::Random(20,20);
    }
    x = A.partialPivLu().solve(b);

    std::cout << "高斯消元解得 x:\n" << res << "\n";
    std::cout << "标准答案 x:\n" << x << "\n";

    std::cout << "\n================================" << "\n\n";

    mitrx::mitrx_LU_solve(A2,res,b2,&L,&U);
    x = A2.partialPivLu().solve(b2);
    std::cout << "要求二：\n";
    std::cout << "LU解得 x:\n" << res << "\n";
    std::cout << "标准答案 x:\n" << x << "\n";
    std::cout << "L:\n" << L << "\n";
    std::cout << "U:\n" << U << "\n";//好舒服啊

    std::cout << "\n================================" << "\n\n";
    A2=Eigen::MatrixXd::Random(20,20);
    b2=Eigen::VectorXd::Random(20);
    while (mitrx::mitrx_main_solve(A2,res,b2)==0)
    {
        A2=Eigen::MatrixXd::Random(20,20);
    }
    mitrx::mitrx_LU_solve(A2,res,b2,&L,&U);
    x = A2.partialPivLu().solve(b2);

    std::cout << "LU解得 x:\n" << res << "\n";
    std::cout << "标准答案 x:\n" << x << "\n";
    std::cout << "L:\n" << L << "\n";
    std::cout << "U:\n" << U << "\n";//好舒服啊
    std::cout << "\n================================" << "\n\n";
}