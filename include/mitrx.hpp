#ifndef __MITRX_HPP__
#define __MITRX_HPP__
#include <Eigen/Dense>
#include <iostream>
namespace mitrx
{
//高斯消元解 Ax=b       Eigen：对喽，千万别用我的求解器
bool mitrx_main_solve(Eigen::MatrixXd A,Eigen::VectorXd& res,Eigen::VectorXd b)
{
    double det=1;
    double col_best;
    int change_row;
    for (int i = 0; i < A.cols()-1; i++)
    {
        change_row=-1;
        col_best=std::abs(A(i,i)); //找最大
        for (int j = i+1; j < A.cols(); j++)
        {
            if (j==i)
            {
                continue;
            }
            
            if ((col_best) < std::abs(A(j,i)))
            {
                col_best= std::abs(A(j,i));
                change_row = j;
            }
        }
        if (change_row !=-1)
        {
            A.row(i).swap(A.row(change_row));//互换 
            std::swap(b(i),b(change_row));//什么？
        }
        if (A(i,i)==0)
        {
            std::cout<<"错误，没有唯一解"<<std::endl;//这里只考虑满秩
            return 0;
        }
        double temp4;
        for (int j = i+1; j < A.cols(); j++)
        {
            temp4=(A(j,i)/A(i,i));
            A.row(j) -= A.row(i)*temp4;
            b(j) -=b(i)*temp4;
        }
    } //这里就让A变为上三角阵了          //solStride是跨越内存数，可能会因为内存对齐而改变

    b(A.rows()-1) /= A(A.rows()-1, A.rows()-1);
    double temp3;
    for (int i =  A.rows()-2 ; i >=0; i--)
    {
        temp3=0;
        for (int j = A.rows()-1; j > i; j--)
        {
            temp3 += A(i,j)*b(j);
        }
        b(i) = (b(i)-temp3)/A(i,i);
    }
    res=b;
    return 1;
}//void
//LU分解
void mitrx_LU_solve(Eigen::MatrixXd& A,Eigen::VectorXd& res,Eigen::VectorXd& b, Eigen::MatrixXd* L_in = nullptr,Eigen::MatrixXd* U_in = nullptr)
{
    if (A.rows()!=A.cols()||A.rows()!=b.rows())
    {
        std::cout<<"错误，非方阵"<<std::endl;//这里只考虑满秩
        return;
    }
    Eigen::MatrixXd L = Eigen::MatrixXd::Identity(A.rows(), A.rows()); // L初始化为主对角线为1的单位阵
    Eigen::MatrixXd U = Eigen::MatrixXd::Zero(A.rows(), A.rows()); 

    //L第一列 U第一行
    for (int i = 0; i < A.rows(); i++)
    {
        U(0,i)=A(0,i);
        if (i==0)
        {
            continue;
        }
        L(i,0)=A(i,0)/U(0,0);
    }

    double temp;
    for (int e = 1; e < A.rows(); e++)
    {
        for (int i = e; i <A.rows(); i++)
        {
            temp=0;
            int j =0;
            for (j = 0; j < e; j++)
            {
                temp+= U(j,i)*L(e,j);
            }
            U(e,i)=A(e,i)-temp;

            temp=0;

            for (j = 0; j < e; j++)
            {
                temp+= U(j,e)*L(i,j);
            }
            L(i,e)=(A(i,e)-temp)/U(e,e);
        }
    }

    Eigen::VectorXd y;
    res.resize(A.rows());
    y.resize(A.rows());
    y(0)=b(0);
    for (int i = 1; i <A.rows(); i++)
    {
        temp=0;
        int j =0;
        for (j = 0; j < i; j++)
        {
            temp+= y(j)*L(i,j);
        }
        y(i)=b(i)-temp;
    }
    res(A.rows()-1)=y(A.rows()-1)/U(A.rows()-1,A.rows()-1);
    for (int i = A.rows()-2; i >=0; i--)
    {
        temp=0;
        int j =0;
        for (j = i+1; j < A.rows(); j++)
        {
            temp+= res(j)*U(i,j);
        }
        res(i)= (y(i)-temp)/U(i,i);
    }

    if(L_in != nullptr)
    {
        *L_in = std::move(L);
    }
    if(U_in != nullptr)
    {
        *U_in = std::move(U);
    }
}
}//namespace 专业
#endif