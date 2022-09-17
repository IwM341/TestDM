#include <iostream>
#include <utils>
using namespace std;

#ifdef TEST_INTERPOL
int main()
{
    Function::UniformGrid UG1(-1.0,1.0,11);
    Function::VectorGrid VG1 = Function::UniformGrid(-1.0,1.0,11);

    Function::LinearInterpolator Int;
    for(size_t i=0;i<10;++i){
        double x = ((double)rand())/RAND_MAX*2-1;
        std::cout << x <<"\t" << UG1.at(UG1.pos(x)) << "\t" <<VG1.at(VG1.pos(x)) << std::endl;
        print(Int(UG1,x),Int(VG1,x));
    }
    return 0;
}

#endif

#ifdef TEST_STRUCT_INTERPOL
int main()
{
    //size_t ind[2] = {1,2};
    //double weight[2] = {1.0,2.0};
    Function::Scheme<2> x{{1,2},{0.5,0.5}};
    print(x);
    print(sizeof(Function::Scheme<2>));
    print(sizeof(Function::Scheme<4>));
    return 0;
}
#endif
#ifdef TEST_GRID_FUNCTION
int main()
{
    Function::UniformGrid UG1(-1.0,1.0,11);
    Function::VectorGrid VG1(std::vector<double>({-100,-20,50,80,100}));
    auto f = [](double x){return x*x*x*x;};
    Function::__GridFunction<1,double,typeof(UG1),Function::LinearInterpolator> Fl(UG1,f);
    Function::__GridFunction<1,double,typeof(UG1),Function::CubicInterpolator> Fc(UG1,f);

    for(size_t i=0;i<10;++i){
        double x = ((double)rand())/RAND_MAX*2-1;
        std::cout << x <<  "\t" << f(x) << "\t" << Fl(x) << "\t" << Fc(x) << std::endl;
    }
    return 0;
}
#endif

#ifdef TEST_M_GRID_FUNCTION
int main()
{
    Function::UniformGrid UG1(-1.0,1.0,11);
    auto f = [](double x,double y){return x*y;};
    Function::GridFunction<2,double,typeof(UG1),Function::LinearInterpolator,
                               typeof(UG1),Function::LinearInterpolator>
            F(UG1,[&UG1,f](double x){return Function::GridFunction<1,double,typeof(UG1),Function::LinearInterpolator>(UG1,[x,f](double y){return f(x,y);});});

    for(size_t i=0;i<10;++i){
        double x = ((double)rand())/RAND_MAX*2-1;
        double y = ((double)rand())/RAND_MAX*2-1;
        print("(",x,y,")",F(x,y),'\t',f(x,y));
    }
    print("\n");
    F.map([](double x,double y){return x+y;});
    for(size_t i=0;i<10;++i){
        double x = ((double)rand())/RAND_MAX*2-1;
        double y = ((double)rand())/RAND_MAX*2-1;
        print("(",x,y,")",F(x,y),'\t',x+y);
    }
    return 0;
}
#endif
#ifdef TEST_M_GRID_FUNCTION_PRINT
int main()
{
    Function::UniformGrid UG1(-1.0,1.0,11);
    auto f = [](double x,double y){return x*y;};
    Function::GridFunction<2,double,typeof(UG1),Function::LinearInterpolator,
                               typeof(UG1),Function::LinearInterpolator>
            F(UG1,[&UG1,f](double x){return Function::GridFunction<1,double,typeof(UG1),Function::LinearInterpolator>(UG1,[x,f](double y){return f(x,y);});});

    std::cout << F.gridStr() <<std::endl;

    std::cout << F.toString() <<std::endl;
}
#endif
#ifdef TEST_M_GRID_FUNCTION_VECTOR
int main()
{
    Function::UniformGrid UG1(-1.0,1.0,11);
    auto f = [](double x,double y){return x*y;};
    Function::GridFunction<double,typeof(UG1),Function::LinearInterpolator,
                               typeof(UG1),Function::LinearInterpolator>
            F(UG1,[&UG1,f](double x){return Function::GridFunction<double,typeof(UG1),Function::LinearInterpolator>(UG1,[x,f](double y){return f(x,y);});});

    auto values = F.AllValues();
    std::cout << std::to_string(values)<<std::endl;
    for(auto & it : values){
        it += 1;
    }
    F.loadIter(values.begin(),values.end());

    std::cout << F.toString() <<std::endl;
}
#endif

#ifdef TEST_HISTO
int main()
{
    Function::UniformGrid UG1(-1.0,1.0,11);
    auto f = [](double x,double y){return x*y;};
    //Function::GridFunction F(UG1,[&UG1,f](double x){return Function::GridFunction<double,typeof(UG1),Function::LinearInterpolator>(UG1,[x,f](double y){return f(x,y);});});
    Function::Histogramm<double,Function::UniformGrid<double>,Function::UniformGrid<double>> H(UG1,UG1);


    size_t M = 10000000;
    size_t N = 24;
    double sm = 0;
    for(size_t i=0;i<M;++i){
        double x = 0;
        double y = 0;
        for(size_t j=0;j<N;++j){
            x += rand()/((double)RAND_MAX)-0.5;
            y += rand()/((double)RAND_MAX)-0.5;
        }
        if(-1<x && x<1 && -1<y &&y<1)
            sm += 1.0/M;
        H.putValue(1.0/M,x,y);
    }

    double sigma_2 = N/12.0;

    std::cout << "-------------------"<<std::endl;
    std::cout << std::to_string(H.AllValues()) <<std::endl;
    //std::cout << H.toFunction().toString() <<std::endl;

    std::cout << "-------------------"<<std::endl;
    auto F =
    Function::GridFunctionCreator2<Function::SchemeHisto>::Create(Function::diffGrid1(H.Grid),[&H,sigma_2](double x){
        return Function::GridFunctionCreator1<Function::SchemeHisto>::Create(Function::diffGrid1(H.Grid),[x,sigma_2](double y){
            return 1.0/(2*M_PI*sigma_2)*exp(-(x*x+y*y)/(2*sigma_2));
        });
    });
    //std::cout << F.toString() <<std::endl;
    std::cout << "-------------------"<<std::endl;
    std::cout << (F-H.toFunction()).toString() <<endl;

    std::cout << H.summ() << std::endl;
    std::cout << sm << std::endl;

    std::cout << 0.04*vector_sum(F.AllValues()) <<std::endl;
}
#endif
#ifdef TEST_HISTO_BOLTSMAN

inline double rand_un(){
    return rand()/(RAND_MAX+1.0);
}
inline double rand_un_incl(){
    return rand()/((double)RAND_MAX);
}
vec3 collision(const vec3 &V,double mu,double sinTheta){
    vec3 Vperp(V.y,-V.x,0);
    double cosTheta = sqrt(1-sinTheta*sinTheta);
    return V*(sinTheta*sinTheta-(1-mu)/(1+mu)*cosTheta*cosTheta) + Vperp*sinTheta*cosTheta*2/(1+mu);
}
vec3 collision(const vec3 &V,const vec3 &V1,double mu,double sinTheta){
    return V1 + collision(V-V1,mu,sinTheta);
}

vec3 collision_mk(const vec3 &V,double V1_disp,double mu){
    double phi = 2*M_PI*rand_un();
    double V1 = V1_disp*sqrt(-2*log(1-rand_un()));
    return collision(V,vec3(V1*cos(phi),V1*sin(phi),0),mu,rand_un_incl()*2-1.0);
}

template <typename T,typename U>
auto dot(const std::vector<T> & X,const std::vector<U> & Y){
    decltype(std::declval<T>()*std::declval<U>()) summ = 0;
    for(size_t i=0;i<X.size();++i){
        summ += X[i]*Y[i];
    }
    return summ;
}

template <typename T,typename U>
auto dot(const std::vector<std::vector<T>> & X,const std::vector<U> & Y){
    std::vector<decltype(std::declval<T>()*std::declval<U>())> summ(X.size(),0);
    for(size_t i=0;i<X.size();++i){
        summ[i] += dot(X[i],Y);
    }
    return summ;
}

template <typename T,typename U>
auto dot(const std::vector<U> & Y,const std::vector<std::vector<T>> & X){
    std::vector<decltype(std::declval<T>()*std::declval<U>())> sum(Y.size(),0);
    for(size_t i=0;i<X.size();++i){
        for(size_t j=0;j<Y.size();++j){
            sum[i] += X[j][i]*Y[j];
        }
    }
    return sum;
}

template <typename  T>
void matrix_dot(std::vector<T> &result,const std::vector<std::vector<T>> &Mat,const std::vector<T> &vec){
    for(size_t i=0;i<result.size();++i){
        result[i] = dot(Mat[i],vec);
    }
}

template <typename T>
std::vector<std::vector<T>> E_Matrix(size_t N){
    return Vector(N,[N](size_t i){
        return Vector(N,[i](size_t j){
            return (T)(i==j);
        });
    });
}

template <typename T>
std::vector<std::vector<T>> MatrixDiagonal(const std::vector<T> &Diag){
    return Vector(Diag.size(),[N = Diag.size(),&Diag](size_t i){
        return Vector(N,[i,Diag](size_t j){
            return (i==j ? Diag[j] : 0);
        });
    });
}

template <typename T>
auto MatrixTranspose(const std::vector<std::vector<T>> &A){
    std::vector<std::vector<T>> Result(A[0].size());
    for(size_t i=0;i<Result.size();++i){
        Result[i].resize(A.size());
        for(size_t j=0;j<A.size();++j){
            Result[i][j] = A[j][i];
        }
    }
    return Result;
}


template <typename T>
void MatrixMult(const std::vector<std::vector<T>> &A,const std::vector<std::vector<T>> &B,std::vector<std::vector<T>> &Result){
    for(size_t i=0;i<A.size();++i){
        for(size_t j=0;j<B[0].size();++j){
            Result[i][j] = 0;
            for(size_t k=0;k<B.size();++k){
                Result[i][j] += A[i][k]*B[k][j];
            }
        }
    }
}

template <typename T>
std::vector<std::vector<T>> MatrixMult(const std::vector<std::vector<T>> &A,const std::vector<std::vector<T>> &B){
    auto BT = TransposeC(B);
    std::vector<std::vector<T>> Ret(A.size());
    for(size_t i=0;i<A.size();++i){
        Ret[i].resize(B[0].size());
        for(size_t j=0;j<B[0].size();++j){
            Ret[i][j] = dot(A[i],BT[j]);
        }
    }
    return Ret;
}


template <typename T>
std::vector<std::vector<T>> MatrixPow(const std::vector<std::vector<T>> &A,size_t N){

    if(N == 0){
        return E_Matrix<T>(A.size());
    }
    std::vector<std::vector<T>> Result = E_Matrix<T>(A.size());
    size_t max2pow = 0;
    size_t pow_of2 = 1;
    while(pow_of2*2<=N){
        max2pow++;
        pow_of2*=2;
    }

    Result = A;
    N %= pow_of2;
    pow_of2 /= 2;
    //max2pow --;

    while(max2pow){
        Result = MatrixMult(Result,Result);
        if(N & pow_of2){
            Result = MatrixMult(Result,A);
        }

        N %=pow_of2;
        max2pow --;
        pow_of2 /= 2;
    }
    return Result;
}


int main()
{

    Function::UniformGrid UG1(-4.0,4.0,11);
    auto f = [](double x,double y){return x*y;};
    //Function::GridFunction F(UG1,[&UG1,f](double x){return Function::GridFunction<double,typeof(UG1),Function::LinearInterpolator>(UG1,[x,f](double y){return f(x,y);});});
    Function::Histogramm<double,Function::UniformGrid<double>,Function::UniformGrid<double>> H(UG1,UG1);

    Function::GridExtractor<decltype (H)>::type_template<decltype (H)> Matrix(UG1,UG1);

    size_t M = 1000;
    double V_disp = 1;
    double mu  = 1;
    for(size_t i=0;i<Matrix.Grid.size()-1;++i){
        double x = 0.5*(Matrix.Grid.at(i+1)+Matrix.Grid.at(i));
        for(size_t j=0;j<Matrix.values[i].Grid.size()-1;++j){
            double y = 0.5*(Matrix.values[i].Grid.at(j+1)+Matrix.values[i].Grid.at(j));
            Matrix.values[i].values[j] =decltype (H)(UG1,UG1);
            for(size_t k = 0;k<M;++k){
                vec3 collision_result = collision_mk(vec3(x,y,0),V_disp,mu);
                Matrix.values[i].values[j].putValue(1.0/M,collision_result.x,collision_result.y);
            }
        }
    }

    std::vector<decltype (H)> pre_vec = Matrix.AllValues();
    std::vector<std::vector<double>> MatT = Vector(pre_vec.size(),[&pre_vec](size_t i){return pre_vec[i].AllValues();});
    auto E = E_Matrix<double> (MatT.size());

    auto Mat = MatrixTranspose(MatT);

    std::vector<double> Ones(MatT.size(),1);
    auto SumA = MatrixDiagonal(Vector(MatT.size(),[&MatT](size_t i){return vector_sum(MatT[i]);}));

    auto S = (Mat-SumA);

    std::cout << dot(Ones,S) << std::endl;
}
#endif
