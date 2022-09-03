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
