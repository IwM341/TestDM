#include <iostream>
#include <utils>
using namespace std;

int main()
{
    Function::VectorGrid<double> G1 = Function::UniformGrid(-1.0,1.0,11);
    for(size_t i=0;i<10;++i){
        double x = ((double)rand())/RAND_MAX*2-1;
        print(G1.pos(x),G1.at(G1.pos(x)),x,G1.at(G1.pos(x)+1));
    }
    return 0;
}
