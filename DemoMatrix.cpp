#include <iostream>
#include "SparseMatrix.h"
#include "RegMatrix.h"
#include <cassert>

int main()
{
    double fill1[] = {1., 0.5, 0.005, 3.};
    double fill2[] = {2., -2.5, 3, 0};
    double fill3[] = {0., -1, 1, 0.5};

    RegMatrix r1(2,2,fill1), r2(2,2,fill2), r3(2,2,fill3);
    SparseMatrix s1(2,2,fill1), s2(2,2,fill2), s3(2,2,fill3);
    assert( r1==s1 );
    assert( r2==s2 );
    assert( r3==s3 );
    
    RegMatrix r4(r1), r5(s1);
    SparseMatrix s4(r1), s5(s1);
    assert( r4==s4 );
    assert( r5==s5 );


    // Only check that pass compilation
    r4 = -r5;
    s4 = -s5;
    r1 += r4;
    s1 += s4;
    r2 *= r5;
    s2 *= s5;
    s3 *= s3;
    r3 *= r3;
    s3 + s4;
    s3 + r4;
    r4 + s3;
    s4 + r3;
    s3 * s4;
    s3 * r4;
    r4 * s3;
    s4 * r3;
    
    bool b;
    double d;
    b = r1.det(d);
    b = s1.det(d);
    r2 = r1.transpose();
    s2 = s1.transpose();
    d= r1.sum();
    d= s1.sum();
    
    return 0;
}
