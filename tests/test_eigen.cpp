#include <iostream>
#include <vector>
#include <Eigen/Dense>

using namespace Eigen;

class Test
{
  public:
    //EIGEN_MAKE_ALIGNED_OPERATOR_NEW
    Eigen::ArrayXf a;
    int size;

    Test(int size, float val)
    {
      this->size = size;
      this->a.setZero(size);
      this->a(0) = val;
      std::cout << "Create obj with first val: " << this->a(0) << std::endl;
    }

    float get(int c) { return this->a[c]; }
    void print()
    {
      std::cout << "The content of obj: ";
      for (int i = 0 ; i < this->size ; i++)
        std::cout << this->a(i) << " " << std::endl;
      std::cout << std::endl;
    }
        
};

Test obj_global(4, 0.2);

void test_obj()
{
  Test obj(5, 0.5);
  std::cout << "Getting only one coef: " << obj.get(0) << std::endl;
  obj.print();

  std::cout << "Getting only one coef: " << obj_global.get(0) << std::endl;
  obj2.print();

}

void test_block()
{
  float arr_f[] = {1.f, 2.f, 3.f, 4.f, 5.f, 6.f, 7.f, 8.f};
  std::vector<float> xx(arr_f, arr_f + 8);
  Eigen::Map<Eigen::ArrayXXf> yy(&xx[0], 4, 2);
  Eigen::ArrayXXf zz = yy.block(1, 0, 2, 2);
  std::cout << "yy = " << yy << std::endl;
  std::cout << "zz = " << zz << std::endl;
  zz *= 2.f;
  std::cout << "yy = " << yy << std::endl;
  std::cout << "zz = " << zz << std::endl;
  std::cout << "yy[0,:] = " << yy.row(0) << std::endl;
  std::cout << "yy[0,:].T = " << yy.row(0).transpose() << std::endl;
  std::cout << std::endl;
}


void test_dot_output()
{
  MatrixXf mat = MatrixXf::Random(3,3);
  ArrayXXf input = ArrayXXf::Random(2,3);
  ArrayXXf output = ArrayXXf::Zero(2,3);

  output.row(0) = mat * input.matrix().row(0).transpose();
  std::cout << "Test do product output:" << std::endl;
  std::cout << "mat : " << mat << std::endl;
  std::cout << "input : " << input << std::endl;
  std::cout << "output : " << output << std::endl;
  VectorXf u = mat * input.matrix().row(0).transpose();
  std::cout << "output : " << u << std::endl;
  std::cout << std::endl;
}

void test_adjoint()
{
  VectorXcf v = VectorXcf::Random(3);

  std::cout << "Test adjoint:" << std::endl;
  std::cout << "v = " << v << std::endl;
  std::cout << " v^H = " << v.adjoint() << std::endl;
  std::cout << "v v^H = " << v * v.adjoint() << std::endl;
  std::cout << "v^H v = " << v.adjoint() * v << std::endl;
  float n2 = (v.adjoint() * v).real()(0,0);
  std::cout << "v^H v (as float) = " << n2 << std::endl;
  std::cout << std::endl;
}


int main()
{
  Matrix2d a;
  a << 1, 2,
    3, 4;
  MatrixXd b(2,2);
  b << 2, 3,
    1, 4;
  std::cout << "a + b =\n" << a + b << std::endl;
  std::cout << "a - b =\n" << a - b << std::endl;
  std::cout << "Doing a += b;" << std::endl;
  a += b;
  std::cout << "Now a =\n" << a << std::endl;
  Vector3d v(1,2,3);
  Vector3d w(1,0,0);
  std::cout << "-v + w - v =\n" << -v + w - v << std::endl;
  std::cout << std::endl;

  ArrayXcf v1 = ArrayXcf::Random(3);
  ArrayXcf v2 = ArrayXcf::Random(3);

  std::cout << "Testing vector operations:" << std::endl;

  std::cout << "v1 = " << v1 << std::endl;
  std::cout << "v2 = " << v2 << std::endl;
  std::cout << "v1 + v2 = " << v1 + v2 << std::endl;
  std::cout << "v1 * v2 = " << v1 * v2 << std::endl;
  std::cout << "v1.abs2() = " << v1.abs2() << std::endl;
  std::cout << std::endl;

  ArrayXXf a1 = ArrayXXf::Random(3, 2);
  //Map<ArrayXf, 1, 2> c2(a1.data, 3);
  std::cout << "Testing slicing operation" << std::endl;
  std::cout << "a1 = " << a1 << std::endl;
  std::cout << "a1[:,1] = " << a1.col(1) << std::endl;
  std::cout << std::endl;

  ArrayXXf all_ones;
  all_ones = ArrayXXf::Zero(4,5);
  all_ones = 1.;
  std::cout << "all_ones = " << all_ones << std::endl;
  std::cout << std::endl;

  ArrayXf x = ArrayXf::Zero(5);
  x.head(2) = 1.;
  std::cout << "x = " << x << std::endl;
  x.segment(2,2) = 2.;
  std::cout << "x = " << x << std::endl;
  std::cout << std::endl;

  test_block();
  test_dot_output();
  test_adjoint();
  test_obj();
}

