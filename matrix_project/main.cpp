#include "Matrix.hpp"

int main()
{
    try
    {
        Matrix<double> matrix5(5, 5);
        for (unsigned i = 0; i < 5; ++i)
        {
            for (unsigned j = 0; j <= i; ++j)
            {
                matrix5(i, j) = static_cast<double>(rand() % 10 + 1);
            }
        }
        
        cout << "5x5 Lower Triangular Matrix:" << endl;
        matrix5.print();
        cout << "Determinant: " << matrix5.determinant() << endl;
    }
    catch (const exception &e)
    {
        cerr << "Error: " << e.what() << endl;
    }
    
    return 0;
}
