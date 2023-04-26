//Kirill Arkhipov
//k.arkhipov@innopolis.university
#include<iostream>
#include<vector>
#include<iomanip>
#include<math.h>
#include<cstdio>
#include <random>

using namespace std;


class colVector {
public:
    vector<double> vec;

    colVector() = default;

    // Method, which swaps two rows in the vector.
    colVector swapRow(int i, int j) {
        colVector result(vec.size());
        for (int i = 0; i < vec.size(); i++) {
            result[i] = vec[i];
        }
        double temp = vec[i];
        result[i] = vec[j];
        result[j] = temp;
        return result;
    }

    // The constructor, which creates a vector of size n and fills it with zeros.
    colVector(int n) {
        vec.resize(n);
        fill(vec.begin(), vec.end(), 0);
    }

    // Overload method, that allows to read the vector from the console.
    friend istream &operator>>(istream &in, colVector &v) {
        for (int i = 0; i < v.vec.size(); i++) {
            in >> v.vec[i];
        }
        return in;
    }

    // Overload method, that allows to print the vector to the console.
    friend ostream &operator<<(ostream &out, colVector &v) {

        for (int i = 0; i < v.vec.size(); i++) {


            out << fixed << setprecision(4) << v.vec[i] << " ";
        }
        return out;
    }

    // Overload method, that allows to access the vector elements by index.
    double &operator[](int i) {
        return vec[i];
    }

    // Overload method, which returns the size of the vector.
    int size() {
        return vec.size();
    }

    // Overload method, which multiplies the vector by a number.
    colVector operator*(double a) {
        colVector res(vec.size());
        for (int i = 0; i < vec.size(); i++) {
            res[i] = vec[i] * a;
        }
        return res;
    }

    // Overload method, which multiplies the vector by another vector.
    colVector operator*(colVector &v) {
        colVector res(vec.size());
        for (int i = 0; i < vec.size(); i++) {
            res[i] = vec[i] * v[i];
        }
        return res;
    }

    // Overload method, which adds the vector to another vector.
    colVector operator+(colVector &v) {
        colVector res(vec.size());
        for (int i = 0; i < vec.size(); i++) {
            res[i] = vec[i] + v[i];
        }
        return res;
    }

    // Overload method, which subtracts the vector from another vector.
    colVector operator-(colVector &v) {
        colVector res(vec.size());
        for (int i = 0; i < vec.size(); i++) {
            res[i] = vec[i] - v[i];
        }
        return res;
    }

    // Overload method, which divides the vector by another vector.
    colVector operator/(colVector &v) {
        colVector res(vec.size());
        for (int i = 0; i < vec.size(); i++) {
            res[i] = vec[i] / v[i];
        }
        return res;
    }

    // Overload method, which divides the vector by a number.
    colVector operator/(double a) {
        colVector res(vec.size());
        for (int i = 0; i < vec.size(); i++) {
            res[i] = vec[i] / a;
        }
        return res;
    }

    // Method, which returns the norm of the vector.
    colVector normalize() {
        colVector res(vec.size());
        double sum = 0;
        for (int i = 0; i < vec.size(); i++) {
            sum += vec[i] * vec[i];
        }
        sum = sqrt(sum);
        for (int i = 0; i < vec.size(); i++) {
            res[i] = vec[i] / sum;
        }
        return res;
    }

    // Method, which check, that the vector contains only 0s.
    bool isZero() {
        for (int i = 0; i < vec.size(); i++) {
            if (vec[i] != 0) {
                return false;
            }
        }
        return true;
    }
};

class Matrix {
public:
    int columns;
    int rows;
    vector<colVector> myMatrix;

    Matrix() {

    }

    int norm(Matrix e1, Matrix e2) {
        Matrix res = e1 - e2;
        Matrix res_abs = res.absolute_value();
        double max = res_abs.myMatrix[0][0];
        for (int i = 0; i < res_abs.rows; i++) {
            res_abs.myMatrix[i][0] = res_abs.myMatrix[i][0] / e1.myMatrix[i][0];
        }
        for (int i = 0; i < res_abs.rows; i++) {
            if (res_abs.myMatrix[i][0] > max) {
                max = res_abs.myMatrix[i][0];
            }

        }
        cout << res_abs << "\n";
        return max;
    }

    Matrix operator=(colVector m) {
        columns = 1;
        rows = m.size();
        for (int i = 0; i < rows; i++) {
            myMatrix[i][0] = m[i];
        }
        return *this;
    }

    Matrix operator=(Matrix m) {
        columns = m.columns;
        rows = m.rows;
        myMatrix = m.myMatrix;
        return *this;
    }

    Matrix(int n, int m) {
        columns = m;
        rows = n;
        for (int i = 0; i < n; i++) {
            myMatrix.push_back(colVector(m));
        }
    }

    Matrix operator*(colVector &v) {
        Matrix res(rows, 1);
        for (int i = 0; i < rows; i++) {
            for (int j = 0; j < columns; j++) {
                res.myMatrix[i][0] += myMatrix[i][j] * v[j];
            }
        }
        return res;
    }

    // Overload method, that allows to read the matrix from the console.
    friend istream &operator>>(istream &in, Matrix &m) {
        for (int i = 0; i < m.myMatrix.size(); i++) {
            in >> m.myMatrix[i];
        }
        return in;
    }

    // Overload method, that allows to print the matrix to the console.
    friend ostream &operator<<(ostream &out, Matrix &m) {
        for (int i = 0; i < m.myMatrix.size(); i++) {
            out << m.myMatrix[i] << "\n";
        }
        return out;
    }

    Matrix operator+(Matrix m) {
        Matrix res(rows, columns);

        for (int i = 0; i < rows; i++) {
            for (int j = 0; j < columns; j++) {
                res.myMatrix[i][j] = myMatrix[i][j] + m.myMatrix[i][j];
            }
        }
        return res;
    }

    Matrix operator-(Matrix &m) {
        Matrix res(rows, columns);

        for (int i = 0; i < rows; i++) {
            for (int j = 0; j < columns; j++) {
                res.myMatrix[i][j] = myMatrix[i][j] - m.myMatrix[i][j];
            }
        }
        return res;
    }

    Matrix operator*(Matrix &m) {
        Matrix res(rows, m.columns);

        for (int i = 0; i < rows; i++) {
            for (int j = 0; j < m.columns; j++) {
                for (int k = 0; k < columns; k++) {
                    res.myMatrix[i][j] += myMatrix[i][k] * m.myMatrix[k][j];
                }
            }
        }
        return res;
    }

    Matrix transpose() {
        Matrix res(columns, rows);
        for (int i = 0; i < rows; i++) {
            for (int j = 0; j < columns; j++) {
                res.myMatrix[j][i] = myMatrix[i][j];
            }
        }
        return res;
    }

    Matrix absolute_value() {
        Matrix res(rows, columns);
        for (int i = 0; i < rows; i++) {
            for (int j = 0; j < columns; j++) {
                res.myMatrix[i][j] = abs(myMatrix[i][j]);
            }
        }
        return res;
    }


};

class SquareMatrix : public Matrix {
public:
    bool check_diagonally_dominant(Matrix &m) {
        for (int i = 0; i < rows; i++) {
            double sum = 0;
            for (int j = 0; j < columns; j++) {
                if (i != j) {
                    sum += abs(m.myMatrix[i][j]);
                }
            }
            if (abs(m.myMatrix[i][i]) < sum) {
                return false;
            }
        }
        return true;
    }

    SquareMatrix make_dioganl_matrix() {
        SquareMatrix res(rows);
        for (int i = 0; i < rows; i++) {
            for (int j = 0; j < columns; j++) {
                if (i != j) {
                    res.myMatrix[i][j] = 0;
                } else {
                    res.myMatrix[i][j] = myMatrix[i][j];
                }
            }
        }
        return res;
    }

    Matrix operator*(colVector &v) {
        Matrix res(rows, 1);
        for (int i = 0; i < rows; i++) {
            for (int j = 0; j < columns; j++) {
                res.myMatrix[i][0] += myMatrix[i][j] * v[j];
            }
        }
        return res;
    }


    SquareMatrix() {
        this->rows = 0;
        this->columns = 0;
        this->myMatrix = {};
    };


    SquareMatrix(int n) {
        this->rows = n;
        this->columns = n;
        for (int i = 0; i < n; i++) {
            myMatrix.push_back(colVector(n));

        }
    }

    friend istream &operator>>(istream &in, SquareMatrix &m) {
        for (int i = 0; i < m.columns; i++) {

            in >> m.myMatrix[i];
        }
        return in;
    }

    // Overload method, that allows to print the matrix to the console.
    friend ostream &operator<<(ostream &out, SquareMatrix &m) {
        for (int i = 0; i < m.columns; i++) {
            out << m.myMatrix[i] << "\n";
        }
        return out;
    }

    SquareMatrix operator+(SquareMatrix &m) {
        Matrix *right;
        Matrix *left;
        right = &m;
        left = this;
        Matrix res(rows, columns);
        res = *left + *right;
        SquareMatrix *res1;
        res1 = (SquareMatrix *) &res;
        return *res1;
    }

    SquareMatrix operator-(SquareMatrix &m) {
        Matrix *right;
        Matrix *left;
        right = &m;
        left = this;
        Matrix res(rows, columns);
        res = *left - *right;
        SquareMatrix *res1;
        res1 = (SquareMatrix *) &res;
        return *res1;
    }

    SquareMatrix operator*(SquareMatrix &m) {
        Matrix *right;
        Matrix *left;
        right = &m;
        left = this;
        Matrix res(rows, columns);
        res = *left * *right;
        SquareMatrix *res1;
        res1 = (SquareMatrix *) &res;
        return *res1;
    }

    SquareMatrix transpose() {
        Matrix *left;
        left = this;
        Matrix res(rows, columns);
        res = left->transpose();
        SquareMatrix *res1;
        res1 = (SquareMatrix *) &res;
        return *res1;
    }


};

class IdentityMatrix : public SquareMatrix {

public:
    IdentityMatrix() {
        this->rows = 0;
        this->columns = 0;
        this->myMatrix = {};
    };

    IdentityMatrix(int n) {
        this->rows = n;
        this->columns = n;
        for (int i = 0; i < n; i++) {
            myMatrix.push_back(colVector(n));

        }

        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                myMatrix[i][j] = 0;
            }
            myMatrix[i][i] = 1;
        }
    }

    IdentityMatrix operator=(SquareMatrix m) {
        for (int i = 0; i < m.rows; i++) {
            for (int j = 0; j < m.columns; j++) {
                this->myMatrix[i][j] = m.myMatrix[i][j];
            }
        }
        return *this;
    }

    IdentityMatrix resize(int n) {
        IdentityMatrix res(n);
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                res.myMatrix[i][j] = myMatrix[i][j];
            }
        }
        return res;
    }
};


class EluminationMatrix : public SquareMatrix {
public:
    int i, j;

    EluminationMatrix(int i, int j, Matrix &m) {
        i = i - 1;
        j = j - 1;
        IdentityMatrix im(m.rows);
        double c = m.myMatrix[i][j] / m.myMatrix[j][j];
        for (int l = 0; l < m.rows; l++) {
            im.myMatrix[i][l] += -c * im.myMatrix[j][l];
        }
        this->myMatrix = im.myMatrix;
        this->rows = im.rows;
        this->columns = im.columns;
    }

};

class PermutationMatrix : public SquareMatrix {
public:
    PermutationMatrix(int i, int j, Matrix &m) {
        i = i - 1;
        j = j - 1;
        IdentityMatrix im(m.rows);
        for (int l = 0; l < m.rows; l++) {
            swap(im.myMatrix[i][l], im.myMatrix[j][l]);
        }
        this->myMatrix = im.myMatrix;
        this->rows = im.rows;
        this->columns = im.columns;
    }
};

SquareMatrix REF(SquareMatrix matrix) {
    int k = 1;
    IdentityMatrix I(matrix.columns);
    for (int i = 0; i < matrix.columns - 1; i++) {
        double m = matrix.myMatrix[i][i];
        int index;
        for (int j = i; j < matrix.rows; j++) {
            if (m < fabs(matrix.myMatrix[j][i])) {
                m = max(fabs(matrix.myMatrix[j][i]), m);
                index = j;
            }
        }
        PermutationMatrix permutationMatrix(i + 1, index + 1, matrix);
        PermutationMatrix permutationMatrixIdentity(i + 1, index + 1, I);

        matrix = permutationMatrix * matrix;
        cout << "step #" << k << ": permutation\n";
        k += 1;
        cout << matrix;

        for (int j = i + 1; j < matrix.rows; j++) {
            EluminationMatrix eluminationMatrix(j + 1, i + 1, matrix);
            matrix = eluminationMatrix * matrix;
            cout << "step #" << k << ": elimination\n";
            k += 1;
            cout << matrix;
        }


    }
    return matrix;
}

double det(SquareMatrix matrix) {
    matrix = REF(matrix);
    double determinant = 1;
    for (int i = 0; i < matrix.columns; i++) {
        determinant = determinant * matrix.myMatrix[i][i];
    }
    cout << "result:\n";
    return determinant;
}

class AugmentedMatrix : public Matrix {
public:
    SquareMatrix mainMatrix;
    IdentityMatrix I;

    AugmentedMatrix(SquareMatrix _matrix) {
        this->mainMatrix = _matrix;
        this->rows = _matrix.rows;
        this->columns = _matrix.columns;

        this->I = IdentityMatrix(_matrix.columns);


    }

    friend ostream &operator<<(ostream &out, AugmentedMatrix &m) {
        for (int i = 0; i < m.mainMatrix.rows; i++) {
            out << fixed << setprecision(4) << m.mainMatrix.myMatrix[i] << m.I.myMatrix[i] << "\n";
        }
        return out;
    }


};

AugmentedMatrix RREF(AugmentedMatrix matrix) {
    int k = 1;
    for (int i = 0; i < matrix.mainMatrix.columns - 1; i++) {
        bool flag = false;
        double m = fabs(matrix.mainMatrix.myMatrix[i][i]);
        int index = i + 1;
        for (int j = i + 1; j < matrix.mainMatrix.rows; j++) {
            if (m < fabs(matrix.mainMatrix.myMatrix[j][i])) {
                m = max(fabs(matrix.mainMatrix.myMatrix[j][i]), m);
                index = j;
                flag = true;
            }
        }
        if (flag) {
            PermutationMatrix permutationMatrix(i + 1, index + 1, matrix.mainMatrix);
            matrix.mainMatrix = permutationMatrix * matrix.mainMatrix;
            PermutationMatrix permutationMatrixIden(i + 1, index + 1, matrix.I);
            matrix.I = permutationMatrixIden * matrix.I;
            k += 1;
        }


        for (int j = i + 1; j < matrix.mainMatrix.rows; j++) {
            EluminationMatrix eluminationMatrix(j + 1, i + 1, matrix.mainMatrix);
            bool equal = true;
            IdentityMatrix flag(eluminationMatrix.rows);

            for (int h = 0; h < matrix.mainMatrix.rows; h++) {
                for (int l = 0; l < matrix.mainMatrix.columns; l++) {
                    if (eluminationMatrix.myMatrix[h][l] != flag.myMatrix[h][l]) {
                        equal = false;
                        break;
                    }
                }
            }
            if (!equal) {
                matrix.mainMatrix = eluminationMatrix * matrix.mainMatrix;
                matrix.I = eluminationMatrix * matrix.I;
                k += 1;
            }

        }

    }
// back
    for (int i = matrix.mainMatrix.columns - 1; i > 0; i--) {
        for (int j = i - 1; j >= 0; j--) {
            EluminationMatrix eluminationMatrix(j + 1, i + 1, matrix.mainMatrix);
            bool equal = true;
            IdentityMatrix flag(eluminationMatrix.rows);

            for (int h = 0; h < matrix.mainMatrix.rows; h++) {
                for (int l = 0; l < matrix.mainMatrix.columns; l++) {
                    if (eluminationMatrix.myMatrix[h][l] != flag.myMatrix[h][l]) {
                        equal = false;
                        break;
                    }
                }
            }
            if (!equal) {
                matrix.mainMatrix = eluminationMatrix * matrix.mainMatrix;
                matrix.I = eluminationMatrix * matrix.I;
                k += 1;
            }

        }

    }


    return matrix;
}


AugmentedMatrix normalization(AugmentedMatrix &matrix) {

    for (int i = 0; i < matrix.mainMatrix.rows; i++) {
        double m = matrix.mainMatrix.myMatrix[i][i];
        for (int j = 0; j < matrix.mainMatrix.columns; j++) {
            if (i == j) {
                for (int k = 0; k < matrix.rows; k++) {
                    matrix.mainMatrix.myMatrix[i][k] = matrix.mainMatrix.myMatrix[i][k] / m;
                    matrix.I.myMatrix[i][k] = matrix.I.myMatrix[i][k] / m;
                    if (matrix.I.myMatrix[i][k] >= 0)
                        matrix.I.myMatrix[i][k] = max(matrix.I.myMatrix[i][k], 0.00);
                    if (matrix.I.myMatrix[i][k] <= 0 and fabs(matrix.I.myMatrix[i][k]) <= 0.001)
                        matrix.I.myMatrix[i][k] = 0.00;
                }


            }
        }

    }
    return matrix;
}


#ifdef WIN32
#define GNUPLOT_NAME "C:\\gnuplot\\bin\\gnuplot -persist"
#else
#define GNUPLOT_NAME "gnuplot -persist"
#endif

int main() {
#ifdef WIN32
    FILE *pipe = _popen(GNUPLOT_NAME, "w");
#else
    FILE* pipe = popen(GNUPLOT_NAME, "w");
#endif
    if (pipe != NULL) {

        int m, temp;
        cin >> m;
        Matrix A_prew = Matrix(m, 2);
        Matrix A_T_prew = Matrix(2, m);


        fprintf(pipe, "%s\n", "plot '-' title 'Data' with points, '-' title 'Least "
                              "Squares Approximation' with lines");




        int index = 0;
        colVector t, b;
        Matrix B(m, 1);
        for (int i = 0; i < 2 * m; i++) {
            if (i % 2 == 0) {
                cin >> temp;
                t.vec.push_back(temp);
            } else {
                cin >> temp;
                b.vec.push_back(temp);
                B.myMatrix[index][0] = temp;
                index += 1;

            }
        }
        for (int j = 0; j < m; j++) {
            A_T_prew.myMatrix[0][j] = 1;
            A_T_prew.myMatrix[1][j] = t.vec[j];

        }
        int degree;
        cin >> degree;


        Matrix A = Matrix(m, degree + 1);
        Matrix A_T = Matrix(degree + 1, m);
        for (int i = 0; i < m; i++) {
            A.myMatrix[i][0] = 1;
            A.myMatrix[i][1] = t.vec[i];
        }
        for (int i = 0; i < m; i++) {
            for (int j = 0; j < degree + 1; j++) {
                if (j > 1)
                    A.myMatrix[i][j] = A.myMatrix[i][j - 1] * t.vec[i];
            }
        }


        cout << "A:\n" << A;
        A_T = A.transpose();
        Matrix A_T_A = A_T * A;
        SquareMatrix A_T_A_Square(A_T_A.rows);
        for (int i = 0; i < A_T_A.rows; i++)
            for (int j = 0; j < A_T_A.columns; j++)
                A_T_A_Square.myMatrix[i][j] = A_T_A.myMatrix[i][j];


        cout << "A_T*A:\n" << A_T_A;

        AugmentedMatrix augmentedMatrix(A_T_A_Square);
        AugmentedMatrix rref = RREF(augmentedMatrix);
        augmentedMatrix = normalization(rref);
        Matrix inverse = augmentedMatrix.I;

        cout << "(A_T*A)^-1:\n" << inverse;

        Matrix A_T_B = A_T * B;
        cout << "A_T*b:\n" << A_T_B;

        Matrix Answer = inverse * A_T_B;
        cout << "x~:\n" << Answer;

        for (int i = 0; i < m; i++) {
            int x, y;
            x = t.vec[i];
            y = b.vec[i];
            fprintf(pipe, "%d\t%d\n", x, y);
        }
        fprintf(pipe, "%s\n", "e");
        for (double x = -10; x < 10; x+=0.1) {

            double y = 0;
            for (int i = 0; i < Answer.rows; i++)
                y += Answer.myMatrix[i][0] * pow(x, i);

            fprintf(pipe, "%f\t%f\n", x, y);
        }

        fprintf(pipe, "%s\n", "e");
        fflush(pipe);
        pclose(pipe);
    }
}