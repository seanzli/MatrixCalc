#include <cmath>
#include <algorithm>
#include <vector>
#include <iostream>
#include <cstdio>
#include <memory>
#include <numeric>
#include <string.h>

class MatrixCalc {
public:
    MatrixCalc() {};
    MatrixCalc(unsigned row, unsigned col) :
    m_row(row), m_col(col) 
    {
        m_size = m_row * m_col;
        m_data.resize(m_size);
    }

    MatrixCalc(double a[], unsigned row, unsigned col) :
    m_row(row), m_col(col) 
    {
        m_size = m_row * m_col;
        m_data.resize(m_size);
        m_data = std::vector<double>(a, a + m_size);
    }

    MatrixCalc(unsigned row, unsigned col, const double a[]) :
    m_row(row), m_col(col) 
    {
        m_size = m_row * m_col;
        m_data.resize(m_size);
        m_data = std::vector<double>(a, a + m_size);
    }

    ~MatrixCalc() {};

    inline unsigned rows() const {return m_row;}
    inline unsigned cols() const {return m_col;}
    inline unsigned size() const {return m_size;}

    double operator()(const unsigned& idx) const
    {
        return idx >= this->size() ? 0.0 : m_data[idx];
    }

    double& operator()(const unsigned& idx)
    {
        return m_data[idx];
    }

    double operator()(const unsigned& row, const unsigned& col) const
    {
        unsigned idx = row * m_col + col;
        return idx >= this->size() ? 0.0 : m_data[idx];
    }

    double& operator()(const unsigned& row, const unsigned& col)
    {
        unsigned idx = row * m_col + col;
        return m_data[idx];
    }

    double operator[](const unsigned& idx) const
    {
        return idx >= this->size() ? 0.0 : m_data[idx];
    }

    double& operator[](const unsigned& idx)
    {
        return m_data[idx];
    }

    bool operator==(const MatrixCalc& op) const
    {
        if (m_col != op.cols() || m_row != op.rows())
            return false;
        for (auto i = 0; i < m_size; i++) {
            if ((*this)(i) != op(i))
                return false;
        }
        return true;
    }

    MatrixCalc T() const
    {
        MatrixCalc out(m_col, m_row);
        for (auto i = 0; i < m_row; i++) {
            for (auto j = 0; j < m_col; j++) {
                out(j, i) = m_data[i * m_col + j];
            }
        }
        return out;
    }

    MatrixCalc operator+(const MatrixCalc& op) const
    {
        MatrixCalc out(m_row, m_col);
        if (op.rows() != m_row || op.cols() != m_col)
            return out;
        for (auto i = 0; i < this->size(); i++) {
            out[i] = this->m_data[i] + op[i];
        }
        return out;
    }

    MatrixCalc operator+(const double& op) const
    {
        MatrixCalc out(m_row, m_col);
        for (auto i = 0; i < this->size(); i++) {
            out[i] = this->m_data[i] + op;
        }
        return out;
    }

    MatrixCalc operator-(const MatrixCalc& op) const
    {
        MatrixCalc out(m_row, m_col);
        if (op.rows() != m_row || op.cols() != m_col)
            return out;
        for (auto i = 0; i < this->size(); i++) {
            out[i] = this->m_data[i] - op[i];
        }
        return out;
    }
    MatrixCalc operator-(const double& op) const
    {
        MatrixCalc out(m_row, m_col);
        for (auto i = 0; i < this->size(); i++) {
            out[i] = this->m_data[i] - op;
        }
        return out;
    }

    MatrixCalc operator*(const MatrixCalc& op) const
    {
        MatrixCalc out(m_row, op.cols());

        if (m_col != op.rows())
            return out;
        for (auto i = 0; i < m_row; i++) {
            for (auto j = 0; j < op.cols(); j++) {
                for (int k = 0; k < m_col; k++) {
                    out(i, j) += m_data[i * m_col + k] * op(k, j);
                }
            }
        }
        return out;
    }

    MatrixCalc operator*(const double& op) const
    {
        MatrixCalc out(m_row, m_col);
        for (auto i = 0; i < this->size(); i++) {
            out[i] = this->m_data[i] * op;
        }
        return out;
    }

    MatrixCalc multiply(const MatrixCalc& op) const
    {
        MatrixCalc out(m_row, m_col);

        if (m_col != op.cols() || m_row != op.rows())
            return out;
        for (auto i = 0; i < m_size; i++)
            out[i] = op[i] * m_data[i];
        return out;
    }

    MatrixCalc inv()
    {
        MatrixCalc out(m_row, m_col);
        double* in = m_data.data();
        double* out_temp = out.m_data.data();
        if (m_row != m_col)
            return out;
        inv(in, out_temp, m_row);
        return out;
    }

    MatrixCalc operator/(MatrixCalc & op) 
    {
        return (*this) * op.inv();
    }

    MatrixCalc operator/(const double& op) const
    {
        if (fabs(op) < 1e-20)
            return MatrixCalc::inf(m_row, m_col);
        return (*this)*(1.0/op);
    }

    static MatrixCalc one(unsigned row, unsigned col)
    {
        MatrixCalc out(row, col);
        out.m_data = std::vector<double>(row*col, 1);
        return out;
    }

    static MatrixCalc zero(unsigned row, unsigned col)
    {
        MatrixCalc out(row, col);
        return out;
    }

    static MatrixCalc eye(unsigned row)
    {
        MatrixCalc out(row, row);
        for (auto i = 0; i < row; i++)
            out(i,i) = 1;
        return out;
    }

    static MatrixCalc eye(unsigned row, unsigned col)
    {
        MatrixCalc out(row, col);
        for (auto i = 0; i < std::min(row, col); i++)
            out(i,i) = 1;
        return out;
    }

    static MatrixCalc inf(unsigned row, unsigned col)
    {
        MatrixCalc out(row, col);
        out.m_data = std::vector<double>(row*col, INT64_MAX);
        return out;
    }

    double sum() const
    {
        return std::accumulate(m_data.begin(), m_data.end(), 0.0);
    }
    
    double mean() const
    {
        return this->sum() / m_size;
    }

    MatrixCalc sum_row() const
    {
        MatrixCalc out(m_row, 1);
        for (auto i = 0; i < m_row; i++) {
            double sum = 0.0;
            for (auto j = 0; j < m_col; j++) {
                sum += (*this)(i,j);
            }
            out[i] = sum;
        }
        return out;
    }

    MatrixCalc sum_col() const
    {
        MatrixCalc out(1, m_col);
        for (auto i = 0; i < m_col; i++) {
            double sum = 0.0;
            for (auto j = 0; j < m_row; j++) {
                sum += (*this)(j,i);
            }
            out[i] = sum;
        }
        return out;
    }

    MatrixCalc mean_row() const
    {
        MatrixCalc out = (*this).sum_row();
        for (auto i = 0; i < m_row; i++)
            out[i] /= m_col;
        return out;
    }

    MatrixCalc mean_col() const
    {
        MatrixCalc out = (*this).sum_col();
        for (auto i = 0; i < m_col; i++)
            out[i] /= m_row;
        return out;
    }

    MatrixCalc getRow(unsigned row) const
    {
        MatrixCalc out(1, m_col);
        for (auto i = 0; i < m_col; i++) {
            out[i] = (*this)(row, i);
        }
        return out;
    }

    MatrixCalc getCol(unsigned col) const
    {
        MatrixCalc out(m_row, 1);
        for (auto i = 0; i < m_row; i++) {
            out[i] = (*this)(i, col);
        }
        return out;
    }

    double std() const
    {
        double _mean = this->mean();
        double sum = 0.0;
        for (auto i = 0; i < m_size; i++) {
            sum += (m_data[i] - _mean);
        }
        return sum / (m_size - 1);
    }

    MatrixCalc std_row() const
    {
        MatrixCalc _mean = this->mean_row();
        int row = _mean.rows();
        MatrixCalc res(row, 1);
        for (auto i = 0; i < m_row; i++) {
            for (auto j = 0; j < m_col; j++) {
                res[i] += ((*this)(i, j) - _mean(i));
            }
            res[i] /= (m_col - 1);
        }
        return res;
    }

    MatrixCalc std_col() const
    {
        MatrixCalc _mean = this->mean_col();
        int col = _mean.cols();
        MatrixCalc res(1, col);
        for (auto j = 0; j < m_col; j++) {
            for (auto i = 0; i < m_row; i++) {
                res[j] += ((*this)(i, j) - _mean(j));
            }
            res[j] /= (m_row - 1);
        }
        return res;
    }

    friend std::ostream &operator<<(std::ostream& out ,const MatrixCalc& op)
    {
        for (auto i = 0; i < op.rows(); i++) {
            for (auto j = 0; j < op.cols(); j++) {
                out << op(i, j) << ',';
            }
            out << std::endl;
        }
        return out;
    }

private:
    unsigned m_row, m_col;
    unsigned m_size;
    std::vector<double> m_data;

    bool inv(double M[], double M_Inv[], int nDim)
    {
        bool bFlag = false;
        int j=0;

        double *LU = new double[m_row*m_col]();
        double *b = new double[m_row]();
        double *Indx = new double[m_row]();
        memcpy(LU, M, sizeof(double)*nDim*nDim);

        // LU decomposition 
        if (LU_Decom( LU, Indx, nDim) == false) {
            delete[] LU;
            delete[] b;
            delete[] Indx;
            return false;
        }

        // Solve Ax=b for  unit vectors b_1..b_n
        for (j=0; j<nDim; j++ ) 
        {
            b[j]= 1.0;                     // Set b to j-th unit vector 
            LU_BackSub( LU, Indx, b, nDim);           // Solve Ax=b
            SetCol(M_Inv, b, nDim, nDim, j+1);
        };

        delete[] LU;
        delete[] b;
        delete[] Indx;
        return true;
    }

    bool LU_Decom(double M[], double Inx[], int nDim)
    {
        // Constants
        const int n = nDim;
        const double tiny = 1.0e-20;       // A small number

        // Variables
        int     i=0,j=0,imax=0,k=0;
        double  aAmax, Sum, Dum;

        // Loop over rows to get scaling information
        for (i=0; i<n; i++)
        {
            aAmax = 0.0;
            for (j=0;j<n;j++)
            {
                if (fabs( M[i * n + j] ) > aAmax )
                    aAmax=fabs( M[i * n + j] );
            }
            if (aAmax==0.0)
            {
                // No nonzero largest element
                //cerr << "ERROR: Singular matrix A in LU_Decomp";
                return false;
            };
            Inx[i] = 1.0/aAmax;           // V stores the implicit scaling of each row
        };

        // Loop over columns of Crout's method
        for ( j=0; j<n; j++ ) 
        {
            if (j > 0) 
            {
                for ( i=0; i<j; i++ )
                {   // This is equation 2.3.12 except for i=j
                    Sum =  M[i * n + j];
                    if (i>0) 
                    {
                        for ( k=0; k<i; k++ )  Sum -=  M[i * n + k] * M[k * n + j] ;
                        M[i * n + j]  = Sum;
                    };
                };
            };

            aAmax=0.0;                  // Initialize for the search of the largest
            // pivot element

            for ( i=j; i<n; i++ ) 
            {     // This is i=j of equation 2.3.12 and 
                Sum =  M[i * n + j] ;             // i=j+1..N of equation 2.3.13
                if (j > 0) 
                {
                    for ( k=0; k<j; k++ ) 
                        Sum -= M[i * n + k] * M[k * n + j];

                    M[i * n + j] = Sum;
                };

                Dum = Inx[i]*fabs(Sum);     // Figure of merit for the pivot

                if (Dum >= aAmax) 
                {       // Is it better than the best so far ?
                    imax  = i;
                    aAmax = Dum;
                };
            };

            if (j != imax) 
            {            
                // Do we need to interchange rows?
                for ( k=0; k<n; k++) 
                {    
                    // Yes, do so ...
                    Dum = M[imax * n + k];
                    M[imax * n + k] =  M[j * n + k];
                    M[j * n + k] = Dum;
                }
                Inx[imax] = Inx[j];           // Also interchange the scale factor 
            };

            Inx[j] = imax;

            if (j != n-1)
            {             
                // Now finally divide by the pivot element
                if (M[j * n + j] == 0.0)
                {      
                    // If the pivot element is zero the matrix 
                    M[j * n + j] = tiny;          // is singular (at least to the precision of
                };                        // the algorithm). For some applications on

                Dum=1.0/M[j * n + j];           // singular matrices, it is desirable to 
                for (i=j+1;i<n;i++)
                {     
                    // substitute tiny for zero. 
                    M[i * n + j]=M[i * n + j]*Dum;
                };
            };

        };   // Go back for the next column in the reduction

        if (M[(n-1) * n + (n-1)]==0.0) M[(n-1) * n + (n-1)]=tiny; 

        return true;
    }

    void LU_BackSub (double M[], double Inx[], double b[], int nDim)
    {
        // Constants
        const int  n = nDim;

        // Local variables
        int     ii=0,i=0,ll=0,j=0;
        double  Sum=0;

        // Start
        // When ii is set to a nonegative value
        // it will become the first nonvanishing element of B. 
        ii = -1;  

        for (i=0; i<n; i++)
        {         
            // We now do the forward substitution.
            ll = (int) Inx[i];         // The only wrinkle is to unscramble the 
            Sum = b[ll];                // permutation as we go.
            b[ll] = b[i];
            if (ii != -1) 
            {
                for (j=ii; j<i; j++) Sum -= M[i * n + j]*b[j];
            }
            else 
            {
                if(fabs(Sum) > 1e-20) 
                    ii = i;   // A nonzero element was encountered, so from 
            };                          // now on we will have to do the sums in the
            b[i] = Sum;                 // loop above.
        };

        for (i=n-1; i>=0; i--) 
        {     
            // Now we do the back substitution, eqn 2.3.7.
            Sum=b[i];
            if (i<n-1) 
            {
                for (j=i+1;j<n;j++)
                {
                    Sum = Sum-M[i * n + j]*b[j];
                };
            };
            b[i] = Sum/M[i * n + i];         // Store a component of the solution vector X.
        }
    }

    void SetCol(double M[], double Cv[], int nRow, int nCol, int k)
    {
        int i=0;
        if ( k>nCol || k<1)
        {
            return;
        }
        for (i=0; i<nRow; i++)
        {
            M[i*nCol+k-1] = Cv[i];
        }
    }

};