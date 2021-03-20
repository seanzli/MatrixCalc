#include <cmath>
#include <algorithm>
#include <memory>

class MatrixCalc {
public:
    MatrixCalc() {m_data == nullptr;};
    MatrixCalc(unsigned row, unsigned col) :
    m_row(row), m_col(col) 
    {
         m_data = new double[m_row * m_col]();
         m_size = m_row * m_col;
    }

    MatrixCalc(unsigned row, unsigned col, const double a[]) :
    m_row(row), m_col(col) 
    {
        m_data = new double[m_row * m_col]();
        memcpy(m_data, a, sizeof(double)*row*col);
        m_size = m_row * m_col;
    }

    ~MatrixCalc() {
        if (m_data != nullptr)
            delete m_data;
    };

    inline unsigned rows() const {return m_row;}
    inline unsigned cols() const {return m_col;}
    inline unsigned size() const {return m_size;}

    double operator()(const unsigned& idx) const
    {
        return idx >= this->size() ? 0.0 : m_data[idx];
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

    MatrixCalc T() const
    {
        MatrixCalc out(m_col, m_row);
        for (int i = 0; i < m_row; i++) {
            for (int j = 0; j < m_col; j++) {
                out(j, i) = m_data[i * m_col + j];
            }
        }
        return out;
    }

    MatrixCalc operator+(const MatrixCalc& op) const
    {
        MatrixCalc out(m_row, m_col);
        if (op.rows() != m_row || op.cols() != m_col)
            return;
        for (int i = 0; i < this->size(); i++) {
            out[i] = this->m_data[i] + op[i];
        }
        return out;
    }

    MatrixCalc operator-(const MatrixCalc& op) const
    {
        MatrixCalc out(m_row, m_col);
        if (op.rows() != m_row || op.cols() != m_col)
            return;
        for (int i = 0; i < this->size(); i++) {
            out[i] = this->m_data[i] - op[i];
        }
        return out;
    }

    MatrixCalc operator*(const MatrixCalc& op) const
    {
        MatrixCalc out(m_row, op.cols());

        if (m_col != op.rows())
            return out;
        for (int i = 0; i < m_row; i++) {
            for (int j = 0; j < op.cols(); j++) {
                for (int k = 0; k < m_col; k++) {
                    out(i, j) += m_data[i * m_col + k] * op(k, j);
                }
            }
        }
        return out;
    }

    MatrixCalc multiply(const MatrixCalc& op) const
    {
        MatrixCalc out(m_row, m_col);

        if (m_col != op.cols() || m_row != op.rows())
            return out;
        for (int i = 0; i < m_size; i++)
            out[i] = op[i] * m_data[i];
        return out;
    }

private:
    unsigned m_row, m_col;
    unsigned m_size;
    double* m_data;

};