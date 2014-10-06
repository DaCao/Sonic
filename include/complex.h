//   complex.h - declaration of class
//   of complex2 number
//
//   The code is property of LIBROW
//   You can use it on your own
//   When utilizing credit LIBROW site



// to be recovered below:

//#ifndef _COMPLEX_H_
//#define _COMPLEX_H_

//class complex2
//{
//protected:
//	//   Internal presentation - real and imaginary parts
//	double m_re;
//	double m_im;
//
//public:
//	//   Imaginary unity
//	static const complex2 i;
//	static const complex2 j;
//
//	//   Constructors
//	complex2(): m_re(0.), m_im(0.) {}
//	complex2(double re, double im): m_re(re), m_im(im) {}
//	complex2(double val): m_re(val), m_im(0.) {}
//
//	//   Assignment
//	complex2& operator= (const double val)
//	{
//		m_re = val;
//		m_im = 0.;
//		return *this;
//	}
//
//	//   Basic operations - taking parts
//	double re() const { return m_re; }
//	double im() const { return m_im; }
//
//	//   Conjugate number
//	complex2 conjugate() const
//	{
//		return complex2(m_re, -m_im);
//	}
//
//	//   Norm   
//	double norm() const
//	{
//		return m_re * m_re + m_im * m_im;
//	}
//
//	//   Arithmetic operations
//	complex2 operator+ (const complex2& other) const
//	{
//		return complex2(m_re + other.m_re, m_im + other.m_im);
//	}
//
//	complex2 operator- (const complex2& other) const
//	{
//		return complex2(m_re - other.m_re, m_im - other.m_im);
//	}
//
//	complex2 operator* (const complex2& other) const
//	{
//		return complex2(m_re * other.m_re - m_im * other.m_im,
//			m_re * other.m_im + m_im * other.m_re);
//	}
//
//	complex2 operator/ (const complex2& other) const
//	{
//		const double denominator = other.m_re * other.m_re + other.m_im * other.m_im;
//		return complex2((m_re * other.m_re + m_im * other.m_im) / denominator,
//			(m_im * other.m_re - m_re * other.m_im) / denominator);
//	}
//
//	complex2& operator+= (const complex2& other)
//	{
//		m_re += other.m_re;
//		m_im += other.m_im;
//		return *this;
//	}
//
//	complex2& operator-= (const complex2& other)
//	{
//		m_re -= other.m_re;
//		m_im -= other.m_im;
//		return *this;
//	}
//
//	complex2& operator*= (const complex2& other)
//	{
//		const double temp = m_re;
//		m_re = m_re * other.m_re - m_im * other.m_im;
//		m_im = m_im * other.m_re + temp * other.m_im;
//		return *this;
//	}
//
//	complex2& operator/= (const complex2& other)
//	{
//		const double denominator = other.m_re * other.m_re + other.m_im * other.m_im;
//		const double temp = m_re;
//		m_re = (m_re * other.m_re + m_im * other.m_im) / denominator;
//		m_im = (m_im * other.m_re - temp * other.m_im) / denominator;
//		return *this;
//	}
//
//	complex2& operator++ ()
//	{
//		++m_re;
//		return *this;
//	}
//
//	complex2 operator++ (int)
//	{
//		complex2 temp(*this);
//		++m_re;
//		return temp;
//	}
//
//	complex2& operator-- ()
//	{
//		--m_re;
//		return *this;
//	}
//
//	complex2 operator-- (int)
//	{
//		complex2 temp(*this);
//		--m_re;
//		return temp;
//	}
//
//	complex2 operator+ (const double val) const
//	{
//		return complex2(m_re + val, m_im);
//	}
//
//	complex2 operator- (const double val) const
//	{
//		return complex2(m_re - val, m_im);
//	}
//
//	complex2 operator* (const double val) const
//	{
//		return complex2(m_re * val, m_im * val);
//	}
//
//	complex2 operator/ (const double val) const
//	{
//		return complex2(m_re / val, m_im / val);
//	}
//
//	complex2& operator+= (const double val)
//	{
//		m_re += val;
//		return *this;
//	}
//
//	complex2& operator-= (const double val)
//	{
//		m_re -= val;
//		return *this;
//	}
//
//	complex2& operator*= (const double val)
//	{
//		m_re *= val;
//		m_im *= val;
//		return *this;
//	}
//
//	complex2& operator/= (const double val)
//	{
//		m_re /= val;
//		m_im /= val;
//		return *this;
//	}
//
//	friend complex2 operator+ (const double left, const complex2& right)
//	{
//		return complex2(left + right.m_re, right.m_im);
//	}
//
//	friend complex2 operator- (const double left, const complex2& right)
//	{
//		return complex2(left - right.m_re, -right.m_im);
//	}
//
//	friend complex2 operator* (const double left, const complex2& right)
//	{
//		return complex2(left * right.m_re, left * right.m_im);
//	}
//
//	friend complex2 operator/ (const double left, const complex2& right)
//	{
//		const double denominator = right.m_re * right.m_re + right.m_im * right.m_im;
//		return complex2(left * right.m_re / denominator,
//			-left * right.m_im / denominator);
//	}
//
//	//   Boolean operators
//	bool operator== (const complex2 &other) const
//	{
//		return m_re == other.m_re && m_im == other.m_im;
//	}
//
//	bool operator!= (const complex2 &other) const
//	{
//		return m_re != other.m_re || m_im != other.m_im;
//	}
//
//	bool operator== (const double val) const
//	{
//		return m_re == val && m_im == 0.;
//	}
//
//	bool operator!= (const double val) const
//	{
//		return m_re != val || m_im != 0.;
//	}
//
//	friend bool operator== (const double left, const complex2& right)
//	{
//		return left == right.m_re && right.m_im == 0.;
//	}
//
//	friend bool operator!= (const double left, const complex2& right)
//	{
//		return left != right.m_re || right.m_im != 0.;
//	}
//};
//#endif



// to be recovered above.

