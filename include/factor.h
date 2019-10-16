/*
 * factor.h
 *
 *  Created on: Feb 8, 2013
 *      Author: radu
 *
 * Copyright (c) 2015, International Business Machines Corporation
 * and University of California Irvine. All rights reserved.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 * DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
 * FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
 * DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
 * SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
 * CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
 * OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
 * OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */


/// \file factor.h
/// \brief A table based factor for graphical models
/// \author Radu Marinescu radu.marinescu@ie.ibm.com

#ifndef IBM_MERLIN_FACTOR_H_
#define IBM_MERLIN_FACTOR_H_


#include <float.h>

#include "enum.h"
#include "util.h"
#include "variable_set.h"
#include "index.h"

namespace merlin {

///
/// \brief Factor for graphical models.
///
/// Table based representation of a factor for graphical models. A 
/// factor encodes a potential (sometimes a probability distribution)
/// defined over a subset of discrete random variables, called a *scope*, and 
/// associates each configuration of the variables in the scope with a 
/// positive real value (sometimes a probability value). The scope is assumed
/// to be sorted lexicogaphically (e.g., [x1,x2,x3]) Also, the indexing of
/// configurations in the factor table is assumed to be based on the BigEndian
/// convention, namely the *first* variable in the ordered scope changes
/// the fastest, then the *second* variable changes its value and so on.
/// For example, consider a factor over binary variables [x1,x2,x3].
/// The corresponding factor table is indexed as follows (internally):
///
/// 0: [0,0,0]    4: [0,0,1]
/// 1: [1,0,0]    5: [1,0,1]
/// 2: [0,1,0]    6: [0,1,1]
/// 3: [1,1,0]    7: [1,1,1]
///
/// (For converting between source and target orders, use the convert_index class)
///
class factor {
public:
	// Typedefs:

	typedef double value;					///< A real value.
	typedef variable_set::vindex vindex;	///< Variable identifiers (0...N-1)
	typedef variable_set::vsize vsize;    	///< Variable values (0...K-1)

	// Constructors and destructor:

	///
	/// \brief Copy-constructor.
	///
	/// Constructs a copy from an object of the same type.
	///
	factor(factor const& f) :
			v_(f.v_), t_(f.t_), c_(f.c_) {
	};

	///
	/// \brief Scalar constructor.
	///
	/// Creates a constant factor over an empty set of variables.
	///	\param s The scalar value used to initialize the table (default 1.0)
	///
	factor(const value s = 1.0) :
			v_(), t_(1, s), c_(-1) {
	};

	///
	/// \brief Constructor.
	///
	/// Creates a constant factor over a given set of variables.
	/// \param vs 	The set of variables defining the scope
	/// \param s 	The scalar value used to initialize the table (default 1.0)
	///
	factor(variable_set const& vs, value s = 1.0) : v_(vs), t_(), c_(-1) {
		t_.resize(vs.num_states());
		set_dims();
		fill(s);
	};

	///
	/// \brief Constructor.
	///
	/// Creates a factor over a given set of variables and table.
	///	\param vs 	The input set of variables
	/// \param t 	The input table
	///
	factor(variable_set const& vs, value* t) :
			v_(vs), t_(), c_(-1) {
		t_.resize(v_.num_states());
		set_dims();
		std::copy(t, t + t_.size(), t_.begin());
	};

	///
	/// \brief Class destructor
	///
	virtual ~factor() {};

	// Assignments & copy constructors:

	///
	/// \brief Assignment operator.
	///
	/// The operator performs a deep copy of the object.
	///	\param rhs	A factor object of the same type
	///	\return a reference to the current factor whose content was copied from 
	/// the factor received as argument.
	///
	factor& operator=(factor const& rhs) {
		if (this != &rhs) {
			std::vector<double> tmp;
			t_.swap(tmp);                    // force vector to release memory
			v_ = rhs.v_;
			t_ = rhs.t_;
			c_ = rhs.c_;
			set_dims();                 // then reassign
		}
		return *this;
	};

	///
	/// \brief Swap the object contents.
	///
	/// Exchange the contents of "this" and the factor received as input.
	///	\param f 	The factor to exchange the content of
	///
	void swap(factor& f) {
		if (&f != this) {
			v_.swap(f.v_);
			t_.swap(f.t_);
			int tmp = c_;
			c_ = f.c_;
			f.c_ = tmp;
		}
	};

	///
	/// \brief Clone from pointer.
	///
	/// \return a pointer to the cloned factor.
	///
	virtual factor* clone() {
		factor *f = new factor(*this);
		return f;
	};

	///
	/// \brief Set the domain sizes of the scope's variables.
 	///
	void set_dims() {
	};

	// Accessor functions:

	///
	/// \brief Size of the factor's scope.
	///
	/// \return the number of variables in the factor's scope.
	///
	size_t nvar() const {
		return v_.nvar();
	};
	
	///
	/// \brief Scope of the factor.
	///
	/// \return the scope (as set of variable ids) of the factor.
	///
	const variable_set& vars() const {
		return v_;
	};

	///
	/// \brief Scope of the factor.
	///
	/// \return the scope (as set of variable ids) of the factor.
	///
	const variable_set& variables() const {
		return v_;
	};
	
	///
	/// \brief Variable dimensions.
	///
	/// \return the domain sizes (i.e., dimensions) of the variables in 
	/// the factor's scope.
	///
	const vsize* dims() const {
		return v_.dims();
	};

	///
	/// \brief Table of factor values.
	///
	/// \return the table containing the factor values. Each value in the table
	/// corresponds to a particular configuration of the scope variables.
	///
	const value* table() const {
		return &t_[0];
	};

	///
	/// \brief Size of the factor's table.
 	///
 	/// \return the size of the table storing the factor values. It is equal to
 	/// the product of the domain sizes of the variables in the factor's scope.
 	///
	size_t num_states() const {
		return t_.size();
	};

	///
	/// \brief Size of the factor's table.
 	///
 	/// \return the size of the table storing the factor values. It is equal to
 	/// the product of the domain sizes of the variables in the factor's scope.
 	///
	size_t numel() const {
		return t_.size();
	};

	// Boolean checks on factor properties:

	///
	/// \brief Empty factor.
	///
	/// Check if the table of factor values is empty.
	/// \return *true* if the table is empty and *false* otherwise.
	///
	bool isempty() const {
		return t_.empty();
	};

	///
	/// \brief Valid factor.
	///
	/// Check if the table of factor values contains valid entries, namely
	/// it does not contain any *NaN* elements.
	/// \return *true* if all table entries are valid and *false* otherwise.
	///
	bool isnan() const {
		bool b = false;
		for (size_t i = 0; i < num_states(); i++)
			b |= isnan(t_[i]);
		return b;
	};

	///
	/// \brief Finite factor.
	///
	/// Check if the table of factor values contains *finite* real valued entries.
	/// \return *true* if all table entries are finite and *false* otherwise.
	///
	bool isfinite() const {
		bool b = true;
		for (size_t i = 0; i < num_states(); i++)
			b &= isfinite(t_[i]);
		return b;
	};

	///
	/// \brief Scalar (constant) factor.
	///
	/// Check if the factor is a *constant* (or *scalar*).
	/// \return *true* if the table size is 1 and *false* otherwise.
	///	
	bool isscalar() const {
		return numel() == 1;
	};

	// Direct value accessor:

	///
	/// \brief Access table elements (const)
	///
	/// Direct access (read-only) to a factor value.
	/// \param v 	Index of the table element.
	/// \return the table element at the given index.
	///		
	value operator[](vsize v) const {
		return t_[v];
	};

	///
	/// \brief Rewrite table elements (non-const)
	///
	/// Direct access (read and write) to a factor value.
	/// \param v 	Index of the table element.
	/// \return the non-const reference to table element at the given index.
	///		
	value& operator[](vsize v) {
		return t_[v];
	};

	///
	/// \brief Access table elements (safe)
	///
	/// Direct access (read-only) to a factor value.
	/// \param i 	Index of the table element.
	/// \return the table element at the given index.
	///			
	value get(vsize i) const {
		return t_.at(i);
	};

	///
	/// \brief Rewrite table elements (safe)
	///
	/// Re-write a factor value.
	/// \param i 	Index of the table element.
	/// \param v 	New value to be written in the table.
	///		
	void set(vsize i, value v) {
		t_.at(i) = v;
	}

	// Filling the factor table:

	///
	/// \brief Fill with constant (in-place).
	///
	/// Set all of the factor table entries to a constant value.
	/// \param v 	The constant value to be used for filling.
	/// \return the reference to the updated factor.
	///
	factor& fill(value v) {
		std::fill(t_.begin(), t_.end(), v);
		return *this;
	};

	///
	/// \brief Fill with random values in [0,1).
	///
	/// Set all of the factor table entries to random values drawn uniformly
	/// at random between 0 and 1.
	/// \return the reference to the updated factor.
	///
	factor& randomize() {
		std::generate(t_.begin(), t_.end(), randu);
		return *this;
	};

	///
	/// \brief Fill a Bayesian network CPT with 1/K values.
	///
	/// Set all the factor table enties to 1/K, where K is the domain size of
	/// the child node.
	/// \return the reference to the updated factor.
	///
	factor& fill_uniform_bayes() {
		assert(c_ >= 0);
		variable_set::const_iterator i = v_.begin();
		for (; i != v_.end(); ++i) {
			if ((*i).label() == (size_t)c_) {
				value v = 1.0/(*i).states();
				std::fill(t_.begin(), t_.end(), v);
			}
		}

		return *this;
	}

	///
	/// \brief Fill a Bayesian network CPT with random values in [0, 1).
	///
	/// Set all the factor table entries to random values drawn uniformly
	/// at random between 0 and 1. The CPT is also normalized.
	/// \return the reference to the updated factor.
	///
	factor& fill_random_bayes() {
		assert(c_ >= 0);
		std::generate(t_.begin(), t_.end(), randu);
		variable_set pa;
		variable ch;
		bool has_child = false;
		for (variable_set::const_iterator i = v_.begin(); i != v_.end(); ++i) {
			if (i->label() != (size_t)c_) {
				pa |= *i;
			} else {
				ch = *i;
				has_child = true;
			}
		}

		// Safety check
		assert(has_child);

		// Normalize the CPT
		if (pa.size() == 0) {
			normalize();
		} else {
			index_config idx1(pa); // parent configurations
			config_index idx2(v_); // factor configurations
			for (size_t i = 0; i < pa.num_states(); ++i) {

				double s = 0.0;
				std::map<size_t, size_t> cfg = idx1.convert(i);
				std::vector<size_t> indeces;
				for (size_t k = 0; k < ch.states(); ++k) {
					cfg[ch.label()] = k;
					size_t j = idx2.convert(cfg);
					indeces.push_back(j);
					s += t_[j];
				}

				for (size_t j = 0; j < indeces.size(); ++j) {
					size_t p = indeces[j];
					t_[p] = (s == 0.0) ? 0.0 : (t_[p] / s);
				}
			}
		}

		return *this;
	}

	///
	/// \brief Fill with 1/N values.
	///
	/// Set all of the factor table entries to 1/N, where N is the size of the
	/// table.
	/// \return the reference to the updated factor.
	///	
	factor& fill_uniform() {
		return fill(1.0 / num_states());
	};

	// Unary transformations (in-place); see outside class def'n for const versions

	///
	/// \brief Absolute value transformation.
	///
	/// Set all table values to their absolute values. This is an in-place
	/// trasformation.
	///	\return the reference to the transformed factor.
	///
	inline factor& abs(void) {
		std::transform(t_.begin(), t_.end(), t_.begin(), unOpAbs());
		return *this;
	};

	///
	/// \brief Exponential value transformation.
	///
	/// Set all table values to their exponential values. This is an in-place
	/// trasformation.
	///	\return the reference to the transformed factor.
	///	
	inline factor& exp(void) {
		std::transform(t_.begin(), t_.end(), t_.begin(), unOpExp());
		return *this;
	};

	///
	/// \brief Natural logarithm value transformation.
	///
	/// Set all table values to their natural log (base e) values. This is an in-place
	/// trasformation.
	///	\return the reference to the transformed factor.
	///		
	inline factor& log(void) {
		std::transform(t_.begin(), t_.end(), t_.begin(), unOpLog());
		return *this;
	};

	///
	/// \brief Logarithm (log 2) value transformation.
	///
	/// Set all table values to their log (base 2) values. This is an in-place
	/// trasformation.
	///	\return the reference to the transformed factor.
	///			
	inline factor& log2(void) {
		std::transform(t_.begin(), t_.end(), t_.begin(), unOpLogL(std::log(2)));
		return *this;
	};

	///
	/// \brief Logarithm (log 10) value transformation.
	///
	/// Set all table values to their log (base 10) values. This is an in-place
	/// trasformation.
	///	\return the reference to the transformed factor.
	///				
	inline factor& log10(void) {
		std::transform(t_.begin(), t_.end(), t_.begin(), unOpLog10());
		return *this;
	}
	;

	// Functors defined for unary operations (transformations) to the factor table:

	///
	/// \brief Functor for absolute value transformation.
	///
	struct unOpAbs {
		value operator()(value a) {
			return std::abs(a);
		}
	};

	///
	/// \brief Functor for exponential value transformation.
	///	
	struct unOpExp {
		value operator()(value a) {
			return std::exp(a);
		}
	};

	///
	/// \brief Functor for natural log value transformation.
	///
	struct unOpLog {
		value operator()(value a) {
			return std::log(a);
		}
	};

	///
	/// \brief Functor for log value transformation.
	///
	struct unOpLogL {
		value l;
		unOpLogL(value L) :
				l(L) {
		}
		;
		value operator()(value a) {
			return std::log(a) / l;
		}
	};

	///
	/// \brief Functor for log base 10 value transformation
	///
	struct unOpLog10 {
		value operator()(value a) {
			return std::log10(a);
		}
	};

	// Basic factor operations (+,-,*,/):

	///
	/// \brief Sum (const) operation for two factors (A + B).
	///
	factor operator+(const factor& B) const {
		return binaryOp(B, binOpPlus());
	};

	///
	/// \brief Sum (non-const) operation for two factors (A += B).
	///
	factor& operator+=(const factor& B) {
		return binaryOpIP(B, binOpPlus());
	};

	///
	/// \brief Minus (const) operation for two factors (A - B).
	///	
	factor operator-(const factor& B) const {
		return binaryOp(B, binOpMinus());
	};

	///
	/// \brief Minus (non-const) operation for two factors (A -= B).
	///	
	factor& operator-=(const factor& B) {
		return binaryOpIP(B, binOpMinus());
	};

	///
	/// \brief Multiplication (const) operation for two factors (A * B).
	///	
	factor operator*(const factor& B) const {
		return binaryOp(B, binOpTimes());
	};

	///
	/// \brief Multiplication (non-const) operation for two factors (A *= B).
	///	
	factor& operator*=(const factor& B) {
		return binaryOpIP(B, binOpTimes());
	};

	///
	/// \brief Division (const) operation for two factors (A / B).
	///	
	factor operator/(const factor& B) const {
		return binaryOp(B, binOpDivide());
	};

	///
	/// \brief Division (non-const) operation for two factors (A /= B).
	///	
	factor& operator/=(const factor& B) {
		return binaryOpIP(B, binOpDivide());
	};

	///
	/// \brief Sum (const) operation for a factor and a scalar (A + c).
	///
	factor operator+(const value B) const {
		return binaryOp(B, binOpPlus());
	};

	///
	/// \brief Sum (non-const) operation for a factor and a scalar (A += c).
	///	
	factor& operator+=(const value B) {
		return binaryOpIP(B, binOpPlus());
	};

	///
	/// \brief Minus (const) operation for a factor and a scalar (A - c).
	///	
	factor operator-(const value B) const {
		return binaryOp(B, binOpMinus());
	};

	///
	/// \brief Minus (non-const) operation for a factor and a scalar (A -= c).
	///	
	factor& operator-=(const value B) {
		return binaryOpIP(B, binOpMinus());
	};

	///
	/// \brief Multiplication (const) operation for a factor and a scalar (A * c).
	///	
	factor operator*(const value B) const {
		return binaryOp(B, binOpTimes());
	};

	///
	/// \brief Multiplication (non-const) operation for a factor and a scalar (A *= c).
	///
	factor& operator*=(const value B) {
		return binaryOpIP(B, binOpTimes());
	};

	///
	/// \brief Division (const) operation for a factor and a scalar (A / c).
	///	
	factor operator/(const value B) const {
		return binaryOp(B, binOpDivide());
	};

	///
	/// \brief Division (non-const) operation for a factor and a scalar (A /= c).
	///	
	factor& operator/=(const value B) {
		return binaryOpIP(B, binOpDivide());
	};

	///
	/// \brief Power (const) operation for a factor and a scalar (A ^ c).
	///	
	factor operator^(const value B) const {
		return binaryOp(B, binOpPower());
	};

	///
	/// \brief Power (non-const) operation for a factor and a scalar (A ^= c).
	///	
	factor& operator^=(const value B) {
		return binaryOpIP(B, binOpPower());
	};

	// Above operators use the following internal definitions:

	// Binary operations (eg A + B); returns new object

	///
	/// \brief Binary operation over factors.
	///
	/// Functor for binary operations between two factors A and B.
	/// \param B 	The factor to be combined with
	/// \param Op 	The binary operation
	/// \return a new factor corresponding to the operation between 
	///		the input factors A and B.
	///
	template<typename Function> factor binaryOp(const factor& B,
			Function Op) const {
		variable_set v = v_ + B.v_;  						// expand scope to union
		factor F(v);             					//  and create target factor
		subindex s1(v, v_), s2(v, B.v_); 	// index over A and B & do the op
		for (size_t i = 0; i < F.num_states(); ++i, ++s1, ++s2)
			F[i] = Op(t_[s1], B[s2]);
		return F; 										// return the new copy
	};

	///
	/// \brief Binary operation between a factor and a scalar.
	///
	/// Functor for binary operations between a factor A and a scalar B.
	/// \param B 	The scalar value to be combined with
	/// \param Op 	The binary operation
	/// \return a new factor corresponding to the operation between 
	///		the input factor A and scalar value B.
	///
	template<typename Function> factor binaryOp(const value B,
			Function Op) const {
		factor F = *this;
		F.binaryOpIP(B, Op);
		return F; // for scalar args, define with an in-place operator
	};


	// Binary in-place operations (eg A += B); returns reference to modified A
	
	///
	/// \brief Binary operation over factors (in-place).
	///
	/// Functor for in-place binary operations between two factors A and B.
	/// \param B 	The factor to be combined with
	/// \param Op 	The binary operation
	/// \return a reference to the modified factor A corresponding to the 
	/// 	operation between the input factors A and B.
	///
	template<typename Function> factor& binaryOpIP(const factor& B,
			Function Op) {
		variable_set v = v_ + B.v_;  							// expand scope to union
		if (v != v_)
			*this = binaryOp(B, Op); // if A's scope is too small, call binary op
		else {
			subindex s2(v_, B.v_);       		// otherwise create index over B
			for (size_t i = 0; i < num_states(); ++i, ++s2)
				Op.IP(t_[i], B[s2]);	// and do the operations
		}
		return *this;
	};

	///
	/// \brief Binary operation between a factor and a scalar (in-place).
	///
	/// Functor for in-place binary operations between a factor A and a scalar B.
	/// \param B 	The scalar value to be combined with
	/// \param Op 	The binary operation
	/// \return a reference to the modified factor A corresponding to the 
	/// 	operation between  the input factor A and scalar value B.
	///
	template<typename Function> factor& binaryOpIP(const value B, Function Op) {
		for (size_t i = 0; i < num_states(); i++)
			Op.IP(t_[i], B);
		return *this;	// simplifies for scalar args
	};

	// Functors defined for binary operations on the factor table : Op(a,b) and Op.IP(a,b) (in-place version)
	
	///
	/// \brief Functor for binary operation + (summation) on the factor table.
	///
	struct binOpPlus {
		value operator()(value a, const value b) {
			return a + b;
		}
		;
		value& IP(value& a, const value b) {
			return a += b;
		}
		;
	};

	///
	/// \brief Functor for binary operation - (substraction) on the factor table.
	///
	struct binOpMinus {
		//value  operator()(value  a, const value b) { return a-b; };
		//value&         IP(value& a, const value b) { return a-=b;};
		value operator()(value a, const value b) {
			return (b != -infty()) ? a - b : a;
		}
		;
		value& IP(value& a, const value b) {
			return (b != -infty()) ? a -= b : a;
		}
		;
	};

	///
	/// \brief Functor for binary operation * (multiplication) on the factor table.
	///	
	struct binOpTimes {
		value operator()(value a, const value b) {
			return a * b;
		}
		;
		value& IP(value& a, const value b) {
			return a *= b;
		}
		;
	};

	///
	/// \brief Functor for binary operation / (division) on the factor table.
	///	
	struct binOpDivide {
		value operator()(value a, const value b) {
			return (b) ? a / b : 0;
		}
		;
		value& IP(value& a, const value b) {
			return (b) ? a /= b : a = 0;
		}
		;
	};

	///
	/// \brief Functor for binary operation ^ (power) on the factor table.
	///		
	struct binOpPower {
		value operator()(value a, const value b) {
			return std::pow(a, b);
		}
		;
		value& IP(value& a, const value b) {
			return a = std::pow(a, b);
		}
		;
	};

	// Partition function, entropy, and normalization:

	///
	/// \brief Normalize the factor (const).
	///
	/// \return a copy of the factor resulted from the normalization.
	///
	factor normalized() const {
		factor F = *this;
		F.normalize();
		return F;
	};

	///
	/// Normalize the factor (non-const)
 	///
 	/// \return a reference to the modified factor as a result of normalization.
 	///
	factor& normalize() {
		double Z = sum();
		if (Z != 0)
			*this /= Z;
		return *this;
	};

	///
	/// \brief Log-partition function.
	///
	/// \return the value that represents the log-partition function of the factor.
	///
	double logpartition() const {
		return std::log(sum());
	};

	//double logpartition() const { return std::log(sum())/std::log(2.0); };

	///
	/// \brief Entropy function.
	///
	/// \return the value that represents the entropy function of the factor.
	///
	value entropy(void) const {
		value H = 0, Z = 0;
		for (size_t i = 0; i < num_states(); i++) {
			Z += t_[i];
			double L = std::log(t_[i]);
			if (!isinf(L))
				H -= t_[i] * L; //std::log(t_[i]);
		}
		H /= Z;
		H += std::log(Z);
		return H;
	}
	;

	// Elimination operators (sum, max, min, ...):

	///
	/// Elimination by summation.
	///
	/// Sum out a subset of variables from the factor's scope.
	/// \param sum_out 	The set of variables to be summed-out
	/// \return a copy of the factor resulting from the variable elimination operation.
	///
	factor sum(variable_set const& sum_out) const {
		variable_set t = v_ - sum_out;
		return marginal(t);
	};

	///
	/// Eliminate all variables by summation.
	///
	/// Sum out all variables in the factor's scope.
	/// \return a constant value resulting from the variable elimination operation.
	///	
	value sum() const {
		return std::accumulate(t_.begin(), t_.end(), 0.0, std::plus<value>());
	};

	///
	/// Elimination by weighted summation.
	///
	/// Weighted-Sum out a subset of variables from the factor's scope.
	/// \param sum_out 	The set of variables to be summed-out
	/// \param pow 		The exponent of the weighted sum operator
	/// \return a copy of the factor resulting from the weighted variable elimination operation.
	///	
	factor sum_power(variable_set const& sum_out, value pow) const {
		if (pow == 1.0)
			return sum(sum_out);
		else if (pow == -infty())
			return min(sum_out);
		else if (pow == infty())
			return max(sum_out);
		else {
			factor F = *this;
			F.log();
			F *= pow;
			F = F.logsumexp(sum_out);
			F /= pow;
			F.exp();
			return F;
		}
	};

	///
	/// Elimination by summation (log-space).
	///
	/// Sum out (log-space) a subset of variables from the factor's scope.
	/// \param sum_out 	The set of variables to be summed-out
	/// \return a copy of the factor resulting from the variable elimination operation.
	///	
	factor logsumexp(const variable_set& sum_out) const {
		variable_set target = v_ - sum_out;
		factor mx = maxmarginal(target);
		factor Scaled = *this - mx;
		Scaled.exp();
		mx += Scaled.marginal(target).log();
		return mx;
	};

	///
	/// \brief Sigma operator.
	///
	/// The special *sigma* operator used by weighed mini-buckets for 
	/// marginal MAP queries.
	///
	factor sigma(size_t n) {
		factor F = *this;
		F /= F.max();
		F ^= (double)n;
		return F;
	}

	///
	/// Elimination by maximization.
	///
	/// Max out a subset of variables from the factor's scope.
	/// \param sum_out 	The set of variables to be maxed-out
	/// \return a copy of the factor resulting from the variable elimination operation.
	///
	factor max(variable_set const& sum_out) const {
		variable_set t = v_ - sum_out;
		return maxmarginal(t);
	};

	///
	/// Maximum value.
	///
	/// Max out all variables in the factor's scope.
	/// \return the maximum value of the factor.
	///		
	value max() const {
		value const & (*max)(value const &, value const &) = std::max<value>;
		return std::accumulate(t_.begin(), t_.end(), -infty(), max);
	};

	///
	/// Elimination by minimization.
	///
	/// Minmize out a subset of variables from the factor's scope.
	/// \param sum_out 	The set of variables to be minmized-out
	/// \return a copy of the factor resulting from the variable elimination operation.
	///	
	factor min(variable_set const& sum_out) const {
		variable_set t = v_ - sum_out;
		return minmarginal(t);
	};

	///
	/// Minimum value.
	///
	/// Minimize out all variables in the factor's scope.
	/// \return the minimum value of the factor.
	///			
	value min() const {
		value const & (*min)(value const &, value const &) = std::min<value>;
		return std::accumulate(t_.begin(), t_.end(), infty(), min);
	};

	///
	/// \brief Argmax operation.
	///
	/// \return the argument of the maximization operation.
	///
	size_t argmax() const {
		return std::distance(t_.begin(), std::max_element(t_.begin(), t_.end()));
	};

	///
	/// \brief Argmin operation.
	///
	/// \return the argument of the minimization operation.
	///	
	size_t argmin() const {
		return std::distance(t_.begin(), std::min_element(t_.begin(), t_.end()));
	};


	///
	/// \brief Condition on a single variable.
	///
	/// \param v_rem 	The variable to be conditioned on
	/// \param v_state	The value of the conditioning variable
	/// \return a copy of the factor resulting from the conditioning operation.
	///
	factor condition(const variable_set& v_rem, vsize v_state) const {
		variable_set v_keep = vars() - v_rem;
		factor F(v_keep, 0.0);
		//superindex sup(src,vKeep,vState); size_t N=vKeep.nrStates();
		//for (size_t i=0;i<N;++i,++sup) F[i]=t_[sup];
		subindex src(vars(), v_rem), dst(vars(), v_keep);
		for (size_t i = 0; i < num_states(); ++i, ++src, ++dst) {
			if (src == v_state) {
				F[dst] = t_[i];  // !!! terrible; needs supindex
			}
		}
		return F;
	};

	///
	/// \brief Condition on an evidence vector
	///
	/// \param evidence	The evidence vector
	/// \return a copy of the factor resulting from the conditioning operation.
	///
	factor condition(std::vector<int>& evidence) const {
		variable_set v_rem;
		for (variable_set::const_iterator i = vars().begin();
				i != vars().end(); ++i) {
			int v = i->label();
			if (evidence[v] >= 0) {
				v_rem |= (*i);
			}
		}

		factor F = *this;
		for (variable_set::const_iterator i = v_rem.begin();
				i != v_rem.end(); ++i) {
			variable_set vs(*i);
			int val = evidence[i->label()];
			F = F.condition(vs, val);
		}

		return F;
	}

	///
	/// \brief Condition on a single variable.
	///
	/// \param v_rem 	The variable to be conditioned on
	/// \param v_state	The value of the conditioning variable
	/// \return a copy of the factor resulting from the conditioning operation.
	///	
	factor slice(const variable_set& v_rem, vsize v_state) const {
		return condition(v_rem, v_state);
	};

	///
	/// \brief Embed extra variables in the factor.
	///
	/// \param v 	The set variables to be embeded
	/// \return a copy of the factor resulting from embedding the new variables.
	///
	factor embed(const variable_set& v) const {
		if (vars() == v)
			return *this;
		else
			return (*this) + factor(v / vars(), 0.0);
	};

	///
	/// \brief Sample the distribution (factor). 
	///
	size_t sample() const {
		double x = 0.0, y = randu() * sum();
		for (size_t i = 0; i < num_states(); ++i)
			if ((x += t_[i]) > y)
				return i;
		return num_states();
	};

	///
	/// \brief Sample the distribution (factor). 
	///	
	vsize draw() const {
		return sample();
	};

	///
	/// \brief Marginal over a set of varibles.
	///
	/// Compute the marginal over a subset of variables by summing out
	/// the remaining variables in the factor's scope. 
	/// \param target 	The scope of the marginal
	/// \return a copy of the factor resulting from marginalization.
	///
	factor marginal(variable_set const& target) const {
		factor F(target & vars(), 0.0);
		subindex s(v_, F.vars());
		for (size_t i = 0; i < num_states(); ++i, ++s)
			F[s] += t_[i];
		return F;
	};

	///
	/// \brief Weighted marginal over a set of varibles.
	///
	/// Compute the weighted marginal over a subset of variables by applying
	/// the weighted sum operator on the remaining variables in the factor's scope.
	/// \param target 	The scope of the marginal
	/// \param w 		The weight used by the weighted sum operator
	/// \return a copy of the factor resulting from marginalization.
	///
	factor marginal(variable_set const& target, const value w) const {
		if (w == infty()) { // max-marginal
			return maxmarginal(target);
		} else { // weighted marginal
			factor FF = *this;
			FF ^= (1.0/w);
			factor F(target & vars(), 0.0);
			subindex s(v_, F.vars());
			for (size_t i = 0; i < num_states(); ++i, ++s)
				F[s] += FF[i];
			return F;
		}
	};

	///
	/// \brief Max-Marginal over a set of varibles.
	///
	/// Compute the max-marginal over a subset of variables by maximizing out
	/// the remaining variables in the factor's scope. 
	/// \param target 	The scope of the marginal
	/// \return a copy of the factor resulting from marginalization.
	///
	factor maxmarginal(variable_set const& target) const {
		factor F(target & vars(), -infty());
		subindex s(v_, F.vars());
		for (size_t i = 0; i < num_states(); ++i, ++s)
			F[s] = (F[s] > t_[i]) ? F[s] : t_[i];
		return F;
	};

	///
	/// \brief Min-Marginal over a set of varibles.
	///
	/// Compute the min-marginal over a subset of variables by minimizing out
	/// the remaining variables in the factor's scope. 
	/// \param target 	The scope of the marginal
	/// \return a copy of the factor resulting from marginalization.
	///	
	factor minmarginal(variable_set const& target) const {
		factor F(target & vars(), infty());
		subindex s(v_, F.vars());
		for (size_t i = 0; i < num_states(); ++i, ++s)
			F[s] = (F[s] > t_[i]) ? t_[i] : F[s];
		return F;
	}
	;

	///
	/// \brief Distance measures.
	///
	MER_ENUM( Distance , L1,L2,LInf,KL,HPM,MAS,OptGap );

	///
	/// \brief Distance between two factors.
	///
	/// Compute the distance between the distributions corresponding to the
	/// two factors (e.g., KL-distance).
	/// \param F2 		The factor to compute the distance from
	/// \param type 	The distance measure
	/// \return a real value representing the distance between the two factors.
	///
	double distance(factor const& F2, Distance type = Distance::L2) const {
		assert(vars() == F2.vars());
		factor F(*this), Ftmp;               // make a copy for manipulation
		double dist = -1.0;                   // local variables
		value Z;
		switch (type) {
		case Distance::L2:                // L2, sum of squared errors
			F -= F2;
			F *= F;
			dist = F.sum();
			break;
		case Distance::L1:               // L1, sum of absolute errors  (=TV !!)
			F -= F2;
			dist = F.abs().sum();
			break;
		case Distance::LInf:              // L-infinity, max absolute error
			F -= F2;
			dist = F.abs().max();
			break;
		case Distance::KL:                // KL-divergence (relative entropy)
			Z = sum();
			F /= F2;
			F *= F2.sum() / Z;
			F.log();
			F *= *this;
			dist = F.sum() / Z;
			break;
		case Distance::HPM:               // Hilbert's projective metric
			F /= F2;
			F.log();
			dist = F.max() - F.min(); //   aka "dynamic range"
			break;
		case Distance::MAS:               // "MAS" error value (not a metric)
			F.log();
			Ftmp = F2;
			F /= Ftmp.log();
			dist = std::max(F.max(), 1.0 / F.min()) - 1.0;
			break;
		case Distance::OptGap:            // "Primal/Dual Gap"-like
			double mx1, mx2, gap1, gap2;
			mx1 = mx2 = gap1 = gap2 = 0.0;
			for (size_t i = 0; i < num_states(); ++i) {
				if (mx1 < F[i]) {
					gap2 = F2[i];
					mx1 = F[i];
				} else if (mx1 == F[i])
					gap2 = std::min(gap2, F2[i]);
				if (mx2 < F2[i]) {
					gap1 = F[i];
					mx2 = F2[i];
				} else if (mx2 == F2[i])
					gap1 = std::min(gap1, F[i]);
			}
			return (mx1 - gap1) + (mx2 - gap2);
			break;
			//case Distance::Hellinger:   (!!)
			//  F^=0.5; F-=F2^0.5; F*=F; dist=(0.5*F.sum())^0.5;  // straightforward computation
			//  F*=F2; F^=0.5; dist=(1-F.sum())^0.5;              // alternate computation if F,F2 normalized
			//	break;
		default:
			throw std::runtime_error("Invalid distance type");
		}
		return dist;
	};

	///
	/// \brief Norm of the factor.
	///
	/// Compute the norm of the factor (e.g., L2 norm).
	/// \param type 	The distance measure
	/// \return a real value representing the norm of the factor.
	///
	double norm(Distance type = Distance::L2) const {
		factor F(*this);                // make a copy for manipulation
		double dist = -1.0;               //
		switch (type) {
		case Distance::L2:                // L2, sum of squared errors
			F *= F;
			dist = F.sum();
			break;
		case Distance::L1:                // L1, sum of absolute errors
			dist = F.abs().sum();
			break;
		case Distance::LInf:                // L-infinity, max absolute error
			dist = F.abs().max();
			break;
		case Distance::KL:       // KL-divergence (relative entropy => entropy?)
			return entropy();
			break;
		case Distance::HPM:               // Hilbert's projective metric
			F.log();
			dist = F.max() - F.min();      //   aka "dynamic range"
			break;
		default:
			throw std::runtime_error("Invalid norm type");
		}
		return dist;
	};

	///
	/// \brief Decomposition methods.
 	///
	MER_ENUM( Decomp , L2,L2_HPM,L2_MAS );

	///
	/// \brief Decompose a factor into a sum of smaller factors.
	///
	/// \param vlist 		The list of scopes to be used in the decomposition
	/// \param method 		The decomposition method
	/// \return a vector of new factors representing the decomposition.
	std::vector<factor> decomp_sum(std::vector<variable_set> vlist,
			factor::Decomp method) const {
		int nF = vlist.size();
		double mx, mn;
		std::vector<factor> Flist(nF);

		factor tmp, F = *this;
		switch (method) {
		case Decomp::L2: //L2
			double Cn, Cd;
			Cd = F.numel();
			Cn = F.sum(); // /Cd*(1-1.0/nF);
			for (int j = 0; j < nF; j++) {
				Flist[j] = F.marginal(vlist[j]);
				double D = Cd / Flist[j].numel();
				Flist[j] /= D;
				Flist[j] -= Cn / Cd * (1.0 - 1.0 / (nF - j));
				F -= Flist[j];
				Cn -= Flist[j].sum() * D;
			}
			break;
		case Decomp::L2_HPM: //L2+HPM
			Flist = decomp_sum(vlist, Decomp::L2);
			for (int j = 0; j < nF; j++)
				F -= Flist[j];
			mx = F.max();
			mn = F.min();
			for (int j = 0; j < nF; j++)
				Flist[j] += (mx + mn) / 2 / nF;
			break;
		case Decomp::L2_MAS: //L2+MAS
			Flist = decomp_sum(vlist, Decomp::L2);
			F = Flist[0];
			for (int j = 1; j < nF; j++)
				F += Flist[j];
			F /= *this;
			F.log();
			mx = F.max();
			mn = F.min();
			for (int j = 0; j < nF; j++)
				Flist[j] *= std::exp(-(mx + mn) / 2 / nF);
			break;
		}
		return Flist;
	}

	///
	/// \brief Decompose a factor into a product of factors.
	///
	/// \param vlist 		The list of scopes to be used in the decomposition
	/// \param method 		The decomposition method
	/// \return a vector of new factors representing the decomposition.
	std::vector<factor> decomp_prod(std::vector<variable_set> vlist,
			factor::Decomp method) const {
		factor F = *this;
		F.log();
		std::vector<factor> Flist = F.decomp_sum(vlist, method);
		for (size_t j = 0; j < vlist.size(); j++)
			Flist[j].exp();
		return Flist;
	}

	///
	/// \brief Output operator (friend).
	///
	/// Write the (formatted) content of the factor to an output stream.
	/// \param out 		The output stream
	/// \param f 		The Factor to be written out
	/// \return a reference to the modified output stream containing the 
	///		content of the factor received as input.
	///
	friend std::ostream& operator<<(std::ostream& out, const factor& f) {
		out << "Factor over " << f.variables() << ":";
		for (size_t j = 0; j < f.t_.size(); j++)
			out << " " << f.t_[j];
		return out;
	};

	void set_child(int c) {
		c_ = c;
	}

	int get_child() const {
		return c_;
	}

protected:

	variable_set v_;					///< Variable list vector (*scope*).
	std::vector<value> t_;				///< Table of values.
	int c_;								///< Index of the child variable (for Bayes nets only), -1 for Markov models

	///
	/// \brief Calculate the factor table size.
 	///
 	/// Compute the actual size of the factor table by multiplying the
 	/// domain sizes of the variables in its scope.
 	/// \return the factor table size.
 	///
	vsize calc_numel() const {
		vsize n = 1;
		vsize const* d = dims();
		for (size_t i = 0; i < nvar(); i++)
			n *= d[i];
		return (n > 1) ? n : 1;
	}
	;

	///
	/// \brief Check if a real value is finite.
	///	
	static inline bool isfinite(double v) {
		return (v <= DBL_MAX && v >= -DBL_MAX);
	};

	///
	/// \brief Check if a real value is not-a-number.
	///
	static inline bool isnan(value v) {
		return (v != v);
	};

	///
	/// \brief Return the infinity numerical limit.
	///
	static inline value infty() {
		return std::numeric_limits<value>::infinity();
	};
};

// "Static" functions that operate on Factor class variables

inline factor abs(const factor& A) {
	factor F = A;
	F.abs();
	return F;
}

inline factor exp(const factor& A) {
	factor F = A;
	F.exp();
	return F;
}

inline factor log(const factor& A) {
	factor F = A;
	F.log();
	return F;
}

inline factor log2(const factor& A) {
	factor F = A;
	F.log();
	F /= log(2.0);
	return F;
}

inline factor log10(const factor& A) {
	factor F = A;
	F.log10();
	return F;
}


template<class InputIterator>
inline factor mean(InputIterator first, InputIterator last) {
	size_t N = 0;
	factor F(0.0);
	for (; first != last; ++first, ++N)
		F += *first;
	if (N)
		F /= N;
	return F;
}
template<class InputIterator>
inline factor geomean(InputIterator first, InputIterator last) {
	size_t N = 0;
	factor F(1.0);
	for (; first != last; ++first, ++N)
		F *= *first;
	if (N)
		F ^= 1.0 / N;
	return F;
}

} // namespace

#endif /* IBM_MERLIN_FACTOR_H_ */
