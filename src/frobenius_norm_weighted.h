/* Routines for linear multiple output using sparse group lasso regression.
 Intended for use with R.
 Copyright (C) 2014 Martin Vincent

 This program is free software: you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation, either version 3 of the License, or
 (at your option) any later version.

 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.

 You should have received a copy of the GNU General Public License
 along with this program.  If not, see <http://www.gnu.org/licenses/>
 */

using namespace sgl;

#ifndef FOBENIUS_NORM_WEIGHTED_HPP_
#define FOBENIUS_NORM_WEIGHTED_HPP_

//type_X : sgl::matrix or sgl::sparse_matrix
//type_Y : sgl::matrix or sgl::sparse_matrix

template < typename type_X, typename type_Y >
class FrobeniusLossWeighted {

public:

	const natural n_variables; // number of parameters

private:

	type_Y const& Y; //response - matrix of size n_samples x n_responses
	matrix const& W; //vector of size n_samples x n_responses

	matrix lp; //linear predictors - matrix of size n_samples x n_responses

public:

	typedef hessian_diagonal<false> hessian_type;

	typedef DataPackage_3<
		MatrixData<type_X>,
		MultiResponse<type_Y, 'Y'>,
		Data<sgl::matrix, 'W'> > data_type;

	FrobeniusLossWeighted() 	:
		n_variables(0),
		Y(null_matrix),
		W(null_matrix),
		lp(null_matrix)	{}

	FrobeniusLossWeighted(data_type const& data) :
		n_variables(data.get_B().n_responses),
		Y(data.get_B().response),
		W(data.get_C().data),
		lp(data.get_A().n_samples, n_variables) {}

	void set_lp(matrix const& lp) {
		this->lp = lp;
	}

	void set_lp_zero() {
		lp.zeros();
	}

	const matrix gradients() const {
		return static_cast<double>(2)*trans(W%(lp-Y));
	}

	void compute_hessians() const	{
		return;
	}

  const hessian_type::representation hessians(natural i) const {
		return static_cast<double>(2)*trans(W.row(i));
	}

	const numeric sum_values() const {
		return accu(W%(lp-Y)%(lp-Y));
	}

};

typedef ObjectiveFunctionType <
	GenralizedLinearLossDense < 
		FrobeniusLossWeighted < matrix, matrix > > > frobenius_w
;

typedef ObjectiveFunctionType <
	GenralizedLinearLossSparse <
		FrobeniusLossWeighted < sparse_matrix, matrix > > > frobenius_w_spx
;

typedef ObjectiveFunctionType <
	GenralizedLinearLossDense <
		FrobeniusLossWeighted < matrix, sparse_matrix > > > frobenius_w_spy
;

typedef ObjectiveFunctionType <
	GenralizedLinearLossSparse <
		FrobeniusLossWeighted < sparse_matrix, sparse_matrix > > > frobenius_w_spx_spy
;

#endif /* FOBENIUS_NORM_WEIGHTED_HPP_ */
