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

#ifndef FOBENIUS_NORM_HPP_
#define FOBENIUS_NORM_HPP_

using namespace sgl;

//type_X : sgl::matrix or sgl::sparse_matrix
//type_Y : sgl::matrix or sgl::sparse_matrix

template < typename type_X, typename type_Y >
class FrobeniusLoss {

public:

	const natural n_samples;
	const natural n_responses;
	const natural n_variables; // number of parameters

private:

	type_Y const& Y; //response - matrix of size n_samples x n_responses
	matrix lp; //linear predictors - matrix of size n_samples x n_responses

public:

	typedef hessian_identity<true> hessian_type; //constant hessians of type double * Id

	typedef DataPackage_2<
		MatrixData<type_X>,
		MultiResponse<type_Y, 'Y'> > data_type
	;

	FrobeniusLoss()	:
	 	n_samples(0),
		n_responses(0),
		n_variables(0),
		Y(null_matrix),
		lp(n_samples, n_responses)	{}

	FrobeniusLoss(data_type const& data) :
		n_samples(data.get_A().n_samples),
		n_responses(data.get_B().n_responses),
		n_variables(n_responses),
		Y(data.get_B().response),
		lp(n_samples, n_responses) {}

	void set_lp(matrix const& lp) {
		this->lp = lp;
	}

	void set_lp_zero() {
		lp.zeros();
	}

	const matrix gradients() const {
		return static_cast<double>(2)/static_cast<double>(n_samples)*trans(lp-Y);
	}

	void compute_hessians() const {
		return;
	}

  const double hessians(natural i) const {
		return static_cast<double>(2)/static_cast<double>(n_samples);
	}

	const numeric sum_values() const	{
		return accu((lp-Y)%(lp-Y))/static_cast<double>(n_samples);
	}

};

// Std design matrix
typedef ObjectiveFunctionType <
	GenralizedLinearLossDense <
		FrobeniusLoss < matrix, matrix > > > frobenius
;

typedef ObjectiveFunctionType <
	GenralizedLinearLossSparse <
		FrobeniusLoss < sparse_matrix, matrix > > > frobenius_spx
;

typedef ObjectiveFunctionType <
	GenralizedLinearLossDense <
		FrobeniusLoss < matrix, sparse_matrix > > > frobenius_spy
;

typedef ObjectiveFunctionType <
	GenralizedLinearLossSparse <
		FrobeniusLoss < sparse_matrix, sparse_matrix > > > frobenius_spx_spy
;

#endif /* FOBENIUS_NORM_HPP_ */
