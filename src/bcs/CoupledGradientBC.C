/****************************************************************/
/*               DO NOT MODIFY THIS HEADER                      */
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*           (c) 2010 Battelle Energy Alliance, LLC             */
/*                   ALL RIGHTS RESERVED                        */
/*                                                              */
/*          Prepared by Battelle Energy Alliance, LLC           */
/*            Under Contract No. DE-AC07-05ID14517              */
/*            With the U. S. Department of Energy               */
/*                                                              */
/*            See COPYRIGHT for full restrictions               */
/****************************************************************/

#include "CoupledGradientBC.h"

template <>
InputParameters
validParams<CoupledGradientBC>()
{
  InputParameters params = validParams<IntegratedBC>();
  params.addRequiredCoupledVar("v", "The coupled variable supplying the gradient.");
  return params;
}

CoupledGradientBC::CoupledGradientBC(const InputParameters & parameters)
  : IntegratedBC(parameters), _grad_v(coupledGradient("v")), _v_id(coupled("v"))
{
}

Real
CoupledGradientBC::computeQpResidual()
{
  return -_test[_i][_qp] * _normals[_qp] * _grad_v[_qp];
}

Real
CoupledGradientBC::computeQpOffDiagJacobian(unsigned jvar)
{
  if (jvar == _v_id)
    return -_test[_i][_qp] * _normals[_qp] * _grad_phi[_j][_qp];

  else
    return 0.;
}
