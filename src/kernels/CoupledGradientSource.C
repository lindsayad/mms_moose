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

#include "CoupledGradientSource.h"

template <>
InputParameters
validParams<CoupledGradientSource>()
{
  InputParameters params = validParams<Kernel>();
  params.addRequiredCoupledVar("v", "the coupled var");
  return params;
}

CoupledGradientSource::CoupledGradientSource(const InputParameters & parameters)
  : Kernel(parameters), _grad_v(coupledGradient("v")), _v_id(coupled("v"))
{
}

Real
CoupledGradientSource::computeQpResidual()
{
  return -_grad_v[_qp] * _grad_v[_qp] * _test[_i][_qp];
}

Real
CoupledGradientSource::computeQpJacobian()
{
  return 0;
}

Real
CoupledGradientSource::computeQpOffDiagJacobian(unsigned jvar)
{
  if (jvar == _v_id)
    return -2. * _grad_v[_qp] * _grad_phi[_j][_qp] * _test[_i][_qp];

  else
    return 0;
}
